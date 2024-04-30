#!/usr/bin/env python3
"""Run GATK on a locally mutated GIAB BAM and score against mutation truth."""

from __future__ import annotations

import argparse
import gzip
import json
import platform
import subprocess
import sys
from pathlib import Path
from typing import Any, Dict, Tuple

from score_giab_bam import (
    BCFTOOLS_IMAGE,
    GATK_IMAGE,
    SAMTOOLS_IMAGE,
    REPO_DIR,
    TRUTH_VCF_URL,
    _build_gatk_flags,
    _docker_base,
    _read_config,
    _run,
    _safe_region,
)


SCORE_DIR = Path(__file__).resolve().parent
if str(REPO_DIR) not in sys.path:
    sys.path.insert(0, str(REPO_DIR))

from utils.scoring import AdvancedScorer, HAPPY_DOCKER_IMAGE, HappyScorer  # noqa: E402


def _latest_mutated_manifest() -> Path:
    manifests = sorted(
        (SCORE_DIR / "giab_mutated").glob("*.manifest.json"),
        key=lambda p: p.stat().st_mtime,
    )
    manifests = [p for p in manifests if not p.name.startswith("HISEQ")]
    if not manifests:
        raise SystemExit("No mutated GIAB manifest found. Run ./mutate_giab_bam.sh first.")
    return manifests[-1]


def _host_path(path_value: str) -> Path:
    if path_value.startswith("/workspace/"):
        return SCORE_DIR.parent / path_value.removeprefix("/workspace/")
    return Path(path_value)


def _count_vcf(path: Path) -> int:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as handle:
        return sum(1 for line in handle if line.strip() and not line.startswith("#"))


def _read_vcf_records(path: Path) -> list[tuple[str, int, str, str]]:
    opener = gzip.open if str(path).endswith(".gz") else open
    records = []
    with opener(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            chrom, pos_s, _id, ref, alts = fields[:5]
            for alt in alts.split(","):
                records.append((chrom, int(pos_s), ref, alt))
    return records


def _docker_cmd(*mounts: Tuple[Path, str], entrypoint: str | None = None) -> list[str]:
    cmd = ["docker", "run", "--rm"]
    if platform.machine() == "arm64":
        cmd.extend(["--platform", "linux/amd64"])
    if entrypoint:
        cmd.extend(["--entrypoint", entrypoint])
    for host, container in mounts:
        cmd.extend(["-v", f"{host.resolve()}:{container}"])
    return cmd


def _prepare_vcf_for_happy(path: Path, out_dir: Path) -> Path:
    prepared = out_dir / f"{path.stem}.happy_input.vcf"
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8") as src, prepared.open("w", encoding="utf-8") as dst:
        has_info = False
        has_type_info = False
        has_format = False
        has_contig = False
        for line in src:
            if line.startswith("##INFO=<ID=SYNTHETIC"):
                has_info = True
            if line.startswith("##INFO=<ID=TYPE"):
                has_type_info = True
            if line.startswith("##FORMAT=<ID=GT"):
                has_format = True
            if line.startswith("##contig="):
                has_contig = True

            if line.startswith("#CHROM"):
                if not has_info:
                    dst.write("##INFO=<ID=SYNTHETIC,Number=0,Type=Flag,Description=\"Local synthetic injected mutation\">\n")
                if not has_type_info:
                    dst.write("##INFO=<ID=TYPE,Number=1,Type=String,Description=\"Local synthetic variant type\">\n")
                if not has_format:
                    dst.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
                if not has_contig:
                    chrom = None
                    with opener(path, "rt", encoding="utf-8") as scan:
                        for scan_line in scan:
                            if scan_line.startswith("#"):
                                continue
                            chrom = scan_line.split("\t", 1)[0]
                            break
                    if chrom:
                        dst.write(f"##contig=<ID={chrom}>\n")

                fields = line.rstrip("\n").split("\t")
                if len(fields) < 10:
                    dst.write(line.rstrip("\n") + "\tFORMAT\tHG002\n")
                else:
                    dst.write(line)
                continue

            if line.startswith("#"):
                dst.write(line)
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) == 8:
                dst.write(line.rstrip("\n") + "\tGT\t1/1\n")
            else:
                dst.write(line)
    return prepared


def _bgzip_and_index_vcf(path: Path, out_dir: Path) -> Path:
    path = _prepare_vcf_for_happy(path, out_dir)
    if str(path).endswith(".vcf.gz"):
        gz_path = path
        if not Path(str(gz_path) + ".tbi").exists():
            _run(_docker_base((gz_path.parent, "/vcf")) + [
                BCFTOOLS_IMAGE,
                "bcftools", "index", "-t", f"/vcf/{gz_path.name}",
            ])
        return gz_path

    gz_path = out_dir / f"{path.stem}.vcf.gz"
    _run(_docker_base((path.parent, "/in"), (out_dir, "/out")) + [
        BCFTOOLS_IMAGE,
        "bcftools", "view",
        f"/in/{path.name}",
        "-Oz",
        "-o", f"/out/{gz_path.name}",
    ])
    _run(_docker_base((out_dir, "/out")) + [
        BCFTOOLS_IMAGE,
        "bcftools", "index", "-t", f"/out/{gz_path.name}",
    ])
    return gz_path


def _slice_giab_truth(region: str, out_dir: Path, truth_vcf_url: str) -> Path:
    truth = out_dir / "giab.truth.vcf.gz"
    _run(_docker_base((out_dir, "/out")) + [
        BCFTOOLS_IMAGE,
        "bcftools", "view",
        "-r", region,
        truth_vcf_url,
        "-Oz",
        "-o", "/out/giab.truth.vcf.gz",
    ])
    _run(_docker_base((out_dir, "/out")) + [
        BCFTOOLS_IMAGE,
        "bcftools", "index", "-t", "/out/giab.truth.vcf.gz",
    ])
    return truth


def _merge_truth_vcfs(giab_truth: Path, synthetic_truth: Path, out_dir: Path) -> Path:
    merged = out_dir / "merged.truth.vcf.gz"
    _run(_docker_base((out_dir, "/out")) + [
        BCFTOOLS_IMAGE,
        "bcftools", "concat",
        "-a",
        "-D",
        "/out/giab.truth.vcf.gz",
        f"/out/{synthetic_truth.name}",
        "-Oz",
        "-o", "/out/merged.truth.vcf.gz",
    ])
    _run(_docker_base((out_dir, "/out")) + [
        BCFTOOLS_IMAGE,
        "bcftools", "index", "-t", "/out/merged.truth.vcf.gz",
    ])
    return merged


def _ensure_reference_sdf(reference: Path) -> Path:
    sdf = reference.parent / f"{reference.stem}.sdf"
    if sdf.exists() and sdf.is_dir():
        return sdf

    print(f"   Creating RTG SDF reference: {sdf}", flush=True)
    _run(_docker_cmd((reference.parent, "/ref"), entrypoint="rtg") + [
        HAPPY_DOCKER_IMAGE,
        "format",
        "-o", f"/ref/{sdf.name}",
        f"/ref/{reference.name}",
    ])
    return sdf


def _score_synthetic_only(output_vcf: Path, truth_vcf: Path) -> dict:
    truth = set(_read_vcf_records(truth_vcf))
    query = set(_read_vcf_records(output_vcf))
    truth_positions = {(chrom, pos) for chrom, pos, _ref, _alt in truth}
    query_at_truth_positions = {
        rec for rec in query if (rec[0], rec[1]) in truth_positions
    }

    tp = len(query_at_truth_positions & truth)
    fp = len(query_at_truth_positions - truth)
    fn = len(truth - query_at_truth_positions)
    precision = tp / (tp + fp) if tp + fp else 0.0
    recall = tp / (tp + fn) if tp + fn else 0.0
    f1 = 2 * precision * recall / (precision + recall) if precision + recall else 0.0
    region_fp = len(query - query_at_truth_positions)

    return {
        "tp_snp": tp,
        "fp_snp": fp,
        "fn_snp": fn,
        "precision_snp": precision,
        "recall_snp": recall,
        "f1_snp": f1,
        "tp_indel": 0,
        "fp_indel": 0,
        "fn_indel": 0,
        "precision_indel": 0.0,
        "recall_indel": 0.0,
        "f1_indel": 0.0,
        "truth_total_snp": len(truth),
        "query_total_snp": len(query_at_truth_positions),
        "region_fp_total": region_fp,
        "simple_score": 100.0 * f1,
        "combined_final": f1,
    }


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Score GATK calls from a locally mutated GIAB BAM against synthetic truth"
    )
    parser.add_argument("--manifest", default=str(_latest_mutated_manifest()))
    parser.add_argument("--config", default=str(REPO_DIR / "configs" / "gatk.conf"))
    parser.add_argument("--reference", default=str(REPO_DIR / "datasets" / "reference" / "chr20" / "chr20.fa"))
    parser.add_argument("--out-dir", default=str(SCORE_DIR / "giab_mutated_score"))
    parser.add_argument("--giab-truth-vcf-url", default=TRUTH_VCF_URL)
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--memory-gb", type=int, default=8)
    args = parser.parse_args()

    manifest_path = Path(args.manifest)
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    region = manifest["region"]
    mutated_bam = _host_path(manifest["mutated_bam"])
    truth_vcf = _host_path(manifest["truth_vcf"])
    mutations_vcf = _host_path(manifest.get("mutations_vcf", manifest["truth_vcf"]))
    reference = Path(args.reference)
    config = _read_config(Path(args.config))

    out_dir = Path(args.out_dir) / _safe_region(region)
    out_dir.mkdir(parents=True, exist_ok=True)
    rg_bam = out_dir / "input.with_rg.bam"
    output_vcf = out_dir / "gatk.output.vcf.gz"

    # Add a local RG/SM tag for GATK, because the public GIAB BAM slice lacks RG.
    _run(_docker_base((mutated_bam.parent, "/bam"), (out_dir, "/out")) + [
        SAMTOOLS_IMAGE,
        "samtools", "addreplacerg",
        "-r", r"@RG\tID:HG002_MUT\tSM:HG002_MUT\tPL:ILLUMINA\tLB:GIAB_MUT\tPU:unit1",
        "-o", "/out/input.with_rg.bam",
        f"/bam/{mutated_bam.name}",
    ])
    _run(_docker_base((out_dir, "/out")) + [
        SAMTOOLS_IMAGE,
        "samtools", "index",
        "/out/input.with_rg.bam",
    ])

    gatk_cmd = _docker_base(
        (reference.parent, "/ref"),
        (out_dir, "/out"),
    ) + [
        GATK_IMAGE,
        "gatk",
        "--java-options", f"-Xmx{args.memory_gb}g",
        "HaplotypeCaller",
        "-R", f"/ref/{reference.name}",
        "-I", "/out/input.with_rg.bam",
        "-O", "/out/gatk.output.vcf.gz",
        "-L", region,
        "--native-pair-hmm-threads", str(args.threads),
    ] + _build_gatk_flags(config)
    _run(gatk_cmd)

    synthetic_truth_vcf_for_happy = _bgzip_and_index_vcf(truth_vcf, out_dir)
    mutations_vcf_for_happy = _bgzip_and_index_vcf(mutations_vcf, out_dir)
    giab_truth_vcf = _slice_giab_truth(region, out_dir, args.giab_truth_vcf_url)
    truth_vcf_for_happy = _merge_truth_vcfs(giab_truth_vcf, synthetic_truth_vcf_for_happy, out_dir)
    ref_sdf = _ensure_reference_sdf(reference)

    happy_scorer = HappyScorer()
    metrics = happy_scorer.score_vcf(
        truth_vcf=str(truth_vcf_for_happy),
        query_vcf=str(output_vcf),
        reference_fasta=str(reference),
        confident_bed=None,
        region=region,
        reference_sdf=str(ref_sdf),
        mutations_vcf=str(mutations_vcf_for_happy),
    )
    if metrics is None:
        raise SystemExit("SN107 hap.py scoring failed; no valid metrics produced.")

    advanced_score = AdvancedScorer.compute_advanced_score(metrics)
    combined_final = advanced_score / 100.0
    exact_debug_metrics = _score_synthetic_only(output_vcf, truth_vcf)

    report = {
        "workflow": "GIAB BAM -> local synthetic injection -> GATK -> compare to synthetic truth",
        "scoring_system": "SN107 HappyScorer + AdvancedScorer",
        "score_schema_version": "0.1.1",
        "region": region,
        "mutated_bam": str(mutated_bam),
        "truth_vcf": str(truth_vcf),
        "merged_truth_vcf": str(truth_vcf_for_happy),
        "giab_truth_vcf": str(giab_truth_vcf),
        "mutations_vcf": str(mutations_vcf),
        "output_vcf": str(output_vcf),
        "happy_vcf": str(out_dir / f"happy_{output_vcf.stem}.vcf.gz"),
        "happy_summary_csv": str(out_dir / f"happy_{output_vcf.stem}.summary.csv"),
        "config": str(Path(args.config).resolve()),
        "reference": str(reference.resolve()),
        "reference_sdf": str(ref_sdf.resolve()),
        "truth_variant_count": _count_vcf(truth_vcf),
        "output_variant_count": _count_vcf(output_vcf),
        "metrics": metrics,
        "advanced_score": advanced_score,
        "combined_final": combined_final,
        "exact_match_debug_metrics": exact_debug_metrics,
        "note": "SN107-like local scoring: hap.py/vcfeval metrics filtered to injected mutations, then AdvancedScorer.",
    }
    report_path = out_dir / "score_report.json"
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    print("\nMutated GIAB BAM GATK score")
    print(f"  region={region}")
    print(f"  score={advanced_score:.2f}/100")
    print(f"  combined_final={combined_final:.4f}")
    print(f"  SNP F1={metrics['f1_snp']:.4f} INDEL F1={metrics['f1_indel']:.4f}")
    print(f"  truth variants={report['truth_variant_count']} output variants={report['output_variant_count']}")
    print(f"  output={output_vcf}")
    print(f"  truth={truth_vcf}")
    print(f"  hap.py summary={report['happy_summary_csv']}")
    print(f"  report={report_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
