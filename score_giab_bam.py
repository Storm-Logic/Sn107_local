#!/usr/bin/env python3
"""Run GATK on a local GIAB BAM slice and score against GIAB truth VCF."""

from __future__ import annotations

import argparse
import gzip
import json
import math
import platform
import subprocess
import sys
from pathlib import Path
from typing import Any, Dict, Iterable, Optional, Tuple


SCORE_DIR = Path(__file__).resolve().parent
REPO_DIR = SCORE_DIR.parent / "mine_107"
GATK_IMAGE = "broadinstitute/gatk:4.5.0.0"
BCFTOOLS_IMAGE = "quay.io/biocontainers/bcftools:1.20--h8b25389_0"
SAMTOOLS_IMAGE = "quay.io/biocontainers/samtools:1.20--h50ea8bc_0"
TRUTH_VCF_URL = (
    "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/"
    "AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/"
    "HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
)


def _safe_region(region: str) -> str:
    return region.replace(":", "_").replace("-", "_")


def _load_manifest(path: Path) -> Dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _latest_giab_manifest() -> Path:
    manifests = sorted((SCORE_DIR / "giab_bam").glob("*.manifest.json"), key=lambda p: p.stat().st_mtime)
    if not manifests:
        raise SystemExit("No GIAB BAM manifest found. Run ./create_giab_bam.sh first.")
    return manifests[-1]


def _read_config(path: Path) -> Dict[str, Any]:
    values: Dict[str, Any] = {}
    for line_num, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        if "=" not in stripped:
            raise ValueError(f"{path}:{line_num}: expected key=value")
        key, _, raw = stripped.partition("=")
        raw = raw.strip()
        if raw.lower() in {"true", "false"}:
            value: Any = raw.lower() == "true"
        else:
            try:
                value = int(raw)
            except ValueError:
                try:
                    value = float(raw)
                except ValueError:
                    value = raw
        values[key.strip()] = value
    return values


def _docker_base(*mounts: Tuple[Path, str]) -> list[str]:
    cmd = ["docker", "run", "--rm"]
    if platform.machine() == "arm64":
        cmd.extend(["--platform", "linux/amd64"])
    for host, container in mounts:
        cmd.extend(["-v", f"{host.resolve()}:{container}"])
    return cmd


def _run(cmd: list[str]) -> None:
    print("+ " + " ".join(cmd), flush=True)
    subprocess.run(cmd, check=True)


def _build_gatk_flags(config: Dict[str, Any]) -> list[str]:
    mapping = {
        "min_base_quality_score": "--min-base-quality-score",
        "min_mapping_quality_score": "--minimum-mapping-quality",
        "base_quality_score_threshold": "--base-quality-score-threshold",
        "standard_min_confidence_threshold_for_calling": "--standard-min-confidence-threshold-for-calling",
        "emit_ref_confidence": "--emit-ref-confidence",
        "pcr_indel_model": "--pcr-indel-model",
        "min_pruning": "--min-pruning",
        "max_alternate_alleles": "--max-alternate-alleles",
        "min_dangling_branch_length": "--min-dangling-branch-length",
        "max_num_haplotypes_in_population": "--max-num-haplotypes-in-population",
        "adaptive_pruning_initial_error_rate": "--adaptive-pruning-initial-error-rate",
        "pruning_lod_threshold": "--pruning-lod-threshold",
        "active_probability_threshold": "--active-probability-threshold",
        "min_assembly_region_size": "--min-assembly-region-size",
        "max_assembly_region_size": "--max-assembly-region-size",
        "assembly_region_padding": "--assembly-region-padding",
        "pair_hmm_gap_continuation_penalty": "--pair-hmm-gap-continuation-penalty",
        "phred_scaled_global_read_mismapping_rate": "--phred-scaled-global-read-mismapping-rate",
        "heterozygosity": "--heterozygosity",
        "indel_heterozygosity": "--indel-heterozygosity",
        "sample_ploidy": "--sample-ploidy",
        "contamination_fraction_to_filter": "--contamination-fraction-to-filter",
        "max_reads_per_alignment_start": "--max-reads-per-alignment-start",
    }
    bool_flags = {
        "recover_all_dangling_branches": "--recover-all-dangling-branches",
        "dont_use_soft_clipped_bases": "--dont-use-soft-clipped-bases",
    }
    flags: list[str] = []
    for key, flag in mapping.items():
        if key in config:
            flags.extend([flag, str(config[key])])
    for key, flag in bool_flags.items():
        if bool(config.get(key, False)):
            flags.append(flag)
    return flags


def _open_vcf(path: Path):
    return gzip.open(path, "rt", encoding="utf-8") if path.suffix == ".gz" else path.open("r", encoding="utf-8")


def _read_vcf_keys(path: Path) -> Dict[str, set[Tuple[str, int, str, str]]]:
    result = {"snp": set(), "indel": set()}
    with _open_vcf(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                continue
            chrom, pos_s, _id, ref, alts = fields[:5]
            for alt in alts.split(","):
                kind = "snp" if len(ref) == len(alt) == 1 else "indel"
                result[kind].add((chrom, int(pos_s), ref, alt))
    return result


def _f1(tp: int, fp: int, fn: int) -> tuple[float, float, float]:
    precision = tp / (tp + fp) if tp + fp else 0.0
    recall = tp / (tp + fn) if tp + fn else 0.0
    f1 = 2 * precision * recall / (precision + recall) if precision + recall else 0.0
    return precision, recall, f1


def _score(query_vcf: Path, truth_vcf: Path) -> Dict[str, Any]:
    query = _read_vcf_keys(query_vcf)
    truth = _read_vcf_keys(truth_vcf)
    metrics: Dict[str, Any] = {}
    for kind in ("snp", "indel"):
        tp = len(query[kind] & truth[kind])
        fp = len(query[kind] - truth[kind])
        fn = len(truth[kind] - query[kind])
        precision, recall, f1 = _f1(tp, fp, fn)
        metrics.update({
            f"tp_{kind}": tp,
            f"fp_{kind}": fp,
            f"fn_{kind}": fn,
            f"precision_{kind}": precision,
            f"recall_{kind}": recall,
            f"f1_{kind}": f1,
            f"truth_total_{kind}": len(truth[kind]),
            f"query_total_{kind}": len(query[kind]),
        })
    total_truth = metrics["truth_total_snp"] + metrics["truth_total_indel"]
    weighted_f1 = (
        (metrics["f1_snp"] * metrics["truth_total_snp"] + metrics["f1_indel"] * metrics["truth_total_indel"])
        / total_truth
        if total_truth
        else 0.0
    )
    metrics["simple_score"] = 100.0 * weighted_f1
    metrics["combined_final"] = weighted_f1
    return metrics


def main() -> int:
    parser = argparse.ArgumentParser(description="Run GATK on GIAB BAM and score against GIAB truth")
    parser.add_argument("--manifest", default=str(_latest_giab_manifest()), help="GIAB BAM manifest from create_giab_bam.py")
    parser.add_argument("--config", default=str(REPO_DIR / "configs" / "gatk.conf"), help="GATK config file")
    parser.add_argument("--reference", default=str(REPO_DIR / "datasets" / "reference" / "chr20" / "chr20.fa"))
    parser.add_argument("--truth-vcf-url", default=TRUTH_VCF_URL)
    parser.add_argument("--out-dir", default=str(SCORE_DIR / "giab_score"))
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--memory-gb", type=int, default=8)
    args = parser.parse_args()

    manifest_path = Path(args.manifest)
    manifest = _load_manifest(manifest_path)
    region = manifest["region"]
    bam = Path(manifest["bam"])
    reference = Path(args.reference)
    config = _read_config(Path(args.config))

    out_dir = Path(args.out_dir) / _safe_region(region)
    out_dir.mkdir(parents=True, exist_ok=True)
    truth_vcf = out_dir / "truth.vcf.gz"
    output_vcf = out_dir / "gatk.output.vcf.gz"
    rg_bam = out_dir / "input.with_rg.bam"

    # Slice GIAB truth to the same region.
    _run(_docker_base((out_dir, "/out")) + [
        BCFTOOLS_IMAGE,
        "bcftools", "view",
        "-r", region,
        args.truth_vcf_url,
        "-Oz",
        "-o", "/out/truth.vcf.gz",
    ])
    _run(_docker_base((out_dir, "/out")) + [
        BCFTOOLS_IMAGE,
        "bcftools", "index",
        "-t",
        "/out/truth.vcf.gz",
    ])

    # GATK HaplotypeCaller requires at least one sample, which is derived from
    # read-group SM tags. The public GIAB novoalign BAM has no @RG lines, so add
    # a local read group to the sliced BAM before calling.
    _run(_docker_base((bam.parent, "/bam"), (out_dir, "/out")) + [
        SAMTOOLS_IMAGE,
        "samtools", "addreplacerg",
        "-r", r"@RG\tID:HG002\tSM:HG002\tPL:ILLUMINA\tLB:GIAB\tPU:unit1",
        "-o", "/out/input.with_rg.bam",
        f"/bam/{bam.name}",
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

    metrics = _score(output_vcf, truth_vcf)
    report = {
        "region": region,
        "bam": str(bam),
        "reference": str(reference),
        "config": str(Path(args.config).resolve()),
        "truth_vcf": str(truth_vcf),
        "output_vcf": str(output_vcf),
        "metrics": metrics,
        "note": "Exact allele matching against sliced GIAB truth; this is not hap.py/vcfeval.",
    }
    report_path = out_dir / "score_report.json"
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    print("\nGIAB BAM GATK score")
    print(f"  region={region}")
    print(f"  score={metrics['simple_score']:.2f}/100")
    print(f"  SNP F1={metrics['f1_snp']:.4f} INDEL F1={metrics['f1_indel']:.4f}")
    print(f"  output={output_vcf}")
    print(f"  truth={truth_vcf}")
    print(f"  report={report_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
