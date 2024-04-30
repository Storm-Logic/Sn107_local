#!/usr/bin/env python3
"""Inject local fake variants into a GIAB BAM slice.

The host entrypoint runs this script inside Docker with pysam installed. The
worker edits reads at selected positions, writes a mutated BAM, indexes it,
and writes mutation/truth VCFs. This is a lightweight local analogue of the
SN107 read-level mutation step; it is not HelixForge.
"""

from __future__ import annotations

import argparse
import json
import os
import random
import re
import secrets
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Dict, List, Set, Tuple


SCORE_DIR = Path(__file__).resolve().parent
IMAGE_NAME = "local-sn107-pysam-mutator:latest"
BASES = "ACGT"
EVENT_TYPES = ("snp", "del", "ins")
CONSUMES_QUERY = {0, 1, 4, 7, 8}
CONSUMES_REF = {0, 2, 3, 7, 8}


def _latest_manifest() -> Path:
    manifests = sorted((SCORE_DIR / "giab_bam").glob("*.manifest.json"), key=lambda p: p.stat().st_mtime)
    if not manifests:
        raise SystemExit("No GIAB BAM manifest found. Run ./create_giab_bam.sh first.")
    return manifests[-1]


def _safe(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value).strip("_")


def _to_container_path(path: Path, workspace: Path) -> str:
    resolved = path.resolve()
    return "/workspace/" + str(resolved.relative_to(workspace.resolve()))


def _run_host(args: argparse.Namespace) -> int:
    workspace = SCORE_DIR.parent
    dockerfile = SCORE_DIR / "Dockerfile.mutate"
    subprocess.run(["docker", "build", "-t", IMAGE_NAME, "-f", str(dockerfile), str(SCORE_DIR)], check=True)

    manifest = Path(args.manifest).resolve()
    out_dir = Path(args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    seed = args.seed if args.seed is not None else secrets.randbits(32)

    manifest_data = json.loads(manifest.read_text(encoding="utf-8"))
    for key in ("bam", "bai"):
        if key in manifest_data:
            manifest_data[key] = _to_container_path(Path(manifest_data[key]), workspace)
    container_manifest_host = out_dir / ".container_input_manifest.json"
    container_manifest_host.write_text(json.dumps(manifest_data, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    cmd = [
        "docker", "run", "--rm",
        "-v", f"{workspace}:/workspace",
        "-w", "/workspace/Score",
        IMAGE_NAME,
        "python3", "mutate_giab_bam.py",
        "--worker",
        "--manifest", _to_container_path(container_manifest_host, workspace),
        "--out-dir", _to_container_path(out_dir, workspace),
        "--reference", _to_container_path(Path(args.reference), workspace),
        "--mutations", str(args.mutations),
        "--allele-fraction", str(args.allele_fraction),
        "--seed", str(seed),
        "--min-coverage", str(args.min_coverage),
        "--variant-types", args.variant_types,
    ]
    return subprocess.call(cmd)


def _parse_region(region: str) -> Tuple[str, int, int]:
    chrom, coords = region.split(":", 1)
    start_s, end_s = coords.split("-", 1)
    return chrom, int(start_s), int(end_s)


def _choose_alt(ref: str, rng: random.Random) -> str:
    return rng.choice([b for b in BASES if b != ref])


def _add_cigar(cigar: List[Tuple[int, int]], op: int, length: int) -> None:
    if length <= 0:
        return
    if cigar and cigar[-1][0] == op:
        cigar[-1] = (op, cigar[-1][1] + length)
    else:
        cigar.append((op, length))


def _parse_variant_types(value: str) -> List[str]:
    if value.lower() in {"all", "*"}:
        return list(EVENT_TYPES)
    parsed = [item.strip().lower() for item in value.split(",") if item.strip()]
    invalid = [item for item in parsed if item not in EVENT_TYPES]
    if invalid:
        raise SystemExit(f"Invalid --variant-types values: {', '.join(invalid)}. Use snp,del,ins or all.")
    if not parsed:
        raise SystemExit("--variant-types must include at least one of snp,del,ins")
    return parsed


def _write_vcf(path: Path, chrom: str, variants: List[Dict[str, Any]]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        handle.write("##fileformat=VCFv4.2\n")
        handle.write("##source=local-read-level-giab-mutator\n")
        handle.write("##INFO=<ID=SYNTHETIC,Number=0,Type=Flag,Description=\"Local synthetic injected mutation\">\n")
        handle.write("##INFO=<ID=TYPE,Number=1,Type=String,Description=\"Local synthetic variant type\">\n")
        handle.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
        handle.write(f"##contig=<ID={chrom}>\n")
        handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002\n")
        for variant in sorted(variants, key=lambda item: (item["pos"], item["ref"], item["alt"])):
            handle.write(
                f"{chrom}\t{variant['pos']}\t.\t{variant['ref']}\t{variant['alt']}"
                f"\t60\tPASS\tSYNTHETIC;TYPE={variant['kind']}\tGT\t1/1\n"
            )


def _vcf_variant_for_event(fasta, chrom: str, pos: int, event: Dict[str, Any]) -> Dict[str, Any]:
    kind = event["kind"]
    if kind == "snp":
        return {"kind": "snp", "pos": pos, "ref": event["ref"], "alt": event["alt"]}
    if kind == "ins":
        ref = event["ref"]
        return {"kind": "ins", "pos": pos, "ref": ref, "alt": ref + event["inserted"]}
    if kind == "del":
        anchor_pos = pos - 1
        anchor = fasta.fetch(chrom, anchor_pos - 1, anchor_pos).upper()
        return {"kind": "del", "pos": anchor_pos, "ref": anchor + event["deleted"], "alt": anchor}
    raise ValueError(f"Unknown event kind: {kind}")


def _apply_events_to_read(read, events_by_pos: Dict[int, Dict[str, Any]]) -> bool:
    original_seq = read.query_sequence
    if not original_seq or not read.cigartuples:
        return False
    original_quals = list(read.query_qualities or [30] * len(original_seq))
    new_seq: List[str] = []
    new_quals: List[int] = []
    new_cigar: List[Tuple[int, int]] = []
    qpos = 0
    rpos = read.reference_start
    changed = False

    for op, length in read.cigartuples:
        if op in {0, 7, 8}:
            for _ in range(length):
                pos1 = rpos + 1
                event = events_by_pos.get(pos1)
                base = original_seq[qpos]
                qual = original_quals[qpos]
                if event and event["kind"] == "del":
                    _add_cigar(new_cigar, 2, 1)
                    changed = True
                else:
                    if event and event["kind"] == "snp":
                        base = event["alt"]
                        changed = True
                    new_seq.append(base)
                    new_quals.append(qual)
                    _add_cigar(new_cigar, op, 1)
                    if event and event["kind"] == "ins":
                        inserted = event["inserted"]
                        new_seq.extend(inserted)
                        new_quals.extend([qual] * len(inserted))
                        _add_cigar(new_cigar, 1, len(inserted))
                        changed = True
                qpos += 1
                rpos += 1
        elif op in {1, 4}:
            new_seq.extend(original_seq[qpos:qpos + length])
            new_quals.extend(original_quals[qpos:qpos + length])
            _add_cigar(new_cigar, op, length)
            qpos += length
        elif op in {2, 3}:
            _add_cigar(new_cigar, op, length)
            rpos += length
        elif op in {5, 6}:
            _add_cigar(new_cigar, op, length)
        else:
            _add_cigar(new_cigar, op, length)

    if not changed:
        return False

    read.query_sequence = "".join(new_seq)
    read.query_qualities = new_quals
    read.cigartuples = new_cigar
    for tag in ("MD", "NM", "UQ"):
        if read.has_tag(tag):
            read.set_tag(tag, None)
    return True


def _run_worker(args: argparse.Namespace) -> int:
    import pysam

    manifest = json.loads(Path(args.manifest).read_text(encoding="utf-8"))
    input_bam = Path(manifest["bam"])
    reference = Path(args.reference)
    chrom, start, end = _parse_region(manifest["region"])
    rng = random.Random(args.seed)
    variant_types = _parse_variant_types(args.variant_types)

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    name = f"{Path(input_bam).stem}.mutated.seed{args.seed}"
    output_bam = out_dir / f"{name}.bam"

    fasta = pysam.FastaFile(str(reference))

    # First pass: collect candidate positions and read names covering them.
    coverage: Dict[int, Set[str]] = defaultdict(set)
    ref_support: Dict[int, Set[str]] = defaultdict(set)
    base_counts: Dict[int, Counter] = defaultdict(Counter)
    with pysam.AlignmentFile(input_bam, "rb") as bam:
        for read in bam.fetch(chrom, start - 1, end):
            if (
                read.is_unmapped
                or read.is_duplicate
                or read.is_secondary
                or read.is_supplementary
                or read.query_sequence is None
            ):
                continue
            for qpos, rpos in read.get_aligned_pairs(matches_only=True):
                if qpos is None or rpos is None:
                    continue
                pos1 = rpos + 1
                if start <= pos1 <= end:
                    base = read.query_sequence[qpos].upper()
                    if base in BASES:
                        coverage[pos1].add(read.query_name)
                        base_counts[pos1][base] += 1
                        ref_base = fasta.fetch(chrom, pos1 - 1, pos1).upper()
                        if base == ref_base:
                            ref_support[pos1].add(read.query_name)

    candidates = [
        pos for pos, names in coverage.items()
        if (
            start + 1 <= pos <= end - 1
            and
            len(names) >= args.min_coverage
            and len(ref_support[pos]) >= args.min_coverage
            and fasta.fetch(chrom, pos - 1, pos).upper() in BASES
            and fasta.fetch(chrom, pos - 2, pos - 1).upper() in BASES
            and not any(base_counts[pos][base] for base in BASES if base != fasta.fetch(chrom, pos - 1, pos).upper())
        )
    ]
    if len(candidates) < args.mutations:
        raise SystemExit(
            f"Only {len(candidates)} candidate positions with coverage >= {args.min_coverage}; "
            f"requested {args.mutations}"
        )

    rng.shuffle(candidates)
    selected = []
    for pos in candidates:
        if all(abs(pos - existing) >= 10 for existing in selected):
            selected.append(pos)
            if len(selected) == args.mutations:
                break
    if len(selected) < args.mutations:
        raise SystemExit(
            f"Only {len(selected)} spaced candidate positions with coverage >= {args.min_coverage}; "
            f"requested {args.mutations}"
        )
    selected = sorted(selected)

    mutation_by_pos: Dict[int, Dict[str, Any]] = {}
    for idx, pos in enumerate(selected):
        kind = variant_types[idx % len(variant_types)]
        ref = fasta.fetch(chrom, pos - 1, pos).upper()
        if kind == "snp":
            mutation_by_pos[pos] = {"kind": "snp", "ref": ref, "alt": _choose_alt(ref, rng)}
        elif kind == "ins":
            mutation_by_pos[pos] = {"kind": "ins", "ref": ref, "inserted": rng.choice(BASES)}
        elif kind == "del":
            mutation_by_pos[pos] = {"kind": "del", "ref": ref, "deleted": ref}

    # Select read names to edit per mutation.
    mutate_read_at_pos: Dict[Tuple[str, int], Dict[str, Any]] = {}
    for pos, event in mutation_by_pos.items():
        names = sorted(ref_support[pos])
        n_to_mutate = max(1, int(round(len(names) * args.allele_fraction)))
        for read_name in rng.sample(names, min(n_to_mutate, len(names))):
            mutate_read_at_pos[(read_name, pos)] = event

    with pysam.AlignmentFile(input_bam, "rb") as src, pysam.AlignmentFile(output_bam, "wb", template=src) as dst:
        for read in src:
            if read.is_unmapped or read.query_sequence is None:
                dst.write(read)
                continue
            events_for_read: Dict[int, Dict[str, Any]] = {}
            for qpos, rpos in read.get_aligned_pairs(matches_only=True):
                if qpos is None or rpos is None:
                    continue
                pos1 = rpos + 1
                event = mutate_read_at_pos.get((read.query_name, pos1))
                if event:
                    events_for_read[pos1] = event
            if events_for_read:
                _apply_events_to_read(read, events_for_read)
            dst.write(read)

    pysam.index(str(output_bam))
    variants = [
        _vcf_variant_for_event(fasta, chrom, pos, event)
        for pos, event in mutation_by_pos.items()
    ]
    variant_counts = Counter(variant["kind"] for variant in variants)
    mutations_vcf = out_dir / f"{name}.mutations.vcf"
    truth_vcf = out_dir / f"{name}.truth.vcf"
    _write_vcf(mutations_vcf, chrom, variants)
    _write_vcf(truth_vcf, chrom, variants)

    mutated_manifest = {
        "source_manifest": str(Path(args.manifest)),
        "source_bam": str(input_bam),
        "reference": str(reference),
        "region": manifest["region"],
        "mutated_bam": str(output_bam),
        "mutated_bai": str(output_bam) + ".bai",
        "mutations_vcf": str(mutations_vcf),
        "truth_vcf": str(truth_vcf),
        "mutation_count": len(variants),
        "variant_counts": dict(variant_counts),
        "variant_types": variant_types,
        "allele_fraction": args.allele_fraction,
        "seed": args.seed,
        "min_coverage": args.min_coverage,
        "note": "Local read-level SNP/INDEL mutation injection for GIAB BAM; not HelixForge.",
    }
    manifest_path = out_dir / f"{name}.manifest.json"
    manifest_path.write_text(json.dumps(mutated_manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    print("\nCreated mutated GIAB BAM")
    print(f"  BAM: {output_bam}")
    print(f"  BAI: {output_bam}.bai")
    print(f"  mutations: {mutations_vcf}")
    print(f"  truth: {truth_vcf}")
    print(f"  manifest: {manifest_path}")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description="Inject fake SNP/INDEL mutations into a GIAB BAM slice")
    parser.add_argument("--worker", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("--manifest", default=str(_latest_manifest()), help="Input GIAB BAM manifest")
    parser.add_argument("--out-dir", default=str(SCORE_DIR / "giab_mutated"), help="Output directory")
    parser.add_argument("--reference", default=str(SCORE_DIR.parent / "mine_107" / "datasets" / "reference" / "chr20" / "chr20.fa"))
    parser.add_argument("--mutations", type=int, default=25)
    parser.add_argument("--allele-fraction", type=float, default=1.0)
    parser.add_argument("--seed", type=int, default=None, help="Random seed. Omit for a new random injection each run.")
    parser.add_argument("--min-coverage", type=int, default=12)
    parser.add_argument("--variant-types", default="all", help="Comma-separated variant types: snp,del,ins, or all")
    args = parser.parse_args()
    if args.worker:
        return _run_worker(args)
    return _run_host(args)


if __name__ == "__main__":
    raise SystemExit(main())
