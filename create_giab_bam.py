#!/usr/bin/env python3
"""Create a local region BAM from public GIAB HG002 data using Docker samtools."""

from __future__ import annotations

import argparse
import json
import platform
import re
import subprocess
import sys
from pathlib import Path


SCORE_DIR = Path(__file__).resolve().parent
DEFAULT_BAM_URL = (
    "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/"
    "AshkenazimTrio/HG002_NA24385_son/NIST_Illumina_2x250bps/"
    "novoalign_bams/HG002.GRCh38.2x250.bam"
)
SAMTOOLS_IMAGE = "quay.io/biocontainers/samtools:1.20--h50ea8bc_0"


def _safe_name(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value).strip("_")


def _default_region() -> str:
    config_path = SCORE_DIR / "miner_config.json"
    if config_path.exists():
        try:
            data = json.loads(config_path.read_text(encoding="utf-8"))
            return data.get("challenge", {}).get("region") or "chr20:1000000-1010000"
        except Exception:
            pass
    return "chr20:1000000-1010000"


def _docker_base(output_dir: Path) -> list[str]:
    cmd = ["docker", "run", "--rm"]
    if platform.machine() == "arm64":
        cmd.extend(["--platform", "linux/amd64"])
    cmd.extend([
        "-v", f"{output_dir.resolve()}:/data",
        SAMTOOLS_IMAGE,
    ])
    return cmd


def _run(cmd: list[str]) -> None:
    print("+ " + " ".join(cmd), flush=True)
    subprocess.run(cmd, check=True)


def main() -> int:
    parser = argparse.ArgumentParser(description="Create a region BAM from public GIAB HG002 GRCh38 data")
    parser.add_argument("--region", default=_default_region(), help="Region to slice, for example chr20:1000000-1010000")
    parser.add_argument("--bam-url", default=DEFAULT_BAM_URL, help="Remote coordinate-sorted BAM URL")
    parser.add_argument("--out-dir", default=str(SCORE_DIR / "giab_bam"), help="Output directory")
    parser.add_argument("--name", help="Output basename without .bam")
    args = parser.parse_args()

    output_dir = Path(args.out_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    name = args.name or f"HG002_GRCh38_{_safe_name(args.region)}"
    output_bam = output_dir / f"{name}.bam"

    view_cmd = _docker_base(output_dir) + [
        "samtools", "view",
        "-b",
        "-h",
        "-o", f"/data/{output_bam.name}",
        args.bam_url,
        args.region,
    ]
    index_cmd = _docker_base(output_dir) + [
        "samtools", "index",
        f"/data/{output_bam.name}",
    ]
    quickcheck_cmd = _docker_base(output_dir) + [
        "samtools", "quickcheck",
        "-v",
        f"/data/{output_bam.name}",
    ]
    count_cmd = _docker_base(output_dir) + [
        "samtools", "view",
        "-c",
        f"/data/{output_bam.name}",
    ]

    try:
        _run(view_cmd)
        _run(index_cmd)
        _run(quickcheck_cmd)
        count = subprocess.check_output(count_cmd, text=True).strip()
    except subprocess.CalledProcessError as exc:
        print(f"GIAB BAM creation failed: {exc}", file=sys.stderr)
        return exc.returncode or 1

    manifest = {
        "source": "GIAB HG002 GRCh38 Illumina 2x250 novoalign",
        "source_bam_url": args.bam_url,
        "region": args.region,
        "bam": str(output_bam),
        "bai": str(output_bam) + ".bai",
        "read_count": int(count),
        "samtools_image": SAMTOOLS_IMAGE,
    }
    manifest_path = output_dir / f"{name}.manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    print("\nCreated GIAB region BAM")
    print(f"  BAM: {output_bam}")
    print(f"  BAI: {output_bam}.bai")
    print(f"  reads: {count}")
    print(f"  manifest: {manifest_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
