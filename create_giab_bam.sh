#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")"
python3 create_giab_bam.py "$@"
