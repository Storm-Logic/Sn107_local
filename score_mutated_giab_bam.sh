#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")"
python3 score_mutated_giab_bam.py "$@"
