#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")"

REGION="${REGION:-chr20:1000000-1500000}"
MUTATIONS="${MUTATIONS:-1000}"
ALLELE_FRACTION="${ALLELE_FRACTION:-1.0}"
GATK_CONFIG="${GATK_CONFIG:-gatk.conf}"

./create_giab_bam.sh --region "$REGION"
mutate_args=(--mutations "$MUTATIONS" --allele-fraction "$ALLELE_FRACTION")
if [[ -n "${SEED:-}" ]]; then
  mutate_args+=(--seed "$SEED")
fi
./mutate_giab_bam.sh "${mutate_args[@]}"
./score_mutated_giab_bam.sh --config "$GATK_CONFIG"
