#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
find NNSProofs -name '*.lean' -exec awk '/^(public )?import Mathlib\./ {print $NF}' {} + |
  sort -u | xargs lake exe cache get
