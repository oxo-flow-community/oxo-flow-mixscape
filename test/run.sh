#!/usr/bin/env bash
# Acceptance test for the oxo-flow-mixscape port.
# Usage: ./test/run.sh            (uses ./main.oxoflow)
set -euo pipefail
cd "$(dirname "$0")/.."
OXO=${OXO:-oxo-flow}

echo "==> validate"
"$OXO" validate main.oxoflow

echo "==> lint (warnings are acceptable, errors are not)"
"$OXO" lint main.oxoflow

echo "==> dry-run with default config"
# Engine headline wording changed in oxo-flow 0.20.1 (Traitome/oxo-flow#432):
# the old "DAG: (dry-run) N rules would execute" became
# "Plan: would run: N | skip: M | completed: K". Accept either so CI passes
# with both the latest release and older pins.
"$OXO" dry-run main.oxoflow --samples first:1 > /tmp/oxo-dryrun-$$.txt 2>&1
grep -qE "would (execute|run)" /tmp/oxo-dryrun-$$.txt

echo "==> debug: expanded commands contain no literal {wildcards}"
"$OXO" debug main.oxoflow 2>&1 | grep -q '{sample}' && { echo "unexpanded wildcards in debug output"; exit 1; } || true

echo "PASS"
