#!/usr/bin/env bash
# Build the user help, and fail if it has more warnings than it had.
#
#   docs/user/tools/build.sh [output dir]     (default: docs/user/_build)
#
# The imported Qt-era pages carry warnings of their own (orphan pages, stale
# references), so -W would fail every build. Instead the count is a ratchet:
# warnings-baseline.txt holds the number the tree has now; more fails, and when
# a change removes some, lower the baseline in the same change.
set -euo pipefail
here="$(cd "$(dirname "$0")" && pwd)"
docs="$(dirname "$here")"
out="${1:-$docs/_build}"
log="$(mktemp)"
sphinx-build -q -b html "$docs/source" "$out" 2>"$log" || { cat "$log"; exit 1; }
count="$(grep -cE 'WARNING|ERROR' "$log" || true)"
baseline="$(cat "$here/warnings-baseline.txt")"
echo "Built $out: $count warnings (baseline $baseline)"
if [ "$count" -gt "$baseline" ]; then
  echo "More warnings than the baseline. New ones are among:"
  grep -E 'WARNING|ERROR' "$log"
  exit 1
fi
if [ "$count" -lt "$baseline" ]; then
  echo "Fewer warnings than the baseline: lower tools/warnings-baseline.txt to $count."
fi
