#!/usr/bin/env bash
# Does the changeset between two commits need the desktop build?
#
#   changes-need-build.sh <base-sha> <head-sha> [merge-base|direct]
#
# Prints "true" or "false". A changeset needs the build unless EVERY changed
# file is documentation: any *.md, anything under docs/, LICENSE, or a test
# baseline under server/.test-baselines/ (per-test JUnit results and their
# summaries: evidence about a build, never an input to one). A mixed
# changeset builds. So does anything that cannot be decided (an unknown base,
# as on the first push of a branch, or a failing git): when in doubt, build.
#
# "merge-base" (the default) diffs head against its merge base with base,
# which is what a pull request changes; "direct" diffs the two commits, which
# is what a push to a branch changes.
set -u
base="${1:-}"; head="${2:-}"; mode="${3:-merge-base}"

if [ -z "$base" ] || [ -z "$head" ] \
   || ! git cat-file -e "${base}^{commit}" 2>/dev/null \
   || ! git cat-file -e "${head}^{commit}" 2>/dev/null; then
  echo "true"; exit 0
fi

if [ "$mode" = "direct" ]; then range="${base}..${head}"; else range="${base}...${head}"; fi
if ! changed="$(git diff --name-only "$range" 2>/dev/null)"; then
  echo "true"; exit 0
fi

# No changed files at all: nothing to build.
if [ -z "$changed" ]; then echo "false"; exit 0; fi

if printf '%s\n' "$changed" | grep -qvE '(\.md$|^docs/|^LICENSE$|^server/\.test-baselines/)'; then
  echo "true"
else
  echo "false"
fi
