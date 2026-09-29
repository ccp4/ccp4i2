#!/usr/bin/env bash
# Which CI validation does a changeset need?
#
#   changes-need.sh <base-sha> <head-sha> [merge-base|direct]
#
# Prints three lines, suitable for $GITHUB_OUTPUT:
#
#   backend=true|false     the Python unit + API suites
#   frontend=true|false    the client's vitest + tsc
#   build=true|false       the mac/win/linux desktop packages
#
# Replaces changes-need-build.sh, which answered only the build question, and
# answered it for every non-documentation change -- so a server-only pull
# request packaged three desktop apps. Those builds run `npm ci` and
# `npm run package-*` in client/ and packages/ and read nothing under server/,
# so they cannot tell you anything about a backend change. On one measured
# server-only PR they were 87% of the billable minutes and none of the signal
# (macOS bills at 10x, Windows at 2x).
#
# "merge-base" (the default) diffs head against its merge base with base, which
# is what a pull request changes; "direct" diffs the two commits, which is what
# a push to a branch changes.
#
# When in doubt, run everything: an unknown base (the first push of a branch), a
# failing git, or a path no rule below claims.
#
# NOTE on the version lockstep. The desktop app pins the EXACT backend version
# it expects (client/main/ccp4i2-server-version.ts) and that must equal
# ccp4i2.__version__. It is tempting to make a change to server/ccp4i2/__init__.py
# force a desktop build for that reason. It is not needed: a version bump always
# touches the client pin too (scripts/cut-alpha.sh edits both files and greps to
# confirm), so the client rule already claims it. The invariant itself is checked
# directly and cheaply by check-version-lockstep.sh, which is a better home for
# it than three package builds that never read the pin.
set -u

base="${1:-}"; head="${2:-}"; mode="${3:-merge-base}"

everything() {
  echo "backend=true"; echo "frontend=true"; echo "build=true"; exit 0
}

nothing() {
  echo "backend=false"; echo "frontend=false"; echo "build=false"; exit 0
}

if [ -z "$base" ] || [ -z "$head" ] \
   || ! git cat-file -e "${base}^{commit}" 2>/dev/null \
   || ! git cat-file -e "${head}^{commit}" 2>/dev/null; then
  everything
fi

if [ "$mode" = "direct" ]; then range="${base}..${head}"; else range="${base}...${head}"; fi
if ! changed="$(git diff --name-only "$range" 2>/dev/null)"; then
  everything
fi

# Nothing changed at all.
[ -n "$changed" ] || nothing

# Documentation only: a *.md anywhere, docs/, LICENSE, or a test baseline under
# server/.test-baselines/ (evidence about a run, never an input to one).
if ! printf '%s\n' "$changed" \
     | grep -qvE '(\.md$|^docs/|^LICENSE$|^server/\.test-baselines/)'; then
  nothing
fi

# Any CI change validates itself with the lot: these files decide what runs, so
# a mistake in them is invisible to a narrowed run.
if printf '%s\n' "$changed" | grep -qE '^\.github/'; then
  everything
fi

backend=false
frontend=false

# packages/ is the shared API contract: consumed by the server (PyPI
# ccp4i2-api) and npm-installed by the client build, so it claims both.
if printf '%s\n' "$changed" | grep -qE '^(server/|packages/)'; then
  backend=true
fi
if printf '%s\n' "$changed" \
     | grep -qE '^(client/|packages/|package(-lock)?\.json$)'; then
  frontend=true
fi

# Anything not claimed by a rule above is unfamiliar; run everything rather
# than guess that it does not matter.
if printf '%s\n' "$changed" \
     | grep -qvE '(\.md$|^docs/|^LICENSE$|^server/|^packages/|^client/|^package(-lock)?\.json$)'; then
  everything
fi

echo "backend=${backend}"
echo "frontend=${frontend}"
# The desktop package is built from client/ and packages/ alone.
echo "build=${frontend}"
