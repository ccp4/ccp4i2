#!/usr/bin/env bash
# Cut the next CCP4i2 alpha.
#
# The `django` branch is protected by a ruleset: EVERY change must go through a
# pull request with all required checks green, and merge COMMITS are disabled
# (squash/rebase only). A direct `git push ... django` is rejected with GH006,
# so the release is a TWO-STEP flow, not one push:
#
#   Step 1 (--pr, the default): bump every version location in lockstep, commit
#   on a release-vX branch, push it, and open the PR against django. Stops there.
#
#   Step 2 (--tag): AFTER that PR has been reviewed, gone green, and been
#   squash-merged, tag the merged commit on django and push ONLY the tag. The
#   tag push is what fires .github/workflows/release.yml (verify -> PyPI wheel ->
#   mac/win/linux installers -> GitHub Release).
#
# The ONE source of truth is server/ccp4i2/__init__.py (MAJOR/MINOR/PATCH +
# PRERELEASE). Step 1 derives the next version and keeps the desktop app's
# exact-pin default (client/main/ccp4i2-server-version.ts) in sync; step 2 reads
# back whatever version actually landed on django and tags that, so the tag can
# never disagree with the merged bump.
#
# Full walkthrough, gotchas (artifact-storage 403s, verifying the real packaged
# build): docs/RELEASING.md.
#
# Usage:
#   scripts/cut-alpha.sh                 # step 1: bump aN->a(N+1), branch, push, PR
#   scripts/cut-alpha.sh --version 3.1.0b1   # step 1 with an explicit version
#   scripts/cut-alpha.sh --tag           # step 2 (after merge): tag django, push tag
#   scripts/cut-alpha.sh --dry-run       # show what the chosen step WOULD do
#   scripts/cut-alpha.sh --no-push       # step 1: commit on the branch locally only
#
# Requires: run from the repo root, on/for the `django` branch, with push access
# to the ccp4 remote (uses `gh auth token` if the plain remote push is unauth'd)
# and the `gh` CLI. Refuses to run with a dirty tree (except the files it edits).

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$REPO_ROOT"

INIT="server/ccp4i2/__init__.py"
CLIENT_VER="client/main/ccp4i2-server-version.ts"
LOCK="server/ccp4i2/requirements-runtime.txt"
PYPROJECT="server/pyproject.toml"
BRANCH="django"
# The upstream remote is found by URL, not by name: it has been called both
# `ccp4` and `origin`, and a rename must not break the release.
REMOTE="$(git remote -v | awk '/github.com[:\/]ccp4\/ccp4i2(\.git)? \(push\)/ {print $1; exit}')"

DRY_RUN=0
NO_PUSH=0
TAG_MODE=0
EXPLICIT_VERSION=""

while [ $# -gt 0 ]; do
  case "$1" in
    --dry-run) DRY_RUN=1 ;;
    --no-push) NO_PUSH=1 ;;
    --tag) TAG_MODE=1 ;;
    --version) EXPLICIT_VERSION="${2:?--version needs an argument}"; shift ;;
    -h|--help) sed -n '2,40p' "$0"; exit 0 ;;
    *) echo "unknown arg: $1" >&2; exit 2 ;;
  esac
  shift
done

say() { printf '\033[1;36m==>\033[0m %s\n' "$*"; }
die() { printf '\033[1;31mERROR:\033[0m %s\n' "$*" >&2; exit 1; }

# --- Preconditions --------------------------------------------------------
command -v gh >/dev/null || die "gh CLI not found"
[ -n "$REMOTE" ] || die "no git remote points at github.com/ccp4/ccp4i2; add one"
[ -f "$INIT" ] || die "run from repo root ($INIT missing)"

CUR_BRANCH="$(git rev-parse --abbrev-ref HEAD)"
[ "$CUR_BRANCH" = "$BRANCH" ] || die "on branch '$CUR_BRANCH'; releases cut from '$BRANCH'. git checkout $BRANCH"

# Modified TRACKED files (ignoring the ones we're about to edit) are a hard
# stop — they'd be swept into the release commit or reset away. Untracked files
# are left alone (they can't accidentally enter a `git add <specific files>`).
DIRTY="$(git status --porcelain --untracked-files=no | grep -vE " ($INIT|$CLIENT_VER)$" || true)"
if [ -n "$DIRTY" ]; then
  printf '\033[1;31mERROR:\033[0m modified tracked files present; commit or stash first:\n%s\n' "$DIRTY" >&2
  exit 1
fi

# --- Sync to upstream ------------------------------------------------------
say "Fetching $REMOTE/$BRANCH"
git fetch "$REMOTE" --tags -q

version_on() {  # read the MAJOR.MINOR.PATCH+PRERELEASE recorded at a git ref
  git show "$1:$INIT" | python3 -c "import re,sys; t=sys.stdin.read(); g=lambda k:re.search(rf'^{k} *= *(.+)', t, re.M).group(1).strip().strip('\"'); print(f\"{g('MAJOR')}.{g('MINOR')}.{g('PATCH')}{g('PRERELEASE')}\")"
}
push_url() {  # authenticated push URL if the plain remote push is unauth'd
  local token; token="$(gh auth token 2>/dev/null || true)"
  [ -n "$token" ] && echo "https://x-access-token:${token}@github.com/ccp4/ccp4i2.git" || echo "$REMOTE"
}

# electron-updater's GitHub provider SKIPS any release whose git tag is not
# valid semver (it calls semver.valid on the tag from the releases feed), so the
# tag must be semver — v3.1.0-alpha.58, not the PEP 440 v3.1.0a58 — or installed
# apps never see the release and report "No published versions on GitHub". The
# Python package version stays PEP 440 (3.1.0a58); only the tag is semver-ised.
# release.yml's verify-version converts the tag back to PEP to match __version__.
pep_to_semver() {
  printf '%s' "$1" | sed -E 's/([0-9])a([0-9]+)/\1-alpha.\2/; s/([0-9])b([0-9]+)/\1-beta.\2/; s/([0-9])rc([0-9]+)/\1-rc.\2/'
}

# --- Step 2 (--tag): tag the merged bump on django, push the tag -----------
if [ "$TAG_MODE" = 1 ]; then
  MERGED_VER="$(version_on "$REMOTE/$BRANCH")"
  TAG="v$(pep_to_semver "$MERGED_VER")"
  say "Version on $REMOTE/$BRANCH: $MERGED_VER   ->   tag $TAG"
  # The bump must already be merged. If the tip of django is not a release
  # commit for this version, the PR from step 1 has not landed yet.
  git log -1 --format='%s' "$REMOTE/$BRANCH" | grep -q "release: ccp4i2 $MERGED_VER" \
    || die "the tip of $REMOTE/$BRANCH is not 'release: ccp4i2 $MERGED_VER' — has the release PR merged? (step 1 opens it; merge it first)"
  if git ls-remote --tags "$REMOTE" "refs/tags/$TAG" | grep -q "$TAG"; then
    die "tag $TAG already exists on $REMOTE — this release was already tagged."
  fi
  if curl -fsS "https://pypi.org/pypi/ccp4i2/${MERGED_VER}/json" >/dev/null 2>&1; then
    die "ccp4i2 $MERGED_VER is ALREADY on PyPI (immutable) — already released."
  fi
  if [ "$DRY_RUN" = 1 ]; then
    say "DRY RUN — would tag $REMOTE/$BRANCH ($(git rev-parse --short "$REMOTE/$BRANCH")) as $TAG and push it (fires release.yml)."
    exit 0
  fi
  git tag -a "$TAG" "$REMOTE/$BRANCH" -m "CCP4i2 $MERGED_VER"
  say "Pushing $TAG to $REMOTE (this triggers the release)"
  git push "$(push_url)" "$TAG"
  say "Released: $TAG pushed. Watch: gh run watch --repo ccp4/ccp4i2 \$(gh run list --repo ccp4/ccp4i2 --workflow release.yml --branch $TAG --limit 1 --json databaseId -q '.[0].databaseId')"
  exit 0
fi

# --- Step 1 continues: sync local django so the release branch forks from it
BEHIND="$(git rev-list --count HEAD.."$REMOTE/$BRANCH" 2>/dev/null || echo 0)"
if [ "$BEHIND" -gt 0 ]; then
  [ "$DRY_RUN" = 1 ] && die "local $BRANCH is $BEHIND behind $REMOTE/$BRANCH; sync first (dry-run won't touch git)"
  say "Local $BRANCH is $BEHIND behind; hard-resetting to $REMOTE/$BRANCH"
  git reset --hard "$REMOTE/$BRANCH" -q
fi

# --- Compute the next version ---------------------------------------------
read_field() { grep -E "^$1 *=" "$INIT" | head -1 | sed -E 's/.*= *//; s/ *#.*//; s/"//g'; }
MAJOR="$(read_field MAJOR)"; MINOR="$(read_field MINOR)"; PATCH="$(read_field PATCH)"
CUR_PRE="$(read_field PRERELEASE)"
CUR_VER="${MAJOR}.${MINOR}.${PATCH}${CUR_PRE}"

if [ -n "$EXPLICIT_VERSION" ]; then
  NEW_VER="$EXPLICIT_VERSION"
  # derive the new PRERELEASE suffix by stripping the numeric prefix
  NEW_PRE="${NEW_VER#${MAJOR}.${MINOR}.${PATCH}}"
  [ "${MAJOR}.${MINOR}.${PATCH}${NEW_PRE}" = "$NEW_VER" ] || \
    die "--version $NEW_VER doesn't share the ${MAJOR}.${MINOR}.${PATCH} base; edit $INIT by hand for a MAJOR/MINOR/PATCH change"
else
  # bump the trailing integer of an alpha/beta/rc suffix: a6 -> a7, b2 -> b3
  if [[ "$CUR_PRE" =~ ^([abr]|rc)([0-9]+)$ ]]; then
    NEW_PRE="${BASH_REMATCH[1]}$(( BASH_REMATCH[2] + 1 ))"
  else
    die "current PRERELEASE '$CUR_PRE' isn't an alpha/beta/rc suffix; use --version for a non-incremental bump"
  fi
  NEW_VER="${MAJOR}.${MINOR}.${PATCH}${NEW_PRE}"
fi

TAG="v$(pep_to_semver "$NEW_VER")"
say "Current: $CUR_VER   ->   New: $NEW_VER   (tag $TAG)"

# Refuse a version already on PyPI (immutable — re-cut would fail publish-pypi).
if curl -fsS "https://pypi.org/pypi/ccp4i2/${NEW_VER}/json" >/dev/null 2>&1; then
  die "ccp4i2 $NEW_VER is ALREADY on PyPI (immutable). Bump further."
fi
# Refuse an existing tag.
if git ls-remote --tags "$REMOTE" "refs/tags/$TAG" | grep -q "$TAG"; then
  die "tag $TAG already exists on $REMOTE. Bump further."
fi

# --- Consistency: lock vs wheel floor for ccp4i2-api ----------------------
FLOOR="$(grep -E 'ccp4i2-api *>=' "$PYPROJECT" | head -1 | sed -E 's/.*>= *([0-9][^",;[:space:]]*).*/\1/')"
LOCKPIN="$(grep -E '^ccp4i2-api==' "$LOCK" | head -1 | sed -E 's/^ccp4i2-api==//; s/[[:space:]#].*//')"
say "ccp4i2-api: wheel floor >=$FLOOR   lock pin ==$LOCKPIN"
verlte() { python3 -c "import re,sys; f=lambda v:[int((re.match(r'\d*',p).group() or 0)) for p in v.split('.')]; sys.exit(0 if f('$1')>=f('$2') else 1)"; }
verlte "$LOCKPIN" "$FLOOR" || die "runtime lock ccp4i2-api==$LOCKPIN < wheel floor >=$FLOOR. Update $LOCK (and re-test) before releasing — the app installs the LOCK version."

# --- Apply the version bumps (step 1) -------------------------------------
TODAY="$(date +'%Y, %-m, %-d')"
RELEASE_BRANCH="release-$TAG"
if [ "$DRY_RUN" = 1 ]; then
  say "DRY RUN — step 1 would:"
  echo "    bump $INIT        PRERELEASE = \"$NEW_PRE\"  (+ date datetime($TODAY))"
  echo "    bump $CLIENT_VER  default pin -> \"$NEW_VER\""
  echo "    commit 'release: ccp4i2 $NEW_VER' on branch $RELEASE_BRANCH"
  echo "    push $RELEASE_BRANCH and open a PR into $BRANCH"
  echo "    (then, after the PR is merged: scripts/cut-alpha.sh --tag)"
  exit 0
fi

say "Bumping $INIT"
sed -i.bak -E "s/^PRERELEASE = \".*\"/PRERELEASE = \"$NEW_PRE\"/" "$INIT"
sed -i.bak -E "s/^__version_date__ = datetime\(.*\)/__version_date__ = datetime($TODAY)/" "$INIT"
rm -f "$INIT.bak"

say "Bumping $CLIENT_VER exact-pin default"
sed -i.bak -E "s/(CCP4I2_SERVER_VERSION_FLOOR \|\| )\"[^\"]*\"/\1\"$NEW_VER\"/" "$CLIENT_VER"
rm -f "$CLIENT_VER.bak"

# Sanity: the file now yields exactly NEW_VER
GOT="$(python3 -c "import re; t=open('$INIT').read(); g=lambda k:re.search(rf'^{k} *= *(.+)', t, re.M).group(1).strip().strip('\"'); print(f\"{g('MAJOR')}.{g('MINOR')}.{g('PATCH')}{g('PRERELEASE')}\")")"
[ "$GOT" = "$NEW_VER" ] || die "post-edit version is '$GOT', expected '$NEW_VER' — check $INIT"
grep -q "|| \"$NEW_VER\"" "$CLIENT_VER" || die "client pin didn't update to $NEW_VER"

# --- Commit on a release branch, push, open the PR ------------------------
# django is PR-only (a direct push is rejected with GH006), so the bump lands
# via a PR; the tag is pushed separately by --tag once it has merged.
git checkout -b "$RELEASE_BRANCH" 2>/dev/null || git checkout "$RELEASE_BRANCH"
git add "$INIT" "$CLIENT_VER"
git commit -q -m "release: ccp4i2 $NEW_VER

Automated alpha cut via scripts/cut-alpha.sh (bump PRERELEASE + exact-pin
default in lockstep). After this PR merges, run scripts/cut-alpha.sh --tag to
tag the merged commit and fire the release workflow.

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"

if [ "$NO_PUSH" = 1 ]; then
  say "Committed on $RELEASE_BRANCH locally. --no-push: NOT pushing / no PR. To continue:"
  echo "    git push $REMOTE $RELEASE_BRANCH && gh pr create --base $BRANCH --head $RELEASE_BRANCH"
  exit 0
fi

say "Pushing $RELEASE_BRANCH to $REMOTE"
git push -u "$(push_url)" "$RELEASE_BRANCH"

say "Opening the release PR into $BRANCH"
gh pr create --repo ccp4/ccp4i2 --base "$BRANCH" --head "$RELEASE_BRANCH" \
  --title "release: ccp4i2 $NEW_VER" \
  --body "Version bump to \`$NEW_VER\` — PRERELEASE + the desktop exact-pin default, in lockstep. Merge this (squash), then run \`scripts/cut-alpha.sh --tag\`: it tags the merged commit on \`$BRANCH\` and pushes the tag, firing the release workflow (PyPI wheel + mac/win/linux installers + GitHub Release)." \
  || die "gh pr create failed — the branch pushed, so open the PR by hand (base $BRANCH, head $RELEASE_BRANCH)."

say "Step 1 done. Next: review + wait for all checks green + squash-merge the PR, then run:"
echo "    scripts/cut-alpha.sh --tag"
