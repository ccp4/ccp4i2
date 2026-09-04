#!/usr/bin/env bash
# Cut the next CCP4i2 alpha in one command.
#
# Encodes the whole pre-flight that was previously hand-driven: sync the branch,
# bump every version location in lockstep, verify they agree, commit, tag, push.
# Pushing the tag triggers .github/workflows/release.yml (verify -> PyPI wheel ->
# mac/win/linux installers -> GitHub Release).
#
# The ONE source of truth is server/ccp4i2/__init__.py (MAJOR/MINOR/PATCH +
# PRERELEASE). This script derives the next version and keeps the desktop app's
# exact-pin default (client/main/ccp4i2-server-version.ts) in sync.
#
# Usage:
#   scripts/cut-alpha.sh              # bump PRERELEASE aN -> a(N+1), tag, push
#   scripts/cut-alpha.sh --version 3.1.0b1   # set an explicit version
#   scripts/cut-alpha.sh --dry-run    # show what it WOULD do, change nothing
#   scripts/cut-alpha.sh --no-push    # commit + tag locally, don't push (no release)
#
# Requires: run from the repo root, on/for the `django` branch, with push access
# to the ccp4 remote (uses `gh auth token` if the plain remote push is unauth'd).
# Refuses to run with a dirty tree (except the version files it edits).

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
EXPLICIT_VERSION=""

while [ $# -gt 0 ]; do
  case "$1" in
    --dry-run) DRY_RUN=1 ;;
    --no-push) NO_PUSH=1 ;;
    --version) EXPLICIT_VERSION="${2:?--version needs an argument}"; shift ;;
    -h|--help) sed -n '2,30p' "$0"; exit 0 ;;
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

TAG="v${NEW_VER}"
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

# --- Apply the version bumps ----------------------------------------------
TODAY="$(date +'%Y, %-m, %-d')"
if [ "$DRY_RUN" = 1 ]; then
  say "DRY RUN — would set:"
  echo "    $INIT        PRERELEASE = \"$NEW_PRE\"  (+ __version_date__ = datetime($TODAY))"
  echo "    $CLIENT_VER  default pin -> \"$NEW_VER\""
  echo "    commit: 'release: ccp4i2 $NEW_VER'"
  echo "    tag+push: $TAG -> $REMOTE   (triggers release.yml)"
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

# --- Commit, tag, push -----------------------------------------------------
git add "$INIT" "$CLIENT_VER"
git commit -q -m "release: ccp4i2 $NEW_VER

Automated alpha cut via scripts/cut-alpha.sh (bump PRERELEASE + exact-pin
default in lockstep). Pushing tag $TAG triggers the release workflow.

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
git tag -a "$TAG" -m "CCP4i2 $NEW_VER"

TOKEN="$(gh auth token 2>/dev/null || true)"
PUSH_URL="$REMOTE"
[ -n "$TOKEN" ] && PUSH_URL="https://x-access-token:${TOKEN}@github.com/ccp4/ccp4i2.git"

if [ "$NO_PUSH" = 1 ]; then
  say "Committed + tagged locally ($TAG). --no-push: NOT pushing. To release:"
  echo "    git push $REMOTE $BRANCH && git push $REMOTE $TAG"
  exit 0
fi

say "Pushing $BRANCH + $TAG to $REMOTE (this triggers the release)"
git push "$PUSH_URL" "$BRANCH"
git push "$PUSH_URL" "$TAG"

say "Released: tag $TAG pushed. Watch: gh run watch --repo ccp4/ccp4i2 \$(gh run list --repo ccp4/ccp4i2 --workflow release.yml --branch $TAG --limit 1 --json databaseId -q '.[0].databaseId')"
