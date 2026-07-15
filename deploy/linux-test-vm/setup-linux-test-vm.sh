#!/bin/bash
# Runs ON the Ubuntu 24.04 test VM. Installs a CCP4 build, exercises the exact
# install mechanism the Electron app uses (pip-into-ccp4-python), and fetches
# BOTH desktop artifacts (AppImage + .deb) so you can compare sandbox behaviour
# on a locked-down kernel.
#
#   bash setup-linux-test-vm.sh <CCP4_BUILD_URL> [BRANCH]
#
# Deliberately does NOT abort on a failed pip install -- reproducing pip's
# complaints on a fresh build is the point. Everything is logged to
# ~/ccp4i2-setup.log so you can read the raw failure.

CCP4_URL="${1:?Usage: setup-linux-test-vm.sh <CCP4_BUILD_URL> [BRANCH]}"
BRANCH="${2:-django}"
HOME_DIR="$HOME"
CCP4_PARENT="$HOME_DIR/ccp4-test"
LOG="$HOME_DIR/ccp4i2-setup.log"
DESKTOP="$HOME_DIR/Desktop"

exec > >(tee -a "$LOG") 2>&1
echo "=== CCP4i2 Linux test VM setup  $(date -u) ==="
echo "Kernel userns restriction (1 = locked down like modern Ubuntu):"
cat /proc/sys/kernel/apparmor_restrict_unprivileged_userns 2>/dev/null || echo "  (knob absent -- not a restricted kernel)"

# --- 1. Download + extract CCP4 -------------------------------------------
mkdir -p "$CCP4_PARENT"
TARBALL="$CCP4_PARENT/$(basename "${CCP4_URL%%\?*}")"
if [[ ! -f "$TARBALL" ]]; then
  echo "--- Downloading CCP4 build ---"
  echo "    $CCP4_URL"
  wget --continue --show-progress -O "$TARBALL" "$CCP4_URL"
else
  echo "[OK] CCP4 tarball already present: $TARBALL"
fi

echo "--- Extracting (this takes a few minutes) ---"
tar -xf "$TARBALL" -C "$CCP4_PARENT"

# Find the CCP4 root (the dir containing bin/ccp4.setup-sh)
CCP4_SETUP="$(find "$CCP4_PARENT" -maxdepth 3 -name ccp4.setup-sh -path '*/bin/*' 2>/dev/null | head -1)"
if [[ -z "$CCP4_SETUP" ]]; then
  echo "ERROR: could not find bin/ccp4.setup-sh under $CCP4_PARENT" >&2
  echo "       Extracted contents:" >&2
  ls -la "$CCP4_PARENT" >&2
  exit 1
fi
CCP4_ROOT="$(cd "$(dirname "$CCP4_SETUP")/.." && pwd)"
echo "[OK] CCP4 root: $CCP4_ROOT"

# Some tarball builds ship a relocation fixup (BINARY.setup / install.sh).
if [[ -x "$CCP4_ROOT/BINARY.setup" ]]; then
  echo "--- Running BINARY.setup (path relocation) ---"
  ( cd "$CCP4_ROOT" && ./BINARY.setup ) || echo "[WARN] BINARY.setup returned non-zero (may be fine)"
fi

# --- 2. Verify ccp4-python + gemmi ----------------------------------------
echo "--- Sourcing CCP4 environment ---"
# shellcheck disable=SC1090
source "$CCP4_SETUP"

echo "--- ccp4-python sanity ---"
ccp4-python --version || echo "[WARN] ccp4-python --version failed"
ccp4-python -c "import gemmi; print('gemmi', gemmi.__version__, 'OK')" \
  || echo "[WARN] gemmi import failed"

# --- 3. Exercise the install mechanism (the flaky bit) --------------------
echo ""
echo "=== ATTEMPTING: ccp4-python -m pip install ccp4i2 (verbose) ==="
echo "    Any failure below is the thing we're here to diagnose."
set -x
ccp4-python -m pip install --upgrade pip 2>&1 | tail -5
ccp4-python -m pip install --verbose ccp4i2
PIP_RC=$?
set +x
echo "=== pip install ccp4i2 exit code: $PIP_RC ==="
if [[ $PIP_RC -ne 0 ]]; then
  echo "[EXPECTED-POSSIBLE] pip failed. Full transcript in $LOG."
  echo "--- pip check (pre-existing breakage in this build) ---"
  ccp4-python -m pip check || true
fi

# --- 4. Fetch BOTH desktop artifacts (AppImage + .deb) from CI ------------
echo ""
echo "--- Fetching CCP4i2 Linux artifacts from CI (if gh is available/authed) ---"
mkdir -p "$DESKTOP"
if command -v gh >/dev/null 2>&1 && gh auth status >/dev/null 2>&1; then
  gh run download --repo ccp4/ccp4i2 --name linux-appimages --dir "$DESKTOP" \
    && chmod +x "$DESKTOP"/*.AppImage 2>/dev/null \
    && echo "[OK] Artifacts in $DESKTOP:" && ls -1 "$DESKTOP"/*.AppImage "$DESKTOP"/*.deb 2>/dev/null \
    || echo "[WARN] gh download failed -- grab the artifacts manually."
else
  echo "gh not installed/authed. Download the Linux artifacts manually from:"
  echo "  https://github.com/ccp4/ccp4i2/actions/workflows/electron-multiplatform-build.yml"
  echo "  (artifact: linux-appimages -- contains BOTH the .AppImage and the .deb)"
fi

# --- 4b. Install the .deb (the sandbox-correct path) ----------------------
DEB="$(ls "$DESKTOP"/*.deb 2>/dev/null | head -1)"
if [[ -n "$DEB" ]]; then
  echo ""
  echo "--- Installing the .deb (postinst setuids chrome-sandbox) ---"
  echo "    sudo apt-get install -y '$DEB'"
  sudo apt-get install -y "$DEB" || sudo dpkg -i "$DEB" || echo "[WARN] .deb install failed -- see above."
  echo "chrome-sandbox perms (want '-rwsr-xr-x root root'):"
  ls -l /opt/ccp4i2-django/chrome-sandbox 2>/dev/null || echo "  (chrome-sandbox not found at expected path)"
fi

# --- 5. Persist the CCP4 env for future logins ----------------------------
if ! grep -q 'ccp4.setup-sh' "$HOME_DIR/.bashrc" 2>/dev/null; then
  echo "source '$CCP4_SETUP'" >> "$HOME_DIR/.bashrc"
  echo "[OK] Added CCP4 source line to ~/.bashrc"
fi

cat <<EOF

=== Setup done ===
CCP4 root:   $CCP4_ROOT
Full log:    $LOG

Compare sandbox behaviour on this locked-down kernel:
  - .deb (should work WITH the sandbox, no flags):  ccp4i2-django
    (launch from an RDP terminal; the app is installed under /opt/ccp4i2-django)
  - AppImage (the format that struggles here):      $DESKTOP/*.AppImage

Headless install path:
  ccp4-python -m pip install ccp4i2

EOF
