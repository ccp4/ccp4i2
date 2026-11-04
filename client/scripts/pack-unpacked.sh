#!/usr/bin/env bash
# Archive the unpacked app trees electron-builder leaves in release/ so the app
# can be launched without going through an installer: what a distribution
# (CCP4's) needs to drop into its own tree. Run after `npm run package-*`.
#
#   client/scripts/pack-unpacked.sh            # all trees found in release/
#
# electron-builder always materialises the unpacked tree before wrapping it in
# a dmg/exe/AppImage/deb: release/mac*/ccp4i2x.app, release/linux*-unpacked,
# release/win*-unpacked. This packs each one it finds into
# release/unpacked/ccp4i2x-<version>-<platform>-<arch>.{tar.gz,zip}:
#
#   macOS    the signed (and, on a release, notarised) ccp4i2x.app, as a
#            tarball so symlinks and extended attributes survive; double-click
#            or `open` it. Keep tar for unpacking too (Finder is fine).
#   Linux    ccp4i2x/ with the binary, ccp4i2x.sh (launcher that sorts out the
#            Chromium sandbox) and setup-sandbox.sh (root step that makes
#            chrome-sandbox setuid, as the .deb's postinst does).
#   Windows  ccp4i2x/ with ccp4i2x.exe; run it directly.
#
# Environment:
#   RELEASE_DIR        electron-builder's output dir (default: client/release)
#   UNPACKED_VERSION   version in the archive name (default: package.json's)
#
# Runs under bash on all three runners (Git Bash on Windows, where 7-Zip makes
# the zip; `zip` is the fallback).
set -euo pipefail

client="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
release="${RELEASE_DIR:-$client/release}"
version="${UNPACKED_VERSION:-$(node -p "require('$client/package.json').version")}"
assets="$client/assets/linux"
out="$release/unpacked"
stage="$out/stage"
rm -rf "$out"
mkdir -p "$out" "$stage"

readme() {  # <platform>  -> README.txt text
    cat <<EOF
ccp4i2x $version, unpacked ($1)

This is the CCP4i2 desktop app as a plain directory tree: no installer, no
auto-update packaging. It is what a distribution puts in its own tree.

The app needs a CCP4 installation whose ccp4-python has the matching
ccp4i2==$version backend installed (pip install --no-deps, see
server/ccp4i2/requirements-runtime.txt in the ccp4i2 repository). It looks
for CCP4 in \$CCP4, on PATH (bin/ccp4-python) and in the usual places, and
asks on first launch if it finds none.

EOF
    case "$1" in
        macos)   echo "Launch: open ccp4i2x.app (or copy it to /Applications)." ;;
        linux)   cat <<'EOF'
Launch: ccp4i2x/ccp4i2x.sh
  The launcher runs the app sandboxed when the kernel allows it, and falls
  back to --no-sandbox (printing why) when it does not. To keep the sandbox
  on kernels that restrict unprivileged user namespaces (Ubuntu 24.04+,
  Debian 13), run once as root after unpacking:
    sudo ccp4i2x/setup-sandbox.sh
  which makes chrome-sandbox setuid root, exactly as the .deb's postinst does.
  The setuid bit is ignored on a filesystem mounted nosuid.
EOF
        ;;
        windows) echo "Launch: ccp4i2x\\ccp4i2x.exe" ;;
    esac
}

arch_of() {  # electron-builder's output dir name -> arch
    case "$1" in
        mac|linux-unpacked|win-unpacked) echo x64 ;;
        *universal*)                     echo universal ;;
        *arm64*)                         echo arm64 ;;
        *ia32*)                          echo ia32 ;;
        *)                               echo unknown ;;
    esac
}

made=0

# macOS: the .app bundle, tarred in place (no copy of a 1 GB tree).
for app in "$release"/mac*/ccp4i2x.app; do
    [ -d "$app" ] || continue
    dir="$(basename "$(dirname "$app")")"
    arch="$(arch_of "$dir")"
    name="ccp4i2x-$version-macos-$arch"
    mkdir -p "$stage/$name"
    readme macos > "$stage/$name/README.txt"
    tar -czf "$out/$name.tar.gz" -C "$(dirname "$app")" ccp4i2x.app -C "$stage/$name" README.txt
    made=$((made + 1))
done

# Linux: stage as ccp4i2x/ with the launcher and the root-time sandbox step.
for tree in "$release"/linux*-unpacked; do
    [ -d "$tree" ] || continue
    arch="$(arch_of "$(basename "$tree")")"
    name="ccp4i2x-$version-linux-$arch"
    mkdir -p "$stage/$name"
    cp -Rp "$tree" "$stage/$name/ccp4i2x"
    cp "$assets/ccp4i2x.sh" "$assets/setup-sandbox.sh" "$stage/$name/ccp4i2x/"
    chmod 755 "$stage/$name/ccp4i2x/ccp4i2x.sh" "$stage/$name/ccp4i2x/setup-sandbox.sh"
    readme linux > "$stage/$name/README.txt"
    tar -czf "$out/$name.tar.gz" -C "$stage/$name" ccp4i2x README.txt
    made=$((made + 1))
done

# Windows: stage as ccp4i2x\ and zip.
for tree in "$release"/win*-unpacked; do
    [ -d "$tree" ] || continue
    arch="$(arch_of "$(basename "$tree")")"
    name="ccp4i2x-$version-windows-$arch"
    mkdir -p "$stage/$name"
    cp -Rp "$tree" "$stage/$name/ccp4i2x"
    readme windows > "$stage/$name/README.txt"
    if command -v 7z >/dev/null 2>&1; then
        (cd "$stage/$name" && 7z a -tzip -bso0 -bsp0 "$out/$name.zip" ccp4i2x README.txt)
    else
        (cd "$stage/$name" && zip -qr "$out/$name.zip" ccp4i2x README.txt)
    fi
    made=$((made + 1))
done

rm -rf "$stage"
if [ "$made" = 0 ]; then
    echo "pack-unpacked: no unpacked tree under $release (run npm run package-* first)" >&2
    exit 1
fi
ls -lh "$out"
