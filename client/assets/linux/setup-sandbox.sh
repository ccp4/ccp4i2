#!/usr/bin/env bash
# Make Chromium's setuid sandbox helper usable in an unpacked ccp4i2x tree.
# Run as root once, after unpacking (a distribution's install step, or an
# admin). It is what the .deb's postinst does; nothing can do it without root,
# and the setuid bit is ignored on a filesystem mounted nosuid.
#
# Without this, ccp4i2x.sh still runs the app, but on kernels that restrict
# unprivileged user namespaces (Ubuntu 24.04+, Debian 13) it has to pass
# --no-sandbox, which gives up the renderer sandbox.
set -eu
here="$(cd "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")" && pwd)"
helper="$here/chrome-sandbox"
[ -f "$helper" ] || { echo "setup-sandbox: no $helper" >&2; exit 1; }
if [ "$(id -u)" != "0" ]; then
    echo "setup-sandbox: must run as root (sudo $0)" >&2
    exit 1
fi
chown root:root "$helper"
chmod 4755 "$helper"
echo "setup-sandbox: $helper is now setuid root"
