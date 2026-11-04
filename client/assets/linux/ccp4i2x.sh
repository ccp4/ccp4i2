#!/usr/bin/env bash
# Launch ccp4i2x from an unpacked tree (no installer, no AppImage, no .deb).
#
# Chromium needs a sandbox for its renderer processes and can get one two ways:
# unprivileged user namespaces (no setup, but Ubuntu 24.04+ and Debian 13
# restrict them), or the setuid-root helper chrome-sandbox beside the binary
# (what the .deb's postinst arranges; nothing can arrange it in a tree a user
# unpacked). When neither is available Chromium aborts rather than run
# unsandboxed, so this launcher checks both and only then adds --no-sandbox.
# The flag must be on the command line: Chromium's zygote reads it before any
# of the app's own code runs.
#
# Any arguments are passed through to the app.

here="$(cd "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")" && pwd)"
app="$here/ccp4i2x"
helper="$here/chrome-sandbox"

sandbox_ok() {
    # The setuid helper: owned by root, mode 4755, on a filesystem that
    # honours setuid (an NFS or home mount with nosuid silently does not).
    if [ -u "$helper" ] && [ "$(stat -c %u "$helper" 2>/dev/null)" = "0" ]; then
        return 0
    fi
    # Otherwise unprivileged user namespaces must be allowed.
    local restrict clone
    restrict="$(cat /proc/sys/kernel/apparmor_restrict_unprivileged_userns 2>/dev/null || echo 0)"
    clone="$(cat /proc/sys/kernel/unprivileged_userns_clone 2>/dev/null || echo 1)"
    [ "$restrict" = "0" ] && [ "$clone" != "0" ]
}

if sandbox_ok; then
    exec "$app" "$@"
else
    echo "ccp4i2x: no usable Chromium sandbox on this kernel (chrome-sandbox is not" >&2
    echo "setuid root and unprivileged user namespaces are restricted); running with" >&2
    echo "--no-sandbox. To keep the sandbox, as root: chown root:root '$helper' &&" >&2
    echo "chmod 4755 '$helper', or sysctl -w kernel.apparmor_restrict_unprivileged_userns=0" >&2
    exec "$app" --no-sandbox "$@"
fi
