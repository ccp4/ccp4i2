#!/bin/bash
# Runs as root on remove/purge. Unload and remove the AppArmor profile the
# postinst installed. (chrome-sandbox lives under /opt and is removed with the
# package, so nothing to undo there.)
#
# NOTE: no braced shell variables — electron-builder macro-substitutes this file
# (comments included) and a braced variable fails the build. Use brace-less $VAR.
set -e

PROFILE=/etc/apparmor.d/ccp4i2-django
if [ -f "$PROFILE" ]; then
  if command -v apparmor_parser >/dev/null 2>&1; then
    apparmor_parser -R "$PROFILE" 2>/dev/null || true
  fi
  rm -f "$PROFILE"
fi

exit 0
