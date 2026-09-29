#!/usr/bin/env bash
# Does the desktop app expect exactly the backend version this tree builds?
#
#   check-version-lockstep.sh
#
# An alpha app and its backend are strictly bound: the app pins the EXACT
# version it requires (CCP4I2_REQUIRED_SERVER_VERSION in
# client/main/ccp4i2-server-version.ts) and refuses anything else, so that pin
# and ccp4i2.__version__ must agree.
#
# Nothing checked this. release.yml's verify-version compares the TAG with
# __version__ and the ccp4i2-api lock with the wheel floor, but never the client
# pin; scripts/cut-alpha.sh edits both files together and greps to confirm, so
# the invariant held only as long as every bump went through that script. A
# hand-edit to either file would have shipped an app that refuses its own
# backend.
#
# Cheap enough to run on every pull request, which is the point: three desktop
# package builds cost ~110 billable minutes and never read this pin.
set -u

root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
init="$root/server/ccp4i2/__init__.py"
pin_file="$root/client/main/ccp4i2-server-version.ts"

for f in "$init" "$pin_file"; do
  [ -f "$f" ] || { echo "::error::missing $f"; exit 1; }
done

field() { grep -E "^$1 *=" "$init" | head -1 | sed -E 's/.*= *//; s/ *#.*//; s/"//g'; }
version="$(field MAJOR).$(field MINOR).$(field PATCH)$(field PRERELEASE)"

# export const CCP4I2_REQUIRED_SERVER_VERSION =
#   process.env.CCP4I2_SERVER_VERSION_FLOOR || "3.1.0a85";
pin="$(grep -oE 'CCP4I2_SERVER_VERSION_FLOOR *\|\| *"[^"]+"' "$pin_file" \
        | head -1 | sed -E 's/.*"([^"]+)".*/\1/')"

echo "ccp4i2.__version__      = ${version}"
echo "desktop required pin    = ${pin}"

if [ -z "$version" ] || [ -z "$pin" ]; then
  echo "::error::could not read one of the versions (version='${version}' pin='${pin}')"
  exit 1
fi

if [ "$version" != "$pin" ]; then
  cat >&2 <<MSG
::error::Version lockstep broken: ccp4i2.__version__ is '${version}' but the
desktop app requires exactly '${pin}'. An alpha app refuses any other backend,
so it would not run. Bump both together -- scripts/cut-alpha.sh does, or edit
server/ccp4i2/__init__.py and client/main/ccp4i2-server-version.ts to match.
MSG
  exit 1
fi

echo "Version lockstep OK: ${version}"
