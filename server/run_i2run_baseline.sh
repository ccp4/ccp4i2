#!/bin/bash
# Full i2run suite runner — captures per-test results (JUnit XML) + a verbose
# log so two CCP4 builds can be diffed test-by-test (e.g. to spot regressions
# when migrating CCP4 versions).
#
# Usage:
#   CCP4_SETUP=/path/to/ccp4-XXXX/bin/ccp4.setup-sh \
#   CCP4_LABEL=ccp4-XXXX \
#   bash server/run_i2run_baseline.sh
#
# Outputs to server/.test-baselines/<CCP4_LABEL>/{results.xml,pytest.log}.
# Diff two runs:
#   diff <(...) — or compare the two results.xml with any JUnit differ.
#
# NOTE: no `set -u` — ccp4.setup-sh references unbound vars (e.g. MANPATH) and
# would abort under nounset.

# Server dir = this script's directory (portable; no hardcoded user paths).
SERVER_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

CCP4_SETUP="${CCP4_SETUP:?Set CCP4_SETUP=/path/to/ccp4-XXXX/bin/ccp4.setup-sh}"
CCP4_LABEL="${CCP4_LABEL:-$(basename "$(dirname "$(dirname "${CCP4_SETUP}")")")}"

OUT_DIR="${SERVER_DIR}/.test-baselines/${CCP4_LABEL}"
mkdir -p "${OUT_DIR}"
LOG="${OUT_DIR}/pytest.log"
XML="${OUT_DIR}/results.xml"

# shellcheck disable=SC1090
source "${CCP4_SETUP}"
cd "${SERVER_DIR}" || exit 99

{
  echo "==== i2run run ===="
  echo "date:        $(date)"
  echo "ccp4 setup:  ${CCP4_SETUP}"
  echo "ccp4-python: $(command -v ccp4-python)"
  echo "git HEAD:    $(git rev-parse --short HEAD 2>/dev/null) ($(git rev-parse --abbrev-ref HEAD 2>/dev/null))"
  echo "==================="
} | tee "${LOG}"

# Serial run (no xdist) for deterministic, attributable per-test results.
ccp4-python -m pytest ccp4i2/tests/i2run/ \
  -v -ra --tb=short \
  --durations=25 \
  --continue-on-collection-errors \
  --junitxml="${XML}" \
  2>&1 | tee -a "${LOG}"

STATUS=${PIPESTATUS[0]}
echo "==== pytest exit status: ${STATUS} ====" | tee -a "${LOG}"
exit "${STATUS}"
