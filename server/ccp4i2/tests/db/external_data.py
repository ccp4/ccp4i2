"""Where this tier's pre-built project zips live, and how to skip without them.

Most of ``tests/db/`` builds its fixtures by importing a CCP4i2 project zip
from ``test101/ProjectZips`` -- a checkout that sits *beside* this repository
and that nothing fetches. CLAUDE.md is explicit that tests must not depend on
it ("Unit tests must not depend on external data (test101, ProjectZips). Use
``demo_data/`` from the repo or ``I2_TOP`` for paths").

``tests/api/unit/`` already skips when the directory is absent, which is why
that tier is green in CI on a machine that has never heard of test101. This
tier never got the same treatment, so on any such machine its tests *failed*
rather than skipped -- 33 of them -- and because one ``setUp`` created its
scratch directory with a bare ``mkdir()`` and never got to its ``tearDown``,
a single missing zip cascaded into every later test in that class failing on
``FileExistsError``.

Skipping is the stopgap, not the destination: these fixtures should be
rebuilt on the repo's own ``demo_data/`` so the tier runs everywhere. Until
then, a missing zip should say so plainly instead of looking like a bug.
"""

import tempfile
from pathlib import Path

import pytest

#: server/ccp4i2/tests/db/external_data.py -> the repo root
_REPO_ROOT = Path(__file__).resolve().parent.parent.parent.parent.parent

#: The sibling checkout holding the pre-built project zips.
TEST_ZIPS_DIR = _REPO_ROOT.parent / "test101" / "ProjectZips"

SKIP_REASON = (
    f"Pre-built project zips not found: {TEST_ZIPS_DIR}. This tier's fixtures "
    "still come from the external test101 checkout; see external_data.py."
)

#: Skip a class or function that needs those zips.
requires_project_zips = pytest.mark.skipif(
    not TEST_ZIPS_DIR.exists(), reason=SKIP_REASON
)

#: Scratch directory the tests import projects into.
#:
#: Deliberately under the system temp directory. It used to be
#: ``<repo>/CCP4I2_TEST_PROJECT_DIRECTORY`` -- inside the checkout, and not
#: gitignored -- so a run that died before its tearDown left an untracked
#: directory in the working tree, ready to be committed by accident.
PROJECTS_SCRATCH_DIR = Path(tempfile.gettempdir()) / "ccp4i2_db_test_projects"
