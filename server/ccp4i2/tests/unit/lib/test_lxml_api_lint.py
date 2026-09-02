"""The lxml-only-API check, and a ratchet holding the tree clean.

`cmapcoeff` passed lxml's `pretty_print` to `xml.etree.ElementTree` for years.
Every such line raises the first time it runs, and none of them is visible to a
type checker: `CData` and the report nodes define `__getattr__`, so any
attribute must be assumed to exist. A grep is therefore the only check
available, and this is the test of the grep.

The last test is the ratchet. Running the script is a CI step, but somebody
running pytest locally should learn about a new one too.
"""

import importlib.util
import sys
from pathlib import Path

import pytest


def _repo_root() -> Path:
    for candidate in Path(__file__).resolve().parents:
        if (candidate / "scripts" / "check_lxml_api.py").is_file():
            return candidate
    return None


REPO_ROOT = _repo_root()

pytestmark = pytest.mark.skipif(
    REPO_ROOT is None,
    reason="scripts/check_lxml_api.py is not in an installed package",
)


@pytest.fixture(scope="module")
def lint():
    spec = importlib.util.spec_from_file_location(
        "check_lxml_api", REPO_ROOT / "scripts" / "check_lxml_api.py"
    )
    module = importlib.util.module_from_spec(spec)
    sys.modules["check_lxml_api"] = module
    spec.loader.exec_module(module)
    return module


def write(tmp_path: Path, body: str) -> Path:
    path = tmp_path / "module_under_test.py"
    path.write_text(body)
    return path


def kinds(findings):
    return [(f[2], f[6]) for f in findings]


class TestTheLxmlOnlyRule:
    def test_flags_pretty_print_in_a_standard_library_module(self, lint, tmp_path):
        path = write(
            tmp_path,
            "import xml.etree.ElementTree as etree\n"
            "def f(root):\n"
            "    return etree.tostring(root, pretty_print=True)\n",
        )
        assert kinds(lint.scan_file(str(path), "m.py")) == [
            ("pretty_print", "no-lxml-import")
        ]

    def test_allows_it_where_lxml_is_actually_imported(self, lint, tmp_path):
        path = write(
            tmp_path,
            "from lxml import etree\n"
            "def f(root):\n"
            "    return etree.tostring(root, pretty_print=True)\n",
        )
        assert lint.scan_file(str(path), "m.py") == []

    def test_import_lxml_dot_something_also_counts(self, lint, tmp_path):
        path = write(
            tmp_path, "import lxml.etree as etree\nx = etree.tostring(1, pretty_print=True)\n"
        )
        assert lint.scan_file(str(path), "m.py") == []

    def test_lxml_named_only_in_a_comment_does_not_exempt(self, lint, tmp_path):
        """The difference between a real finding and a missed one."""
        path = write(
            tmp_path,
            "# this file used to use lxml\n"
            '"""Ported away from lxml."""\n'
            "import xml.etree.ElementTree as etree\n"
            "x = etree.tostring(1, pretty_print=True)\n",
        )
        assert kinds(lint.scan_file(str(path), "m.py")) == [
            ("pretty_print", "no-lxml-import")
        ]

    @pytest.mark.parametrize(
        "line,expected",
        [
            ("node.xpath('a/b')", "xpath"),
            ("node.getparent()", "getparent"),
            ("node.iterancestors()", "iterancestors"),
            ("line = node.sourceline", "sourceline"),
            ("node.text = CDATA('x')", "CDATA"),
        ],
    )
    def test_the_rest_of_the_banned_surface(self, lint, tmp_path, line, expected):
        path = write(tmp_path, f"import xml.etree.ElementTree as etree\n{line}\n")
        assert kinds(lint.scan_file(str(path), "m.py")) == [
            (expected, "no-lxml-import")
        ]


class TestAlwaysBanned:
    def test_xpath0_is_flagged_even_in_an_lxml_module(self, lint, tmp_path):
        """lxml has no xpath0 either -- it was a Qt-era report-parser helper."""
        path = write(tmp_path, "from lxml import etree\nx = node.xpath0('a')\n")
        assert kinds(lint.scan_file(str(path), "m.py")) == [("xpath0", "always")]

    def test_lxml_html_clean_is_flagged_in_an_lxml_module(self, lint, tmp_path):
        path = write(
            tmp_path, "from lxml import etree\nfrom lxml.html.clean import Cleaner\n"
        )
        assert kinds(lint.scan_file(str(path), "m.py")) == [
            ("lxml.html.clean", "always")
        ]


class TestEscapeHatches:
    def test_an_annotated_line_is_allowed(self, lint, tmp_path):
        path = write(
            tmp_path,
            "import xml.etree.ElementTree as etree\n"
            "x = etree.tostring(1, pretty_print=True)  # lxml-ok: root is lxml here\n",
        )
        assert lint.scan_file(str(path), "m.py") == []

    def test_commented_out_code_is_reported_but_inert(self, lint, tmp_path):
        path = write(
            tmp_path,
            "import xml.etree.ElementTree as etree\n"
            "# x = etree.tostring(1, pretty_print=True)\n",
        )
        findings = lint.scan_file(str(path), "m.py")
        assert len(findings) == 1
        assert findings[0][4] is True, "a commented-out call cannot raise"

    def test_a_file_that_will_not_parse_is_skipped(self, lint, tmp_path):
        path = write(tmp_path, "print 'python 2'\nx.xpath('a')\n")
        assert lint.scan_file(str(path), "m.py") == []


def test_the_tree_is_clean(lint):
    """The ratchet: no live finding anywhere in the shipped package.

    Fixing one is usually a one-line port -- `.xpath(p)` to `.findall(p)`,
    `pretty_print=True` to an `ET.indent()` beforehand. Where the object really
    is an lxml one, import lxml deliberately or annotate the line with
    `# lxml-ok: why`.
    """
    findings = []
    for path, relpath in lint.walk(lint.SCAN_ROOTS):
        findings.extend(f for f in lint.scan_file(path, relpath) if not f[4])

    assert findings == [], "\n".join(
        f"{f[0]}:{f[1]}: {f[2]} -- {f[3]}" for f in findings
    )
