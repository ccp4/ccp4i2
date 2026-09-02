#!/usr/bin/env python3
"""Fail on lxml-only APIs in modules that do not import lxml.

CCP4i2 was written against lxml and has been moving to the standard library:
109 modules import lxml, 134 import ``xml.etree.ElementTree``, 8 import both.
Most report modules say ``import xml.etree.ElementTree as etree``, so a reader
who knows the history sees ``etree.`` and writes ``pretty_print=True`` under a
reasonable assumption that nothing in the file contradicts.

The result raises ``TypeError`` the first time that line runs, and until
2026-08-25 the framework printed the exception to a console nobody reads and
marked the job Finished. `cmapcoeff` carried such a line for years -- no
cmapcoeff job ever wrote a ``program.xml`` -- and it was found only when C1
started failing jobs on it. Four more instances came out of one grep. That grep
is this file.

The rule: **an lxml-only API in a module that does not import lxml is an
error.** It is a grep rather than a type check because it has to be. `CData`
and the report nodes are dynamic -- both define ``__getattr__``, so a type
checker must assume any attribute may exist, and this whole family is
structurally invisible to static analysis.

Stdlib only, no CCP4 environment needed:

    python3 scripts/check_lxml_api.py           # the gate; exit 1 on a finding
    python3 scripts/check_lxml_api.py --all     # also list inert (commented) hits
    python3 scripts/check_lxml_api.py --json

To allow a line the check cannot see is fine -- an lxml object handed in from
another module, say -- put ``# lxml-ok: why`` on it. The reason is required.
"""

import argparse
import ast
import json
import os
import re
import sys

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

SCAN_ROOTS = [os.path.join("server", "ccp4i2")]

# Directories whose contents are not shipped code.
SKIP_DIRS = {"__pycache__", ".git", "node_modules", "tests"}

# --------------------------------------------------------------------------
# The banned surface: lxml APIs with no standard-library equivalent, each of
# which raises rather than misbehaving. Measured on lxml 4.9.4 in
# ccp4-20260702, every cross-library operation raises TypeError at once, so
# what is being searched for is an API mismatch, not a mixed tree.
# --------------------------------------------------------------------------

LXML_ONLY = {
    "pretty_print": (
        re.compile(r"\bpretty_print\s*="),
        "lxml keyword; stdlib formats with ET.indent(tree) before serialising",
    ),
    "xpath": (
        re.compile(r"\.xpath\s*\("),
        "lxml method; stdlib has .find()/.findall() with ElementPath syntax",
    ),
    "getparent": (
        re.compile(r"\.getparent\s*\("),
        "lxml method; stdlib elements do not know their parent",
    ),
    "iterancestors": (
        re.compile(r"\.iterancestors\s*\("),
        "lxml method; stdlib elements do not know their parent",
    ),
    "sourceline": (
        re.compile(r"\.sourceline\b"),
        "lxml attribute; stdlib elements do not record their source line",
    ),
    "CDATA": (
        re.compile(r"\bCDATA\s*\("),
        "lxml constructor; stdlib has no CDATA support",
    ),
}

# Banned everywhere, whatever a module imports, because no library provides
# them any more. This is where port residue goes: when a port removes an API,
# add its names here so the next caller is caught by CI rather than by a user.
ALWAYS_BANNED = {
    "lxml.html.clean": (
        re.compile(r"\blxml\.html\.clean\b|\bfrom\s+lxml\.html\s+import\s+clean\b"),
        "moved out of lxml in 5.x; import lxml_html_clean with a fallback",
    ),
    "xpath0": (
        re.compile(r"\.xpath0\s*\("),
        "Qt-era CCP4ReportParser helper; nothing defines it now. Use .find()",
    ),
}

ALLOW_RE = re.compile(r"#\s*lxml-ok:\s*(\S.*)")


def imports_lxml(tree: ast.AST) -> bool:
    """True if the module imports lxml under any name.

    Parsed rather than grepped: a comment or a docstring mentioning lxml must
    not count as importing it, and that is the difference between a real
    finding and a missed one.
    """
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            if any(a.name == "lxml" or a.name.startswith("lxml.") for a in node.names):
                return True
        elif isinstance(node, ast.ImportFrom):
            module = node.module or ""
            if module == "lxml" or module.startswith("lxml."):
                return True
    return False


def scan_file(path: str, relpath: str):
    """Findings for one file: (relpath, line, kind, why, inert, source, rule)."""
    with open(path, encoding="utf-8", errors="replace") as handle:
        source = handle.read()

    try:
        tree = ast.parse(source)
    except SyntaxError:
        # Python 2 leftovers exist in this tree; nothing to say about them.
        return []

    uses_lxml = imports_lxml(tree)
    findings = []

    for number, line in enumerate(source.splitlines(), start=1):
        stripped = line.strip()
        allowed = ALLOW_RE.search(line)
        # Commented-out code cannot raise, but it is how the next copy-paste
        # starts, so it is reported rather than ignored -- just not fatally.
        inert = stripped.startswith("#")

        if allowed:
            continue

        # The lxml-only set is a finding only where the module never imports
        # lxml. Note that one deliberate import exempts the whole file: the
        # rule is a grep, so prefer the narrowest fix that keeps a module
        # standard-library throughout.
        candidates = [] if uses_lxml else list(LXML_ONLY.items())
        for kind, (pattern, why) in candidates:
            if pattern.search(line):
                findings.append(
                    (relpath, number, kind, why, inert, stripped, "no-lxml-import")
                )

        for kind, (pattern, why) in ALWAYS_BANNED.items():
            if pattern.search(line):
                findings.append(
                    (relpath, number, kind, why, inert, stripped, "always")
                )

    return findings


def walk(roots):
    for root in roots:
        absolute = os.path.join(REPO_ROOT, root)
        for directory, subdirs, files in os.walk(absolute):
            subdirs[:] = [d for d in subdirs if d not in SKIP_DIRS]
            for name in sorted(files):
                if not name.endswith(".py"):
                    continue
                path = os.path.join(directory, name)
                yield path, os.path.relpath(path, REPO_ROOT)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--all",
        action="store_true",
        help="also list inert (commented-out) occurrences, which never fail the gate",
    )
    parser.add_argument("--json", action="store_true", help="machine-readable output")
    parser.add_argument(
        "paths",
        nargs="*",
        help="files or directories to scan (default: the whole package)",
    )
    args = parser.parse_args(argv)

    roots = args.paths or SCAN_ROOTS
    findings = []
    for path, relpath in walk(roots):
        findings.extend(scan_file(path, relpath))

    live = [f for f in findings if not f[4]]
    inert = [f for f in findings if f[4]]

    if args.json:
        print(
            json.dumps(
                {
                    "live": [
                        {"file": f[0], "line": f[1], "api": f[2], "why": f[3]}
                        for f in live
                    ],
                    "inert": [
                        {"file": f[0], "line": f[1], "api": f[2]} for f in inert
                    ],
                },
                indent=2,
            )
        )
        return 1 if live else 0

    for relpath, number, kind, why, _, source, rule in live:
        context = (
            "in a module that does not import lxml"
            if rule == "no-lxml-import"
            else "which no library provides any more"
        )
        print(f"{relpath}:{number}: '{kind}' {context}")
        print(f"    {source}")
        print(f"    {why}")

    if args.all and inert:
        print()
        print("Commented out, so inert -- but the next copy-paste starts here:")
        for relpath, number, kind, _, _, source, _rule in inert:
            print(f"  {relpath}:{number}: {kind}  |  {source[:90]}")

    if live:
        print()
        print(f"{len(live)} call(s) to an API that is not there.")
        print("Each raises the first time its line runs. Port it to the standard")
        print("library, import lxml deliberately, or annotate with '# lxml-ok: why'.")
        return 1

    print(f"No missing-API calls ({len(inert)} inert occurrence(s) in comments).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
