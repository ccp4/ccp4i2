#!/usr/bin/env python3
"""Measure error-handling anti-patterns across wrappers/ and pipelines/.

Regenerates the counts quoted in docs/error-handling-remediation.md so that
document's burn-down is measured rather than asserted.  Stdlib only - runs
under any Python 3.9+, no CCP4 environment needed.

    python3 scripts/scan_error_handling.py                    # summary + per-pipeline table
    python3 scripts/scan_error_handling.py --files            # also list offending files
    python3 scripts/scan_error_handling.py --predict-red-list # i2run tests C1 may flip
    python3 scripts/scan_error_handling.py --json             # machine-readable

The counts are proxies, not defects.  A bare `except:` around a cleanup path is
harmless; one wrapping a parse of the program's primary output is not.  Use the
numbers to pick a reading order, not to declare a file guilty.
"""

import argparse
import json
import os
import re
import sys
from collections import Counter, defaultdict

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCAN_ROOTS = [
    os.path.join("server", "ccp4i2", "wrappers"),
    os.path.join("server", "ccp4i2", "pipelines"),
]

# --------------------------------------------------------------------------
# Patterns.  Keys are the metric names used in the tracker document.
# --------------------------------------------------------------------------

LINE_PATTERNS = {
    "bare_except": re.compile(r"^\s*except\s*:"),
    "print_call": re.compile(r"^\s*print\s*\("),
    "return_failed": re.compile(r"return\s+(?:self|CPluginScript)\.FAILED"),
}

SRC_PATTERNS = {
    # except ...: followed immediately by pass/continue - a swallowed failure
    "swallowed": re.compile(r"except[^\n:]*:\s*\n\s*(?:pass|continue)\s*\n"),
    # Unguarded indexing into an XML query result -> IndexError on absent nodes
    "unguarded_index": re.compile(r"\.(?:findall|xpath)\([^)]*\)\s*\[\s*0\s*\]"),
    "append_error_report": re.compile(r"appendErrorReport\s*\("),
    "logger_call": re.compile(r"\blogger\.\w+\("),
}

# appendErrorReport(...) with its argument list, for the severity/code analysis
CALL_RE = re.compile(r"appendErrorReport\s*\(([^()]*(?:\([^()]*\)[^()]*)*)\)")

# def <hook>(...): plus its indented body, for per-hook return analysis
def _hook_re(name):
    return re.compile(
        r"^(\s*)def " + name + r"\s*\([^)]*\):\n(.*?)(?=\n\1def |\n\S|\Z)",
        re.S | re.M,
    )

LIFECYCLE_HOOKS = [
    "validity",
    "runTimeValidity",
    "processInputFiles",
    "makeCommandAndScript",
    "startProcess",
    "process",
    "postProcessCheck",
    "processOutputFiles",
    "checkOutputData",
]


def iter_python_files(root):
    for dirpath, dirnames, filenames in os.walk(root):
        dirnames[:] = [d for d in dirnames if d != "__pycache__"]
        for name in sorted(filenames):
            if name.endswith(".py"):
                yield os.path.join(dirpath, name)


def scan_file(path):
    """Return (counts, detail) for one file."""
    with open(path, encoding="utf-8", errors="replace") as handle:
        src = handle.read()
    lines = src.split("\n")

    counts = Counter()
    counts["loc"] = len(lines)

    for line in lines:
        for key, pattern in LINE_PATTERNS.items():
            if pattern.search(line):
                counts[key] += 1

    for key, pattern in SRC_PATTERNS.items():
        counts[key] += len(pattern.findall(src))

    # C1: processOutputFiles whose return value process() discards
    body_match = _hook_re("processOutputFiles").search(src)
    if body_match:
        body = body_match.group(2)
        if re.search(r"return\s+(?:self|CPluginScript)\.FAILED", body):
            counts["pof_returns_failed"] = 1
        if re.search(r"return\s+(?:self|CPluginScript)\.UNSATISFACTORY", body):
            counts["pof_returns_unsatisfactory"] = 1

    # C2: appendErrorReport calls that leave severity to the substring heuristic
    for match in CALL_RE.finditer(src):
        args = match.group(1).strip()
        counts["aer_total"] += 1
        if "severity" in args:
            counts["aer_explicit_severity"] += 1
        parts = [a.strip() for a in re.split(r",(?![^\[\(]*[\]\)])", args) if a.strip()]
        if len(parts) == 1:
            counts["aer_code_only"] += 1

    # Failure returns with no accompanying error report in the preceding lines
    for index, line in enumerate(lines):
        if LINE_PATTERNS["return_failed"].search(line):
            context = "\n".join(lines[max(0, index - 6): index + 1])
            if "appendErrorReport" not in context and "errorReport" not in context:
                counts["silent_failed_return"] += 1

    hooks = [h for h in LIFECYCLE_HOOKS if re.search(r"^\s*def " + h + r"\b", src, re.M)]
    return counts, hooks


def scan():
    totals = Counter()
    files_with = defaultdict(list)
    hook_files = Counter()
    per_group = defaultdict(Counter)

    for root in SCAN_ROOTS:
        abs_root = os.path.join(REPO_ROOT, root)
        if not os.path.isdir(abs_root):
            sys.exit("not found: %s (run from a CCP4i2 checkout)" % abs_root)
        is_pipelines = root.endswith("pipelines")
        for path in iter_python_files(abs_root):
            rel = os.path.relpath(path, REPO_ROOT)
            counts, hooks = scan_file(path)
            totals.update(counts)
            totals["files"] += 1
            for key, value in counts.items():
                if value and key != "loc":
                    files_with[key].append(rel)
            for hook in hooks:
                hook_files[hook] += 1
            if is_pipelines:
                # group by the directory directly under pipelines/; skip loose
                # files at the top level (pipelines/__init__.py and friends)
                tail = os.path.relpath(path, abs_root).split(os.sep)
                if len(tail) > 1:
                    per_group[tail[0]].update(counts)

    return totals, files_with, hook_files, per_group


HEADLINE = [
    ("bare_except", "bare `except:` clauses"),
    ("swallowed", "`except ...: pass` / `continue`"),
    ("unguarded_index", "`.findall(...)[0]` / `.xpath(...)[0]`"),
    ("print_call", "`print()` calls"),
    ("logger_call", "`logger.*()` calls"),
    ("silent_failed_return", "`return FAILED` with no error report nearby"),
    ("pof_returns_failed", "processOutputFiles returning FAILED (C1)"),
    ("pof_returns_unsatisfactory", "processOutputFiles returning UNSATISFACTORY (C1)"),
    ("aer_total", "appendErrorReport calls"),
    ("aer_explicit_severity", "  ... of which pass an explicit severity (C2)"),
    ("aer_code_only", "  ... of which pass a code and nothing else (C2)"),
]


def build_subjob_graph():
    """module name -> set of task names it spawns via makePluginObject('...')."""
    graph = defaultdict(set)
    spawn_re = re.compile(r"makePluginObject\(\s*['\"](\w+)['\"]")
    for root in SCAN_ROOTS:
        for path in iter_python_files(os.path.join(REPO_ROOT, root)):
            with open(path, encoding="utf-8", errors="replace") as handle:
                src = handle.read()
            module = os.path.splitext(os.path.basename(path))[0]
            graph[module].update(spawn_re.findall(src))
    return graph


def predict_red_list(files_with):
    """Which i2run tests run a task that reaches a C1-affected wrapper?

    C1 is `processOutputFiles()` returning FAILED into a `process()` that
    discards it.  Once that return is honoured, these tests can flip red.
    "Can", not "will": the wrapper's failure branch only fires on a specific
    condition, which a happy-path test usually does not hit.  A test that DOES
    flip was passing while its own wrapper said the run had failed - which is
    exactly the latent bug worth having found.
    """
    affected = {os.path.splitext(os.path.basename(p))[0]
                for p in files_with.get("pof_returns_failed", ())}
    graph = build_subjob_graph()

    def reaches(task, seen=None):
        seen = set() if seen is None else seen
        if task in seen:
            return None
        seen.add(task)
        if task in affected:
            return task
        for child in sorted(graph.get(task, ())):
            found = reaches(child, seen)
            if found:
                return found
        return None

    test_dir = os.path.join(REPO_ROOT, "server", "ccp4i2", "tests", "i2run")
    if not os.path.isdir(test_dir):
        return [], [], 0

    task_res = [re.compile(r"args\s*=\s*\[\s*['\"](\w+)['\"]"),
                re.compile(r"i2run\(\s*\[\s*['\"](\w+)['\"]")]
    at_risk, already_skipped, clear = [], [], 0
    for name in sorted(os.listdir(test_dir)):
        if not (name.startswith("test_") and name.endswith(".py")):
            continue
        with open(os.path.join(test_dir, name), encoding="utf-8", errors="replace") as handle:
            src = handle.read()
        tasks = set()
        for pattern in task_res:
            tasks.update(pattern.findall(src))
        culprits = {t: reaches(t) for t in sorted(tasks) if reaches(t)}
        if culprits:
            target = already_skipped if "mark.skip" in src else at_risk
            target.append((name, culprits))
        elif tasks:
            clear += 1
    return at_risk, already_skipped, clear


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--files", action="store_true",
                        help="list the files contributing to each metric")
    parser.add_argument("--json", action="store_true",
                        help="emit machine-readable JSON instead of a table")
    parser.add_argument("--predict-red-list", action="store_true",
                        help="i2run tests that can flip red when C1 is fixed")
    args = parser.parse_args()

    totals, files_with, hook_files, per_group = scan()

    if args.predict_red_list and not args.json:
        at_risk, already_skipped, clear = predict_red_list(files_with)
        print("Predicted C1 red list - i2run tests running a task that reaches")
        print("a wrapper whose processOutputFiles() returns FAILED.")
        print()
        for name, culprits in at_risk:
            trail = ", ".join(t if t == r else "%s -> %s" % (t, r)
                              for t, r in sorted(culprits.items()))
            print("  %-34s %s" % (name, trail))
        print()
        print("  %d active test files at risk" % len(at_risk))
        print("  %d already-skipped files would also be affected if re-enabled"
              % len(already_skipped))
        print("  %d test files run only unaffected tasks" % clear)
        print()
        print("  At risk != will fail: the wrapper's failure branch only fires on a")
        print("  specific condition a happy-path test usually misses.  A test that")
        print("  does flip was green while its own wrapper said the run had failed.")
        return

    if args.json:
        payload = {
            "totals": dict(totals),
            "files_per_metric": {k: sorted(set(v)) for k, v in files_with.items()},
            "lifecycle_hook_files": dict(hook_files),
            "per_pipeline": {k: dict(v) for k, v in sorted(per_group.items())},
            "predicted_red_list": [
                {"test": name, "reaches": culprits}
                for name, culprits in predict_red_list(files_with)[0]
            ],
        }
        print(json.dumps(payload, indent=2, sort_keys=True))
        return

    print("CCP4i2 error-handling scan")
    print("%d Python files, %d lines, under wrappers/ and pipelines/"
          % (totals["files"], totals["loc"]))
    print()
    print("%-52s %7s %7s" % ("metric", "count", "files"))
    print("-" * 68)
    for key, label in HEADLINE:
        print("%-52s %7d %7d"
              % (label, totals[key], len(set(files_with.get(key, ())))))

    print()
    print("Lifecycle hooks implemented (files defining each):")
    print("-" * 68)
    for hook in LIFECYCLE_HOOKS:
        print("  %-24s %4d" % (hook, hook_files[hook]))

    print()
    print("Per-pipeline (server/ccp4i2/pipelines/*), ranked by lines of code:")
    header = ("%-34s %7s %7s %7s %7s %7s %7s"
              % ("pipeline", "LOC", "bare", "swall", "print", "[0]idx", "retFAIL"))
    print(header)
    print("-" * len(header))
    for name, counts in sorted(per_group.items(), key=lambda kv: -kv[1]["loc"]):
        print("%-34s %7d %7d %7d %7d %7d %7d"
              % (name, counts["loc"], counts["bare_except"], counts["swallowed"],
                 counts["print_call"], counts["unguarded_index"],
                 counts["return_failed"]))

    if args.files:
        print()
        for key, label in HEADLINE:
            paths = sorted(set(files_with.get(key, ())))
            if not paths:
                continue
            print()
            print("%s (%d files):" % (label, len(paths)))
            for path in paths:
                print("    %s" % path)


if __name__ == "__main__":
    main()
