#!/usr/bin/env python3
"""Compare two i2run baselines test-by-test.

    python3 scripts/diff_i2run_baselines.py pre-C1 post-C1

Takes either baseline labels (resolved under server/.test-baselines/<label>/
of this checkout, or of the checkout given by --root) or paths to results.xml.
Reports tests that changed outcome, and for each newly-failing test the first
line of its failure message -- which is what triage starts from.

Stdlib only; no CCP4 environment needed.
"""
import argparse
import sys
import xml.etree.ElementTree as ET
from pathlib import Path


def resolve(label_or_path: str, root: Path) -> Path:
    candidate = Path(label_or_path)
    if candidate.is_file():
        return candidate
    for base in (root, Path.cwd()):
        guess = base / "server" / ".test-baselines" / label_or_path / "results.xml"
        if guess.is_file():
            return guess
    raise SystemExit(f"No results.xml for {label_or_path!r}")


def outcomes(results_xml: Path) -> dict:
    """{test id: (outcome, first line of any message)} from a JUnit file."""
    out = {}
    for case in ET.parse(results_xml).iter("testcase"):
        name = f"{case.get('classname', '')}::{case.get('name', '')}"
        outcome, message = "passed", ""
        for child in case:
            tag = child.tag
            if tag in ("failure", "error", "skipped"):
                outcome = {"failure": "failed", "error": "error",
                           "skipped": "skipped"}[tag]
                text = (child.get("message") or child.text or "").strip()
                message = text.splitlines()[0] if text else ""
                break
        out[name] = (outcome, message)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("before")
    ap.add_argument("after")
    ap.add_argument("--root", default=None,
                    help="checkout holding server/.test-baselines (default: "
                         "this script's repository)")
    args = ap.parse_args()

    root = Path(args.root) if args.root else Path(__file__).resolve().parent.parent
    before = outcomes(resolve(args.before, root))
    after = outcomes(resolve(args.after, root))

    def tally(d):
        counts = {}
        for outcome, _ in d.values():
            counts[outcome] = counts.get(outcome, 0) + 1
        return ", ".join(f"{v} {k}" for k, v in sorted(counts.items()))

    print(f"{args.before:>10}: {len(before):3d} tests -- {tally(before)}")
    print(f"{args.after:>10}: {len(after):3d} tests -- {tally(after)}")

    newly_failing, newly_passing, other, vanished, appeared = [], [], [], [], []
    for name in sorted(set(before) | set(after)):
        was = before.get(name)
        now = after.get(name)
        if was is None:
            appeared.append(name)
            continue
        if now is None:
            vanished.append(name)
            continue
        if was[0] == now[0]:
            continue
        if now[0] in ("failed", "error") and was[0] == "passed":
            newly_failing.append((name, now[1]))
        elif was[0] in ("failed", "error") and now[0] == "passed":
            newly_passing.append(name)
        else:
            other.append((name, was[0], now[0]))

    def section(title, rows, render):
        print(f"\n{title}: {len(rows)}")
        for row in rows:
            print(render(row))

    section("NEWLY FAILING", newly_failing,
            lambda r: f"  {r[0]}\n      {r[1][:160]}")
    section("NEWLY PASSING", newly_passing, lambda r: f"  {r}")
    section("OTHER CHANGES", other, lambda r: f"  {r[0]}: {r[1]} -> {r[2]}")
    if appeared:
        section("ONLY IN AFTER", appeared, lambda r: f"  {r}")
    if vanished:
        section("ONLY IN BEFORE", vanished, lambda r: f"  {r}")

    return 1 if newly_failing else 0


if __name__ == "__main__":
    sys.exit(main())
