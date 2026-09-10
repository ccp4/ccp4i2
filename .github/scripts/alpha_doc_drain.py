#!/usr/bin/env python3
"""Drain the alpha-testing Google Doc's Feedback section into GitHub issues.

Fetches the doc's HTML export (link-shared, no auth), parses the `Feedback`
section's per-tester `<h2>` blocks, skips the placeholders, and opens one issue
per new block. Deduplicates by a hidden fingerprint marker embedded in each
issue it creates, so it needs no state file and can run from any branch.

Env: GITHUB_REPOSITORY (default ccp4/ccp4i2), DRY_RUN (true to print only).
Needs `gh` authenticated (GITHUB_TOKEN in Actions).
"""
import hashlib
import os
import re
import subprocess
import sys
import urllib.request
from html.parser import HTMLParser

DOC_ID = "1mLbsfvJV0JHdHbwOogGfytc4fWh3M-YKtGq9ulXv81M"
EXPORT = f"https://docs.google.com/document/d/{DOC_ID}/export?format=html"
DOC_URL = f"https://docs.google.com/document/d/{DOC_ID}/edit"
REPO = os.environ.get("GITHUB_REPOSITORY", "ccp4/ccp4i2")
DRY = os.environ.get("DRY_RUN", "").lower() in ("1", "true", "yes")
LABELS = ["alpha-feedback", "from-tester-doc"]
FP_RE = re.compile(r"<!--\s*alpha-doc-fp:([0-9a-f]+)\s*-->")

# Blocks that are the doc's own scaffolding, never real feedback.
PLACEHOLDER_NAMES = {"your name here"}
PLACEHOLDER_BODIES = {
    "", "anonymous feedback here…", "anonymous feedback here...",
    "your feedback here…", "your feedback here...",
}


def fetch(url):
    req = urllib.request.Request(url, headers={"User-Agent": "ccp4i2-alpha-drain"})
    with urllib.request.urlopen(req, timeout=45) as r:
        return r.read().decode("utf-8", "replace")


class Doc(HTMLParser):
    """Flatten the doc to an ordered list of (tag, text) for h1/h2/h3/p."""

    def __init__(self):
        super().__init__()
        self.items = []
        self._tag = None
        self._buf = []

    def handle_starttag(self, tag, attrs):
        if tag in ("h1", "h2", "h3", "p"):
            self._flush()
            self._tag = tag
            self._buf = []

    def handle_endtag(self, tag):
        if tag == self._tag:
            self._flush()

    def handle_data(self, data):
        if self._tag is not None:
            self._buf.append(data)

    def _flush(self):
        if self._tag is not None:
            text = re.sub(r"\s+", " ", "".join(self._buf)).strip()
            self.items.append((self._tag, text))
            self._tag = None
            self._buf = []


def feedback_entries(items):
    """(name, body) for each h2 block under the Feedback h1, before Appendix."""
    start = next((i for i, (t, x) in enumerate(items)
                  if t == "h1" and x.strip().lower() == "feedback"), None)
    if start is None:
        return []
    out = []
    i, n = start + 1, len(items)
    while i < n:
        tag, txt = items[i]
        if tag == "h1":
            break
        if tag in ("h2", "h3") and txt.strip().lower().startswith("appendix"):
            break
        if tag == "h2":
            name = txt.strip()
            body, j = [], i + 1
            while j < n and items[j][0] not in ("h1", "h2"):
                if items[j][0] in ("p", "h3") and items[j][1].strip():
                    body.append(items[j][1].strip())
                j += 1
            out.append((name, " ".join(body).strip()))
            i = j
            continue
        i += 1
    return out


def is_placeholder(name, body):
    return (name.strip().lower() in PLACEHOLDER_NAMES
            or body.strip().lower() in PLACEHOLDER_BODIES
            or not body.strip())


def fingerprint(name, body):
    h = hashlib.sha256((name + "\x1e" + body).encode("utf-8")).hexdigest()
    return h[:12]


def existing_fingerprints():
    out = subprocess.run(
        ["gh", "issue", "list", "--repo", REPO, "--label", "from-tester-doc",
         "--state", "all", "--limit", "800", "--json", "body"],
        capture_output=True, text=True)
    if out.returncode != 0:
        return set()
    import json
    seen = set()
    for issue in json.loads(out.stdout or "[]"):
        seen.update(FP_RE.findall(issue.get("body", "") or ""))
    return seen


def ensure_label():
    subprocess.run(["gh", "label", "create", "from-tester-doc",
                    "--repo", REPO, "--color", "0e8a16",
                    "--description", "Drained from the alpha-testing Google Doc"],
                   capture_output=True, text=True)


def summary(body, words=10):
    s = re.split(r"(?<=[.!?])\s", body.strip())[0]
    parts = s.split()
    return " ".join(parts[:words]) + ("…" if len(parts) > words else "")


def create_issue(name, body, fp):
    title = f"[alpha feedback] {name}: {summary(body)}"[:250]
    who = "anonymous" if name.strip().lower() == "anonymous" else name.strip()
    issue_body = (
        f"{body}\n\n"
        f"— reported by **{who}** in the "
        f"[alpha-testing feedback doc]({DOC_URL}).\n\n"
        f"<!-- alpha-doc-fp:{fp} -->"
    )
    args = ["gh", "issue", "create", "--repo", REPO, "--title", title,
            "--body-file", "-"]
    for lbl in LABELS:
        args += ["--label", lbl]
    if DRY:
        print(f"[DRY] would create: {title}")
        return
    r = subprocess.run(args, input=issue_body, capture_output=True, text=True)
    print(("created: " if r.returncode == 0 else "FAILED: ")
          + (r.stdout or r.stderr).strip().splitlines()[-1])


def main():
    doc = Doc()
    doc.feed(fetch(EXPORT))
    entries = feedback_entries(doc.items)
    real = [(n, b) for (n, b) in entries if not is_placeholder(n, b)]
    print(f"Feedback blocks: {len(entries)} total, "
          f"{len(real)} real (placeholders skipped){' [DRY RUN]' if DRY else ''}")
    if not real:
        return
    if not DRY:
        ensure_label()
    seen = set() if DRY else existing_fingerprints()
    new = 0
    for name, body in real:
        fp = fingerprint(name, body)
        if fp in seen:
            print(f"skip (already imported): {name}")
            continue
        create_issue(name, body, fp)
        new += 1
    print(f"{'would create' if DRY else 'created'} {new} new issue(s)")


if __name__ == "__main__":
    main()
