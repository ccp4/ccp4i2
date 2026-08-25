"""Make a third-party HTML report self-contained inside its job directory.

Programs like MrParse emit HTML that points at their own installation:
``<script src="/abs/path/to/site-packages/mrparse/html/vue.min.js">`` and
similar. That works when the file is opened straight from disk on the
machine that produced it — what the Qt-era CCP4i2 did — and fails every
other way:

* Served over HTTP the absolute path is interpreted as a URL path on the
  server origin, so every asset 404s and the report renders as unstyled
  headings.
* The path names *this machine's* CCP4 installation, so the report breaks
  as soon as the project is exported, imported elsewhere, or relocated.

Localising means copying the referenced directory into the job directory
and rewriting the references to be relative to the HTML file. Relative is
the load-bearing choice — it encodes neither a filesystem location nor a
project or job identity, so the report survives being zipped up by
project export, re-imported under a different project id, moved to a new
project root, or opened directly from disk. (An earlier rewrite in
``mrparse_report`` pointed assets at ``/database/projectid/N/jobnumber/M/
file/...``, which is both a Qt-era route that no longer exists and an
identity that changes on import.)

It also keeps every request inside the report's own directory, which is
the subtree a file grant covers — see ``ccp4i2_api.file_grants``.

The cost is duplication: MrParse's asset directory is ~1.1 MB per job.
That buys a job directory that renders anywhere, and it is small beside
the coordinate and reflection files sitting next to it.
"""

import logging
import re
import shutil
from pathlib import Path
from typing import Mapping, Optional

logger = logging.getLogger(f"ccp4i2:{__name__}")


def vendored_asset(name: str) -> Path:
    """Path to an asset shipped with CCP4i2 under ``docs/report_files``.

    Used for libraries a third-party report loads from a CDN. Vendoring
    them keeps a report readable on an offline beamline machine and years
    after the CDN stopped serving that version — a job report is a record,
    not a live page.
    """
    from ccp4i2.report.core import CURRENT_CSS_VERSION

    return (
        Path(__file__).parent.parent
        / "docs"
        / "report_files"
        / CURRENT_CSS_VERSION
        / name
    )


def discover_install_dir(html_text: str, marker: str) -> Optional[Path]:
    """Find the absolute install directory an HTML report references.

    ``marker`` is the tail the directory ends with (e.g. ``mrparse/html``).
    Returns None when the report carries no such absolute reference, which
    is also how a report that has already been localised looks.
    """
    pattern = re.compile(r'["\'](/(?:[^"\']*?/)?' + re.escape(marker) + r')(?:/[^"\']*)?["\']')
    for match in pattern.finditer(html_text):
        candidate = Path(match.group(1))
        if candidate.is_dir():
            return candidate
    return None


def localise_report_assets(
    html_path: Path,
    *,
    marker: str,
    subdirectory: str,
    extra_assets: Optional[Mapping[str, Path]] = None,
    output_path: Optional[Path] = None,
) -> Optional[Path]:
    """Copy an HTML report's asset directory in beside it and relativise it.

    Args:
        html_path: the report as the program wrote it. Left untouched.
        marker: tail of the install directory to look for, e.g. "mrparse/html".
        subdirectory: name to copy that directory to, beside the report.
        extra_assets: URLs appearing in the report mapped to local files to
            serve in their place — CDN references, which are localised the
            same way and for the same reasons as the install-directory ones.
        output_path: where to write the rewritten HTML (default: the input
            name with an ``_i2`` suffix, matching the convention the xia2
            reports already use).

    Returns the rewritten file's path, or None if there was nothing to
    localise — no absolute references, or the install directory is gone
    (a project exported from another machine, say, where the copied assets
    are already in place and the rewritten HTML came along with them).
    """
    html_path = Path(html_path)
    if output_path is None:
        output_path = html_path.with_name(f"{html_path.stem}_i2{html_path.suffix}")

    html_text = html_path.read_text(encoding="utf-8", errors="replace")
    install_dir = discover_install_dir(html_text, marker)
    if install_dir is None:
        return output_path if output_path.exists() else None

    # Report classes re-run on every view of a finished job, so only pay the
    # copy once — the assets are a fixed part of the installation, not job output.
    destination = html_path.parent / subdirectory
    already_copied = destination.is_dir() and any(destination.iterdir())
    if not already_copied:
        try:
            shutil.copytree(install_dir, destination, dirs_exist_ok=True)
        except OSError as err:
            logger.warning("Could not copy %s to %s: %s", install_dir, destination, err)
            return None

    # One textual substitution covers both the markup references and any
    # JavaScript that builds paths from the same string at runtime (MrParse
    # sets a `const mrparse_html_dir` the Vue bundle then interpolates).
    localised = html_text.replace(str(install_dir), subdirectory)

    for url, source in (extra_assets or {}).items():
        if url not in localised:
            continue
        source = Path(source)
        if not source.is_file():
            logger.warning("Vendored asset %s is missing; leaving %s alone", source, url)
            continue
        target = destination / source.name
        if not target.exists():
            shutil.copyfile(source, target)
        localised = localised.replace(url, f"{subdirectory}/{source.name}")

    output_path.write_text(localised, encoding="utf-8")
    return output_path
