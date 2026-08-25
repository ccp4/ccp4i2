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
    subdirectory: str,
    marker: Optional[str] = None,
    install_dir: Optional[Path] = None,
    extra_assets: Optional[Mapping[str, Path]] = None,
    url_rewrites: Optional[Mapping[str, str]] = None,
    output_path: Optional[Path] = None,
) -> Optional[Path]:
    """Copy an HTML report's asset directory in beside it and relativise it.

    Args:
        html_path: the report as the program wrote it. Left untouched.
        subdirectory: name to copy that directory to, beside the report.
        marker: tail of the install directory to look for, e.g. "mrparse/html".
            Ignored when install_dir is given.
        install_dir: the directory to copy, when the caller already knows it.
            Reports that reference their assets only by URL carry no absolute
            path to discover, so the source has to be named — xia2's report
            loads its libraries from CDNs, but DIALS ships those same
            libraries, at those same versions, in dials/static.
        extra_assets: URLs appearing in the report mapped to local files to
            serve in their place — CDN references, which are localised the
            same way and for the same reasons as the install-directory ones.
            Each file is copied to the top of ``subdirectory``.
        url_rewrites: URLs mapped to paths *within* the copied directory, for
            assets that are already there. Preferred over extra_assets when
            an asset resolves siblings relatively — KaTeX's stylesheet asks
            for ``fonts/…``, so flattening it to the top would strand its
            fonts.
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
    if install_dir is None and marker is not None:
        install_dir = discover_install_dir(html_text, marker)
    if install_dir is None or not Path(install_dir).is_dir():
        return output_path if output_path.exists() else None
    install_dir = Path(install_dir)

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

    for url, relative in (url_rewrites or {}).items():
        if not (destination / relative).is_file():
            logger.warning(
                "%s is not in %s; leaving %s alone", relative, destination, url
            )
            continue
        localised = localised.replace(url, f"{subdirectory}/{relative}")

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


# The libraries a DIALS-templated report loads from CDNs, mapped to where the
# same file sits inside ``dials/static``. DIALS ships all of them, at these
# exact versions, because its own report templates offer a local-dependency
# mode (``report_local_dep.html``) — but xia2 renders from ``report_base.html``,
# which hardcodes the CDN block with no way to switch. Taking them from the
# DIALS install rather than vendoring copies here keeps them matched to the
# DIALS that produced the report, and brings KaTeX's 60 font files along.
DIALS_CDN_ASSETS = {
    "https://code.jquery.com/jquery-1.12.0.min.js": "js/jquery-1.12.0.min.js",
    "https://cdn.plot.ly/plotly-latest.min.js": "js/plotly-latest.min.js",
    "https://maxcdn.bootstrapcdn.com/bootstrap/3.3.5/js/bootstrap.min.js": "js/bootstrap.min.js",
    "https://maxcdn.bootstrapcdn.com/bootstrap/3.3.5/css/bootstrap.min.css": "css/bootstrap.min.css",
    "https://cdnjs.cloudflare.com/ajax/libs/KaTeX/0.5.1/katex.min.js": "katex/katex.min.js",
    "https://cdnjs.cloudflare.com/ajax/libs/KaTeX/0.5.1/contrib/auto-render.min.js": "katex/contrib/auto-render.min.js",
    "https://cdnjs.cloudflare.com/ajax/libs/KaTeX/0.5.1/katex.min.css": "katex/katex.min.css",
}


def dials_static_dir() -> Optional[Path]:
    """Where the installed DIALS keeps the libraries its reports use."""
    try:
        import dials
    except ImportError:
        logger.debug("DIALS is not importable; cannot localise its report assets")
        return None
    static = Path(dials.__file__).parent / "static"
    return static if static.is_dir() else None


def localise_dials_report(
    html_path: Path, subdirectory: str = "dials_static"
) -> Optional[Path]:
    """Make a DIALS-templated HTML report load without reaching the network.

    A report that fetches jQuery, Plotly, Bootstrap and KaTeX from CDNs is
    refused by the app's content security policy, so it arrives unstyled and
    with no plots at all. Beyond the policy, a job report is a record: it
    should still render on an offline beamline machine, and years after those
    CDNs stopped serving these versions.

    Copies ``dials/static`` in beside the report (~3.7 MB, of which KaTeX's
    fonts are 2.2 MB) and points the report at it. Returns the rewritten
    file's path, or None if DIALS is not importable.
    """
    static = dials_static_dir()
    if static is None:
        return None
    return localise_report_assets(
        html_path,
        install_dir=static,
        subdirectory=subdirectory,
        url_rewrites=DIALS_CDN_ASSETS,
    )
