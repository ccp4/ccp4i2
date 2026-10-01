# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
# import os
# import sys
# sys.path.insert(0, os.path.abspath('.'))


# -- Project information -----------------------------------------------------

project = 'CCP4i2'
copyright = '2020, CCP4'
author = 'CCP4'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.

extensions = [
    'sphinx_tabs.tabs',
    'sphinx.ext.doctest',
    'sphinx.ext.todo',
    'sphinx.ext.coverage',
    'sphinx.ext.mathjax',
    'sphinxcontrib.contentui'
]

master_doc = 'index'

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#

# MN edited to use the 'read the docs' theme
#import sphinx_rtd_theme
import sphinx_material

#extensions.append("sphinx_rtd_theme")
extensions.append("sphinx_material")

#html_theme = "sphinx_rtd_theme"
#html_theme = 'alabaster'
html_theme = "sphinx_material"

html_theme_options = {
'repo_url': "https://gitlab.com/ccp4i2/rstdocs",
'repo name': "rstdocs",
'repo_type': 'gitlab',
'theme_color': '#00020a',
'color_primary': '#005BBB',
'color_accent': 'yellow',
'logo_icon': '',
'globaltoc_depth': 6
}

html_context = {
#    "display_gitlab": True, # Integrate Gitlab
#    "gitlab_user": "ccp4i2", # Username
#    "gitlab_repo": "rstdocs", # Repo name
#    "gitlab_version": "master", # Version
    "conf_py_path": "/source/", # Path in the checkout to the docs root
}

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

html_css_files = ['https://github.com/bwithd/sphinx-material/blob/main/sphinx_material/sphinx_material/static/stylesheets/application.css']

# The documentation status page, regenerated from the task chooser on every
# build so it cannot go stale (tools/status.py).
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent / "tools"))
import status as _status
_status.write_status(_Path(__file__).resolve().parent / "status.rst")


# ---- Published help: draft banners, and one URL per task -------------------
# The app's Help button opens tasks/<task name>/index.html. A task's page can
# live elsewhere (status.ALIASES, a page shared by several tasks, a document
# not called index), and some tasks have none yet: build_finished writes a
# redirect at that URL for each, so the app needs no map of its own.
import json as _json

_TASK_PAGES = _Path(__file__).resolve().parent / "tasks"

DRAFT_BANNER = """
.. admonition:: Draft

   This page was written for the new interface and has not yet been
   reviewed by a crystallographer who knows the program. The pictures and
   numbers come from a real run; the advice may still change.

"""


def _is_draft(docname):
    parts = docname.split("/")
    if len(parts) < 3 or parts[0] != "tasks":
        return False
    shots = _TASK_PAGES / parts[1] / "shots.json"
    try:
        return bool(_json.loads(shots.read_text(encoding="utf-8")).get("draft"))
    except (OSError, ValueError):
        return False


def _draft_banner(app, docname, source):
    """Put the banner under the page's title, before its first paragraph."""
    if not _is_draft(docname):
        return
    lines = source[0].split("\n")
    # The title is the first line of text followed by an underline (an
    # overline before it is skipped too).
    for i in range(len(lines) - 1):
        under = lines[i + 1].strip()
        if lines[i].strip() and under and set(under) <= set("#=*-~^\"'`") and len(under) >= 3 \
                and not set(lines[i].strip()) <= set("#=*-~^\"'`"):
            source[0] = "\n".join(lines[:i + 2] + [DRAFT_BANNER] + lines[i + 2:])
            return


def _task_redirects(app, exception):
    if exception is not None or app.builder.format != "html":
        return
    out = _Path(app.outdir)
    seen = set()
    for _title, tasks in _status.chooser_categories():
        for task in tasks:
            if task in seen:
                continue
            seen.add(task)
            target = out / "tasks" / task / "index.html"
            state, doc = _status.state(task)
            if doc == f"{task}/index" and target.exists():
                continue
            if target.exists():
                continue  # a page of that name exists already
            to = f"../{doc}.html" if doc else "../index.html"
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_text(
                '<!doctype html><meta charset="utf-8">'
                f'<meta http-equiv="refresh" content="0; url={to}">'
                f'<link rel="canonical" href="{to}">'
                f'<title>{task}</title>'
                + (f'<p>The help for {task} is <a href="{to}">here</a>.</p>' if doc else
                   f'<p>There is no help page for {task} yet: see the '
                   f'<a href="{to}">task list</a>.</p>'),
                encoding="utf-8")


def setup(app):
    app.connect("source-read", _draft_banner)
    app.connect("build-finished", _task_redirects)
