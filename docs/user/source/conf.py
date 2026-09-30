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
