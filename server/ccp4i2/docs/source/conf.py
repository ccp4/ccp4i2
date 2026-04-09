# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
project = 'CCP4i2'
copyright = '2025, CCP4'
author = 'CCP4 Development Team'
release = '2.0'

# -- General configuration ---------------------------------------------------
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.napoleon',  # Support Google/NumPy style docstrings
    'sphinx.ext.viewcode',  # Add links to source code
    'sphinx.ext.intersphinx',  # Link to other projects' documentation
]

# Autosummary settings
autosummary_generate = True
autosummary_imported_members = True

# Napoleon settings for docstring parsing
napoleon_google_docstring = True
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = True
napoleon_include_private_with_doc = False
napoleon_include_special_with_doc = True
napoleon_use_admonition_for_examples = True
napoleon_use_admonition_for_notes = True
napoleon_use_admonition_for_references = True
napoleon_use_ivar = False
napoleon_use_param = True
napoleon_use_rtype = True
napoleon_type_aliases = None

# Autodoc settings
autodoc_default_options = {
    'members': True,
    'member-order': 'bysource',
    'special-members': '__init__',
    'undoc-members': True,
    'exclude-members': '__weakref__'
}
autodoc_typehints = 'description'
autodoc_class_signature = 'separated'

# Intersphinx mapping
intersphinx_mapping = {
    'python': ('https://docs.python.org/3', None),
}

templates_path = ['_templates']
exclude_patterns = []

# -- Options for HTML output -------------------------------------------------
html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']

# Theme options
html_theme_options = {
    'navigation_depth': 4,
    'collapse_navigation': False,
    'sticky_navigation': True,
    'includehidden': True,
    'titles_only': False
}

# -- Options for autodoc -----------------------------------------------------
# Mock imports that might not be available during doc build
autodoc_mock_imports = [
    'django',
    'PySide2',
    'gemmi',
    'mmdb2',
    'mmut',
    'ccp4mg',
]
