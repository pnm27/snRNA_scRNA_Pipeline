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
import os
import sys

sys.path.insert(0, os.path.abspath('../../helper_py_scripts/'))
print(sys.path)

# -- Project information -----------------------------------------------------

project = 'snRNAseq scRNAseq Pipeline'
copyright = '2022, Prashant N M'
author = 'Prashant N M'

# The full version, including alpha/beta/rc tags
release = '0.1'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "myst_parser",
    "sphinx.ext.duration",
    "sphinx.ext.autosectionlabel",
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.todo",
    "sphinxarg.ext",
    "sphinxcontrib.mermaid",
    "hoverxref.extension",
    "sphinxcontrib.bibtex",
    'sphinx.ext.autosummary',
    # "numpydoc",
    # "sphinx_mdinclude",
]

# For todo extension
todo_include_todos = True
source_suffix = [
    ".rst",
    ".md",
]

# For having myst to generate anchors until ### (h3 level) headings
myst_heading_anchors = 3

# Myst extend the figure directive
myst_enable_extensions = [
    "colon_fence",
    "tasklist",
]

# To make available sections that don't have a unique name as a hyperlink target
autosectionlabel_prefix_document = True

# Napoleon settings
napoleon_google_docstring = True
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = False
napoleon_include_private_with_doc = False
napoleon_include_special_with_doc = True
napoleon_use_admonition_for_examples = False
napoleon_use_admonition_for_notes = False
napoleon_use_admonition_for_references = False
napoleon_use_ivar = False
napoleon_use_param = True
napoleon_use_rtype = True
napoleon_preprocess_types = False
napoleon_type_aliases = None
napoleon_attr_annotations = True

# For hover tooltips
hoverxref_roles = [
    'term',
]

# To hoverxref role types
hoverxref_role_types = {
    'hoverxref': 'tooltip',
    'p': 'modal',
    'term': 'tooltip'
}

# For BibTeX
bibtex_bibfiles = ['bibliography.bib']

# For hover on BibTeX
hoverxref_domains =[
    'cite',
]

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
html_theme = 'sphinx_rtd_theme'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = []

# These paths are either relative to html_static_path
# or fully qualified paths (eg. https://...)
html_css_files = [
    'css/tables.css',
    'css/page_style.css',
]