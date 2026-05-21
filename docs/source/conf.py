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
import os, subprocess
import sys

sys.path.insert(0, os.path.abspath('../../helper_py_scripts/'))
# HERE = Path(__file__).parent
# sys.path[:0] = [str(HERE.parent), str(HERE / '../../helper_py_scripts/')]
# print(sys.path)
# print(os.listdir(sys.path[0]))

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
    # "myst_parser",
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
    "sphinxcontrib.jquery",
    "sphinx_design",
    "sphinx_autodoc_typehints",
    "sphinx_copybutton",
    # "sphinx_autoapi",
    "autoapi.extension",
    "myst_nb",
    # "numpydoc",
    # "sphinx_mdinclude",
]

# To prevent section label (autosectionlabel) collisions
autosectionlabel_prefix_document = True
autosectionlabel_maxdepth = 3 # max depth for labelling

# For todo extension
todo_include_todos = True
source_suffix = {
    '.rst': 'restructuredtext',
    '.txt': 'restructuredtext',
    '.md': 'myst-nb',
    '.ipynb': 'myst-nb',
    '.myst': 'myst-nb',
}

# For having myst to generate anchors until ### (h3 level) headings
myst_heading_anchors = 3

# Myst extend the figure directive
myst_enable_extensions = [
    "colon_fence",
    "tasklist",
]

# For mermaid runs
mermaid_params = ['-p' 'puppeteer-config.json']

# To make available sections that don't have a unique name as a hyperlink target
autosectionlabel_prefix_document = True

# Napoleon settings
napoleon_google_docstring = False
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = False
napoleon_include_private_with_doc = False
napoleon_include_special_with_doc = True
napoleon_use_admonition_for_examples = False
napoleon_use_admonition_for_notes = False
napoleon_use_admonition_for_references = False
napoleon_use_ivar = False
napoleon_use_param = True
napoleon_use_rtype = False # sphinx.ext.napoleon setting
napoleon_preprocess_types = False
napoleon_type_aliases = None
napoleon_attr_annotations = True

# reads Python type hints and automatically adds type information
typehints_use_rtype = False

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
exclude_patterns = [
    'build',
    ]

# For autodoc setup to not show on toctree
autosummary_generate = False
# Turn off prepending module names
add_module_names = False
# Sort members by type
autodoc_member_order = 'groupwise'
autodoc_mock_imports = [
    "scanpy",
    "anndata",
    "solo-sc",
    "pegasuspy",
    "matplotlib",
]

# moves types into the description section
autodoc_typehints = "description"

# AutoAPI search dir for source
autoapi_dirs = ['../../helper_py_scripts/']

# DEBUG - Create files in docs/source/autoapi/
autoapi_keep_files = True 

# Skip running API
# def run_apidoc(_):
#     subprocess.call([
#         "sphinx-apidoc",
#         "-e", # each script in separate page
#         "-f", # force build - removes lingering older docs
#         "-o",
#         "api",
#         "../../helper_py_scripts",
#         "../../helper_py_scripts/__pycache__", # exclude
#     ])

# To copy all lines in a multiline command prompt
copybutton_line_continuation_character = "\\"
# To copy HERE-doc
copybutton_here_doc_delimiter = "EOF"

# Document __init__, __repr__, and __str__ methods
# def skip(app, what, name, obj, would_skip, options):
#     if name in ("__init__", "__repr__", "__str__"):
#         return False
#     return would_skip
def setup(app):
    # app.connect("builder-inited", run_apidoc)
    print(app)
#     app.connect("autodoc-skip-member", skip)

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'furo' # 'sphinx_rtd_theme'

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