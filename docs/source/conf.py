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
import sphinx_rtd_theme

sys.path.insert(0, os.path.abspath('../../pycasino'))

# -- Project information -----------------------------------------------------

project = 'pycasino'
copyright = '2024, Salvatore De Angelis'
author = 'Salvatore De Angelis'

# The full version, including alpha/beta/rc tags
release = '04/06/2024'


# -- General configuration ---------------------------------------------------

# We require sphinx >=2 because of sphinxcontrib.bibtex,
needs_sphinx = '2.0'

extensions = [
    'IPython.sphinxext.ipython_console_highlighting',
    'IPython.sphinxext.ipython_directive',
    'sphinx.ext.autodoc',
    'sphinx.ext.mathjax',
    'sphinx.ext.todo',
    'sphinx.ext.autosummary',
    'sphinx.ext.napoleon',
    'sphinx.ext.viewcode',
]

exclude_patterns = ['_build', '**.ipynb_checkpoints']

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = 'sphinx'

# -- Options for bibtex


# -- Options for Napoleon -----------------------------------------------------

napoleon_google_docstring = True
napoleon_numpy_docstring = True
napoleon_include_private_with_doc = False
napoleon_include_special_with_doc = False
napoleon_use_admonition_for_examples = False
napoleon_use_admonition_for_notes = False
napoleon_use_admonition_for_references = False
napoleon_use_ivar = False
napoleon_use_param = False
napoleon_use_rtype = False

# -- Options for TODO ---------------------------------------------------------

todo_include_todos = True



# -- Options for HTML output ----------------------------------------------

html_theme = 'sphinx_rtd_theme'
html_theme_path = [sphinx_rtd_theme.get_html_theme_path()]

html_theme_options = {
    'style_nav_header_background': '#4f8fb8ff',
    'collapse_navigation': False,
    'logo_only': True,
}

# html_logo = 'img/tomopy-logo-wide-mono.svg'
# html_favicon = 'img/tomopy-logo.svg'


# -- Options for Texinfo output -------------------------------------------
# http://www.sphinx-doc.org/en/master/usage/extensions/autodoc.html#confval-autodoc_mock_imports

autodoc_mock_imports = [
    'concurrent',
    'DM3lib',
    'libtomo',
    'tomopy.util.extern',
    'matplotlib',
    'numexpr',
    'numpy',
    'pyfftw',
    'pywt',
    'scipy',
    'skimage',
    'tifffile',
]