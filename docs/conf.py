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

project = 'simos'
copyright = '2025, Konstantin Herb and Laura Völker'
author = 'Konstantin Herb and Laura Völker'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
import sphinx_rtd_theme

extensions = ['sphinx.ext.autodoc','sphinx_rtd_theme','sphinx.ext.napoleon','sphinx.ext.todo','sphinx.ext.viewcode', 'sphinx.ext.autosectionlabel', 'sphinx_exec_code'
]



# Add any paths that contain templates here, relative to this directory.
#templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_rtd_theme'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

# -- Configure Python Autodoc -----------------------------------------------------
# Autodoc actually imports (i.e. runs) a module to discover the docstrings. The machine building the docs shouldn't have
# to build parament itself. So instead, we mock an installed parament by pointing to the source code instead.
import sys
import os
import numpy as np
print("KOKO:", os.path.abspath('..'))
sys.path.insert(0, os.path.abspath('..'))
print("KOKO:", os.path.abspath('../simos'))
sys.path.insert(0, os.path.abspath('../simos'))


# If true, the current module name will be prepended to all description
# unit titles (such as .. function::).
add_module_names = False


napoleon_google_docstring = True
todo_include_todos = True


autodoc_member_order = 'alphabetical'