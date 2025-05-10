# Configuration file for the Sphinx documentation builder.

# -- Project information -----------------------------------------------------

project = 'Benchmarking CCI'
copyright = '2025, Li-Ting Ku'
author = 'Li-Ting Ku'

# The full version, including alpha/beta/rc tags
release = '1.0.0'

# -- General configuration

extensions = [
    "nbsphinx",
    'sphinx.ext.duration',
    'sphinx.ext.doctest',
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.intersphinx',
]

intersphinx_mapping = {
    'python': ('https://docs.python.org/3/', None),
    'sphinx': ('https://www.sphinx-doc.org/en/master/', None),
}
intersphinx_disabled_domains = ['std']

templates_path = ['_templates']

# -- Options for HTML output

html_theme = 'sphinx_rtd_theme'
html_theme_options = {
    'style_nav_header_background': '#489fb5'
}
html_static_path = ['_static'] 
html_css_files = [
    'custom.css'
]
html_context = {
    "display_github": True,
    "github_user": "LitingKu",
    "github_repo": "Benchmarking-CCI",
    "github_version": "main",
    "conf_py_path": "/docs/source/",
}
# -- Options for EPUB output
epub_show_urls = 'footnote'

nbsphinx_allow_errors = True
nbsphinx_execute = 'never'
