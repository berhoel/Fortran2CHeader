"""Configuration file for the Sphinx documentation builder."""

# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
from __future__ import annotations

project = "Fortran2CHeader"
copyright = (  # noqa:A001
    "1999, 2000, 2020, 2023—2023 Berthold Höllmann <berhoel@gmail.com>"
)

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# -- Options for HTML output -------------------------------------------------

from berhoel.sphinx_settings import (  # noqa:E402 ; isort:skip
    ProjectTypes,
    setup,
    favicons,
    language,
    extensions,
    latex_engine,
    configuration,
    latex_elements,
    latex_show_urls,
    exclude_patterns,
    html_static_path,
    napoleon_use_ivar,
    napoleon_use_param,
    napoleon_use_rtype,
    sitemap_url_scheme,
    todo_include_todos,
    intersphinx_mapping,
    autodoc_default_options,
    napoleon_numpy_docstring,
    napoleon_google_docstring,
    napoleon_include_init_with_doc,
    napoleon_include_private_with_doc,
    napoleon_include_special_with_doc,
    napoleon_use_admonition_for_notes,
    napoleon_use_admonition_for_examples,
    napoleon_use_admonition_for_references,
)
from berhoel import sphinx_settings  # noqa:E402 ; isort:skip

globals().update(configuration().configuration())

html_theme = "berhoel_sphinx_theme"
html_theme_path = sphinx_settings.get_html_theme_path()

__all__ = (
    "ProjectTypes",
    "autodoc_default_options",
    "configuration",
    "exclude_patterns",
    "extensions",
    "favicons",
    "html_static_path",
    "intersphinx_mapping",
    "language",
    "latex_elements",
    "latex_engine",
    "latex_show_urls",
    "napoleon_google_docstring",
    "napoleon_include_init_with_doc",
    "napoleon_include_private_with_doc",
    "napoleon_include_special_with_doc",
    "napoleon_numpy_docstring",
    "napoleon_use_admonition_for_examples",
    "napoleon_use_admonition_for_notes",
    "napoleon_use_admonition_for_references",
    "napoleon_use_ivar",
    "napoleon_use_param",
    "napoleon_use_rtype",
    "setup",
    "sitemap_url_scheme",
    "todo_include_todos",
)

extensions.extend(
    [
        "sphinx_argparse_cli",
    ],
)

# (Optional) Logo. Should be small enough to fit the navbar (ideally 24x24).
# Path should be relative to the ``_static`` files directory.
# t m l_logo = "_static/Fortran2CHeader_logo.svg"

html_theme_options = {
    "navbar_links": [("GitLab", "https://github.com/berhoel/Fortran2CHeader", True)]
}
