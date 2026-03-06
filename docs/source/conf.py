import os
import sys

project = "wind3d_field_lines"
author = "wind3d_field_lines contributors"
copyright = "2026, wind3d_field_lines contributors"
release = "0.2.0"

extensions = [
    "sphinx.ext.autodoc",  # Generate API docs from docstrings
    "sphinx.ext.napoleon",  # Enable Google/NumPy style docstrings
    "sphinx.ext.viewcode",  # Add links to highlighted source code
    "sphinx.ext.doctest",  # Enable doctest snippets
]

templates_path = ["_templates"]
exclude_patterns = []

language = "en"
autodoc_typehints = "description"

html_theme = "bizstyle"
html_static_path = ["_static"]
html_context = {
    "display_github": True,
    "github_user": "iijimahr",
    "github_repo": "wind3d_field_lines",
    "github_version": "main",
    "conf_py_path": "/docs/source/",
}

sys.path.insert(0, os.path.abspath("../../src"))
