"""Sphinx configuration for the GLISS documentation."""

import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join("..", "python")))

import gliss  # noqa: E402

project = "GLISS"
author = "GLISS contributors"
release = gliss.__version__
version = release
extensions = ["myst_parser", "sphinx.ext.autodoc", "sphinx.ext.napoleon"]
source_suffix = {".md": "markdown"}
exclude_patterns = ["_build"]
autodoc_member_order = "bysource"
html_title = f"GLISS {release}"
