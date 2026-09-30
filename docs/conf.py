"""Sphinx configuration for the GLISS documentation."""

import os
import sys

# The published site documents the installed package; a local build from a
# checkout without an installed gliss falls back to the source tree.
try:
    import gliss  # noqa: E402
except ImportError:
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
html_static_path = ["_static"]
# versions.js reads ../versions.json, written by ci/publish_docs.py next to
# the version directories of the published site, and adds a version menu.
html_js_files = ["versions.js"]
