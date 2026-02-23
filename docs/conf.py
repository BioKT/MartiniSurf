from __future__ import annotations

import os
import sys
from datetime import datetime

# Allow autodoc/imports from project root.
sys.path.insert(0, os.path.abspath(".."))

project = "MartiniSurf"
author = "MartiniSurf Contributors"
copyright = f"{datetime.now():%Y}, {author}"
release = "1.0.0"

extensions = [
    "myst_parser",
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
]

source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
}

master_doc = "index"

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]
html_title = "MartiniSurf Documentation"
html_logo = "../logo.png"

myst_enable_extensions = [
    "colon_fence",
    "deflist",
]

# Keep docs clean on RTD when markdown headings are reused.
suppress_warnings = ["myst.header"]
