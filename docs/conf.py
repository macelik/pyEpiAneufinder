from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.abspath(".."))

project = "pyEpiAneufinder"
author = "pyEpiAneufinder contributors"
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.autosummary",
]
autosummary_generate = True
templates_path = ["_templates"]
exclude_patterns = ["_build"]
html_theme = "furo"
html_static_path = ["_static"]
