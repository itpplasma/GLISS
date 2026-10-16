"""Build historical documentation from a released checkout and its installed wheel.

The first two releases predate docs/. Their original Python guides are
rendered as archives, with a separately identified, executable version and
Mercier example. Current API pages are never copied into an old release.
"""

import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile


def build(source, tag, commit, output, fixture, switcher):
    if not re.fullmatch(r"v\d+\.\d+\.\d+", tag):
        raise ValueError("a release tag vX.Y.Z is required")
    if not re.fullmatch(r"[0-9a-f]{40}", commit):
        raise ValueError("a full release commit is required")
    actual = subprocess.check_output(
        ["git", "-C", str(source), "rev-parse", "HEAD"], text=True
    ).strip()
    if actual != commit:
        raise ValueError(f"{tag} checkout is {actual}, expected {commit}")
    guide = source / "python" / "README.md"
    guide_text = guide.read_text(encoding="utf-8")

    import gliss

    package = Path(gliss.__file__).resolve().parent
    if not package.is_relative_to(Path(sys.prefix).resolve()):
        raise RuntimeError(f"release docs imported outside their environment: {package}")
    library = Path(gliss._load_library()._name).resolve()
    if not library.is_relative_to(package):
        raise RuntimeError(f"release docs loaded an external native library: {library}")
    version = tag[1:]
    if gliss.__version__ != version or gliss.version() != version:
        raise RuntimeError(f"{tag} does not match the installed Python/native package")

    environment = os.environ.copy()
    environment.pop("GLISS_LIB", None)
    environment.pop("PYTHONPATH", None)
    output.mkdir(parents=True, exist_ok=True)
    provenance = {
        "tag": tag,
        "git_commit": commit,
        "package_version": version,
        "guide_sha256": hashlib.sha256(guide.read_bytes()).hexdigest(),
        "fixture_sha256": hashlib.sha256(fixture.read_bytes()).hexdigest(),
    }
    with tempfile.TemporaryDirectory(prefix=f"gliss-docs-{tag}-") as temporary:
        directory = Path(temporary)
        docs = directory / "source"
        docs.mkdir()
        static = docs / "_static"
        static.mkdir()
        shutil.copy2(switcher, static / "versions.js")
        shutil.copy2(fixture, static / "solovev_q1.035.nc")
        (docs / "guide.md").write_text(guide_text, encoding="utf-8")
        (docs / "conf.py").write_text(
            'project = "GLISS"\n'
            f"release = version = {version!r}\n"
            f"html_title = {'GLISS ' + version + ' archive'!r}\n"
            'extensions = ["myst_parser", "sphinx.ext.doctest"]\n'
            'source_suffix = {".md": "markdown"}\n'
            'html_static_path = ["_static"]\n'
            'html_extra_path = ["source.json"]\n'
            'html_js_files = ["versions.js"]\n', encoding="utf-8",
        )
        (docs / "index.md").write_text(
            f"# GLISS {version}\n\n"
            f"This archive preserves the Python guide shipped in tag `{tag}` "
            f"at commit `{commit}`. Its APIs and conventions describe that "
            "historical package; subsequent releases contain corrections. "
            "For current work, use the "
            "[latest documentation](https://itpplasma.github.io/GLISS/latest/).\n\n"
            f"[Tagged source](https://github.com/itpplasma/GLISS/tree/{commit}) "
            "and [archive provenance](source.json).\n\n"
            "```{toctree}\n:maxdepth: 2\n\nguide\n```\n\n"
            "## Executed release example\n\n"
            "CI builds and installs this tag's wheel, checks its Python and native "
            "versions, and runs the guide's Mercier call with the "
            "[public Solov'ev input](_static/solovev_q1.035.nc). "
            "The original guide remains unchanged; its other examples are "
            "historical user-input workflows and require the stated input files.\n\n"
            "```{doctest}\n"
            ">>> import gliss\n"
            f">>> gliss.__version__ == gliss.version() == {version!r}\nTrue\n"
            ">>> import numpy as np\n"
            ">>> from pathlib import Path\n"
            ">>> equilibrium_export = Path('_static/solovev_q1.035.nc')\n"
            ">>> s, d_mercier = gliss.mercier_profile(equilibrium_export)\n"
            ">>> s.ndim == d_mercier.ndim == 1 and s.size == d_mercier.size > 0\nTrue\n"
            ">>> bool(np.all(np.isfinite(s)) and np.all(np.isfinite(d_mercier)))\nTrue\n"
            "```\n", encoding="utf-8",
        )
        (docs / "source.json").write_text(
            json.dumps(provenance, indent=2) + "\n", encoding="utf-8"
        )
        for builder in ("doctest", "html", "linkcheck"):
            destination = output if builder == "html" else directory / builder
            subprocess.run(
                [sys.executable, "-m", "sphinx", "-W", "-b", builder,
                 str(docs), str(destination)],
                cwd=docs, env=environment, check=True,
            )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--tag", required=True)
    parser.add_argument("--commit", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--fixture", type=Path, required=True)
    parser.add_argument("--switcher", type=Path, required=True)
    arguments = parser.parse_args()
    build(arguments.source.resolve(), arguments.tag, arguments.commit,
          arguments.output.resolve(), arguments.fixture.resolve(),
          arguments.switcher.resolve())
