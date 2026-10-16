"""Place a built HTML documentation tree into the versioned site layout.

usage: python ci/publish_docs.py HTML_DIRECTORY VERSION SITE_DIRECTORY

The site has one directory per version: ``latest`` for the default branch and
``vX.Y.Z`` for release tags. The version directory is replaced, the other
versions are kept, ``versions.json`` lists them (latest first, then releases
newest first), and ``index.html`` redirects to ``latest`` or, without it, to
the newest release. A release directory is written once: republishing an
existing tag is refused so that published release documentation never
changes.
"""

import argparse
import json
import re
import shutil
from pathlib import Path

_RELEASE = re.compile(r"^v(\d+)\.(\d+)\.(\d+)$")


def _ordered(names):
    releases = sorted(
        (name for name in names if _RELEASE.match(name)),
        key=lambda name: tuple(int(part) for part in _RELEASE.match(name).groups()),
        reverse=True,
    )
    return (["latest"] if "latest" in names else []) + releases


def publish(html: Path, version: str, site: Path) -> list:
    if version != "latest" and not _RELEASE.match(version):
        raise ValueError(f"version must be 'latest' or vX.Y.Z, got {version!r}")
    if not (html / "index.html").is_file():
        raise ValueError(f"{html} is not a built HTML tree")
    site.mkdir(parents=True, exist_ok=True)
    target = site / version
    if target.exists():
        if version != "latest":
            raise ValueError(f"release documentation {version} is already published")
        shutil.rmtree(target)
    shutil.copytree(html, target)
    names = [item.name for item in site.iterdir() if item.is_dir()]
    versions = _ordered(names)
    (site / "versions.json").write_text(json.dumps(versions) + "\n", encoding="utf-8")
    home = versions[0]
    (site / "index.html").write_text(
        "<!DOCTYPE html>\n<meta charset=\"utf-8\">\n<title>GLISS documentation</title>\n"
        f"<meta http-equiv=\"refresh\" content=\"0; url={home}/index.html\">\n"
        f"<a href=\"{home}/index.html\">GLISS documentation</a>\n",
        encoding="utf-8",
    )
    (site / ".nojekyll").write_text("", encoding="utf-8")
    return versions


def publish_missing_releases(archives: Path, site: Path) -> None:
    """Add tag-built archives once, rejecting a different source for that tag."""
    releases = _ordered([item.name for item in archives.iterdir() if item.is_dir()])
    for version in releases:
        if version == "latest":
            raise ValueError("release archives cannot replace latest")
        archive = archives / version
        provenance = json.loads((archive / "source.json").read_text(encoding="utf-8"))
        if (provenance.get("tag") != version
                or provenance.get("package_version") != version[1:]
                or not re.fullmatch(r"[0-9a-f]{40}", provenance.get("git_commit", ""))):
            raise ValueError(f"invalid release provenance for {version}")
        target = site / version
        if target.exists():
            previous = json.loads((target / "source.json").read_text(encoding="utf-8"))
            for key in ("tag", "git_commit", "package_version", "guide_sha256"):
                if previous.get(key) != provenance.get(key):
                    raise ValueError(f"published release {version} has different {key}")
        else:
            publish(archive, version, site)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("html", type=Path)
    parser.add_argument("version")
    parser.add_argument("site", type=Path)
    parser.add_argument("--releases", type=Path)
    arguments = parser.parse_args()
    if arguments.releases is not None:
        publish_missing_releases(arguments.releases, arguments.site)
    print(publish(arguments.html, arguments.version, arguments.site))
