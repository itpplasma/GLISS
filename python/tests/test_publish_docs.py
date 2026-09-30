"""Versioned documentation site layout written by ci/publish_docs.py."""

import importlib.util
import json
from pathlib import Path

import pytest

_SCRIPT = Path(__file__).resolve().parents[2] / "ci" / "publish_docs.py"


def _publisher():
    specification = importlib.util.spec_from_file_location("publish_docs", _SCRIPT)
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module.publish


def _html(root: Path, text: str) -> Path:
    root.mkdir(parents=True)
    (root / "index.html").write_text(text, encoding="utf-8")
    return root


def test_publish_orders_versions_and_redirects_to_latest(tmp_path):
    publish = _publisher()
    site = tmp_path / "site"
    assert publish(_html(tmp_path / "a", "0.0.2"), "v0.0.2", site) == ["v0.0.2"]
    assert (site / "index.html").read_text().count("v0.0.2/index.html") == 2
    publish(_html(tmp_path / "b", "0.0.10"), "v0.0.10", site)
    versions = publish(_html(tmp_path / "c", "main"), "latest", site)
    assert versions == ["latest", "v0.0.10", "v0.0.2"]
    assert json.loads((site / "versions.json").read_text()) == versions
    assert "latest/index.html" in (site / "index.html").read_text()
    assert (site / ".nojekyll").is_file()
    publish(_html(tmp_path / "d", "main 2"), "latest", site)
    assert (site / "latest" / "index.html").read_text() == "main 2"
    assert (site / "v0.0.2" / "index.html").read_text() == "0.0.2"


def test_publish_refuses_to_change_a_release(tmp_path):
    publish = _publisher()
    site = tmp_path / "site"
    publish(_html(tmp_path / "a", "first"), "v1.2.3", site)
    with pytest.raises(ValueError, match="already published"):
        publish(_html(tmp_path / "b", "second"), "v1.2.3", site)
    assert (site / "v1.2.3" / "index.html").read_text() == "first"
    with pytest.raises(ValueError, match="latest"):
        publish(_html(tmp_path / "c", "x"), "main", site)
    with pytest.raises(ValueError, match="built HTML"):
        publish(tmp_path / "missing", "latest", site)
