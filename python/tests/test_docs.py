"""Execute every Python block of the documentation quickstart, in order."""

import os
import re
from pathlib import Path

import pytest

pytestmark = pytest.mark.native

_REPOSITORY = Path(__file__).resolve().parents[2]


def test_quickstart_examples_run(native_library, test_data, monkeypatch):
    docs = Path(os.environ.get("GLISS_DOCS", _REPOSITORY / "docs"))
    quickstart = docs / "quickstart.md"
    if not quickstart.is_file():
        pytest.fail(f"the documentation quickstart is missing: {quickstart}")
    blocks = re.findall(
        r"^```python\n(.*?)^```", quickstart.read_text(encoding="utf-8"),
        re.MULTILINE | re.DOTALL,
    )
    assert len(blocks) >= 5
    monkeypatch.setenv("GLISS_TEST_DATA", str(test_data))
    namespace = {"__name__": "quickstart"}
    for index, block in enumerate(blocks):
        code = compile(block, f"{quickstart}:block{index}", "exec")
        exec(code, namespace)
