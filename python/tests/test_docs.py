"""Execute every Python block of the quickstart and VMEC guide, in order."""

import os
import re
from pathlib import Path

import pytest

pytestmark = pytest.mark.native

_REPOSITORY = Path(__file__).resolve().parents[2]


@pytest.mark.parametrize(("name", "minimum_blocks"), [("quickstart", 5), ("vmec", 2)])
def test_documentation_examples_run(native_library, test_data, monkeypatch, name, minimum_blocks):
    if name == "vmec":
        pytest.importorskip("booz_xform", reason="VMEC examples require gliss[vmec]")
    docs = Path(os.environ.get("GLISS_DOCS", _REPOSITORY / "docs"))
    document = docs / f"{name}.md"
    if not document.is_file():
        pytest.fail(f"the documentation page is missing: {document}")
    blocks = re.findall(
        r"^```python\n(.*?)^```", document.read_text(encoding="utf-8"),
        re.MULTILINE | re.DOTALL,
    )
    assert len(blocks) >= minimum_blocks
    monkeypatch.setenv("GLISS_TEST_DATA", str(test_data))
    namespace = {"__name__": name}
    for index, block in enumerate(blocks):
        code = compile(block, f"{document}:block{index}", "exec")
        exec(code, namespace)
