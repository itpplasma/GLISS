import ast
import json
import re
from pathlib import Path

from gliss._closure_registry import CATEGORIES, CLOSURES, PROFILES, manifest_document


ROOT = Path(__file__).parents[2]
FROZEN_MANIFEST = ROOT / "python" / "gliss" / "closure_manifest.json"


def test_frozen_manifest_matches_registry():
    expected = json.dumps(manifest_document(), indent=2, sort_keys=True) + "\n"
    assert FROZEN_MANIFEST.read_text(encoding="utf-8") == expected


def _registered_ctest(path):
    # The evidence file must build a test executable that ctest runs.
    build = (ROOT / "CMakeLists.txt").read_text(encoding="utf-8")
    target = re.search(
        r"add_executable\(\s*(\S+)\s+" + re.escape(path) + r"\s*\)", build
    )
    assert target, f"{path} is not built as a test executable"
    command = r"add_test\(\s*NAME\s+\S+\s+COMMAND\s+" + re.escape(target[1])
    assert re.search(command + r"\b", build), f"{target[1]} is not run by ctest"


def _fortran_check(source, token):
    # A token names a check that runs: a procedure the test calls, or the
    # complete failure message of one of its assertions. A word that merely
    # occurs in the file is not evidence.
    called = re.search(r"\bcall\s+" + re.escape(token) + r"\s*\(", source)
    message = f'"{token}"' in source or f"'{token}'" in source
    return bool(called) or message


def _python_test(path, token):
    # A token names a collected test function. Tests of native
    # implementations must live in a native module, which loads the real
    # library instead of a fake.
    module = ast.parse((ROOT / path).read_text(encoding="utf-8"))
    functions = {
        node.name: node for node in module.body if isinstance(node, ast.FunctionDef)
    }
    assert token.startswith("test_") and token in functions, (path, token)
    arguments = {argument.arg for argument in functions[token].args.args}
    assert not arguments & {"monkeypatch", "contexts"}, (path, token)
    native = any(
        isinstance(node, ast.Assign)
        and any(getattr(target, "id", "") == "pytestmark" for target in node.targets)
        and "native" in ast.unparse(node.value)
        for node in module.body
    )
    return native


def test_every_closure_has_behavioral_evidence():
    for closure in CLOSURES:
        assert closure.category in CATEGORIES
        assert closure.implementations
        assert closure.evidence
        for relative in closure.implementations:
            assert (ROOT / relative).is_file(), closure.identifier
        # Implementations in Fortran, or Python that calls into the library,
        # need native evidence.
        uses_native = any(
            not relative.endswith(".py")
            or "_require_symbols" in (ROOT / relative).read_text(encoding="utf-8")
            for relative in closure.implementations
        )
        for item in closure.evidence:
            if item.path.endswith(".f90"):
                _registered_ctest(item.path)
                source = (ROOT / item.path).read_text(encoding="utf-8")
                assert _fortran_check(source, item.token), (closure.identifier, item)
            else:
                assert item.path.startswith("python/tests/"), item
                native = _python_test(item.path, item.token)
                assert native or not uses_native, (closure.identifier, item)


def test_profiles_cover_registry_without_unknown_edges():
    closures = {closure.identifier: closure for closure in CLOSURES}
    assert len(closures) == len(CLOSURES)
    assert len({profile.identifier for profile in PROFILES}) == len(PROFILES)
    selected = set()
    for profile in PROFILES:
        assert profile.interface in {"native", "python"}
        assert len(profile.selection) == len(profile.closures)
        assert set(profile.closures).issubset(CATEGORIES)
        assert set(CATEGORIES) - {"derivative"} <= set(profile.closures)
        assert set(profile.closures.values()).issubset(closures)
        for category, identifier in profile.closures.items():
            assert closures[identifier].category == category
        source = (ROOT / profile.entrypoint.path).read_text(encoding="utf-8")
        assert profile.entrypoint.token in source
        selected.update(profile.closures.values())
    assert selected == set(closures)
