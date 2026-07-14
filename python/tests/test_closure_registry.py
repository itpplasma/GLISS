import json
from pathlib import Path

from gliss._closure_registry import CATEGORIES, CLOSURES, PROFILES, manifest_document


ROOT = Path(__file__).parents[2]
FROZEN_MANIFEST = ROOT / "python" / "gliss" / "closure_manifest.json"


def test_frozen_manifest_matches_registry():
    expected = json.dumps(manifest_document(), indent=2, sort_keys=True) + "\n"
    assert FROZEN_MANIFEST.read_text(encoding="utf-8") == expected


def test_every_closure_has_behavioral_evidence():
    for closure in CLOSURES:
        assert closure.category in CATEGORIES
        assert closure.implementations
        assert closure.evidence
        for relative in closure.implementations:
            assert (ROOT / relative).is_file(), closure.identifier
        for item in closure.evidence:
            source = (ROOT / item.path).read_text(encoding="utf-8")
            assert item.token in source, (closure.identifier, item)


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
