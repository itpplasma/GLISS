#!/usr/bin/env python3
"""Record source coverage and static references; this is not an equation audit."""

import argparse
import ast
from collections import Counter, defaultdict
import hashlib
import json
from pathlib import Path
import re
import subprocess


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "source_inventory.jsonl"
EVIDENCE = ROOT / "ci/source_inventory_evidence.json"
EXTENSIONS = {".f90", ".py", ".sh", ".c", ".h", ".cmake", ".yml", ".js", ".lock"}
EXPLICIT = {"CMakeLists.txt", "fpm.toml", ".github/manylinux.Dockerfile",
            "docs/Makefile", "pyproject.toml", "python/gliss/closure_manifest.json",
            "ci/source_inventory_evidence.json"}
PREFIXES = ("src/", "app/", "include/", "python/gliss/", "python/tests/",
            "benchmarks/", "ci/", ".github/", "cmake/", "test/", "docs/")


def git(*args):
    return subprocess.check_output(["git", *args], cwd=ROOT, text=True).strip()


def paths():
    # Untracked build products are excluded by Git; the new inventory files are
    # included before their first commit. Binary fixtures and prose are omitted.
    names = git("ls-files", "--cached", "--others", "--exclude-standard").splitlines()
    return sorted(set(name for name in names if name in EXPLICIT or (
        name.startswith(PREFIXES) and Path(name).suffix in EXTENSIONS)))


def role(path):
    if path == "ci/source_inventory_evidence.json":
        return "audit_metadata"
    if path.startswith("src/"):
        return "c_abi_implementation" if path.rsplit("/", 1)[-1].startswith("gliss_") else "native_library"
    if path.startswith("include/"):
        return "c_abi_declaration"
    if path.startswith("python/gliss/"):
        return "python_interface"
    if path.startswith("app/"):
        return "application"
    if path.startswith("benchmarks/"):
        return "benchmark_runner"
    if path.startswith(("test/", "python/tests/")) or path == "cmake/test_gliss_header.c":
        return "verification_source"
    if path.startswith("ci/"):
        return "ci_verifier_or_runner"
    if path.startswith("docs/"):
        return "documentation_runtime"
    return "build_or_ci_configuration"


def module_name(path):
    parts = Path(path).with_suffix("").parts
    if len(parts) > 1 and parts[:2] == ("python", "gliss"):
        return ".".join(parts[1:-1] if parts[-1] == "__init__" else parts[1:])
    return None


def inspect_source(path, source):
    data = {"modules": [], "procedures": [], "c_abi_symbols": [],
            "c_binding_definitions": [],
            "imports": [], "parse_findings": []}
    if path.endswith(".f90"):
        data["modules"] = re.findall(r"^\s*module\s+(?!procedure\b|function\b|subroutine\b)(\w+)",
                                     source, re.M | re.I)
        interface_depth = 0
        interface_lines = set()
        for number, line in enumerate(source.splitlines(), 1):
            if line.lstrip().startswith("!"):
                continue
            if re.match(r"\s*(abstract\s+)?interface\b", line, re.I):
                interface_depth += 1
            elif re.match(r"\s*end\s*interface\b", line, re.I):
                interface_depth = max(0, interface_depth - 1)
            if interface_depth:
                interface_lines.add(number)
            match = re.match(r"\s*use\b\s*(?:,\s*(?:non_)?intrinsic\s*)?(?:::)?\s*(\w+)", line, re.I)
            if match:
                data["imports"].append({"name": match[1].lower(), "line": number, "kind": "fortran_use"})
            if not re.match(r"\s*(end|module|interface)\b", line, re.I):
                match = re.search(r"\b(subroutine|function)\s+(\w+)\s*\(", line, re.I)
                if match:
                    kind = "interface_declaration" if interface_depth else "procedure_definition_candidate"
                    data["procedures"].append({"name": match[2].lower(), "line": number, "kind": kind})
                    if interface_depth:
                        data["imports"].append({"name": match[2].lower(), "line": number,
                                                "kind": "fortran_interface_declaration"})
        folded = re.sub(r"&[ \t]*\n[ \t]*&?", "\n", source)
        bindings = list(re.finditer(r"bind\s*\(\s*c\s*,\s*name\s*=\s*['\"](\w+)['\"]", folded, re.I))
        for match in bindings:
            number = folded[:match.start()].count("\n") + 1
            if number in interface_lines:
                data["imports"].append({"name": match[1], "line": number,
                                        "kind": "c_external_binding_declaration"})
            else:
                data["c_binding_definitions"].append({"name": match[1], "line": number})
        data["c_abi_symbols"] = sorted(set(binding["name"] for binding in data["c_binding_definitions"]
                                           if binding["name"].startswith("gliss_")))
    elif path.endswith(".py"):
        try:
            tree = ast.parse(source, filename=path)
        except SyntaxError as error:
            data["parse_findings"].append(str(error))
            return data
        own_module = module_name(path)
        package = own_module if path.endswith("/__init__.py") else (
            own_module.rsplit(".", 1)[0] if own_module and "." in own_module else "")
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                data["imports"].extend({"name": alias.name, "line": node.lineno,
                                        "kind": "python_import"} for alias in node.names)
            elif isinstance(node, ast.ImportFrom):
                name = node.module or ""
                if node.level:
                    stem = package.split(".")[:len(package.split(".")) - node.level + 1]
                    name = ".".join(stem + ([name] if name else []))
                data["imports"].append({"name": name, "line": node.lineno,
                                        "kind": "python_from"})
                # from . import child can import a submodule or a package attribute.
                data["imports"].extend({"name": name + "." + alias.name,
                                        "line": node.lineno, "kind": "python_from_candidate"}
                                       for alias in node.names if alias.name != "*")
            elif isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
                data["procedures"].append({"name": node.name, "line": node.lineno})
    elif path.endswith((".c", ".h")):
        if path.endswith(".h"):
            data["c_abi_symbols"] = sorted(set(re.findall(r"\b(gliss_\w+)\s*\(", source)))
        for number, line in enumerate(source.splitlines(), 1):
            match = re.match(r"\s*#\s*include\s*[<\"]([^>\"]+)[>\"]", line)
            if match:
                data["imports"].append({"name": match[1], "line": number, "kind": "c_include"})
    return data


def records(base):
    names = paths()
    sources = {name: (ROOT / name).read_text(encoding="utf-8") for name in names}
    inspected = {name: inspect_source(name, source) for name, source in sources.items()}
    providers = defaultdict(list)
    abi_providers = defaultdict(list)
    for name, data in inspected.items():
        for module in data["modules"]:
            providers[module.lower()].append(name)
        python_module = module_name(name)
        if python_module:
            providers[python_module].append(name)
        if name.endswith(".h"):
            providers[Path(name).name].append(name)
        if name.startswith("src/"):
            for binding in data["c_binding_definitions"]:
                abi_providers[binding["name"]].append(name)
    dependencies = {name: [] for name in names}
    external = {name: [] for name in names}
    callers = defaultdict(list)

    def link(consumer, provider, kind, symbol, line):
        if provider == consumer:
            return
        edge = {"path": provider, "kind": kind, "symbol": symbol, "line": line}
        if edge not in dependencies[consumer]:
            dependencies[consumer].append(edge)
            callers[provider].append({"path": consumer, "kind": kind,
                                      "symbol": symbol, "line": line})

    for name, data in inspected.items():
        for imported in data["imports"]:
            named_providers = abi_providers if imported["kind"] == "c_external_binding_declaration" else providers
            matches = named_providers.get(imported["name"], [])
            if matches:
                for provider in matches:
                    link(name, provider, imported["kind"], imported["name"], imported["line"])
            elif imported["kind"] != "python_from_candidate":
                external[name].append(imported)
        for number, line in enumerate(sources[name].splitlines(), 1):
            for symbol in sorted(set(re.findall(r"\b\w+\b", line)) & abi_providers.keys()):
                for provider in abi_providers.get(symbol, []):
                    link(name, provider, "c_abi_token_reference", symbol, number)
        if not name.endswith((".f90", ".h")):
            # This intentionally reports candidates: comments and strings can
            # mention paths without executing or importing them.
            for candidate in names:
                if candidate in sources[name] and candidate != name:
                    position = sources[name].find(candidate)
                    line = sources[name][:position].count("\n") + 1
                    link(name, candidate, "literal_path_reference", candidate, line)
    reviews = json.loads(EVIDENCE.read_text(encoding="utf-8"))
    review_map = defaultdict(list)
    abi_review_ids = defaultdict(list)
    for review in reviews:
        for symbol in review.get("c_abi_symbols", []):
            abi_review_ids[symbol].append(review["id"])
        for name in review["paths"]:
            if name not in sources:
                raise ValueError("Evidence names an excluded or absent source: " + name)
            review_map[name].append({key: value for key, value in review.items() if key != "paths"})
    for name in names:
        data = inspected[name]
        candidates = sorted(set(edge["path"] for edge in callers[name]
                                if role(edge["path"]) == "verification_source"))
        targeted = review_map[name]
        yield {
            "schema_version": 1,
            "source_snapshot_base": base,
            "path": name,
            "source_sha256": hashlib.sha256((ROOT / name).read_bytes()).hexdigest(),
            "role": role(name),
            "language_or_format": Path(name).suffix.lstrip(".") or "build_configuration",
            "symbols": {key: data[key] for key in (
                "modules", "procedures", "c_abi_symbols", "c_binding_definitions")},
            "dependencies": sorted(dependencies[name], key=lambda edge: (edge["path"], edge["line"], edge["kind"])),
            "external_or_unresolved_imports": sorted(external[name], key=lambda edge: (edge["name"], edge["line"])),
            "caller_references": sorted(callers[name], key=lambda edge: (edge["path"], edge["line"], edge["kind"])),
            "reference_coverage": "static imports and token/literal references; dynamic calls and build selection not resolved",
            "read_coverage": "structural metadata scanned; whole implementation and equations unreviewed",
            "mathematical_contract": {"status": "not yet mapped", "scope": "whole component"},
            "units": {"status": "not yet mapped"},
            "assumptions": {"status": "not yet mapped"},
            "independent_oracle": {"status": "not yet qualified for the whole component",
                                   "candidate_verification_sources": candidates,
                                   "candidate_warning": "static reference is not evidence of oracle independence or execution"},
            "derivative_coverage": {"status": "not yet mapped for the whole component"},
            "reviewer": {"whole_component": "unreviewed", "targeted_evidence_owner": "controller record" if targeted else "unassigned"},
            "targeted_contract_evidence": targeted,
            "c_abi_entrypoint_contracts": [{"symbol": symbol,
                                           "status": "bounded evidence only; full contract unreviewed" if abi_review_ids[symbol] else "not yet mapped",
                                           "bounded_evidence_ids": abi_review_ids[symbol]}
                                           for symbol in data["c_abi_symbols"]],
            "unresolved_findings": data["parse_findings"] + [
                "Equation-to-code derivation, units, full derivative paths and independent review remain unaccounted for at component level."],
        }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="check metadata consistency only; no physics validation")
    parser.add_argument("--base", help="pin source snapshot base (default HEAD when writing)")
    args = parser.parse_args()
    existing = []
    if args.check:
        with OUTPUT.open(encoding="utf-8") as stream:
            existing = [json.loads(line) for line in stream]
        if not existing:
            raise SystemExit("Inventory is empty")
    base = args.base or (existing[0]["source_snapshot_base"] if existing else git("rev-parse", "HEAD"))
    base = git("rev-parse", "--verify", base + "^{commit}")
    generated = list(records(base))
    if args.check:
        if generated != existing:
            old = {record["path"]: record for record in existing}
            new = {record["path"]: record for record in generated}
            changed = sorted(name for name in old.keys() | new.keys() if old.get(name) != new.get(name))
            print("Stale inventory (first 12 paths): " + ", ".join(changed[:12]))
            raise SystemExit(1)
    else:
        with OUTPUT.open("w", encoding="utf-8") as stream:
            for count, record in enumerate(generated, 1):
                stream.write(json.dumps(record, sort_keys=True, separators=(",", ":")) + "\n")
                if count % 32 == 0:
                    stream.flush()
    roles = Counter(record["role"] for record in generated)
    edges = sum(len(record["dependencies"]) for record in generated)
    symbols = sum(len(record["symbols"]["c_abi_symbols"]) for record in generated
                  if record["path"].startswith("src/"))
    header_symbols = sum(len(record["symbols"]["c_abi_symbols"]) for record in generated
                         if record["role"] == "c_abi_declaration")
    print(json.dumps({"operation": "consistency_check" if args.check else "write",
                      "source_snapshot_base": base, "records": len(generated),
                      "roles": dict(sorted(roles.items())), "static_reference_edges": edges,
                      "native_gliss_c_symbols": symbols,
                      "public_header_symbols": header_symbols,
                      "whole_component_reviews": 0}, sort_keys=True))


if __name__ == "__main__":
    main()
