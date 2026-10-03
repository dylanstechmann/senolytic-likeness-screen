#!/usr/bin/env python3
"""Senolytic-likeness descriptor screen (ChEMBL + RDKit, no paid ADMET).

This project computes structural descriptors for known senolytic reference
compounds (navitoclax, dasatinib, quercetin, ...) and a candidate pool from
ChEMBL, then ranks candidates by descriptor-space proximity to the nearest
reference compound.

IMPORTANT: descriptor proximity is NOT senolytic activity. This screen is a
chemistry triage aid for hypothesis generation only. No biological, ADMET,
or efficacy claim is made or implied by any score.

Usage:
  python scripts/screen.py --email you@university.edu
  python scripts/screen.py --email you@university.edu --pool-queries pool_queries.txt --top 20
"""
from __future__ import annotations

import argparse
import json
import math
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request, urlopen

CHEMBL_SEARCH_URL = "https://www.ebi.ac.uk/chembl/api/data/molecule/search.json"
DESCRIPTOR_KEYS = ["mw", "logp", "tpsa", "hbd", "hba", "rotb", "rings", "arom_prop", "fsp3"]
DEFAULT_POOL_QUERIES = ["flavonol", "flavone", "BCL-2 inhibitor"]


def http_json(url: str, user_agent: str, timeout: int = 60) -> dict:
    req = Request(url, headers={"User-Agent": user_agent, "Accept": "application/json"})
    with urlopen(req, timeout=timeout) as response:
        return json.loads(response.read().decode())


def now_utc() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def chembl_search(query: str, limit: int, user_agent: str) -> list[dict]:
    """Return ChEMBL molecules matching a text search, with canonical SMILES."""
    params = urlencode({"q": query, "limit": limit})
    data = http_json(f"{CHEMBL_SEARCH_URL}?{params}", user_agent)
    molecules = data.get("molecules", [])
    rows = []
    for mol in molecules:
        smiles = (mol.get("molecule_structures") or {}).get("canonical_smiles")
        chembl_id = mol.get("molecule_chembl_id")
        if smiles and chembl_id:
            rows.append(
                {
                    "chembl_id": chembl_id,
                    "name": mol.get("pref_name") or chembl_id,
                    "smiles": smiles,
                }
            )
    return rows


def compute_descriptors(smiles: str) -> dict | None:
    """RDKit descriptors on the largest fragment of a SMILES (salts stripped)."""
    try:
        from rdkit import Chem
        from rdkit.Chem import Crippen, Descriptors, Lipinski, rdMolDescriptors
    except ImportError:
        raise RuntimeError(
            "RDKit is required: pip install rdkit (see requirements.txt)"
        )
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    fragments = Chem.GetMolFrags(mol, asMols=True)
    if len(fragments) > 1:
        mol = max(fragments, key=lambda m: m.GetNumAtoms())
    heavy = mol.GetNumHeavyAtoms()
    aromatic = sum(1 for atom in mol.GetAtoms() if atom.GetIsAromatic())
    return {
        "mw": round(Descriptors.MolWt(mol), 2),
        "logp": round(Crippen.MolLogP(mol), 2),
        "tpsa": round(Descriptors.TPSA(mol), 2),
        "hbd": Lipinski.NumHDonors(mol),
        "hba": Lipinski.NumHAcceptors(mol),
        "rotb": Lipinski.NumRotatableBonds(mol),
        "rings": rdMolDescriptors.CalcNumRings(mol),
        "arom_prop": round(aromatic / heavy, 3) if heavy else 0.0,
        "fsp3": round(rdMolDescriptors.CalcFractionCSP3(mol), 3),
    }


def standardize(rows: list[dict]) -> None:
    """Add z-score columns (<key>_z) computed over the provided rows."""
    for key in DESCRIPTOR_KEYS:
        values = [row["descriptors"][key] for row in rows if row.get("descriptors")]
        if len(values) < 2:
            mean, std = (values[0] if values else 0.0), 1.0
        else:
            mean = sum(values) / len(values)
            variance = sum((value - mean) ** 2 for value in values) / len(values)
            std = math.sqrt(variance)
        if std == 0:
            std = 1.0
        for row in rows:
            if row.get("descriptors"):
                row.setdefault("z", {})[key] = (row["descriptors"][key] - mean) / std


def z_distance(a: dict, b: dict) -> float:
    """Euclidean distance in standardized descriptor space."""
    return math.sqrt(sum((a[key] - b[key]) ** 2 for key in DESCRIPTOR_KEYS))


def proximity_score(distance: float) -> float:
    """Monotone score in (0, 1]: identical vectors score 1.0."""
    return round(1.0 / (1.0 + distance), 4)


def score_candidates(candidates: list[dict], references: list[dict]) -> None:
    """Attach nearest reference, distance, and proximity score to each candidate."""
    for candidate in candidates:
        best_name, best_distance = None, None
        for reference in references:
            distance = z_distance(candidate["z"], reference["z"])
            if best_distance is None or distance < best_distance:
                best_name, best_distance = reference["name"], distance
        candidate["nearest_reference"] = best_name
        candidate["distance"] = round(best_distance, 3) if best_distance is not None else None
        candidate["score"] = (
            proximity_score(best_distance) if best_distance is not None else None
        )


def read_tsv(path: Path, columns: int = 2) -> list[tuple[str, ...]]:
    rows = []
    for line in path.read_text().splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        parts = tuple(part.strip() for part in stripped.split("\t"))
        if parts[0].lower() in {c.lower() for c in columns_headings(columns)}:
            continue
        if parts:
            rows.append(parts[:columns])
    return rows


def columns_headings(columns: int) -> tuple[str, ...]:
    return ("name", "query")[:columns]


def render_report(run_receipt: dict, references: list[dict], ranked: list[dict], top: int) -> str:
    header = (
        "# Senolytic-likeness descriptor screen\n\n"
        f"Generated by `senolytic-likeness-screen` on {run_receipt['started']}"
        " (ChEMBL + RDKit, local descriptors).\n\n"
        "Reference compounds (known senolytics / senolytic-combo partners): "
        + ", ".join(sorted({r['name'] for r in references}))
        + ".\n\n"
        "## Top candidates by descriptor-space proximity\n\n"
        "| rank | ChEMBL ID | name | score | nearest reference | MW | LogP | TPSA | rings |\n"
        "|---|---|---|---|---|---|---|---|---|\n"
    )
    lines = []
    for position, row in enumerate(ranked[:top], start=1):
        d = row["descriptors"]
        lines.append(
            f"| {position} | {row['chembl_id']} | {row['name']} | {row['score']} "
            f"| {row['nearest_reference']} | {d['mw']} | {d['logp']} | {d['tpsa']} "
            f"| {d['rings']} |"
        )
    caveats = (
        "\n## How to read this (and how not to)\n\n"
        "- The score is a **descriptor-space proximity** to the nearest reference\n"
        "  compound (MW, LogP, TPSA, H-bond donors/acceptors, rotatable bonds,\n"
        "  ring count, aromatic proportion, fraction sp3 — z-standardized over\n"
        "  the combined pool). It is **not** a senolytic activity prediction and\n"
        "  carries **no** ADMET or efficacy information.\n"
        "- Descriptor neighbors of known senolytics are simply chemically similar\n"
        "  molecules; many will be entirely inert in senescence biology.\n"
        "- Salts are stripped (largest fragment) before descriptor calculation.\n"
        "- Use this ranking to triage what to look at first in a real assay or\n"
        "  a structure-based workflow; do not treat it as evidence.\n"
    )
    return header + "\n".join(lines) + caveats


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="screen")
    parser.add_argument("--email", default="senolytic-likeness@local", help="polite-pool contact")
    parser.add_argument("--reference-file", default="reference_compounds.tsv")
    parser.add_argument("--pool-queries-file", default=None, help="one ChEMBL query per line")
    parser.add_argument("--pool-size", type=int, default=15, help="candidates per pool query")
    parser.add_argument("--top", type=int, default=20, help="rows in the report table")
    parser.add_argument("--out", default=".", help="project root")
    args = parser.parse_args(argv)

    root = Path(args.out).resolve()
    user_agent = f"senolytic-likeness-screen/1.0 (mailto:{args.email})"
    reference_rows = read_tsv(root / args.reference_file)
    if not reference_rows:
        print(f"no references in {args.reference_file}", file=sys.stderr)
        return 2

    pool_queries = DEFAULT_POOL_QUERIES
    if args.pool_queries_file:
        lines = [
            line.strip()
            for line in (root / args.pool_queries_file).read_text().splitlines()
            if line.strip() and not line.strip().startswith("#")
        ]
        if lines:
            pool_queries = lines

    run_receipt = {
        "tool": "senolytic-likeness-screen",
        "user_agent": user_agent,
        "started": now_utc(),
        "pool_queries": pool_queries,
        "pool_size": args.pool_size,
        "references": [],
        "candidates": [],
    }

    # 1. Reference set: fetch known senolytics from ChEMBL by name query.
    references: list[dict] = []
    seen_ids: set[str] = set()
    for name, query in reference_rows:
        time.sleep(0.5)
        matches = chembl_search(query, 5, user_agent)
        chosen = None
        for match in matches:
            if match["name"].upper() == name.upper() and match["chembl_id"] not in seen_ids:
                chosen = match
                break
            if chosen is None and query.lower() in (match["name"] or "").lower():
                chosen = match
        if chosen is None and matches:
            chosen = matches[0]
        if not chosen:
            print(f"reference {name}: no ChEMBL match; skipped", file=sys.stderr)
            run_receipt["references"].append({"name": name, "status": "no match"})
            continue
        descriptors = compute_descriptors(chosen["smiles"])
        if descriptors is None:
            print(f"reference {name}: unparseable SMILES; skipped", file=sys.stderr)
            run_receipt["references"].append({"name": name, "status": "bad smiles"})
            continue
        seen_ids.add(chosen["chembl_id"])
        references.append({**chosen, "name": name, "role": "reference", "descriptors": descriptors})
        run_receipt["references"].append(
            {"name": name, "chembl_id": chosen["chembl_id"], "status": "ok"}
        )
    if not references:
        print("no reference compounds resolved; aborting", file=sys.stderr)
        return 1

    # 2. Candidate pool from ChEMBL text searches.
    candidates: list[dict] = []
    for query in pool_queries:
        time.sleep(0.5)
        matches = chembl_search(query, args.pool_size, user_agent)
        kept = 0
        for match in matches:
            if match["chembl_id"] in seen_ids:
                continue
            seen_ids.add(match["chembl_id"])
            descriptors = compute_descriptors(match["smiles"])
            if descriptors is None:
                run_receipt["candidates"].append(
                    {"chembl_id": match["chembl_id"], "status": "bad smiles"}
                )
                continue
            candidates.append({**match, "role": "candidate", "descriptors": descriptors})
            run_receipt["candidates"].append(
                {"chembl_id": match["chembl_id"], "query": query, "status": "ok"}
            )
            kept += 1
        print(f"pool query '{query}': kept {kept} candidates")
    return finish(root, run_receipt, references, candidates, args.top)


def finish(root: Path, run_receipt: dict, references: list[dict], candidates: list[dict], top: int) -> int:
    if not candidates:
        print("no candidates resolved; aborting", file=sys.stderr)
        return 1
    combined = references + candidates
    standardize(combined)
    score_candidates(candidates, references)
    ranked = sorted(candidates, key=lambda row: (row["score"] is None, -(row["score"] or 0)))

    result_columns = [
        "rank", "chembl_id", "name", "score", "distance", "nearest_reference",
        "mw", "logp", "tpsa", "hbd", "hba", "rotb", "rings", "arom_prop", "fsp3", "smiles",
    ]
    lines = ["\t".join(result_columns)]
    for position, row in enumerate(ranked, start=1):
        d = row["descriptors"]
        lines.append(
            "\t".join(
                str(value)
                for value in [
                    position, row["chembl_id"], row["name"], row["score"], row["distance"],
                    row["nearest_reference"], d["mw"], d["logp"], d["tpsa"], d["hbd"],
                    d["hba"], d["rotb"], d["rings"], d["arom_prop"], d["fsp3"], row["smiles"],
                ]
            )
        )
    results_path = root / "data" / "results.tsv"
    results_path.parent.mkdir(parents=True, exist_ok=True)
    results_path.write_text("\n".join(lines) + "\n", encoding="utf-8")

    run_receipt["finished"] = now_utc()
    receipt_path = root / "data" / "provenance" / f"{time.time_ns()}_run.json"
    receipt_path.parent.mkdir(parents=True, exist_ok=True)
    receipt_path.write_text(json.dumps(run_receipt, indent=2), encoding="utf-8")

    report = render_report(run_receipt, references, ranked, top)
    (root / "REPORT.md").write_text(report, encoding="utf-8")

    print(
        f"scored {len(candidates)} candidates against {len(references)} references; "
        f"results: {results_path.name}; report: REPORT.md; provenance: {receipt_path.name}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
