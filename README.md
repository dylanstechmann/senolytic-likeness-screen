# Senolytic-likeness descriptor screen (ChEMBL + RDKit)

Ranks a ChEMBL candidate pool by **descriptor-space proximity** to known
senolytic reference compounds (navitoclax, dasatinib, quercetin, fisetin,
piperlongumine, ruxolitinib, venetoclax) using locally computed RDKit
descriptors. Free data, free software, no paid ADMET.

This is a personal hobby and learning project, developed with substantial assistance from AI coding tools.

**What this is not:** an activity prediction. Descriptor neighbors of
senolytics are simply chemically similar molecules. The score carries no
biological, ADMET, or efficacy information — it is a triage aid for deciding
what to look at first in a real assay or structure-based workflow.

## Pipeline

1. Resolve reference compounds in `reference_compounds.tsv` via the ChEMBL
   molecule search API (`/chembl/api/data/molecule/search.json`).
2. Build a candidate pool from ChEMBL text queries (default: `flavonol`,
   `flavone`, `BCL-2 inhibitor`; override with `pool_queries.txt`).
3. Compute RDKit descriptors on the largest fragment (salts stripped):
   MW, MolLogP, TPSA, HBD, HBA, rotatable bonds, ring count,
   aromatic proportion, fraction sp3.
4. Z-standardize over the combined pool; score each candidate by
   `1/(1 + distance)` to its **nearest** reference in descriptor space.
5. Write `data/results.tsv` (full ranking), `REPORT.md` (top table +
   caveats), and a per-run provenance receipt.

## Usage

```bash
pip install -r requirements.txt
python scripts/screen.py --email you@university.edu
python -m unittest discover -s tests -v
```

Everything runs on a plain laptop; the only third-party dependency is RDKit.

## Reproducing / extending

- Swap the pool queries to focus the screen (e.g., `kinase inhibitor`,
  `flavonol glycoside`) — every run logs its exact queries, retrieval times,
  and per-compound ChEML IDs in `data/provenance/`.
- The reference set is literature senolytics; extend it the same way.
- To turn a top-ranked candidate into a structure picture, run the
  regen-workbench: `regen pubchem <name>` then `regen pymol-png`.

## License

MIT. ChEMBL data are used under the ChEMBL terms of use; cite the ChEMBL
paper (Mendez et al./Zdrazil et al., Nucleic Acids Research).
