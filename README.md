# The TDC Lipophilicity benchmark is one assay — and a predictor that fits nothing retrieves 834 of 840 held-out labels exactly

Supplementary code, pre-registration and evidence ledger.

**Authors:** Enrique Zueco, David Scott Lewis, Haley Yi
**Affiliation:** AIXC Bio, AIXC Research, Zaragoza, Spain
**Contact:** reports@aixcresearch.com
**Venue:** AI for Drug Discovery (AI4DD) workshop, NeurIPS 2026

## What this is

Analysis code, the pre-registration document, and an evidence ledger that backs every quantitative claim in the paper.

## Directory layout

| Path | Contents |
|------|----------|
| `A2_PREDECLARATION_2026-09-04.md` | Pre-registration (claim, refutation criteria, analysis plan), written before T1 and T2 were run |
| `analysis/make_figure.py` | Figure 1: recomputes the T1 aggregates from the benchmark file and ChEMBL 36, asserts them against the evidence, then draws |
| `analysis/bootstrap_ci.py` | Bootstrap 95% CIs (B = 10,000, seed 42) for the Arm B / Arm C MAEs in Table 1 |
| `analysis/guard.py` | Integrity guard: checks every stamped number in the manuscript against the stamp manifest and ledger |
| `evidence/testfold_identity_and_misses.py` | Checks that the scaffold test fold is the fixed test set of TDC's ADMET benchmark group and decomposes the six lookup misses (paper Appendix A) → `evidence/testfold_identity_and_misses.json` |
| `evidence/` | Raw JSON outputs, evidence ledger (`*_ledger.jsonl`) and stamp manifest (`*_stamp_manifest.json`) |
| `evidence/README.md` | Ledger and manifest conventions |
| `analysis/README.md` | How the scripts fit together and how to rerun them |
| `SUPPLEMENTARY_NOTE.md` | Provenance notes and reproduction instructions |
| `MANIFEST.txt` | Size and SHA-256 of every file except `README.md` and `MANIFEST.txt` |
| `refs.bib` | Bibliography |

Every number in the manuscript resolves from the evidence ledger through the stamp manifest. See `SUPPLEMENTARY_NOTE.md` for provenance and reproduction.
