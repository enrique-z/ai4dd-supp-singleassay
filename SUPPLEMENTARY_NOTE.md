# Supplementary Note — single-assay benchmark

## Evidence files

The files in `evidence/` (single_assay_benchmark_ledger.jsonl,
single_assay_benchmark_stamp_manifest.json, recompute_2026_09_05.json,
recheck_2026_09_05.json, t1_canonical.json, t1_arm_c_knn.json,
t1_bootstrap_cis.json, t1_lookup_baseline.json, t1_verification.json,
t2_truncation_provenance.json, interassay_floor_2026_09_05.json,
figure_claims.json) are the outputs of the original run. No value in them has
been changed. In the ledger, the manifest and `t1_verification.json`, a few
free-text provenance annotations (`verdict_reason`, `source` labels, two key
names) were reworded for this release. In `interassay_floor_2026_09_05.json`,
two key names and one note were reworded.

On 2026-10-24, `figure_claims.json` and `interassay_floor_2026_09_05.json`
were replaced by the final versions that the paper was built from. The first
release carried earlier copies of both.

The data paths recorded in the evidence files (`db`, `benchmark` and similar
fields) are written relative to the repository root, e.g.
`data/chembl_36/chembl_36_sqlite/chembl_36.db`. The original run used absolute
local paths. Nothing else in those fields was changed.

## Reproducing the analysis

Requirements: Python 3.10+, NumPy, RDKit, matplotlib, a local SQLite copy of
ChEMBL 36, and the TDC `Lipophilicity_AstraZeneca` file
(`lipophilicity_astrazeneca.tab`, PyTDC 1.1.15).

Run the scripts from the repository root. They read their inputs from `data/`
(`BENCHMARK`, `DB_PATH` in `make_figure.py`; `BENCHMARK`, `DB_URI` in
`bootstrap_ci.py`):

```bash
# data/lipophilicity_astrazeneca.tab
# data/chembl_36/chembl_36_sqlite/chembl_36.db
python3 analysis/make_figure.py    # Figure 1 -> fig1.pdf, evidence/figure_claims.json
python3 analysis/bootstrap_ci.py   # Table 1 CIs -> evidence/t1_bootstrap_cis.json (+ ledger/manifest rows)
python3 analysis/guard.py --selftest
```

`make_figure.py` rebuilds the TDC scaffold split (seed 42), recomputes the
Arm A / Arm B / Arm C errors from the benchmark file and ChEMBL 36, and aborts
if any aggregate disagrees with `evidence/t1_canonical.json`. `bootstrap_ci.py`
recomputes the same per-compound errors, aborts if the MAEs disagree with the
ledger to printed precision, and then runs the seeded bootstrap. Both scripts
overwrite files in `evidence/`, so run them in a copy of the repository if you
want to compare against the shipped evidence.

The census, T1 and T2 evidence files were produced by read-only queries
against ChEMBL 36. Their results are recorded in the JSON files listed above.

To check a manuscript source against the evidence, place `main.tex` in the
directory that contains `evidence/` and run
`python3 analysis/guard.py --tex main.tex --strict`.

## Pre-registration

`A2_PREDECLARATION_2026-09-04.md` is the pre-registration cited in the paper.
For this release, internal workflow notes were removed from it. The claim,
refutation criteria (T1, T2), rival mechanisms, analysis plan and commitment
are unchanged.
