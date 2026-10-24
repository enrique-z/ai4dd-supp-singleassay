# Analysis scripts

The scripts locate `evidence/` as `../evidence`, relative to this directory,
and read their inputs from `data/` at the repository root. The exact commands
are in `SUPPLEMENTARY_NOTE.md`.

## guard.py — the anti-fabrication check

`guard.py --tex main.tex` resolves every stamped literal in the manuscript
against `evidence/*_stamp_manifest.json` and `evidence/*_ledger.jsonl`, and
fails on any of:

- `E1` an unresolved `\ResultCell` placeholder
- `E2` a `\StampedCell` key missing from the manifest
- `E3` value drift: the manuscript value differs from the manifest value
- `E4` the backing ledger row is missing or its verdict is not `SUPPORTED`
- `E5` a manifest row flagged `must_print` (a pre-registered outcome) that the
  manuscript does not print
- `E6` a literal `[TBD]` in the manuscript
- `E7` a bare prose number. This is always an error when it duplicates a
  stamped manifest value, and otherwise a warning unless `--strict` is set.

Run `python3 analysis/guard.py --selftest` to check the guard itself on temporary
fixtures. Passing the guard shows that the manuscript agrees with the
manifest. Whether the manifest agrees with the data is recorded in each ledger
row's `source` and `verdict_reason`.

## make_figure.py — recompute, assert, draw

1. **Load** the benchmark file, ChEMBL 36 and the T1 evidence.
2. **Recompute** every aggregate the figure draws: the TDC scaffold split, the
   Arm A / B / C errors and the MAEs.
3. **Assert** each value against `evidence/t1_canonical.json` and
   `evidence/t1_arm_c_knn.json`, and abort on any disagreement.
4. **Draw** `fig1.pdf`. `pdf.fonttype` is set to 42, so no Type-3 fonts.
5. **Write** `evidence/figure_claims.json`, which declares the claim_ids the
   figure draws so that the guard counts them as reported.

## bootstrap_ci.py — Table 1 confidence intervals

The script recomputes per-compound absolute errors for Arm B (full fold and
covered subset) and Arm C (k = 5). It aborts unless the resulting MAEs match
the stored ledger values to printed precision. It then runs a seeded
bootstrap (B = 10,000, seed 42) over compounds and writes
`evidence/t1_bootstrap_cis.json`, plus the matching ledger rows and manifest
entries. Timestamps are pinned (`SOURCE_DATE_EPOCH` or a fixed sentinel), so
two consecutive runs produce byte-identical output.

## Conventions every script follows

- Read real data only: ChEMBL 36 opened read-only and the published TDC file.
  Nothing is generated with a random-number generator except the seeded
  bootstrap resamples.
- Write ledger rows using the schema in `evidence/README.md`. A row whose
  verdict is not `SUPPORTED` is kept for inspection, never silently dropped.
- Be deterministic: the same inputs give the same ledger rows byte for byte.
