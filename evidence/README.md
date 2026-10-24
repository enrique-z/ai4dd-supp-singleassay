# Evidence ledger — conventions

Every numeric literal in `main.tex` resolves to a row in this folder's ledger,
except the `\SafeNum{...}` architecture constants and nominal parameters.
`guard.py` fails if it cannot resolve one.

## Ledger row schema

One JSONL file per ledger, one row per claim. The full row shape is:

```json
{
  "claim_id": "t1::arm_a_full_mae",
  "metric": "mae",
  "value": 0.024,
  "source": "/abs/path/to/evidence/t1_canonical.json",
  "verdict": "SUPPORTED",
  "verdict_reason": "t1_canonical.json task3.arm_a.full_fold_metrics.mae=0.02399; rounded to 0.0240.",
  "ci_low": null,
  "ci_high": null,
  "must_print": true,
  "timestamp": "2026-09-05T00:00:00Z"
}
```

Early rows in `single_assay_benchmark_ledger.jsonl` carry only `claim_id`,
`verdict` and `verdict_reason`. Their values live in the stamp manifest.

Verdicts:
- `SUPPORTED` — the value traces to a real measurement, the real database or
  a cited source.
- `UNSUPPORTED` — blocks. The literal in `main.tex` must revert to
  `\ResultCell` or be removed.
- `NOT_APPLICABLE` — informational only.
- `CONTESTED` — needs review, and blocks release until resolved.

When a `claim_id` appears more than once, the later row supersedes the
earlier one.

## The `must_print` field

`must_print` controls which ledger rows must appear in the manuscript.

```
must_print: true   -- the guard FAILS if this row is omitted from main.tex.
                      Used for every pre-registered outcome and headline claim.
must_print: false  -- the row is recorded for audit but need not be typeset.
                      Used for decomposition/robustness rows that back a printed
                      range: the range is what gets typeset, not the row itself.
```

A manifest with no `must_print` flags at all falls back to requiring every row
to be typeset.

## The stamp manifest

`<paper>_stamp_manifest.json` maps each `\StampedCell{k1}{k2}{v}` key to the
value it must equal:

```json
{
  "ledger_rows": {
    "f1::n_rows": {
      "value": 4200,
      "claim_id": "f1::n_rows",
      "source": "t2_truncation_provenance.json:q1.total_rows_assay",
      "must_print": false
    }
  }
}
```

`guard.py` walks `main.tex` for every `\StampedCell{k1}{k2}{v}` and asserts
that `manifest.ledger_rows["k1::k2"].value == v`, with numeric tolerance. It
also asserts that the backing ledger row's verdict is `SUPPORTED`. Claim ids
that a figure script lists in `figure_claims.json` count as reported.
