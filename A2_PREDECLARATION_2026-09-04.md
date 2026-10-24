# A2 — Pre-declaration

    paper_slug: single-assay-benchmark
    date:       2026-09-04 (written 2026-09-04T22:00Z; all times UTC)
    venue:      NeurIPS.cc/2026/Workshop/AI4DD
    deadline:   2026-09-06 11:59:00 UTC (re-probed live 2026-09-04T21:50Z, unchanged)
    hypothesis: H15

---

## exploratory_passes_preceding

**This paper has substantial exploratory work preceding pre-declaration, and
none of the findings below were predicted. They were measured first.** Stating
this plainly is the point of the field.

Findings **F1–F5** (identifier retention, single-group concentration, exact value
reproduction, bounded support, and the unmatched-compound contrast) were all
measured before this document was written — first during exploratory analysis,
then re-derived independently on 2026-09-04 between 21:50Z and 22:00Z. They
are **exploratory**, they are **descriptive**, and the paper must present them
as an audit of what is in the file, never as confirmed predictions.

Hypothesis H15 additionally carries a **recorded correction**: an earlier
version of the claim asserted that the three TDC benchmark files "carry no source
identifier of any kind." That was **false** for the file this paper audits —
100.00% of its `Drug_ID` values are literal ChEMBL accessions. The error was
caught by internal review, the claim was inverted, and the inverted claim is
what this paper argues. The false version must not survive anywhere in the draft.

Two things below are **genuinely confirmatory** — they have not been run, their
outcomes are not known to anyone at the time of this commit, and each has a
pre-registered criterion that can lose:

- **T1**, the learning-free baseline (the paper's load-bearing experiment);
- **T2**, the truncation-provenance test (which of two parties applied the range
  bound — this changes which sentence is true, and I do not currently know which).

This file is committed to git **before either result is read.**

---

## claim

The TDC / MoleculeNet "Lipophilicity (AstraZeneca)" benchmark — 4,200 compounds,
a standard ADMET regression task — is in substance a verbatim extract of a
**single ChEMBL assay** (`assay_id 1369278` = `CHEMBL3301363`), reproduced to
floating-point exactness, still carrying its source accession in every row, and
bounded exactly at that assay's declared operating range.

Therefore: **every reported "ADMET generalisation" number on this benchmark is
generalisation across chemical structures within one measurement protocol** — not
across protocol, laboratory, instrument, or value range. The published dataset
description gives size and split type and says nothing about group structure, and
releases no per-row source keys, so a user cannot discover this without
re-deriving it from ChEMBL.

**The claim is about the scope of a licensed inference, not about data quality.**
Single-protocol data is internally consistent and that is a virtue, not a defect.
Any draft sentence that reads as "this dataset is bad" is a defect in the draft.

---

## refutation_criterion

### T1 — the learning-free baseline (load-bearing; the paper depends on it)

A lookup-only predictor with **no fitting of any kind**, evaluated on the
benchmark's own scaffold split, is scored against the public TDC leaderboard
(reported MAE ≈ 0.47–0.535, to be verified against the cited sources before it
is printed).

- **PASS (claim supported):** the learning-free lookup achieves **test-fold
  MAE ≤ 0.10**, i.e. it beats the entire published leaderboard by a factor of
  several while learning nothing. This establishes that the benchmark cannot
  distinguish a model that learned lipophilicity from one that retrieved a public
  label, and that no split protocol on this dataset detects it.
- **FAIL (claim refuted, paper reduced):** the lookup does **not** beat the
  leaderboard's best reported MAE. If so, the contamination argument collapses,
  the "pretraining-corpus contamination" angle is withdrawn, and the paper ships
  as a **shorter, honest provenance audit** of F1–F4 with the baseline reported as
  a negative result — or is dropped. **It will not be padded to fill five pages.**
- **AMBIGUOUS (0.10 < MAE ≤ 0.47):** report the number as measured, state that it
  is competitive-but-not-dominant, and weaken the headline accordingly.

Two coverage figures are pre-registered as mandatory reporting **whatever the
MAE**, because a lookup that covers few compounds is not a baseline:

1. **coverage-with-source** — the fraction of test-fold compounds for which a
   value is retrievable from ChEMBL *including* assay 1369278. Expected high; if
   it is not, T1 is uninterpretable and must be reported as such.
2. **coverage-without-source** — the fraction retrievable from ChEMBL *excluding*
   assay 1369278, i.e. from a genuinely independent protocol. **I predict this
   will be very low**, and if it is, that is itself a finding and must be
   reported as one: it means the benchmark's compounds cannot be validated
   against any other protocol, which is the affirmative evidence for the
   central claim. If coverage-without-source is instead HIGH, the "one protocol
   and no other" framing is wrong and must be softened.

MAE, RMSE and Pearson R are all reported for every arm. The scaffold split must be
the benchmark's own; if it cannot be reproduced deterministically, that fact is
reported and the arm is labelled an approximation, not a replication.

### T2 — truncation provenance (determines which sentence is true)

The benchmark's labels are bounded at exactly [−1.5, 4.5] with zero rows outside.
The source assay's ChEMBL description declares that same range. **I do not know
which party applied the bound**, and the paper must not guess.

- **Outcome A** — ChEMBL's `assay_id 1369278` holds **no** rows outside
  [−1.5, 4.5]: the window is the assay's own, and it propagated into the
  benchmark silently, losing its declaration on the way.
- **Outcome B** — ChEMBL holds rows outside the window that the benchmark
  **dropped**: the truncation was applied downstream, by MoleculeNet or TDC, and
  the paper must attribute it there.
- **Outcome C** — ChEMBL holds out-of-range rows flagged by `standard_relation`
  (`<`/`>`) or by `data_validity_comment`, and a standard `standard_relation='='`
  filter removed them: the truncation has an identifiable mechanism, and the paper
  reports a demonstrated mechanism rather than an inference.

All three are publishable; they are **not** interchangeable, and the draft may
contain only the one that is measured. Note for the draft: the correct
statistical term is **truncation** (out-of-range observations removed), not
**censoring** (retained at the bound). The absence of a pile-up at either bound is
the evidence, and the distinction must be stated rather than left for a reviewer.

---

## rival_mechanism

1. **"This is already known."** That the benchmark derives from an AstraZeneca
   ChEMBL deposition is cited by TDC (reference [11], never called from the body
   text), by PharmaBench (ref. 43), and stated by MoleculeNet. **The paper must
   never claim concealment.** The defensible novelty is the *quantification*
   (99.60% one assay; 99.93% exact value identity; 100% accession retention; 2,889
   available assays and one used, 11.7× the runner-up) and the *consequence*
   (single-protocol generalisation; pretraining contamination). This is stated in
   the introduction, not buried, or a reviewer states it first.
2. **"The concentration is a query artifact."** Ruled out by counting how many
   benchmark molecules the *next* assays hold; if the next-largest holds single
   digits, the concentration is total rather than gradual.
3. **"Bounds could arise from the chemistry, or from curator clipping."** Clipping
   is ruled out by the absence of mass at the bound (F4). The chemistry
   explanation is ruled out by the bounds falling on exactly the two round numbers
   the assay description itself declares. T2 discriminates the remainder.
4. **"Single-assay data is high quality — you are attacking good data."** Correct,
   and conceded **in the abstract**, not the discussion. The concession costs
   nothing and its absence is fatal.
5. **"The lookup baseline is a tautology — you looked up the answer."** The
   pre-registered response is coverage-without-source (above) plus the framing:
   the point is not that a lookup works, it is that *nothing in the benchmark
   protocol detects or prevents it*, and that a ChEMBL-pretrained model has the
   same access by construction.

---

## analysis_plan

1. **F1–F5 re-derivation** (exploratory, already measured; an independent
   re-derivation is made without sight of the target numbers; any disagreement
   is adjudicated before drafting). Read-only against ChEMBL 36 and the TDC file.
2. **T2** — truncation provenance, as specified above.
3. **T1** — learning-free baseline on the benchmark's own scaffold split, with
   both coverage arms.
4. Every printed number becomes a `\StampedCell` backed by a ledger row.
   No number enters `main.tex` without one; the guard is satisfied by producing
   the ledger, never by weakening the guard.

**Scope constraint.** This paper audits **one** TDC file's construction. It does
**not** analyse the three TDC benchmark files as a set or reprint statistics
derived from them. F5 is supporting-only and is the **first thing cut** for space.

This paper concerns one benchmark's construction and the inference it licenses,
not the semantics of ChEMBL's fields (e.g. inequality-flagged potency labels);
the framing must keep that distinction visible.

---

## ledger_rows_expected

Approximately **35–55** rows. Anticipated groups:

- F1 identifier retention: ~6 rows (row count, accession fraction, uniqueness,
  resolution count, InChIKey agreement, mismatch count)
- F2 concentration: ~10 rows (in-assay count and %, next-assay counts, per-type
  assay and record counts, pooled counts, top-5 sizes, largest/runner-up ratio)
- F3 value identity: ~8 rows (denominator, agreement count and %, median |Δ|,
  max |Δ|, and the 3 discrepant rows individually)
- F4 bounded support: ~8 rows (min, max, counts outside each side, counts at each
  bound, median, mean, bin counts supporting exclusion-not-clipping)
- T1 baseline: ~8–12 rows (split sizes, both coverage figures, MAE/RMSE/R per arm,
  leaderboard comparison values)
- T2 provenance: ~4 rows (out-of-range count in source assay, relation
  distribution, validity-comment distribution)

A materially different final count is expected and is not a defect; this is a
declared estimate, not a target to be hit.

---

## Commitment

Committed to git before T1 or T2 is read. If T1 fails its criterion, the failure
is reported in the paper rather than reframed away, and the paper is shortened or
dropped.
