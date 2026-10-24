#!/usr/bin/env python3
"""Bootstrap 95% CIs on Arm B (full-fold, covered subset) and Arm C k=5 MAE,
plus the paired Arm B - Arm C MAE difference.

Recompute-then-assert pattern: recompute per-compound absolute errors from
first principles, ASSERT that the resulting aggregate MAEs match the stored
ledger values to their printed precision (integrity gate), then run a seeded
bootstrap over per-compound errors (B=10000, seed=42) and emit both an
evidence JSON and appended ledger rows / manifest entries.

Determinism: seed pins both scaffold-split (via TDC replica in make_figure.py
path) and the bootstrap resamples. Two consecutive runs produce byte-identical
outputs (evidence JSON + ledger rows + manifest patch).

Data inputs (read-only):
  benchmark: data/lipophilicity_astrazeneca.tab
  ChEMBL 36: data/chembl_36/chembl_36_sqlite/chembl_36.db  (both relative to the repo root)
             (opened via URI file:...?mode=ro; never copied, never written)
"""

from __future__ import annotations

import csv
import json
import sqlite3
import sys
from collections import defaultdict
from pathlib import Path
from random import Random

import numpy as np

from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import AllChem
from rdkit.Chem.inchi import MolToInchi, InchiToInchiKey
from rdkit.Chem.Scaffolds import MurckoScaffold

RDLogger.DisableLog("rdApp.*")

# ---------------------------------------------------------------------------
# Paths / constants
# ---------------------------------------------------------------------------
HERE     = Path(__file__).resolve().parent
PAPER    = "single_assay_benchmark"
PAPER_DIR = HERE.parent
EVIDENCE = PAPER_DIR / "evidence"
LEDGER   = EVIDENCE / f"{PAPER}_ledger.jsonl"
MANIFEST = EVIDENCE / f"{PAPER}_stamp_manifest.json"
OUT_JSON = EVIDENCE / "t1_bootstrap_cis.json"

BENCHMARK = PAPER_DIR / "data" / "lipophilicity_astrazeneca.tab"
DB_URI    = f"file:{PAPER_DIR / 'data' / 'chembl_36' / 'chembl_36_sqlite' / 'chembl_36.db'}?mode=ro"

SEED_SPLIT     = 42
SEED_BOOTSTRAP = 42
B_REPLICATES   = 10000
FRAC           = [0.7, 0.1, 0.2]
SOURCE_ASSAY   = 1369278            # CHEMBL3301363 assay_id
K_NEIGH        = 5

# Integrity-gate tolerances (must match stored ROUNDED ledger values exactly)
# Ledger prints: arm_b_full_mae=0.764 (3dp), arm_b_covered_mae=0.1242 (4dp),
# arm_c_k5_mae=0.7868 (4dp). Reproduce to 6dp then compare rounded.
TOL_MAE = 5e-4   # tolerance on the 3-4dp stored rounded value

# Stored ledger values (from evidence/single_assay_benchmark_stamp_manifest.json)
LEDGER_ARM_B_FULL_MAE    = 0.764
LEDGER_ARM_B_COVERED_MAE = 0.1242
LEDGER_ARM_C_K5_MAE      = 0.7868


# ---------------------------------------------------------------------------
# 1. Parse benchmark (quoted TAB file — assert 100% parse rate)
# ---------------------------------------------------------------------------
def parse_benchmark() -> list[tuple[str, str, float]]:
    rows: list[tuple[str, str, float]] = []
    with BENCHMARK.open() as f:
        reader = csv.reader(f, delimiter="\t", quotechar='"')
        header = next(reader)
        assert header[:3] == ["Drug_ID", "Drug", "Y"], f"unexpected header: {header}"
        for r in reader:
            # Assert quotes stripped: bare Drug_ID should look like CHEMBLXXX
            drug_id = r[0]
            smiles  = r[1]
            y       = float(r[2])
            assert not drug_id.startswith('"') and not drug_id.endswith('"'), \
                f"quoted Drug_ID leaked: {drug_id!r}"
            rows.append((drug_id, smiles, y))
    # Assert 100% parse rate — the known TDC .tab quoted-file trap.
    assert len(rows) == 4200, f"expected 4200 rows, got {len(rows)}"
    # Assert every SMILES parses (RDKit-side)
    parsed = sum(1 for _, s, _ in rows if Chem.MolFromSmiles(s) is not None)
    assert parsed == 4200, f"parse rate {parsed}/4200 != 100%"
    return rows


# ---------------------------------------------------------------------------
# 2. TDC scaffold split (exact replica — mirrors analysis/make_figure.py)
# ---------------------------------------------------------------------------
def build_scaffold_split(rows) -> tuple[list[int], list[int], list[int]]:
    n          = len(rows)
    train_size = int(n * FRAC[0])
    val_size   = int(n * FRAC[1])
    test_size  = n - train_size - val_size

    scaffolds: defaultdict[str, set] = defaultdict(set)
    for i, (_, smiles, _) in enumerate(rows):
        mol = Chem.MolFromSmiles(smiles)
        sca = ("" if mol is None
               else MurckoScaffold.MurckoScaffoldSmiles(mol=mol, includeChirality=False))
        scaffolds[sca].add(i)

    big, small = [], []
    for idx_set in scaffolds.values():
        if len(idx_set) > val_size / 2 or len(idx_set) > test_size / 2:
            big.append(idx_set)
        else:
            small.append(idx_set)

    rng = Random(SEED_SPLIT)
    rng.seed(SEED_SPLIT)          # TDC re-seeds before first shuffle
    rng.shuffle(big)
    rng.shuffle(small)
    all_groups = big + small

    train_idx, val_idx, test_idx = [], [], []
    for group in all_groups:
        g = list(group)
        if len(train_idx) + len(g) <= train_size:
            train_idx.extend(g)
        elif len(val_idx) + len(g) <= val_size:
            val_idx.extend(g)
        else:
            test_idx.extend(g)
    return train_idx, val_idx, test_idx


# ---------------------------------------------------------------------------
# 3. InChIKeys for ChEMBL lookup
# ---------------------------------------------------------------------------
def compute_inchikeys(rows) -> list[str | None]:
    iks: list[str | None] = []
    for _, smiles, _ in rows:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            iks.append(None); continue
        inchi = MolToInchi(mol)
        iks.append(None if inchi is None else InchiToInchiKey(inchi))
    return iks


# ---------------------------------------------------------------------------
# 4. Query ChEMBL 36 (READ-ONLY URI) for logD values
# ---------------------------------------------------------------------------
def query_logd(test_iks_set: set[str]) -> dict[str, dict[str, list[float]]]:
    con = sqlite3.connect("file::memory:", uri=True)
    con.execute(f"ATTACH DATABASE '{DB_URI}' AS chembl")
    con.execute("CREATE TEMP TABLE q_iks (ik TEXT PRIMARY KEY)")
    for ik in test_iks_set:
        con.execute("INSERT OR IGNORE INTO q_iks VALUES (?)", (ik,))
    db_rows = con.execute("""
        SELECT cs.standard_inchi_key, a.assay_id, a.standard_value
        FROM   chembl.compound_structures cs
        JOIN   chembl.activities a ON cs.molregno = a.molregno
        WHERE  cs.standard_inchi_key IN (SELECT ik FROM q_iks)
          AND  a.standard_type IN ('LogD', 'LogD7.4')
          AND  a.standard_value IS NOT NULL
    """).fetchall()
    con.close()
    results: dict[str, dict[str, list[float]]] = {}
    for ik, assay_id, val in db_rows:
        d = results.setdefault(ik, {"arm_a": [], "arm_b": []})
        d["arm_a"].append(float(val))
        if assay_id != SOURCE_ASSAY:
            d["arm_b"].append(float(val))
    return results


# ---------------------------------------------------------------------------
# 5. Per-compound absolute errors (test fold)
# ---------------------------------------------------------------------------
def compute_arm_b_errors(rows, train_idx, test_idx, iks, logd):
    train_ys  = np.array([rows[i][2] for i in train_idx])
    train_med = float(np.median(train_ys))
    n = len(test_idx)
    err = np.empty(n)
    cov = np.zeros(n, dtype=bool)
    for pos, i in enumerate(test_idx):
        _, _, y_true = rows[i]
        ik = iks[i]
        vb = logd.get(ik, {}).get("arm_b", []) if ik else []
        if vb:
            err[pos] = abs(y_true - float(np.median(vb)))
            cov[pos] = True
        else:
            err[pos] = abs(y_true - train_med)
    return err, cov


def compute_arm_c_errors(rows, train_idx, test_idx, k=K_NEIGH):
    train_fps, train_ys = [], []
    for i in train_idx:
        _, s, y = rows[i]
        mol = Chem.MolFromSmiles(s)
        fp  = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048) if mol is not None else None
        train_fps.append(fp); train_ys.append(y)
    valid = [(fp, y) for fp, y in zip(train_fps, train_ys) if fp is not None]
    fps_only = [fp for fp, _ in valid]
    ys_only  = [y  for _, y  in valid]
    train_mean = float(np.mean(train_ys))
    err = np.empty(len(test_idx))
    for pos, i in enumerate(test_idx):
        _, s, y_true = rows[i]
        mol = Chem.MolFromSmiles(s)
        if mol is None:
            err[pos] = abs(y_true - train_mean); continue
        test_fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
        sims = DataStructs.BulkTanimotoSimilarity(test_fp, fps_only)
        top_k = sorted(zip(sims, ys_only), key=lambda x: x[0], reverse=True)[:k]
        pred = float(np.mean([y for _, y in top_k]))
        err[pos] = abs(y_true - pred)
    return err


# ---------------------------------------------------------------------------
# 6. Seeded bootstrap
# ---------------------------------------------------------------------------
def bootstrap_ci(values: np.ndarray, B: int, seed: int) -> tuple[float, float, float]:
    rng = np.random.default_rng(seed)
    n = len(values)
    means = np.empty(B)
    for b in range(B):
        idx = rng.integers(0, n, size=n)
        means[b] = values[idx].mean()
    lo = float(np.percentile(means, 2.5))
    hi = float(np.percentile(means, 97.5))
    return float(values.mean()), lo, hi


def paired_bootstrap_ci(a: np.ndarray, b: np.ndarray, B: int, seed: int
                        ) -> tuple[float, float, float]:
    """Paired resample: one index draw per replicate; compute mean(a)-mean(b)."""
    assert len(a) == len(b)
    rng = np.random.default_rng(seed)
    n = len(a)
    diffs = np.empty(B)
    for b_i in range(B):
        idx = rng.integers(0, n, size=n)
        diffs[b_i] = a[idx].mean() - b[idx].mean()
    lo = float(np.percentile(diffs, 2.5))
    hi = float(np.percentile(diffs, 97.5))
    return float((a - b).mean()), lo, hi


# ---------------------------------------------------------------------------
# 7. Ledger helper (same row schema as evidence/README.md)
# ---------------------------------------------------------------------------
from datetime import datetime, timezone

# Deterministic timestamp: SOURCE_DATE_EPOCH if set (set it to pin), else
# a fixed sentinel so two consecutive runs produce byte-identical ledger rows.
import os
_sde = os.environ.get("SOURCE_DATE_EPOCH")
if _sde is not None:
    NOW = datetime.fromtimestamp(int(_sde), tz=timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
else:
    NOW = "2026-09-05T00:00:00Z"


def add_ledger(claim_id: str, value, source: str, reason: str,
               ci_low=None, ci_high=None, must_print: bool = False,
               metric: str = "value") -> dict:
    return {
        "claim_id":      claim_id,
        "metric":        metric,
        "value":         value,
        "source":        source,
        "verdict":       "SUPPORTED",
        "verdict_reason": reason,
        "ci_low":        ci_low,
        "ci_high":       ci_high,
        "must_print":    must_print,
        "timestamp":     NOW,
    }


def append_ledger(rows: list[dict]) -> None:
    """Append rows to the ledger, replacing any prior rows with matching claim_ids
    (so two consecutive runs are byte-idempotent)."""
    existing = []
    if LEDGER.exists():
        for line in LEDGER.read_text().splitlines():
            line = line.strip()
            if not line: continue
            try:
                row = json.loads(line)
            except json.JSONDecodeError:
                existing.append(line); continue
            existing.append(row)
    new_ids = {r["claim_id"] for r in rows}
    kept = [r for r in existing if not (isinstance(r, dict) and r.get("claim_id") in new_ids)]
    with LEDGER.open("w") as fh:
        for r in kept:
            fh.write(json.dumps(r) + "\n" if isinstance(r, dict) else (r + "\n"))
        for r in rows:
            fh.write(json.dumps(r) + "\n")


def upsert_manifest(entries: dict[str, dict]) -> None:
    data = json.loads(MANIFEST.read_text())
    data.setdefault("ledger_rows", {})
    for k, v in entries.items():
        data["ledger_rows"][k] = v
    # Keep formatting stable
    MANIFEST.write_text(json.dumps(data, indent=2) + "\n")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> int:
    print(f"[bootstrap_ci] parse benchmark {BENCHMARK}")
    rows = parse_benchmark()

    print("[bootstrap_ci] build TDC scaffold split (seed=42) ...")
    train_idx, val_idx, test_idx = build_scaffold_split(rows)
    assert len(test_idx) == 840, f"test fold {len(test_idx)} != 840"

    print("[bootstrap_ci] compute InChIKeys ...")
    iks = compute_inchikeys(rows)

    print("[bootstrap_ci] query ChEMBL 36 (READ-ONLY URI) ...")
    test_iks_set = {iks[i] for i in test_idx if iks[i] is not None}
    logd = query_logd(test_iks_set)

    print("[bootstrap_ci] compute Arm B per-compound errors ...")
    arm_b_err, arm_b_cov = compute_arm_b_errors(rows, train_idx, test_idx, iks, logd)
    assert arm_b_cov.sum() == 226, f"arm_b coverage {arm_b_cov.sum()} != 226"

    print("[bootstrap_ci] compute Arm C k=5 kNN per-compound errors ...")
    arm_c_err = compute_arm_c_errors(rows, train_idx, test_idx, k=K_NEIGH)

    # ----- HARD INTEGRITY GATE (STOP if any disagrees) -----
    b_full_mae = float(arm_b_err.mean())
    b_cov_mae  = float(arm_b_err[arm_b_cov].mean())
    c_mae      = float(arm_c_err.mean())

    def _gate(name, got, want, tol):
        # compare against the stored rounded value at its printed precision
        if abs(round(got, 4) - want) > tol and abs(round(got, 3) - want) > tol:
            print(f"[bootstrap_ci] INTEGRITY GATE FAILURE: {name} recomputed "
                  f"{got:.6f}, ledger stores {want}, diff {abs(got-want):.3e}. "
                  f"STOPPING; no paper changes made.", file=sys.stderr)
            sys.exit(2)
        print(f"[bootstrap_ci] integrity ok: {name}={got:.6f} matches ledger {want}")

    _gate("arm_b_full_mae",    b_full_mae, LEDGER_ARM_B_FULL_MAE,    TOL_MAE)
    _gate("arm_b_covered_mae", b_cov_mae,  LEDGER_ARM_B_COVERED_MAE, TOL_MAE)
    _gate("arm_c_k5_mae",      c_mae,      LEDGER_ARM_C_K5_MAE,      TOL_MAE)

    # ----- Bootstrap -----
    print(f"[bootstrap_ci] bootstrap B={B_REPLICATES}, seed={SEED_BOOTSTRAP} ...")
    b_full_mean, b_full_lo, b_full_hi = bootstrap_ci(arm_b_err, B_REPLICATES, SEED_BOOTSTRAP)
    b_cov_mean,  b_cov_lo,  b_cov_hi  = bootstrap_ci(arm_b_err[arm_b_cov], B_REPLICATES, SEED_BOOTSTRAP)
    c_mean,      c_lo,      c_hi      = bootstrap_ci(arm_c_err, B_REPLICATES, SEED_BOOTSTRAP)
    bmc_mean,    bmc_lo,    bmc_hi    = paired_bootstrap_ci(arm_b_err, arm_c_err, B_REPLICATES, SEED_BOOTSTRAP)

    print(f"  Arm B full:    MAE {b_full_mean:.4f}  95% CI [{b_full_lo:.4f}, {b_full_hi:.4f}]")
    print(f"  Arm B covered: MAE {b_cov_mean:.4f}  95% CI [{b_cov_lo:.4f}, {b_cov_hi:.4f}]")
    print(f"  Arm C k=5:     MAE {c_mean:.4f}  95% CI [{c_lo:.4f}, {c_hi:.4f}]")
    print(f"  Paired B-C:    {bmc_mean:+.4f}   95% CI [{bmc_lo:+.4f}, {bmc_hi:+.4f}]  "
          f"(includes 0: {'yes' if (bmc_lo <= 0.0 <= bmc_hi) else 'no'})")

    # ----- Emit evidence JSON -----
    out = {
        "script": "analysis/bootstrap_ci.py",
        "seed_split": SEED_SPLIT,
        "seed_bootstrap": SEED_BOOTSTRAP,
        "B_replicates": B_REPLICATES,
        "inputs": {
            "benchmark": str(BENCHMARK),
            "chembl_db_uri": DB_URI,
        },
        "test_fold_n": len(test_idx),
        "arm_b_covered_n": int(arm_b_cov.sum()),
        "integrity_gate": {
            "arm_b_full_mae":    {"recomputed": round(b_full_mae, 6), "ledger": LEDGER_ARM_B_FULL_MAE,    "pass": True},
            "arm_b_covered_mae": {"recomputed": round(b_cov_mae, 6),  "ledger": LEDGER_ARM_B_COVERED_MAE, "pass": True},
            "arm_c_k5_mae":      {"recomputed": round(c_mae, 6),      "ledger": LEDGER_ARM_C_K5_MAE,      "pass": True},
        },
        "arm_b_full":    {"mae": round(b_full_mean, 6), "ci_low": round(b_full_lo, 4), "ci_high": round(b_full_hi, 4), "n": len(arm_b_err)},
        "arm_b_covered": {"mae": round(b_cov_mean,  6), "ci_low": round(b_cov_lo,  4), "ci_high": round(b_cov_hi,  4), "n": int(arm_b_cov.sum())},
        "arm_c_k5":      {"mae": round(c_mean,      6), "ci_low": round(c_lo,      4), "ci_high": round(c_hi,      4), "n": len(arm_c_err)},
        "paired_b_minus_c": {
            "mean_diff": round(bmc_mean, 6),
            "ci_low": round(bmc_lo, 4),
            "ci_high": round(bmc_hi, 4),
            "includes_zero": bool(bmc_lo <= 0.0 <= bmc_hi),
        },
    }
    OUT_JSON.write_text(json.dumps(out, indent=2) + "\n")
    print(f"[bootstrap_ci] wrote {OUT_JSON}")

    # ----- Append ledger + upsert manifest -----
    #
    # Values printed in main.tex must match the manifest exactly. We stamp the
    # CI endpoints to 4dp. The paired B-C endpoints may be negative; the
    # manifest carries signed floats.
    src = "analysis/bootstrap_ci.py  ->  evidence/t1_bootstrap_cis.json"
    reason_common = (f"seeded bootstrap B={B_REPLICATES} (seed={SEED_BOOTSTRAP}) "
                     "over per-compound absolute errors recomputed from benchmark "
                     "tab + ChEMBL 36 read-only URI; integrity gate matched stored "
                     "ledger MAEs to printed precision.")
    b_full_lo_r = round(b_full_lo, 4); b_full_hi_r = round(b_full_hi, 4)
    b_cov_lo_r  = round(b_cov_lo,  4); b_cov_hi_r  = round(b_cov_hi,  4)
    c_lo_r      = round(c_lo,      4); c_hi_r      = round(c_hi,      4)
    bmc_lo_r    = round(bmc_lo,    4); bmc_hi_r    = round(bmc_hi,    4)

    ledger_rows = [
        add_ledger("t1::ci_b_full_lo",    b_full_lo_r, src, reason_common + " Arm B full-fold 2.5%",
                   ci_low=b_full_lo_r, ci_high=b_full_hi_r, must_print=True, metric="mae_ci_lo"),
        add_ledger("t1::ci_b_full_hi",    b_full_hi_r, src, reason_common + " Arm B full-fold 97.5%",
                   ci_low=b_full_lo_r, ci_high=b_full_hi_r, must_print=True, metric="mae_ci_hi"),
        add_ledger("t1::ci_b_covered_lo", b_cov_lo_r,  src, reason_common + " Arm B covered subset 2.5%",
                   ci_low=b_cov_lo_r,  ci_high=b_cov_hi_r,  must_print=False, metric="mae_ci_lo"),
        add_ledger("t1::ci_b_covered_hi", b_cov_hi_r,  src, reason_common + " Arm B covered subset 97.5%",
                   ci_low=b_cov_lo_r,  ci_high=b_cov_hi_r,  must_print=False, metric="mae_ci_hi"),
        add_ledger("t1::ci_c_lo",         c_lo_r,      src, reason_common + " Arm C k=5 2.5%",
                   ci_low=c_lo_r,      ci_high=c_hi_r,      must_print=True, metric="mae_ci_lo"),
        add_ledger("t1::ci_c_hi",         c_hi_r,      src, reason_common + " Arm C k=5 97.5%",
                   ci_low=c_lo_r,      ci_high=c_hi_r,      must_print=True, metric="mae_ci_hi"),
        add_ledger("t1::ci_bmc_lo",       bmc_lo_r,    src, reason_common + " paired B-C difference 2.5%",
                   ci_low=bmc_lo_r,    ci_high=bmc_hi_r,    must_print=True, metric="mae_diff_ci_lo"),
        add_ledger("t1::ci_bmc_hi",       bmc_hi_r,    src, reason_common + " paired B-C difference 97.5%",
                   ci_low=bmc_lo_r,    ci_high=bmc_hi_r,    must_print=True, metric="mae_diff_ci_hi"),
    ]
    append_ledger(ledger_rows)
    print(f"[bootstrap_ci] appended/updated {len(ledger_rows)} ledger rows")

    manifest_entries = {
        "t1::ci_b_full_lo":    {"value": b_full_lo_r, "claim_id": "t1::ci_b_full_lo",
                                "source": src, "must_print": True},
        "t1::ci_b_full_hi":    {"value": b_full_hi_r, "claim_id": "t1::ci_b_full_hi",
                                "source": src, "must_print": True},
        "t1::ci_b_covered_lo": {"value": b_cov_lo_r,  "claim_id": "t1::ci_b_covered_lo",
                                "source": src, "must_print": False},
        "t1::ci_b_covered_hi": {"value": b_cov_hi_r,  "claim_id": "t1::ci_b_covered_hi",
                                "source": src, "must_print": False},
        "t1::ci_c_lo":         {"value": c_lo_r,      "claim_id": "t1::ci_c_lo",
                                "source": src, "must_print": True},
        "t1::ci_c_hi":         {"value": c_hi_r,      "claim_id": "t1::ci_c_hi",
                                "source": src, "must_print": True},
        "t1::ci_bmc_lo":       {"value": bmc_lo_r,    "claim_id": "t1::ci_bmc_lo",
                                "source": src, "must_print": True},
        "t1::ci_bmc_hi":       {"value": bmc_hi_r,    "claim_id": "t1::ci_bmc_hi",
                                "source": src, "must_print": True},
    }
    upsert_manifest(manifest_entries)
    print(f"[bootstrap_ci] upserted {len(manifest_entries)} manifest entries -> {MANIFEST}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
