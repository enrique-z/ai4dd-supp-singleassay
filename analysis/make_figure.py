#!/usr/bin/env python3
"""Figure 1 — label distribution and lookup performance (two panels).

Panel (a): Histogram of all 4,200 logD values (bin width 0.05) with
vertical rules at the instrument-declared bounds −1.5 and +4.5.  The
visual argument: the distribution runs smoothly up to both bounds and
then stops dead.  A clipping process would pile mass ON the bound; the
data show no spike — compounds outside range were excluded, not clipped.
An inset zooms the upper tail where the argument is most critical.

Panel (b): ECDF of |error| on the 840-compound test fold for two arms:
  ARM A — lookup including the source assay (CHEMBL3301363 / assay 1369278)
  ARM B — lookup excluding the source assay

Recompute-then-assert pattern:
Every aggregate this figure draws is recomputed from first principles and
asserted against evidence/t1_canonical.json before anything is rendered.
The build fails loudly rather than producing a figure that disagrees with
the prose.  pdf.fonttype=42 prevents Type-3 font desk-rejections.

Usage:
    python3 analysis/make_figure.py          # normal build → fig1.pdf
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

import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["pdf.fonttype"] = 42   # never Type 3 — NeurIPS desk-reject
matplotlib.rcParams["ps.fonttype"]  = 42
import matplotlib.pyplot as plt

from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import AllChem
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem.inchi import MolToInchi, InchiToInchiKey

RDLogger.DisableLog("rdApp.*")

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
HERE       = Path(__file__).resolve().parent
PAPER_DIR  = HERE.parent
EVIDENCE   = PAPER_DIR / "evidence"
T1_JSON        = EVIDENCE / "t1_canonical.json"
T1_ARM_C_JSON  = EVIDENCE / "t1_arm_c_knn.json"
OUT            = PAPER_DIR / "fig1.pdf"
CLAIMS_OUT = EVIDENCE / "figure_claims.json"

BENCHMARK = PAPER_DIR / "data" / "lipophilicity_astrazeneca.tab"
DB_PATH   = PAPER_DIR / "data" / "chembl_36" / "chembl_36_sqlite" / "chembl_36.db"

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
SEED         = 42
FRAC         = [0.7, 0.1, 0.2]
SOURCE_ASSAY = 1369278           # CHEMBL3301363 assay_id

# Assertion tolerances
MAE_TOL      = 1e-6   # tight: MAE values must reproduce to six decimal places
ZERO_TOL     = 0.0    # integer counts: exact equality (stored as float, compared with 0.0)
ARM_C_MAE_TOL = 1e-3  # kNN has minor tie-handling variation; 1e-3 is justified

# Colourblind-safe palette (IBM):  blue / magenta / orange / charcoal hist fill
C_ARM_A = "#648FFF"   # ARM A (including source assay)
C_ARM_B = "#DC267F"   # ARM B (excluding source assay)
C_WALL  = "#FE6100"   # instrument bounds annotation
C_HIST  = "#3D3D3D"   # histogram fill
C_ZERO  = "#18A999"   # ECDF label accent (teal — complements the two main hues)

# ---------------------------------------------------------------------------
# Step 1: parse benchmark
# ---------------------------------------------------------------------------

def parse_benchmark() -> list[tuple[str, str, float]]:
    """Return [(drug_id, smiles, Y)] for all rows.

    The file uses TAB separator and wraps Drug_ID and Drug in literal
    double-quotes; parse with csv.reader to strip them correctly.
    """
    rows: list[tuple[str, str, float]] = []
    with BENCHMARK.open() as f:
        reader = csv.reader(f, delimiter="\t", quotechar='"')
        next(reader)                # skip header
        for r in reader:
            rows.append((r[0], r[1], float(r[2])))
    return rows


# ---------------------------------------------------------------------------
# Step 2: TDC scaffold split (exact replica of TDC 1.1.15)
# ---------------------------------------------------------------------------

def build_scaffold_split(
    rows: list[tuple[str, str, float]],
) -> tuple[list[int], list[int], list[int]]:
    """Return (train_idx, val_idx, test_idx) for the canonical TDC split.

    Replicates create_scaffold_split from tdc/utils/split.py verbatim:
      • MurckoScaffoldSmiles(mol=MolFromSmiles(s), includeChirality=False)
      • Random(seed) instance; re-seeded with same seed before first shuffle
      • Big threshold: size > val_size/2 OR > test_size/2
      • Greedy fill: train → val → remainder to test
    """
    n          = len(rows)
    train_size = int(n * FRAC[0])           # 2940
    val_size   = int(n * FRAC[1])           # 420
    test_size  = n - train_size - val_size  # 840

    scaffolds: defaultdict[str, set] = defaultdict(set)
    for i, (_, smiles, _) in enumerate(rows):
        mol = Chem.MolFromSmiles(smiles)
        sca = (
            ""
            if mol is None
            else MurckoScaffold.MurckoScaffoldSmiles(mol=mol, includeChirality=False)
        )
        scaffolds[sca].add(i)

    big, small = [], []
    for idx_set in scaffolds.values():
        if len(idx_set) > val_size / 2 or len(idx_set) > test_size / 2:
            big.append(idx_set)
        else:
            small.append(idx_set)

    rng = Random(SEED)
    rng.seed(SEED)          # TDC re-seeds the same instance before first shuffle
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
# Step 3: compute InChIKeys
# ---------------------------------------------------------------------------

def compute_inchikeys(rows: list[tuple[str, str, float]]) -> list[str | None]:
    """Return standard InChIKey for each SMILES row (None on parse failure)."""
    iks: list[str | None] = []
    for _, smiles, _ in rows:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            iks.append(None)
            continue
        inchi = MolToInchi(mol)
        iks.append(None if inchi is None else InchiToInchiKey(inchi))
    return iks


# ---------------------------------------------------------------------------
# Step 4: query ChEMBL
# ---------------------------------------------------------------------------

def query_logd(
    test_iks_set: set[str],
) -> dict[str, dict[str, list[float]]]:
    """Return {ik: {'arm_a': [vals], 'arm_b': [vals]}} for matched IKs.

    Uses :memory: + ATTACH so the 96%-full disk is never touched beyond
    the read-only open.  standard_relation is not filtered (as per spec).
    """
    con = sqlite3.connect("file::memory:", uri=True)
    con.execute(f"ATTACH DATABASE 'file:{DB_PATH}?mode=ro' AS chembl")
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
        if ik not in results:
            results[ik] = {"arm_a": [], "arm_b": []}
        results[ik]["arm_a"].append(float(val))
        if assay_id != SOURCE_ASSAY:
            results[ik]["arm_b"].append(float(val))

    return results


# ---------------------------------------------------------------------------
# Step 5: compute per-compound errors
# ---------------------------------------------------------------------------

def compute_errors(
    rows:       list[tuple[str, str, float]],
    train_idx:  list[int],
    test_idx:   list[int],
    iks:        list[str | None],
    logd:       dict[str, dict[str, list[float]]],
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return (arm_a_err, arm_a_cov, arm_b_err, arm_b_cov) for test fold.

    Prediction = median of retrieved logD values;
    fallback    = exact train-fold median (computed fresh, not rounded).
    """
    train_ys  = np.array([rows[i][2] for i in train_idx])
    train_med = float(np.median(train_ys))

    n = len(test_idx)
    arm_a_err = np.empty(n)
    arm_b_err = np.empty(n)
    arm_a_cov = np.zeros(n, dtype=bool)
    arm_b_cov = np.zeros(n, dtype=bool)

    for pos, i in enumerate(test_idx):
        _, _, y_true = rows[i]
        ik   = iks[i]
        lkup = logd.get(ik, {}) if ik is not None else {}

        va = lkup.get("arm_a", [])
        if va:
            arm_a_err[pos]  = abs(y_true - float(np.median(va)))
            arm_a_cov[pos]  = True
        else:
            arm_a_err[pos]  = abs(y_true - train_med)

        vb = lkup.get("arm_b", [])
        if vb:
            arm_b_err[pos]  = abs(y_true - float(np.median(vb)))
            arm_b_cov[pos]  = True
        else:
            arm_b_err[pos]  = abs(y_true - train_med)

    return arm_a_err, arm_a_cov, arm_b_err, arm_b_cov


# ---------------------------------------------------------------------------
# Step 5b: ARM C — k-NN within training fold only (no DB access)
# ---------------------------------------------------------------------------

def compute_arm_c_errors(
    rows:       list[tuple[str, str, float]],
    train_idx:  list[int],
    test_idx:   list[int],
    k:          int = 5,
) -> np.ndarray:
    """Predict each test compound as the mean label of its k nearest train neighbours.

    Fingerprint: ECFP4 (radius=2, nBits=2048).
    Similarity:  Tanimoto via DataStructs.BulkTanimotoSimilarity.
    Coverage:    840/840 by construction (fallback to training mean on parse error).
    """
    # Build training fingerprints and labels in one pass
    train_fps: list = []
    train_ys:  list[float] = []
    for i in train_idx:
        _, smiles, y = rows[i]
        mol = Chem.MolFromSmiles(smiles)
        fp  = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048) if mol is not None else None
        train_fps.append(fp)
        train_ys.append(y)

    valid_train = [(fp, y) for fp, y in zip(train_fps, train_ys) if fp is not None]
    fps_only    = [fp for fp, _ in valid_train]
    ys_only     = [y  for _, y  in valid_train]
    train_mean  = float(np.mean(train_ys))   # fallback for unparse-able test mols

    arm_c_err = np.empty(len(test_idx))
    for pos, i in enumerate(test_idx):
        _, smiles, y_true = rows[i]
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            arm_c_err[pos] = abs(y_true - train_mean)
            continue
        test_fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
        sims    = DataStructs.BulkTanimotoSimilarity(test_fp, fps_only)
        # k nearest by descending similarity
        top_k   = sorted(zip(sims, ys_only), key=lambda x: x[0], reverse=True)[:k]
        pred    = float(np.mean([y for _, y in top_k]))
        arm_c_err[pos] = abs(y_true - pred)

    return arm_c_err


# ---------------------------------------------------------------------------
# Step 6: assert against evidence
# ---------------------------------------------------------------------------

def _chk(name: str, got: float, want: float, tol: float) -> None:
    """Hard-abort if recomputed value disagrees with evidence."""
    if abs(got - want) > tol:
        sys.exit(
            f"[assert] {name}: recomputed {got:.10g}, "
            f"evidence {want:.10g}, diff {abs(got-want):.3e} > tol {tol:.1e}. "
            f"Evidence and figure are out of sync — aborting build."
        )
    print(f"[assert] {name}: {got:.8g} == {want:.8g}  OK")


def run_assertions(
    rows:       list[tuple[str, str, float]],
    arm_a_err:  np.ndarray,
    arm_a_cov:  np.ndarray,
    arm_b_err:  np.ndarray,
    arm_b_cov:  np.ndarray,
    arm_c_err:  np.ndarray,
    t1:         dict,
    t1_arm_c:   dict,
) -> list[str]:
    """Assert all required values; return list of asserted claim_ids."""
    task3 = t1["task3_arm_results"]
    arm_a = task3["arm_a_including_1369278"]
    arm_b = task3["arm_b_excluding_1369278"]

    # --- ARM A ---
    _chk("arm_a_coverage_n",
         float(arm_a_cov.sum()),
         float(arm_a["covered_n"]),
         ZERO_TOL)
    _chk("arm_a_full_fold_mae",
         float(arm_a_err.mean()),
         arm_a["full_fold_metrics"]["mae"],
         MAE_TOL)

    # --- ARM B ---
    _chk("arm_b_coverage_n",
         float(arm_b_cov.sum()),
         float(arm_b["covered_n"]),
         ZERO_TOL)
    _chk("arm_b_full_fold_mae",
         float(arm_b_err.mean()),
         arm_b["full_fold_metrics"]["mae"],
         MAE_TOL)

    # --- ARM C (kNN, k=5) ---
    want_c_mae = t1_arm_c["arm_c_results_all_k"]["5"]["mae"]
    _chk("arm_c_k5_full_fold_mae",
         float(arm_c_err.mean()),
         want_c_mae,
         ARM_C_MAE_TOL)

    # --- Benchmark distribution ---
    ys   = np.array([r[2] for r in rows])
    meta = t1["meta"]
    _chk("benchmark_n_rows",    float(len(rows)),           float(meta["n_rows"]),  ZERO_TOL)
    _chk("benchmark_min",       float(ys.min()),            -1.5,                   ZERO_TOL)
    _chk("benchmark_max",       float(ys.max()),             4.5,                   ZERO_TOL)
    _chk("benchmark_n_below_lo", float((ys < -1.5).sum()),  0.0,                   ZERO_TOL)
    _chk("benchmark_n_above_hi", float((ys >  4.5).sum()),  0.0,                   ZERO_TOL)
    _chk("benchmark_n_at_lo",   float((ys == -1.5).sum()),  1.0,                   ZERO_TOL)
    _chk("benchmark_n_at_hi",   float((ys ==  4.5).sum()),  4.0,                   ZERO_TOL)

    return [
        "t1::arm_a_coverage_n", "t1::arm_a_full_mae",
        "t1::arm_b_coverage_n", "t1::arm_b_full_mae",
        "t1::arm_c_k5_full_mae",
        "f4::val_min", "f4::val_max",
        "f4::n_below_lo", "f4::n_above_hi",
        "f4::n_at_lo", "f4::n_at_hi",
    ]


# ---------------------------------------------------------------------------
# Step 7: draw
# ---------------------------------------------------------------------------

def draw(
    rows:      list[tuple[str, str, float]],
    arm_a_err: np.ndarray,
    arm_a_cov: np.ndarray,
    arm_b_err: np.ndarray,
    arm_b_cov: np.ndarray,
    arm_c_err: np.ndarray,
) -> None:
    """Render the two-panel figure to fig1.pdf (5.5 × 2.6 in)."""
    import matplotlib.gridspec as gridspec
    from matplotlib.patches import Rectangle

    ys     = np.array([r[2] for r in rows])
    n_test = len(arm_a_err)

    # Pre-compute annotation counts from the data (exact arithmetic)
    n_at_lo  = int((ys == -1.5).sum())   # 1
    n_at_hi  = int((ys ==  4.5).sum())   # 4

    # ARM A/B coverage counts (used in legend labels)
    n_a_cov  = int(arm_a_cov.sum())             # 834
    n_b_cov  = int(arm_b_cov.sum())             # 226

    mae_a = float(arm_a_err.mean())
    mae_b = float(arm_b_err.mean())
    mae_c = float(arm_c_err.mean())

    # --- Figure skeleton: panel (a) gets 60 % of inner width ---
    fig = plt.figure(figsize=(5.5, 2.6))
    gs = gridspec.GridSpec(
        1, 2, figure=fig,
        width_ratios=[1.5, 1.0],
        left=0.085, right=0.985, bottom=0.145, top=0.965, wspace=0.38,
    )
    axA = fig.add_subplot(gs[0])
    axB = fig.add_subplot(gs[1])

    # ================================================================
    # PANEL (a): histogram
    # ================================================================
    # Use exact bin edges anchored at -1.60 to avoid float accumulation
    # drift that can shift where exactly-representable values land.
    lo_edge = round(-1.60, 2)
    hi_edge = round( 4.60, 2)
    n_bins  = round((hi_edge - lo_edge) / 0.05)
    bins    = np.linspace(lo_edge, hi_edge, n_bins + 1)

    counts, edges = np.histogram(ys, bins=bins)
    bar_centers   = 0.5 * (edges[:-1] + edges[1:])

    # Main histogram — charcoal bars
    axA.bar(bar_centers, counts,
            width=0.048, color=C_HIST, linewidth=0, zorder=2)

    # Honest axis limits — just enough margin past each wall, no fake space
    X_MIN, X_MAX = -1.85, 4.80
    y_ceil = float(counts.max()) * 1.08   # just above the tallest bar
    axA.set_xlim(X_MIN, X_MAX)
    axA.set_ylim(0, y_ceil)

    # Instrument bound walls — identically styled, full-height, both visible
    for x in (-1.5, 4.5):
        axA.axvline(x, color=C_WALL, lw=1.3, ls="--", zorder=4)

    # ----------------------------------------------------------------
    # Inset: upper tail — placed in the upper-LEFT where the low-logD
    # tail is sparse.  Left edge starts right of the −1.5 wall so both
    # walls remain fully unoccluded in the main axes.
    # ----------------------------------------------------------------
    # Inset bounding box x coverage in data space: [INS_X_LEFT_D, INS_X_RIGHT_D]
    INS_X_LEFT_D  = -1.30   # safely right of the −1.5 wall
    INS_X_RIGHT_D =  0.60
    ul_mask   = (bar_centers >= INS_X_LEFT_D) & (bar_centers <= INS_X_RIGHT_D)
    ul_max    = float(counts[ul_mask].max()) if ul_mask.any() else 0.0
    INS_Y_BOT_D = ul_max + 5.0   # buffer above the tallest bar in this region

    # RUN-TIME GUARD: abort if any bar in the inset footprint would be hidden
    bad = counts[ul_mask] >= INS_Y_BOT_D
    if bad.any():
        sys.exit(
            f"[draw] inset bar-overlap guard failed: bars at "
            f"logD {bar_centers[ul_mask][bad]} have height "
            f"{counts[ul_mask][bad]} >= y_bot {INS_Y_BOT_D:.1f}. "
            f"Adjust INS_X_LEFT_D / INS_X_RIGHT_D or INS_Y_BOT_D."
        )

    x_rng      = X_MAX - X_MIN                           # 6.65
    x_ins_l_ax = (INS_X_LEFT_D  - X_MIN) / x_rng        # ≈ 0.083
    x_ins_r_ax = (INS_X_RIGHT_D - X_MIN) / x_rng        # ≈ 0.368
    y_ins_b_ax = INS_Y_BOT_D / y_ceil                    # ≈ 0.318
    w_ins_ax   = x_ins_r_ax - x_ins_l_ax                 # ≈ 0.285
    h_ins_ax   = 0.95 - y_ins_b_ax                       # ≈ 0.632

    ins = axA.inset_axes([x_ins_l_ax, y_ins_b_ax, w_ins_ax, h_ins_ax])

    # Upper-tail bars: reuse the global bins, filter to [4.28, 4.62]
    tail_mask = (bar_centers >= 4.28) & (bar_centers <= 4.62)
    t_centers = bar_centers[tail_mask]
    t_counts  = counts[tail_mask]
    ins.bar(t_centers, t_counts, width=0.046,
            color=C_HIST, linewidth=0, zorder=2)

    # Upper wall inside inset — clip_on prevents it escaping the inset box
    vl = ins.axvline(4.5, color=C_WALL, lw=1.0, ls="--", zorder=4)
    vl.set_clip_on(True)

    # Count labels on every non-empty bar in the tail.
    for cx, ct in zip(t_centers, t_counts):
        if ct > 0:
            lbl_color = C_WALL if cx > 4.48 else "0.35"
            ins.text(cx, ct + 0.35, str(int(ct)),
                     ha="center", va="bottom", fontsize=7.0,
                     color=lbl_color, zorder=6)

    ins.set_xlim(4.28, 4.62)
    ins.set_ylim(0, float(t_counts.max()) * 1.50)
    # The bar at 4.5 (orange count label "4") previously sat directly above
    # the "4.5" x-tick text, creating visual clutter at print size. Drop the
    # numeric tick labels — the "upper tail" title and the orange wall at 4.5
    # (labelled "n=4" outside the inset) already fix the axis reference.
    ins.set_xticks([4.4, 4.5])
    ins.set_xticklabels([])
    ins.tick_params(axis="y", labelsize=6.0, length=2)
    ins.tick_params(axis="x", length=2, pad=1)
    # Full rectangular frame — unambiguously a separate panel
    for sp in ins.spines.values():
        sp.set_visible(True)
        sp.set_linewidth(0.5)
        sp.set_color("0.55")
    ins.set_title("upper tail", fontsize=7.0, pad=2, color="0.45")

    # Dotted rectangle in the main axes marking the zoomed upper-tail region
    zoom_max = float(counts[(bar_centers >= 4.295) & (bar_centers <= 4.555)].max())
    zoom_rect = Rectangle(
        (4.295, 0), 0.265, zoom_max + 2.5,
        linewidth=0.55, edgecolor="0.62", facecolor="none",
        linestyle=":", zorder=5,
    )
    axA.add_patch(zoom_rect)

    # Wall count annotations — legible height, one per wall, no arrowprops
    axA.text(-1.40, 24.0, f"n={n_at_lo}",
             fontsize=7.5, color=C_WALL, ha="left", va="bottom",
             fontweight="semibold", zorder=5)
    axA.text( 4.38, 24.0, f"n={n_at_hi}",
             fontsize=7.5, color=C_WALL, ha="right", va="bottom",
             fontweight="semibold", zorder=5)

    # Axes dressing
    axA.set_xlabel("Benchmark label (logD7.4)", fontsize=7.5, labelpad=3)
    axA.set_ylabel("Compounds", fontsize=7.5, labelpad=3)
    axA.tick_params(labelsize=6.5)
    axA.text(0.02, 0.97, "(a)", transform=axA.transAxes,
             fontsize=7.5, va="top", fontweight="bold", color="0.25")
    for s in ("top", "right"):
        axA.spines[s].set_visible(False)

    # ================================================================
    # PANEL (b): ECDF of |error|
    # ================================================================
    def ecdf(arr: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        xs   = np.sort(arr)
        ys_e = np.arange(1, len(xs) + 1) / len(xs)
        return xs, ys_e

    xs_a, ys_a = ecdf(arm_a_err)
    xs_b, ys_b = ecdf(arm_b_err)
    xs_c, ys_c = ecdf(arm_c_err)

    # Draw ECDFs — prepend (0, 0) for a clean left edge.
    # MAE values are embedded in the legend labels; no floating labels or
    # dotted MAE lines (three of each collide in the same x-region).
    axB.step(np.r_[0.0, xs_a], np.r_[0.0, ys_a],
             where="post", color=C_ARM_A, lw=1.5, zorder=4,
             label=f"ARM A-pool  (n={n_a_cov}/840, MAE {mae_a:.3f})")
    axB.step(np.r_[0.0, xs_b], np.r_[0.0, ys_b],
             where="post", color=C_ARM_B, lw=1.5, zorder=4,
             label=f"ARM B  (n={n_b_cov}/840, MAE {mae_b:.3f})")
    axB.step(np.r_[0.0, xs_c], np.r_[0.0, ys_c],
             where="post", color=C_ZERO, lw=1.5, zorder=4,
             label=f"ARM C  (n=840/840, MAE {mae_c:.3f})")

    # "689 exact" annotation dropped: the near-vertical rise of the ARM A
    # curve at |error|=0 already makes the point, and the count is stated
    # in the prose.  With three curves the annotation sits on the lines.

    # Legend — MAE values live here; loc=upper left keeps it in the open
    # region above x≈1 where all curves have reached 1.0
    leg = axB.legend(
        fontsize=5.0, loc="lower right", framealpha=0.0,
        handlelength=1.4, labelspacing=0.35,
    )
    for text in leg.get_texts():
        text.set_color("0.25")

    axB.set_xlabel("|error|  (logD units)", fontsize=7.5, labelpad=3)
    axB.set_ylabel("Cumulative fraction", fontsize=7.5, labelpad=3)
    axB.tick_params(labelsize=6.5)
    axB.set_xlim(left=-0.03)
    axB.set_ylim(0, 1.04)
    axB.text(0.02, 0.97, "(b)", transform=axB.transAxes,
             fontsize=7.5, va="top", fontweight="bold", color="0.25")
    for s in ("top", "right"):
        axB.spines[s].set_visible(False)

    # --- save ---
    fig.savefig(OUT, bbox_inches="tight", dpi=400)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> int:
    print("[figure] loading evidence …")
    if not T1_JSON.exists():
        sys.exit(f"[figure] evidence not found: {T1_JSON}")
    t1 = json.loads(T1_JSON.read_text())

    if not T1_ARM_C_JSON.exists():
        sys.exit(f"[figure] ARM C evidence not found: {T1_ARM_C_JSON}")
    t1_arm_c = json.loads(T1_ARM_C_JSON.read_text())

    print("[figure] parsing benchmark …")
    rows = parse_benchmark()
    if len(rows) != 4200:
        sys.exit(f"[figure] expected 4200 rows, got {len(rows)}")

    print("[figure] building TDC scaffold split …")
    train_idx, val_idx, test_idx = build_scaffold_split(rows)
    print(f"[figure]   train={len(train_idx)}  val={len(val_idx)}  test={len(test_idx)}")

    print("[figure] computing InChIKeys …")
    iks = compute_inchikeys(rows)

    print("[figure] querying ChEMBL for test-fold logD values …")
    test_iks_set = {iks[i] for i in test_idx if iks[i] is not None}
    logd = query_logd(test_iks_set)

    print("[figure] computing per-compound prediction errors …")
    arm_a_err, arm_a_cov, arm_b_err, arm_b_cov = compute_errors(
        rows, train_idx, test_idx, iks, logd
    )

    print("[figure] computing ARM C (kNN, k=5) errors …")
    arm_c_err = compute_arm_c_errors(rows, train_idx, test_idx, k=5)
    print(f"[figure]   ARM C MAE = {arm_c_err.mean():.6f}  (n={len(arm_c_err)}/840)")

    print("[figure] running assertions …")
    claim_ids = run_assertions(
        rows, arm_a_err, arm_a_cov, arm_b_err, arm_b_cov, arm_c_err, t1, t1_arm_c
    )
    print("[figure] all assertions passed ✓")

    # Report exact-zero count (used in prose)
    n_zero_a = int((arm_a_err == 0.0).sum())
    print(f"[figure] ARM A |error|==0: {n_zero_a}/{len(test_idx)} "
          f"({100*n_zero_a/len(test_idx):.1f}% of test fold)")

    print(f"[figure] drawing → {OUT} …")
    draw(rows, arm_a_err, arm_a_cov, arm_b_err, arm_b_cov, arm_c_err)
    print(f"[figure] wrote {OUT}  (5.5 × 2.6 in, vector PDF)")

    # Declare drawn claim_ids so the guard counts them as reported
    CLAIMS_OUT.write_text(
        json.dumps(
            {"figure": OUT.name, "claim_ids": sorted(set(claim_ids))},
            indent=1,
        ) + "\n"
    )
    print(f"[figure] declared {len(set(claim_ids))} claim_ids → {CLAIMS_OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
