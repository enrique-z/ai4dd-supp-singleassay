"""Test-fold identity and the six Arm-A-source misses (camera-ready additions).

Checks that the TDC scaffold split (seed 42) test fold equals the TDC ADMET
benchmark group's fixed test set, reproduces Arm A-source on that fold against
ChEMBL 36 assay CHEMBL3301363, and characterises the uncovered test compounds.

Usage (PyTDC 1.1.15, RDKit, pandas, numpy):
  python3 testfold_identity_and_misses.py --tab-dir DIR_WITH_lipophilicity_astrazeneca.tab \
      --group-dir ADMET_GROUP/lipophilicity_astrazeneca --db path/to/chembl_36.db --out OUT_DIR
"""
import argparse, json, os, sqlite3, statistics, sys, inspect
from pathlib import Path
import pandas as pd
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("--tab-dir", required=True)
ap.add_argument("--group-dir", required=True)
ap.add_argument("--db", required=True)
ap.add_argument("--out", required=True)
args = ap.parse_args()
W = Path(args.out); W.mkdir(parents=True, exist_ok=True)
TAB_DIR = Path(args.tab_dir)
GROUP = Path(args.group_dir)
DB = "file:" + os.path.abspath(args.db) + "?mode=ro&immutable=1"
ASSAY_CHEMBL = "CHEMBL3301363"

import tdc
from tdc.single_pred import ADME
from tdc.utils import split as tdc_split

out = {}
try:
    from importlib.metadata import version
    out["pytdc_version"] = version("PyTDC")
except Exception as e:
    out["pytdc_version"] = str(e)

# ---- 1. Real TDC loader + get_split(scaffold, seed=42, frac=[0.7,0.1,0.2])
data = ADME(name="Lipophilicity_AstraZeneca", path=str(TAB_DIR))
df_all = data.get_data()
split = data.get_split(method="scaffold", seed=42, frac=[0.7, 0.1, 0.2])
tr, va, te = split["train"], split["valid"], split["test"]
out["n_all"] = len(df_all)
out["fold_sizes"] = {"train": len(tr), "valid": len(va), "test": len(te)}
out["train_median_y"] = float(np.median(tr["Y"]))
out["trainval_median_y"] = float(np.median(pd.concat([tr, va])["Y"]))

grp_test = pd.read_csv(GROUP / "test.csv")
grp_tv = pd.read_csv(GROUP / "train_val.csv")
s_split = set(te["Drug_ID"]); s_grp = set(grp_test["Drug_ID"])
out["q1"] = {
    "group_test_n": len(grp_test), "group_train_val_n": len(grp_tv),
    "split_test_unique_ids": len(s_split), "group_test_unique_ids": len(s_grp),
    "intersection": len(s_split & s_grp),
    "only_in_split_test": sorted(s_split - s_grp)[:20],
    "only_in_group_test": sorted(s_grp - s_split)[:20],
    "group_trainval_eq_split_train_plus_valid":
        set(grp_tv["Drug_ID"]) == (set(tr["Drug_ID"]) | set(va["Drug_ID"])),
}
# labels identical too?
m = te.merge(grp_test, on="Drug_ID", suffixes=("_s", "_g"))
out["q1"]["labels_identical"] = bool(np.allclose(m["Y_s"], m["Y_g"])) if len(m) else None
out["q1"]["smiles_identical"] = bool((m["Drug_s"] == m["Drug_g"]).all()) if len(m) else None

# seed effect: does get_train_valid_split with different seeds touch test? test is read from file.
from tdc.benchmark_group import base_group
src = inspect.getsource(base_group.BenchmarkGroup.get_train_valid_split)
out["q1"]["get_train_valid_split_reads"] = [l.strip() for l in src.splitlines() if "read_csv" in l or "frac =" in l]

# also: do test sets vary across seeds in a plain get_split? (sanity: seed matters for get_split)
te43 = data.get_split(method="scaffold", seed=43, frac=[0.7, 0.1, 0.2])["test"]
out["q1"]["get_split_seed43_test_overlap_with_seed42"] = len(set(te43["Drug_ID"]) & s_split)

# ---- 2. Arm A-source
con = sqlite3.connect(DB, uri=True)
cur = con.cursor()
assay_id = cur.execute("select assay_id from assays where chembl_id=?", (ASSAY_CHEMBL,)).fetchone()[0]
out["assay_id"] = assay_id
rows = cur.execute("""
  select md.chembl_id, act.molregno, act.standard_value, act.value, act.standard_type
  from activities act join molecule_dictionary md on md.molregno = act.molregno
  where act.assay_id = ?""", (assay_id,)).fetchall()
out["assay_n_records"] = len(rows)
out["assay_n_molecules"] = len({r[1] for r in rows})
byid = {}
for cid, mr, sv, v, st in rows:
    byid.setdefault(cid, []).append(float(sv))
fallback = out["train_median_y"]
preds, cov = [], []
for did in te["Drug_ID"]:
    if did in byid:
        preds.append(statistics.median(byid[did])); cov.append(True)
    else:
        preds.append(fallback); cov.append(False)
te = te.copy()
te["pred"] = preds; te["covered"] = cov
te["abs_err"] = (te["Y"] - te["pred"]).abs()
c = te[te.covered]
out["q2"] = {
    "covered_n": int(te.covered.sum()), "test_n": len(te),
    "covered_pct": round(100 * te.covered.mean(), 2),
    "full_mae": float(te.abs_err.mean()),
    "full_rmse": float(np.sqrt((te.abs_err ** 2).mean())),
    "full_pearson_r": float(np.corrcoef(te.Y, te.pred)[0, 1]),
    "covered_mae": float(c.abs_err.mean()),
    "covered_exact_n_abs_lt_1e-9": int((c.abs_err < 1e-9).sum()),
    "covered_exact_n_abs_le_0.005": int((c.abs_err <= 0.005).sum()),
    "covered_max_abs_err": float(c.abs_err.max()),
    "fallback_train_median": fallback,
}
unc = te[~te.covered]
out["q2"]["uncovered"] = [
    {"Drug_ID": r.Drug_ID, "Y": float(r.Y), "smiles": r.Drug, "abs_err_vs_fallback": float(r.abs_err)}
    for r in unc.itertuples()]

# ---- 3. Relation analysis for the uncovered compounds
from rdkit import Chem
def ik(smi):
    m = Chem.MolFromSmiles(smi)
    return Chem.MolToInchiKey(m) if m else None

assay_rows = cur.execute("""
  select md.chembl_id, act.molregno, act.standard_value, cs.standard_inchi_key, cs.canonical_smiles
  from activities act join molecule_dictionary md on md.molregno = act.molregno
  left join compound_structures cs on cs.molregno = act.molregno
  where act.assay_id = ?""", (assay_id,)).fetchall()
q3 = []
for u in out["q2"]["uncovered"]:
    did, y = u["Drug_ID"], u["Y"]
    rec = {"Drug_ID": did, "Y": y}
    look = cur.execute("select chembl_id, entity_type, entity_id, status, last_active from chembl_id_lookup where chembl_id=?", (did,)).fetchone()
    rec["chembl_id_lookup"] = look
    md = cur.execute("select molregno, pref_name from molecule_dictionary where chembl_id=?", (did,)).fetchone()
    rec["in_molecule_dictionary"] = md
    tdc_ik = ik(u["smiles"]); rec["tdc_smiles_inchikey"] = tdc_ik
    db_ik = None
    if md:
        r = cur.execute("select standard_inchi_key, canonical_smiles from compound_structures where molregno=?", (md[0],)).fetchone()
        db_ik = r[0] if r else None
        rec["db_inchikey"] = db_ik
    same_val = [r for r in assay_rows if abs(r[2] - y) < 1e-9]
    rec["assay_records_with_identical_value"] = []
    for cid, mr, sv, aik, asmi in same_val:
        e = {"chembl_id": cid, "molregno": mr, "value": sv, "inchikey": aik, "smiles": asmi}
        h = cur.execute("""select p.chembl_id, c.chembl_id from molecule_hierarchy mh
              join molecule_dictionary p on p.molregno = mh.parent_molregno
              join molecule_dictionary c on c.molregno = mh.molregno
              where mh.molregno = ? or mh.parent_molregno = ?""", (mr, mr)).fetchall()
        # parent of the assay molecule, and children of the assay molecule
        par = cur.execute("""select p.chembl_id from molecule_hierarchy mh join molecule_dictionary p
              on p.molregno = mh.parent_molregno where mh.molregno = ?""", (mr,)).fetchone()
        e["assay_mol_parent"] = par[0] if par else None
        e["test_id_is_parent_of_assay_mol"] = bool(par and par[0] == did and cid != did)
        if md:
            par_t = cur.execute("""select p.chembl_id from molecule_hierarchy mh join molecule_dictionary p
                  on p.molregno = mh.parent_molregno where mh.molregno = ?""", (md[0],)).fetchone()
            e["test_mol_parent"] = par_t[0] if par_t else None
        ref = db_ik or tdc_ik
        e["inchikey_full_match_test"] = bool(aik and ref and aik == ref)
        e["inchikey_block1_match_test"] = bool(aik and ref and aik[:14] == ref[:14])
        e["inchikey_block1_match_tdc_smiles"] = bool(aik and tdc_ik and aik[:14] == tdc_ik[:14])
        rec["assay_records_with_identical_value"].append(e)
    q3.append(rec)
out["q3"] = q3

# ---- 4. Decomposition
s6 = float(unc.abs_err.sum())
out["q4"] = {"sum_abs_err_uncovered": s6, "sum_over_840": s6 / len(te),
             "covered_sum_abs_err": float(c.abs_err.sum()),
             "full_mae": out["q2"]["full_mae"],
             "uncovered_share_of_total_error": s6 / float(te.abs_err.sum()) if te.abs_err.sum() else None}

te[["Drug_ID", "Y", "pred", "covered", "abs_err"]].to_csv(W / "arm_a_source_test_predictions.csv", index=False)
json.dump(out, open(W / "testfold_identity_and_misses.json", "w"), indent=2, default=str)
print(json.dumps(out, indent=1, default=str))
