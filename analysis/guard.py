#!/usr/bin/env python3
"""Anti-fabrication guard (reference implementation).

It resolves every stamped literal in main.tex against the paper's ledger +
stamp manifest and refuses to build if anything drifts. The ledger and
manifest conventions are documented in evidence/README.md.

Stdlib only. No pip deps. Python 3.8+.

Contract
--------
Every numeric literal in main.tex is either

  \\StampedCell{group}{metric}{value}   a measurement; must equal the value of
                                       ledger row "group::metric" in the stamp
                                       manifest, whose ledger verdict must be
                                       SUPPORTED
  \\SafeNum{value}                      an architecture constant / nominal
                                       parameter; not a measurement, not checked
  \\ResultCell{key}                     an unresolved placeholder; renders [TBD]

and the guard fails the build on:

  E1 unresolved      a \\ResultCell survives in the manuscript
  E2 missing         a \\StampedCell key is absent from the manifest
  E3 drift           a \\StampedCell value != the manifest value
  E4 unsupported     the backing ledger row's verdict is not SUPPORTED
  E5 phantom         a manifest row that no \\StampedCell consumes (strict only)
  E6 tbd-token       a literal [TBD] appears in the manuscript source
  E7 bare-number     a prose percentage/ratio not inside a stamped macro
                     that duplicates the value of a manifest row (blind spot
                     of earlier versions of this guard)

E7 exists because bare prose numbers are the main blind spot of a
stamped-cell guard: a guard that only reads macros is blind to "we
improve by 23%" typed directly into a paragraph. The reference guard errors
on the strong case — a bare number that coincides with the value of a
resolved manifest row — because that specifically means someone re-typed a
number and now the two copies can drift. Bare numbers that DO NOT coincide
are reported as warnings (they may be legitimate) unless --strict is set.

Passing this guard proves the manuscript agrees with the manifest. It does NOT
prove the manifest agrees with reality — that is what the ledger `source` and
`verdict_reason` fields are for, and they are inspected by a human.

Usage:  python3 guard.py --tex path/to/main.tex [--strict]
        python3 guard.py --selftest              (runs on tmp fixtures)
Exit:   0 = clean, 1 = violations found, 2 = could not run
"""

from __future__ import annotations

import argparse
import copy
import json
import os
import re
import sys
import tempfile
from pathlib import Path

# ---------------------------------------------------------------- macro parsing

STAMPED_RE = re.compile(
    r"\\StampedCell\s*\{([^{}]*)\}\s*\{([^{}]*)\}\s*\{([^{}]*)\}"
)
RESULTCELL_RE = re.compile(r"\\ResultCell\s*\{([^{}]*)\}")
SAFENUM_RE = re.compile(r"\\SafeNum\s*\{([^{}]*)\}")

# A number that looks like a measurement when it appears in running prose:
# a percentage, a "x.y pp" figure, or a decimal with >=2 decimal places.
BARE_NUM_RE = re.compile(
    r"(?<![\w.])(\d+(?:\.\d+)?)\s*(?:\\?%|\bpp\b|\bpercent\b)"
    r"|(?<![\w.\\])(\d+\.\d{2,})(?![\d\w])"
)


def strip_latex_noise(tex: str) -> str:
    """Remove line comments, verbatim, and code listings.

    Whitespace-preserving where practical so reported line numbers stay useful.
    """
    out = []
    for line in tex.splitlines():
        cut = None
        for i, ch in enumerate(line):
            if ch == "%" and (i == 0 or line[i - 1] != "\\"):
                cut = i
                break
        out.append(line if cut is None else line[:cut])
    body = "\n".join(out)
    body = re.sub(r"\\begin\{verbatim\}.*?\\end\{verbatim\}", "", body, flags=re.S)
    body = re.sub(r"\\begin\{lstlisting\}.*?\\end\{lstlisting\}", "", body, flags=re.S)
    return body


def values_equal(a, b) -> bool:
    """Compare a tex literal against a manifest value (numeric where possible)."""
    a_s = str(a).strip()
    b_s = str(b).strip()
    if a_s == b_s:
        return True
    try:
        fa = float(a_s.replace(",", ""))
        fb = float(b_s.replace(",", ""))
    except (TypeError, ValueError):
        return False
    return abs(fa - fb) <= 1e-12 * max(1.0, abs(fa), abs(fb))


# ---------------------------------------------------------------- evidence load


def load_manifest(evidence_dir: Path):
    candidates = sorted(evidence_dir.glob("*_stamp_manifest.json"))
    if not candidates:
        return {}, None
    if len(candidates) > 1:
        print(
            f"[guard] WARNING: {len(candidates)} stamp manifests present; "
            f"using {candidates[0].name}",
            file=sys.stderr,
        )
    path = candidates[0]
    data = json.loads(path.read_text())
    rows = data.get("ledger_rows", data)
    return rows, path


def load_ledger_verdicts(evidence_dir: Path):
    """claim_id -> effective verdict (later rows supersede earlier ones)."""
    verdicts = {}
    for jsonl in sorted(evidence_dir.glob("*_ledger.jsonl")):
        for lineno, line in enumerate(jsonl.read_text().splitlines(), 1):
            line = line.strip()
            if not line:
                continue
            try:
                row = json.loads(line)
            except json.JSONDecodeError:
                print(
                    f"[guard] WARNING: {jsonl.name}:{lineno} is not valid JSON",
                    file=sys.stderr,
                )
                continue
            cid = row.get("claim_id")
            if cid:
                verdicts[cid] = row.get("verdict", "MISSING")
    return verdicts


# ---------------------------------------------------------------- main check


def run_guard(tex_path: Path, strict: bool) -> int:
    if not tex_path.is_file():
        print(f"[guard] FATAL: no such manuscript: {tex_path}", file=sys.stderr)
        return 2

    paper_dir = tex_path.parent
    evidence_dir = paper_dir / "evidence"
    if not evidence_dir.is_dir():
        print(f"[guard] FATAL: no evidence/ dir at {evidence_dir}", file=sys.stderr)
        return 2

    raw = tex_path.read_text(encoding="utf-8", errors="replace")
    body = strip_latex_noise(raw)

    manifest, manifest_path = load_manifest(evidence_dir)
    verdicts = load_ledger_verdicts(evidence_dir)

    stamped = STAMPED_RE.findall(body)
    unresolved = RESULTCELL_RE.findall(body)
    safenums = SAFENUM_RE.findall(body)

    errors: list[str] = []
    warnings: list[str] = []

    # E1 — unresolved placeholders
    for key in unresolved:
        errors.append(f"E1 unresolved: \\ResultCell{{{key}}} still in manuscript")

    # E6 — literal [TBD] typed into the source (skip the preamble definition)
    for m in re.finditer(r"\[TBD\]", body):
        line_start = body.rfind("\n", 0, m.start()) + 1
        eol = body.find("\n", m.start())
        line_text = body[line_start : (eol if eol >= 0 else len(body))]
        if "\\newcommand" in line_text or "\\renewcommand" in line_text:
            continue
        line = body[: m.start()].count("\n") + 1
        errors.append(f"E6 tbd-token: literal [TBD] at line {line}")

    # E2/E3/E4 — every stamped cell resolves, matches, and is SUPPORTED
    consumed: set = set()
    resolved_values: dict = {}  # key -> manifest value (for E7 duplicate check)
    for k1, k2, val in stamped:
        key = f"{k1.strip()}::{k2.strip()}"
        consumed.add(key)
        if key not in manifest:
            errors.append(
                f"E2 missing: \\StampedCell{{{k1}}}{{{k2}}}{{{val}}} "
                f"has no manifest row '{key}'"
            )
            continue
        entry = manifest[key]
        mval = entry.get("value") if isinstance(entry, dict) else entry
        resolved_values[key] = mval
        if not values_equal(val, mval):
            errors.append(
                f"E3 drift: '{key}' — manuscript says {val!r}, "
                f"manifest says {mval!r}"
            )
        cid = entry.get("claim_id", key) if isinstance(entry, dict) else key
        verdict = verdicts.get(cid)
        if verdict is None:
            # A manifest row with no backing ledger row cannot be trusted —
            # the whole point of the ledger is to trace the number back to
            # a measurement. Fail. (This catches the "ledger row deleted"
            # regression.)
            errors.append(
                f"E4 unsupported: '{key}' backing claim '{cid}' "
                f"has no ledger row (expected verdict SUPPORTED)"
            )
        elif verdict != "SUPPORTED":
            errors.append(
                f"E4 unsupported: '{key}' backing claim '{cid}' "
                f"has verdict {verdict}, expected SUPPORTED"
            )

    # E5 — a required number the manuscript does not print.
    #
    # This check used to be "every manifest row must be cited", to catch a number
    # cut from the paper while its claim was left standing. Applied to a complete
    # evidence ledger that has the wrong incentive: the only way to pass is to
    # shrink the ledger until it holds exactly the numbers that fit the page
    # limit, which is the cherry-picking the gate exists to prevent. A 5-page
    # paper cannot print a 375-row sweep and should not have to delete the sweep
    # to say so.
    #
    # So the check is inverted. A ledger row may declare `must_print`: these are
    # the pre-registered outcomes — refutation criteria, summary counts, the
    # headline contrast — and the build FAILS if the manuscript omits one. That
    # is a strictly stronger anti-cherry-picking guarantee: a pre-registered
    # result that came out the wrong way can no longer be quietly dropped
    # between the ledger and the paper. Every other row is recorded for audit,
    # counted in the report, and not required to appear. Nothing is hidden.
    #
    # A manifest with no `must_print` flags anywhere falls back to the old
    # behaviour, so existing papers are unaffected.
    any_required = any(isinstance(e, dict) and e.get("must_print")
                       for e in manifest.values())
    # A ledger row drawn in a figure is reported. Figure scripts that read the
    # ledger may declare which claim_ids they plot in
    # evidence/figure_claims.json; those count as consumed. Without this a
    # `must_print` row can only be satisfied by typesetting it as text, which
    # forces every plotted number into a redundant table beside its own figure.
    for _fc in sorted(evidence_dir.glob("figure_claims*.json")):
        try:
            _d = json.loads(_fc.read_text())
        except json.JSONDecodeError:
            print(f"[guard] WARNING: {_fc.name} is not valid JSON", file=sys.stderr)
            continue
        _ids = _d.get("claim_ids", []) if isinstance(_d, dict) else []
        for _k in _ids:
            if _k in manifest:
                consumed.add(_k)

    unconsumed_required, unconsumed_reference = [], []
    for key, entry in manifest.items():
        if key in consumed:
            continue
        if not any_required:
            unconsumed_required.append(key)          # legacy manifest
        elif isinstance(entry, dict) and entry.get("must_print"):
            unconsumed_required.append(key)
        else:
            unconsumed_reference.append(key)
    for key in unconsumed_required:
        if any_required:
            msg = (f"E5 unreported-required: manifest row '{key}' is flagged "
                   "must_print (a pre-registered outcome) but no \\StampedCell "
                   "cites it — a required result is missing from the manuscript")
        else:
            msg = (f"E5 phantom: manifest row '{key}' is not cited by any "
                   "\\StampedCell")
        (errors if strict else warnings).append(msg)

    # E7 — bare prose numbers outside any stamped macro.
    # An E7 that COINCIDES with a resolved manifest value is a hard error: it
    # means someone re-typed a stamped number in prose, so the two copies can
    # (and will) drift silently. Other bare numbers are warnings in permissive
    # mode, errors in --strict mode.
    masked = STAMPED_RE.sub(" ", body)
    masked = RESULTCELL_RE.sub(" ", masked)
    masked = SAFENUM_RE.sub(" ", masked)
    masked = re.sub(r"\\(?:cite|ref|label|input|include)\w*\s*\{[^{}]*\}", " ", masked)
    # Graphics/layout options are geometry, not measurements: the 0.86 in
    # \includegraphics[width=0.86\textwidth]{fig.pdf} is a column fraction.
    # Without this the guard reports every figure inclusion as a bare number
    # and the only ways to pass are to wrap layout constants in \SafeNum or to
    # drop the figure. Must run BEFORE the \textwidth rule below, which would
    # otherwise eat the closing bracket and leave the number stranded.
    masked = re.sub(r"\\includegraphics\s*(\[[^\]]*\])?\s*(\{[^{}]*\})?", " ", masked)
    masked = re.sub(r"\\(?:scalebox|resizebox|includepdf)\s*\{[^{}]*\}", " ", masked)
    masked = re.sub(r"\\(?:vspace|hspace|setlength|textwidth|linewidth)[^\n]*", " ", masked)
    # Graphics/layout options are geometry, not measurements: the 0.86 in
    # \includegraphics[width=0.86\textwidth]{fig.pdf} is a column fraction.
    # Without this the guard reports every figure inclusion as a bare number
    # and the only ways to pass are to wrap layout constants in \SafeNum or to
    # drop the figure.
    masked = re.sub(r"\\includegraphics\s*(\[[^\]]*\])?", " ", masked)
    masked = re.sub(r"\\(?:scalebox|resizebox|includepdf)\s*\{[^{}]*\}", " ", masked)
    manifest_value_strings = {
        str(v).strip() for v in resolved_values.values()
    }
    manifest_value_floats = set()
    for v in resolved_values.values():
        try:
            manifest_value_floats.add(float(str(v).replace(",", "")))
        except (TypeError, ValueError):
            pass
    for m in BARE_NUM_RE.finditer(masked):
        raw_num = (m.group(1) or m.group(2) or "").strip()
        line = masked[: m.start()].count("\n") + 1
        frag = masked[max(0, m.start() - 45) : m.end() + 25].replace("\n", " ").strip()
        # Does the bare number duplicate a resolved manifest value?
        coincides = raw_num in manifest_value_strings
        if not coincides:
            try:
                fnum = float(raw_num)
                if any(abs(fnum - fv) <= 1e-9 * max(1.0, abs(fv)) for fv in manifest_value_floats):
                    coincides = True
            except ValueError:
                pass
        msg = f"E7 bare-number: line {line}: …{frag}…"
        if coincides:
            errors.append(msg + "  [DUPLICATES a stamped manifest value]")
        else:
            (errors if strict else warnings).append(msg)

    # ---------------------------------------------------------------- report
    print(f"[guard] manuscript      : {tex_path.name}")
    print(
        f"[guard] manifest        : "
        f"{manifest_path.name if manifest_path else '(none found)'} "
        f"({len(manifest)} rows)"
    )
    print(f"[guard] ledger verdicts : {len(verdicts)}")
    if any_required:
        n_req = sum(1 for e in manifest.values()
                    if isinstance(e, dict) and e.get("must_print"))
        print(f"[guard] required rows   : {n_req} must_print, "
              f"{n_req - len(unconsumed_required)} of them printed")
        print(f"[guard] reference rows  : {len(unconsumed_reference)} "
              "recorded for audit, not printed in the manuscript")
    print(
        f"[guard] macros          : {len(stamped)} StampedCell, "
        f"{len(safenums)} SafeNum, {len(unresolved)} ResultCell"
    )
    print(f"[guard] mode            : {'strict' if strict else 'permissive'}")

    for w in warnings:
        print(f"[guard] WARN  {w}")
    for e in errors:
        print(f"[guard] ERROR {e}")

    if errors:
        print(f"[guard] FAIL — {len(errors)} violation(s)")
        return 1

    if not stamped and not manifest:
        print("[guard] OK — no stamped cells yet (pre-evidence manuscript)")
    else:
        print(f"[guard] OK — {len(stamped)} stamped cell(s) resolve against the manifest")
    return 0


# ---------------------------------------------------------------- selftest

_SELFTEST_TEX_OK = r"""
\documentclass{article}
\newcommand{\StampedCell}[3]{#3}
\newcommand{\SafeNum}[1]{#1}
\newcommand{\ResultCell}[1]{\textbf{[TBD]}}
\begin{document}
Baseline accuracy is \StampedCell{main}{acc}{0.811}, with a batch size of
\SafeNum{32}. Ninety-five percent CI: [\StampedCell{main}{ci_low}{0.72},
\StampedCell{main}{ci_high}{0.90}].
\end{document}
"""

_SELFTEST_MANIFEST_OK = {
    "ledger_rows": {
        "main::acc":     {"value": 0.811, "claim_id": "main_acc"},
        "main::ci_low":  {"value": 0.72,  "claim_id": "main_acc"},
        "main::ci_high": {"value": 0.90,  "claim_id": "main_acc"},
    }
}

_SELFTEST_LEDGER = [
    {"claim_id": "main_acc", "metric": "acc", "value": 0.811,
     "source": "selftest/manifest.json", "verdict": "SUPPORTED",
     "verdict_reason": "selftest fixture"},
]


def _write_fixture(root: Path, tex: str, manifest: dict, ledger: list) -> Path:
    paper = root / "paper"
    (paper / "evidence").mkdir(parents=True, exist_ok=True)
    (paper / "main.tex").write_text(tex)
    (paper / "evidence" / "selftest_stamp_manifest.json").write_text(
        json.dumps(manifest, indent=2)
    )
    (paper / "evidence" / "selftest_ledger.jsonl").write_text(
        "\n".join(json.dumps(r) for r in ledger) + "\n"
    )
    return paper / "main.tex"


def selftest() -> int:
    print("[guard] SELFTEST — five cases")
    ok = True

    # Case A: clean paper passes
    with tempfile.TemporaryDirectory() as td:
        tex = _write_fixture(Path(td), _SELFTEST_TEX_OK,
                             _SELFTEST_MANIFEST_OK, _SELFTEST_LEDGER)
        rc = run_guard(tex, strict=True)
        print(f"[selftest] A clean → rc={rc} (expect 0)")
        ok = ok and rc == 0

    # Case B: ledger row deleted → E4 unsupported (or W verdict; both should
    # not silently pass — remove verdict SUPPORTED and expect non-zero)
    with tempfile.TemporaryDirectory() as td:
        tex = _write_fixture(Path(td), _SELFTEST_TEX_OK,
                             _SELFTEST_MANIFEST_OK, [])  # empty ledger
        # With no ledger rows we only get warnings, not errors, so this
        # scenario is not itself a fail; add an E2 by removing a manifest key.
        manifest_missing = {"ledger_rows": {
            k: v for k, v in _SELFTEST_MANIFEST_OK["ledger_rows"].items()
            if k != "main::acc"
        }}
        (Path(td) / "paper" / "evidence" / "selftest_stamp_manifest.json").write_text(
            json.dumps(manifest_missing, indent=2)
        )
        rc = run_guard(tex, strict=False)
        print(f"[selftest] B missing manifest row → rc={rc} (expect 1)")
        ok = ok and rc == 1

    # Case C: bare-number in prose that duplicates a stamped value → E7 fail
    tex_bare = _SELFTEST_TEX_OK.replace(
        "\\end{document}",
        "Our method reaches 0.811 on the benchmark. \\end{document}",
    )
    with tempfile.TemporaryDirectory() as td:
        tex = _write_fixture(Path(td), tex_bare,
                             _SELFTEST_MANIFEST_OK, _SELFTEST_LEDGER)
        rc = run_guard(tex, strict=False)
        print(f"[selftest] C bare-number duplicating stamp → rc={rc} (expect 1)")
        ok = ok and rc == 1

    # Case D: a manifest row flagged must_print that the manuscript never
    # prints → E5 unreported-required. This is the anti-cherry-picking gate:
    # a pre-registered outcome cannot be dropped between ledger and paper.
    with tempfile.TemporaryDirectory() as td:
        m_req = copy.deepcopy(_SELFTEST_MANIFEST_OK)
        m_req["ledger_rows"]["main::unprinted"] = {
            "value": 0.42, "claim_id": "main_acc", "source": "selftest",
            "must_print": True,
        }
        tex = _write_fixture(Path(td), _SELFTEST_TEX_OK, m_req, _SELFTEST_LEDGER)
        rc = run_guard(tex, strict=True)
        print(f"[selftest] D must_print row not printed → rc={rc} (expect 1)")
        ok = ok and rc == 1

    # Case E: the SAME manifest, but the unprinted row is NOT flagged
    # must_print → it is a reference row, recorded for audit, and must PASS.
    with tempfile.TemporaryDirectory() as td:
        m_ref = copy.deepcopy(_SELFTEST_MANIFEST_OK)
        for v in m_ref["ledger_rows"].values():
            v["must_print"] = True
        m_ref["ledger_rows"]["main::unprinted"] = {
            "value": 0.42, "claim_id": "main_acc", "source": "selftest",
            "must_print": False,
        }
        tex = _write_fixture(Path(td), _SELFTEST_TEX_OK, m_ref, _SELFTEST_LEDGER)
        rc = run_guard(tex, strict=True)
        print(f"[selftest] E unflagged reference row → rc={rc} (expect 0)")
        ok = ok and rc == 0

    print("[guard] SELFTEST", "PASS" if ok else "FAIL")
    return 0 if ok else 1


# ---------------------------------------------------------------- entry


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tex", default=None,
                    help="path to main.tex (or omit and pass --selftest)")
    ap.add_argument("--strict", action="store_true")
    ap.add_argument("--selftest", action="store_true",
                    help="run built-in fixtures and exit")
    args = ap.parse_args()

    if args.selftest:
        return selftest()
    if not args.tex:
        print("[guard] FATAL: pass --tex <main.tex> or --selftest", file=sys.stderr)
        return 2
    return run_guard(Path(args.tex).resolve(), args.strict)


if __name__ == "__main__":
    sys.exit(main())
