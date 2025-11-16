
# stats_pipeline.py
#
# v2 — Auto-pairs
# - Outputs to ./stat_results
# - Default input: ./results  (so you can just run `python stats_pipeline.py`)
# - For every (scenario, cpu), test:
#     * B0 → EVERY other present variant
#     * PLUS special extras (e.g., I2+I3+I4 → I2+I3+I4+I5), if both present
#
# Usage:
#   1) py -m venv .venv
#   2) .venv\Scripts\activate
#   3) pip install -r requirements.txt
#   4) python stats_pipeline.py            # uses ./results → ./stat_results
#      (optional) python stats_pipeline.py --root "D:\\path\\to\\results" --out "D:\\path\\to\\stat_results"
#
# Outputs in ./stat_results:
#   - raw_interactions.csv
#   - summary_by_group.csv
#   - charts/boxplot_*.png, charts/ecdf_*.png
#   - mw_results.csv

import argparse
import json
import math
import re
from pathlib import Path
from typing import List, Tuple, Dict, Iterable

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu

def find_files(root: Path, name: str) -> List[Path]:
    return [p for p in root.rglob(name) if p.is_file()]

def try_load_json(p: Path):
    try:
        return json.loads(p.read_text(encoding="utf-8"))
    except Exception:
        return None

_CPU_RX = re.compile(r"CPUx(\d+)", re.IGNORECASE)

def infer_cpu_from_path(p: Path):
    m = _CPU_RX.search(str(p))
    if m:
        try:
            return float(m.group(1))
        except Exception:
            return None
    return None

def extract_cpu_from_meta(dirpath: Path):
    for name in ("lab-results.by-variant.json", "lab-results.aggregates.json"):
        cand = dirpath / name
        data = try_load_json(cand)
        if isinstance(data, dict):
            meta = data.get("meta")
            if isinstance(meta, dict):
                cpu = meta.get("cpuThrottlingRate")
                if isinstance(cpu, (int, float)):
                    return float(cpu)
    return infer_cpu_from_path(dirpath)

def ecdf(y: np.ndarray):
    y = np.asarray(y)
    y = y[np.isfinite(y)]
    y.sort()
    n = len(y)
    if n == 0:
        return np.array([]), np.array([])
    x = y
    f = np.arange(1, n+1) / n
    return x, f

def cliffs_delta(x: np.ndarray, y: np.ndarray) -> float:
    x = np.asarray(x); y = np.asarray(y)
    x = x[np.isfinite(x)]; y = y[np.isfinite(y)]
    if len(x) == 0 or len(y) == 0: return math.nan
    xs = np.sort(x); ys = np.sort(y); ny = len(ys)
    lt = 0; j = 0
    for xv in xs:
        while j < ny and ys[j] < xv: j += 1
        lt += j
    gt = 0; j = ny - 1
    for xv in xs:
        while j >= 0 and ys[j] > xv: j -= 1
        gt += (ny - 1 - j)
    n = len(xs) * ny
    return (gt - lt) / n if n else math.nan

def holm_bonferroni(pvals: List[float]) -> List[float]:
    m = len(pvals)
    order = sorted(range(m), key=lambda i: (math.inf if pvals[i] is None else pvals[i]))
    adj = [1.0] * m; prev = 0.0
    for k, i in enumerate(order, start=1):
        p = 1.0 if pvals[i] is None else pvals[i]
        adj_p = (m - k + 1) * p
        adj_p = max(adj_p, prev)
        adj[i] = min(adj_p, 1.0)
        prev = adj[i]
    return adj

SPECIAL_EXTRAS: List[Tuple[str, str]] = [
    ("I2+I3+I4", "I2+I3+I4+I5"),
]

def build_pairs_for_group(variants):
    variants = [v for v in variants if isinstance(v, str) and len(v.strip())]
    uniq = sorted(set(variants))
    pairs = []
    if "B0" in uniq:
        for v in uniq:
            if v != "B0":
                pairs.append(("B0", v))
    for a, b in SPECIAL_EXTRAS:
        if a in uniq and b in uniq and (a, b) not in pairs:
            pairs.append((a, b))
    return pairs

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=None, help="Input root. Default: ./results")
    ap.add_argument("--out", default=None, help="Output dir. Default: ./stat_results")
    args = ap.parse_args()

    cwd = Path.cwd()
    root = (Path(args.root).expanduser().resolve() if args.root else (cwd / "results").resolve())
    out_dir = (Path(args.out).expanduser().resolve() if args.out else (cwd / "stat_results").resolve())
    charts_dir = out_dir / "charts"
    out_dir.mkdir(parents=True, exist_ok=True)
    charts_dir.mkdir(parents=True, exist_ok=True)

    # 1) Collect raw interactions
    rows = []
    for f in find_files(root, "lab-results.json"):
        data = try_load_json(f)
        if not isinstance(data, list):
            continue
        cpu = extract_cpu_from_meta(f.parent) or infer_cpu_from_path(f)
        for rec in data:
            rows.append({
                "scenario": rec.get("scenario"),
                "variant": rec.get("variant"),
                "rep": rec.get("rep"),
                "cpu": cpu,
                "INP": rec.get("INP"),
                "INP_W": rec.get("INP_W"),
                "INP_H": rec.get("INP_H"),
                "INP_R": rec.get("INP_R"),
                "LoAFsum": rec.get("LoAFsum"),
                "LoAFcount": rec.get("LoAFcount"),
                "LoAFsum_500ms": rec.get("LoAFsum_500ms"),
                "LoAFcount_500ms": rec.get("LoAFcount_500ms"),
                "LoAFany_500ms": rec.get("LoAFany_500ms"),
                "source_dir": str(f.parent),
            })
    if not rows:
        print("No lab-results.json files found under:", root); return

    df = pd.DataFrame(rows)
    df.to_csv(out_dir / "raw_interactions.csv", index=False)
    print(f"[OK] Saved {out_dir / 'raw_interactions.csv'} with {len(df)} rows")

    # 2) Summary by (scenario, variant, cpu)
    def q75(x):
        x = np.asarray(x, dtype=float); x = x[np.isfinite(x)]
        return np.quantile(x, 0.75) if len(x) else math.nan
    def q50(x):
        x = np.asarray(x, dtype=float); x = x[np.isfinite(x)]
        return np.quantile(x, 0.5) if len(x) else math.nan

    grp = df.groupby(["scenario", "variant", "cpu"], dropna=False).agg(
        n=("INP", "count"),
        p50_INP_ms=("INP", q50),
        p75_INP_ms=("INP", q75),
        loaf_sum_median=("LoAFsum", q50),
        loaf_any_pct=("LoAFany_500ms", lambda x: 100 * float(np.mean(np.array(x, dtype=float) > 0)) if len(x) else math.nan),
    ).reset_index()
    grp.to_csv(out_dir / "summary_by_group.csv", index=False)
    print(f"[OK] Saved {out_dir / 'summary_by_group.csv'}")

    # 3) Charts per (scenario, cpu)
    for (scenario, cpu), dsub in df.groupby(["scenario", "cpu"]):
        order = sorted([v for v in dsub["variant"].dropna().unique()])
        if not order: continue

        # Boxplot
        fig = plt.figure(figsize=(8, 5))
        data = [dsub.loc[dsub["variant"] == v, "INP"].dropna().values for v in order]
        plt.boxplot(data, labels=order, showfliers=False)
        plt.ylabel("INP (ms)")
        plt.title(f"INP by variant — scenario={scenario}, CPUx{int(cpu) if pd.notna(cpu) else 'NA'}")
        fig.tight_layout()
        fig.savefig(charts_dir / f"boxplot_{scenario}_CPUx{int(cpu) if pd.notna(cpu) else 'NA'}.png", dpi=150)
        plt.close(fig)

        # ECDF
        fig = plt.figure(figsize=(8, 5))
        for v in order:
            y = dsub.loc[dsub["variant"] == v, "INP"].dropna().values
            x, f = ecdf(y)
            if len(x): plt.step(x, f, where="post", label=str(v))
        plt.xlabel("INP (ms)"); plt.ylabel("F(x)")
        plt.title(f"ECDF — scenario={scenario}, CPUx{int(cpu) if pd.notna(cpu) else 'NA'}")
        plt.legend(); fig.tight_layout()
        fig.savefig(charts_dir / f"ecdf_{scenario}_CPUx{int(cpu) if pd.notna(cpu) else 'NA'}.png", dpi=150)
        plt.close(fig)

    print(f"[OK] Saved charts to {charts_dir}")

    # 4) Mann–Whitney U (auto pairs)
    results = []
    for (scenario, cpu), dsub in df.groupby(["scenario", "cpu"]):
        variants = dsub["variant"].dropna().unique().tolist()
        pairs = build_pairs_for_group(variants)
        if not pairs: continue

        pvals = []; tmp = []
        for (a, b) in pairs:
            opt = dsub.loc[dsub["variant"] == b, "INP"].dropna().values
            base = dsub.loc[dsub["variant"] == a, "INP"].dropna().values
            if len(opt) == 0 or len(base) == 0:
                pvals.append(None)
                tmp.append((a, b, math.nan, math.nan, math.nan, math.nan, len(base), len(opt)))
                continue
            res = mannwhitneyu(opt, base, alternative="less")  # optimized < baseline
            U = float(res.statistic); p = float(res.pvalue)
            n = len(opt); m = len(base)
            auc = U / (n * m) if n and m else math.nan
            cd = cliffs_delta(opt, base)
            tmp.append((a, b, U, p, auc, cd, m, n))
            pvals.append(p)

        p_corr = holm_bonferroni(pvals) if pvals else []
        for (row, pc) in zip(tmp, p_corr):
            a, b, U, p, auc, cd, m, n = row
            results.append({
                "scenario": scenario,
                "cpu": cpu,
                "baseline": a,
                "optimized": b,
                "n_baseline": m,
                "n_optimized": n,
                "U": U,
                "p_value": p,
                "p_value_holm": pc,
                "AUC_prob_opt_lt_base": auc,
                "cliffs_delta": cd,
            })

    if results:
        pd.DataFrame(results).to_csv(out_dir / "mw_results.csv", index=False)
        print(f"[OK] Saved {out_dir / 'mw_results.csv'}")
    else:
        print("No comparable pairs found for Mann-Whitney U.")

if __name__ == "__main__":
    main()
