# stats_pipeline.py
#
# v5 — Auto-pairs (two-sided MWU + direction via effect size) + publication-ready tables
# - Outputs to ./stat_results
# - Default input: ./results  (so you can just run `python stats_pipeline.py`)
# - For every (scenario, cpu), test:
#     * B0 → EVERY other present variant
#     * PLUS special extras (e.g., I2+I3+I4 → I2+I3+I4+I5), if both present
#
# Statistical notes:
# - Mann–Whitney U: TWO-SIDED to detect both improvements and degradations.
# - Ties are common in INP (quantization); we force method="asymptotic" for stable behavior.
# - Direction via effect size:
#     PS(A) = P(INP_opt < INP_base) + 0.5*P(ties)    (aka "probability of superiority")
#     delta = 2*PS(A) - 1   (delta > 0 => optimized tends to have smaller INP)
#
# Usage:
#   1) py -m venv .venv
#   2) .venv\Scripts\activate
#   3) pip install -r requirements.txt
#   4) python stats_pipeline.py
#      (optional) python stats_pipeline.py --root "D:\\path\\to\\results" --out "D:\\path\\to\\stat_results"
#
# Outputs in ./stat_results:
#   - raw_interactions.csv
#   - summary_by_group.csv
#   - charts/boxplot_*.png, charts/ecdf_*.png
#   - mw_results.csv
#   - mw_results_enriched.csv            (adds p75(B0), p75(var), Δp75%)
#   - table_1_4_counts.csv               (counts of CPU regimes with p_Holm<0.05 & delta>0 for B0 comparisons)
#   - table_regressions_counts.csv       (counts of CPU regimes with p_Holm<0.05 & delta<0 for B0 comparisons)
#   - regressions_list.csv               (all significant regressions vs B0)
#   - table_1_5_candidates.csv           (top improvements + top regressions per scenario)
#   - run_info.json                      (versions + run parameters)

import argparse
import json
import math
import re
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import List, Tuple

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
    f = np.arange(1, n + 1) / n
    return x, f


def holm_bonferroni(pvals: List[float]) -> List[float]:
    """
    Holm step-down adjustment (controls FWER).
    """
    m = len(pvals)
    order = sorted(range(m), key=lambda i: (math.inf if pvals[i] is None else pvals[i]))
    adj = [1.0] * m
    prev = 0.0
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


def safe_boxplot(data, labels):
    """
    Matplotlib 3.9+: `labels` renamed to `tick_labels`.
    Keep backward compatibility.
    """
    try:
        plt.boxplot(data, tick_labels=labels, showfliers=False)
    except TypeError:
        plt.boxplot(data, labels=labels, showfliers=False)


def _pkg_version(obj) -> str:
    v = getattr(obj, "__version__", None)
    return str(v) if v is not None else "unknown"


def write_run_info(out_dir: Path, root: Path):
    info = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "python": sys.version.replace("\n", " "),
        "numpy": _pkg_version(np),
        "pandas": _pkg_version(pd),
        "matplotlib": _pkg_version(plt.matplotlib),
        "scipy": "unknown",
        "params": {"root": str(root), "out": str(out_dir)},
    }
    try:
        import scipy as _scipy
        info["scipy"] = _pkg_version(_scipy)
    except Exception:
        pass

    (out_dir / "run_info.json").write_text(json.dumps(info, ensure_ascii=False, indent=2), encoding="utf-8")


def safe_pct(opt, base):
    if opt is None or base is None:
        return math.nan
    try:
        opt = float(opt)
        base = float(base)
        if not (math.isfinite(opt) and math.isfinite(base)) or base == 0:
            return math.nan
        return 100.0 * (opt - base) / base
    except Exception:
        return math.nan


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

    write_run_info(out_dir, root)

    # 1) Collect raw interactions
    rows = []
    for f in find_files(root, "lab-results.json"):
        data = try_load_json(f)
        if not isinstance(data, list):
            continue
        cpu = extract_cpu_from_meta(f.parent) or infer_cpu_from_path(f)
        for rec in data:
            rows.append(
                {
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
                }
            )

    if not rows:
        print("No lab-results.json files found under:", root)
        return

    df = pd.DataFrame(rows)
    df.to_csv(out_dir / "raw_interactions.csv", index=False)
    print(f"[OK] Saved {out_dir / 'raw_interactions.csv'} with {len(df)} rows")

    # 2) Summary by (scenario, variant, cpu)
    def q75(x):
        x = np.asarray(x, dtype=float)
        x = x[np.isfinite(x)]
        return np.quantile(x, 0.75) if len(x) else math.nan

    def q50(x):
        x = np.asarray(x, dtype=float)
        x = x[np.isfinite(x)]
        return np.quantile(x, 0.5) if len(x) else math.nan

    grp = (
        df.groupby(["scenario", "variant", "cpu"], dropna=False)
        .agg(
            n=("INP", "count"),
            p50_INP_ms=("INP", q50),
            p75_INP_ms=("INP", q75),
            loaf_sum_median=("LoAFsum", q50),
            loaf_any_pct=("LoAFany_500ms", lambda x: 100 * float(np.mean(np.array(x, dtype=float) > 0)) if len(x) else math.nan),
        )
        .reset_index()
    )

    grp.to_csv(out_dir / "summary_by_group.csv", index=False)
    print(f"[OK] Saved {out_dir / 'summary_by_group.csv'}")

    # 3) Charts per (scenario, cpu)
    for (scenario, cpu), dsub in df.groupby(["scenario", "cpu"]):
        order = sorted([v for v in dsub["variant"].dropna().unique()])
        if not order:
            continue

        # Boxplot
        fig = plt.figure(figsize=(8, 5))
        data = [dsub.loc[dsub["variant"] == v, "INP"].dropna().values for v in order]
        safe_boxplot(data, labels=order)
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
            if len(x):
                plt.step(x, f, where="post", label=str(v))
        plt.xlabel("INP (ms)")
        plt.ylabel("F(x)")
        plt.title(f"ECDF — scenario={scenario}, CPUx{int(cpu) if pd.notna(cpu) else 'NA'}")
        plt.legend()
        fig.tight_layout()
        fig.savefig(charts_dir / f"ecdf_{scenario}_CPUx{int(cpu) if pd.notna(cpu) else 'NA'}.png", dpi=150)
        plt.close(fig)

    print(f"[OK] Saved charts to {charts_dir}")

    # 4) Mann–Whitney U (auto pairs) + effect sizes
    results = []
    for (scenario, cpu), dsub in df.groupby(["scenario", "cpu"]):
        variants = dsub["variant"].dropna().unique().tolist()
        pairs = build_pairs_for_group(variants)
        if not pairs:
            continue

        pvals: List[float] = []
        tmp = []
        for (a, b) in pairs:
            opt = dsub.loc[dsub["variant"] == b, "INP"].dropna().values
            base = dsub.loc[dsub["variant"] == a, "INP"].dropna().values
            if len(opt) == 0 or len(base) == 0:
                pvals.append(None)
                tmp.append((a, b, math.nan, math.nan, math.nan, math.nan, len(base), len(opt)))
                continue

            # Force asymptotic method (ties are common in INP due to quantization)
            res = mannwhitneyu(opt, base, alternative="two-sided", method="asymptotic")
            U = float(res.statistic)
            p = float(res.pvalue)

            n = len(opt)
            m = len(base)
            nm = n * m

            # SciPy's U for x=opt corresponds to P(opt > base) (+ 0.5 ties).
            # We want PS(A) as P(opt < base) (+ 0.5 ties):
            ps_opt_lt_base = 1.0 - (U / nm) if nm else math.nan
            delta = 2.0 * ps_opt_lt_base - 1.0 if math.isfinite(ps_opt_lt_base) else math.nan

            tmp.append((a, b, U, p, ps_opt_lt_base, delta, m, n))
            pvals.append(p)

        p_corr = holm_bonferroni(pvals) if pvals else []
        for (row, pc) in zip(tmp, p_corr):
            a, b, U, p, ps_opt_lt_base, delta, m, n = row
            results.append(
                {
                    "scenario": scenario,
                    "cpu": cpu,
                    "baseline": a,
                    "optimized": b,
                    "n_baseline": m,
                    "n_optimized": n,
                    "U": U,
                    "p_value": p,
                    "p_value_holm": pc,
                    "PS_prob_opt_lt_base": ps_opt_lt_base,  # renamed from "AUC" to avoid ROC-AUC confusion
                    "cliffs_delta": delta,
                }
            )

    if not results:
        print("No comparable pairs found for Mann-Whitney U.")
        return

    mw = pd.DataFrame(results)
    mw.to_csv(out_dir / "mw_results.csv", index=False)
    print(f"[OK] Saved {out_dir / 'mw_results.csv'}")

    # 5) Enrich MW results with p75(B0), p75(variant), Δp75%
    p75_tbl = grp[["scenario", "cpu", "variant", "p75_INP_ms"]].copy()

    mw_en = mw.merge(
        p75_tbl.rename(columns={"variant": "baseline", "p75_INP_ms": "p75_baseline_ms"}),
        on=["scenario", "cpu", "baseline"],
        how="left",
    )
    mw_en = mw_en.merge(
        p75_tbl.rename(columns={"variant": "optimized", "p75_INP_ms": "p75_optimized_ms"}),
        on=["scenario", "cpu", "optimized"],
        how="left",
    )

    mw_en["delta_p75_pct"] = [
        safe_pct(o, b) for o, b in zip(mw_en["p75_optimized_ms"], mw_en["p75_baseline_ms"])
    ]

    mw_en.to_csv(out_dir / "mw_results_enriched.csv", index=False)
    print(f"[OK] Saved {out_dir / 'mw_results_enriched.csv'}")

    # 6) Table like "1.4": number of CPU regimes (out of 4) with significant improvement vs B0
    b0_only = mw_en[mw_en["baseline"] == "B0"].copy()
    b0_only["sig_improve"] = (b0_only["p_value_holm"] < 0.05) & (b0_only["cliffs_delta"] > 0)

    counts = (
        b0_only[b0_only["sig_improve"]]
        .groupby(["scenario", "optimized"])["cpu"]
        .nunique()
        .reset_index()
        .rename(columns={"cpu": "cpu_regimes_with_improvement"})
    )

    table_1_4 = (
        counts.pivot(index="optimized", columns="scenario", values="cpu_regimes_with_improvement")
        .fillna(0)
        .astype(int)
        .reset_index()
        .rename(columns={"optimized": "variant"})
    )

    table_1_4.to_csv(out_dir / "table_1_4_counts.csv", index=False)
    print(f"[OK] Saved {out_dir / 'table_1_4_counts.csv'}")

    # 7) Regression summaries (p_Holm < 0.05 AND delta < 0)
    b0_only["sig_regress"] = (b0_only["p_value_holm"] < 0.05) & (b0_only["cliffs_delta"] < 0)

    reg_counts = (
        b0_only[b0_only["sig_regress"]]
        .groupby(["scenario", "optimized"])["cpu"]
        .nunique()
        .reset_index()
        .rename(columns={"cpu": "cpu_regimes_with_regression"})
    )

    table_reg = (
        reg_counts.pivot(index="optimized", columns="scenario", values="cpu_regimes_with_regression")
        .fillna(0)
        .astype(int)
        .reset_index()
        .rename(columns={"optimized": "variant"})
    )
    table_reg.to_csv(out_dir / "table_regressions_counts.csv", index=False)
    print(f"[OK] Saved {out_dir / 'table_regressions_counts.csv'}")

    reg_list = b0_only[b0_only["sig_regress"]].copy()
    reg_list = reg_list[
        [
            "scenario",
            "cpu",
            "optimized",
            "p75_baseline_ms",
            "p75_optimized_ms",
            "delta_p75_pct",
            "p_value",
            "p_value_holm",
            "PS_prob_opt_lt_base",
            "cliffs_delta",
        ]
    ].sort_values(["scenario", "cpu", "p_value_holm"])
    reg_list.to_csv(out_dir / "regressions_list.csv", index=False)
    print(f"[OK] Saved {out_dir / 'regressions_list.csv'}")

    # 8) Candidates for Table 1.5 (top improvements + top regressions per scenario)
    def pick_top(df_s, flag_col, top_n=5):
        df_s = df_s[df_s[flag_col]].copy()
        return df_s.sort_values("p_value_holm").head(top_n)

    cand = []
    for scen, df_s in b0_only.groupby("scenario"):
        best = pick_top(df_s, "sig_improve", top_n=5)
        worst = pick_top(df_s, "sig_regress", top_n=5)
        best["kind"] = "improvement"
        worst["kind"] = "regression"
        cand.append(best)
        cand.append(worst)

    cand_df = pd.concat(cand, ignore_index=True) if cand else pd.DataFrame()
    if len(cand_df):
        cand_df = cand_df[
            [
                "kind",
                "scenario",
                "cpu",
                "optimized",
                "p75_baseline_ms",
                "p75_optimized_ms",
                "delta_p75_pct",
                "p_value_holm",
                "cliffs_delta",
                "PS_prob_opt_lt_base",
            ]
        ]
        cand_df.to_csv(out_dir / "table_1_5_candidates.csv", index=False)
        print(f"[OK] Saved {out_dir / 'table_1_5_candidates.csv'}")


if __name__ == "__main__":
    main()
