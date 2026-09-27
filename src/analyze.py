"""Reproduce a balanced 2 x 3 x 2 factorial analysis of catalyst activity.

Inference is conditional on the 12 tested catalyst batches. The two observations
per cell are separate electrode tests, NOT independently synthesized batches.
Run from any directory: python path/to/src/analyze.py
"""
from __future__ import annotations

import argparse
import json
from itertools import combinations, product
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
from scipy import stats
import statsmodels
from statsmodels.formula.api import ols
from statsmodels.stats.anova import anova_lm
from statsmodels.stats.multitest import multipletests
from patsy import build_design_matrices

ROOT = Path(__file__).resolve().parents[1]
FACTORS = ["Loading_level", "Carbon", "Locality"]
LEVELS = [["15", "40"], ["Disordered", "Graphitized", "FCX-B"], ["IN", "OUT"]]
FORMULA = "SA ~ C(Loading_level, Sum) * C(Carbon, Sum) * C(Locality, Sum)"
TERMS = {
    "C(Loading_level, Sum)": "Loading",
    "C(Carbon, Sum)": "Carbon",
    "C(Locality, Sum)": "Locality",
    "C(Loading_level, Sum):C(Carbon, Sum)": "Loading x Carbon",
    "C(Loading_level, Sum):C(Locality, Sum)": "Loading x Locality",
    "C(Carbon, Sum):C(Locality, Sum)": "Carbon x Locality",
    "C(Loading_level, Sum):C(Carbon, Sum):C(Locality, Sum)": "Loading x Carbon x Locality",
    "Residual": "Residual",
}
SCOPE = "Conditional electrode-level analysis of these batches; no synthesis-level replication."
UNIT = r"Specific activity ($\mu$A cm$^{-2}_{Pt}$)"


def load_data(path: Path = ROOT / "data/catalyst_activity.csv") -> pd.DataFrame:
    """Preserve supplied measurements; extract nominal factors from sample IDs."""
    df = pd.read_csv(path)
    required = ["Sample", "Loading", "Locality", "Replicate", "ECSA_CO", "ECSA_TEM", "MA", "SA"]
    if set(df.columns) != set(required) or df[required].isna().any().any():
        raise ValueError("The supplied data must have all eight columns with no missing values.")
    numeric = ["Loading", "Replicate", "ECSA_CO", "ECSA_TEM", "MA", "SA"]
    if not np.isfinite(df[numeric].to_numpy(dtype=float)).all():
        raise ValueError("All numeric observations must be finite.")
    if (df[["Loading", "ECSA_CO", "ECSA_TEM", "MA", "SA"]] <= 0).any().any():
        raise ValueError("This dataset expects positive physical measurements.")
    parts = df["Sample"].str.extract(r"^Pt_(IN|OUT)-(15|40)/FCX-(Disordered|Graphitized|B)$")
    if parts.isna().any().any() or not parts[0].eq(df["Locality"]).all():
        raise ValueError("Sample identifiers and locality labels are inconsistent.")
    df["Loading_level"] = parts[1]
    df["Carbon"] = parts[2].replace({"B": "FCX-B"})
    for name, levels in zip(FACTORS, LEVELS):
        df[name] = pd.Categorical(df[name], categories=levels, ordered=True)
    if df.duplicated(["Sample", "Replicate"]).any():
        raise ValueError("Sample / replicate combinations must be unique.")
    counts = df.groupby(FACTORS, observed=False).size()
    if len(df) != 24 or len(counts) != 12 or not counts.eq(2).all():
        raise ValueError("This teaching analysis requires all 12 cells, with two tests per cell.")
    if not df.groupby("Sample")["Replicate"].apply(lambda x: set(x) == {1, 2}).all():
        raise ValueError("Each sample must have replicate labels 1 and 2.")
    if not df.groupby("Sample")["Loading"].nunique().eq(1).all():
        raise ValueError("Measured loading is expected to be fixed within each catalyst sample.")
    # Audit only: do not replace the provided response with rounded-input ratios.
    df["SA_ratio_audit"] = 100 * df["MA"] / df["ECSA_CO"]
    df["SA_audit_difference"] = df["SA"] - df["SA_ratio_audit"]
    return df


def fit_model(df: pd.DataFrame):
    return ols(FORMULA, data=df, missing="raise").fit()


def anova_table(model, response: np.ndarray) -> pd.DataFrame:
    """Type III with sum contrasts: equal-weight marginal hypotheses here."""
    table = anova_lm(model, typ=3).drop(index="Intercept")
    table = table.rename(index=TERMS, columns={"sum_sq": "SS", "df": "df", "PR(>F)": "p_raw"})
    table["MS"] = table["SS"] / table["df"]
    sst = np.sum((response - np.mean(response)) ** 2)
    sse = table.loc["Residual", "SS"]
    table["eta_squared"] = table["SS"] / sst
    table["partial_eta_squared"] = table["SS"] / (table["SS"] + sse)
    table.loc["Residual", "partial_eta_squared"] = np.nan
    terms = table.index != "Residual"
    table.loc[terms, "p_holm_7_terms"] = multipletests(table.loc[terms, "p_raw"], method="holm")[1]
    return table[["SS", "df", "MS", "F", "p_raw", "p_holm_7_terms", "eta_squared", "partial_eta_squared"]]


def manual_decomposition(df: pd.DataFrame) -> pd.DataFrame:
    """Independently compute the orthogonal SS from marginal cell means.

    No regression or statsmodels calculation is used in this verification.
    Each effect is a marginal mean minus its lower-order components.
    """
    means = df.groupby(FACTORS, observed=False)["SA"].mean().to_numpy().reshape(2, 3, 2)
    components = {(): np.full((1, 1, 1), means.mean())}
    labels = ["Loading", "Carbon", "Locality"]
    rows = []
    for size in [1, 2, 3]:
        for subset in combinations(range(3), size):
            average_axes = tuple(i for i in range(3) if i not in subset)
            effect = means.mean(axis=average_axes, keepdims=True) if average_axes else means.copy()
            for previous, component in components.items():
                if set(previous).issubset(subset):
                    effect = effect - component
            components[subset] = effect
            ss = 2 * np.sum(np.broadcast_to(effect, means.shape) ** 2)
            df_term = np.prod([means.shape[i] - 1 for i in subset])
            rows.append({"term": " x ".join(labels[i] for i in subset), "SS_manual": ss, "df_manual": df_term})
    fitted_means = df.groupby(FACTORS, observed=False)["SA"].transform("mean")
    rows.append({"term": "Residual", "SS_manual": np.sum((df["SA"] - fitted_means) ** 2), "df_manual": 12})
    return pd.DataFrame(rows).set_index("term")


def cell_summary(df: pd.DataFrame, model) -> pd.DataFrame:
    table = df.groupby(FACTORS, observed=False)["SA"].agg(n="size", mean="mean", sd="std").reset_index()
    # Pooled model uncertainty, NOT a CI computed separately from n=2 in each cell.
    table["se_pooled"] = np.sqrt(model.mse_resid / table["n"])
    half_width = stats.t.ppf(0.975, model.df_resid) * table["se_pooled"]
    table["ci95_low_pointwise"] = table["mean"] - half_width
    table["ci95_high_pointwise"] = table["mean"] + half_width
    return table


def contrast_result(model, rows: list[dict], weights: list[float]) -> dict:
    """A contrast is a weighted difference of fitted cell means, L beta."""
    design = np.asarray(build_design_matrices([model.model.data.design_info], pd.DataFrame(rows))[0])
    contrast = np.asarray(weights) @ design
    test = model.t_test(contrast)
    interval = np.asarray(test.conf_int()).ravel()
    scalar = lambda x: float(np.asarray(x).item())
    return {"estimate": scalar(test.effect), "se": scalar(test.sd), "t": scalar(test.tvalue),
            "df": model.df_resid, "p_raw": scalar(test.pvalue),
            "ci95_low_pointwise": interval[0], "ci95_high_pointwise": interval[1]}


def simple_locality_effects(model) -> pd.DataFrame:
    rows = []
    for carbon, loading in product(LEVELS[1], LEVELS[0]):
        cells = [{"Carbon": carbon, "Loading_level": loading, "Locality": loc} for loc in ["OUT", "IN"]]
        rows.append({"Carbon": carbon, "Loading_level": loading, "contrast": "OUT - IN",
                     **contrast_result(model, cells, [1, -1])})
    table = pd.DataFrame(rows)
    table["p_holm_6_comparisons"] = multipletests(table["p_raw"], method="holm")[1]
    return table


def interaction_contrasts(model) -> pd.DataFrame:
    rows = []
    for carbon in LEVELS[1]:
        cells = [{"Carbon": carbon, "Loading_level": load, "Locality": loc}
                 for load, loc in [("40", "OUT"), ("40", "IN"), ("15", "OUT"), ("15", "IN")]]
        rows.append({"Carbon": carbon, "contrast": "(OUT-IN)40 - (OUT-IN)15",
                     **contrast_result(model, cells, [1, -1, -1, 1])})
    table = pd.DataFrame(rows)
    table["p_holm_3_comparisons"] = multipletests(table["p_raw"], method="holm")[1]
    return table


def average_effects(model) -> pd.DataFrame:
    rows = []
    for factor, positive, negative in [("Loading_level", "40", "15"), ("Locality", "OUT", "IN")]:
        cells = [dict(zip(FACTORS, values)) for values in product(*LEVELS)]
        n_positive = sum(c[factor] == positive for c in cells)
        weights = [1 / n_positive if c[factor] == positive else -1 / n_positive for c in cells]
        rows.append({"factor": factor, "contrast": f"{positive} - {negative}",
                     **contrast_result(model, cells, weights)})
    return pd.DataFrame(rows)


def figures(df, model, cells, effects, out: Path) -> None:
    out.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                         "axes.spines.top": False, "axes.spines.right": False, "savefig.dpi": 180})
    colors = {"IN": "#007f83", "OUT": "#c06022"}
    fig, axes = plt.subplots(1, 3, figsize=(12, 4.8), sharey=True)
    for ax, carbon in zip(axes, LEVELS[1]):
        for loc in LEVELS[2]:
            subset = cells[(cells.Carbon == carbon) & (cells.Locality == loc)]
            shift = -0.025 if loc == "IN" else 0.025
            x = np.arange(2) + shift
            ax.errorbar(x, subset["mean"], yerr=subset["mean"] - subset["ci95_low_pointwise"],
                        fmt="o-", color=colors[loc], lw=2, capsize=4, label=loc)
            for k, loading in enumerate(LEVELS[0]):
                raw = df[(df.Carbon == carbon) & (df.Locality == loc) & (df.Loading_level == loading)]
                ax.scatter(k + shift + np.array([-0.035, 0.035]), raw.SA, s=26,
                           color=colors[loc], alpha=0.6, marker="x", zorder=4)
        ax.set(title=carbon, xticks=[0, 1], xticklabels=["15 wt%", "40 wt%"], xlabel="Nominal Pt loading", ylim=(10, 95))
        ax.grid(axis="y", alpha=0.18)
    axes[0].set_ylabel(UNIT)
    axes[-1].legend(frameon=False, title="Pt locality", loc="upper left")
    fig.suptitle("The OUT–IN activity gap changes with loading and carbon", fontsize=15, y=0.99)
    fig.text(0.5, 0.015, "Crosses: electrode tests. Circles: means. Bars: pointwise 95% CIs using pooled electrode error (12 df).\nConditional on the tested catalyst batches; n = 2 electrodes per cell.", ha="center", fontsize=9)
    fig.tight_layout(rect=[0, 0.11, 1, 0.94])
    for extension in ["png", "svg"]:
        fig.savefig(out / f"interaction_plot.{extension}", bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.5, 5.4))
    for i, row in effects.iterrows():
        color = "#365ba5" if row.Loading_level == "15" else "#a53b60"
        ax.errorbar(row.estimate, i, xerr=[[row.estimate-row.ci95_low_pointwise], [row.ci95_high_pointwise-row.estimate]],
                    fmt="o", color=color, capsize=4, markersize=7)
    ax.axvline(0, color="#777777", linestyle="--", linewidth=1)
    ax.set(yticks=range(len(effects)), yticklabels=[f"{r.Carbon} | {r.Loading_level} wt%" for r in effects.itertuples()],
           xlabel=r"OUT minus IN specific activity ($\mu$A cm$^{-2}_{Pt}$)", title="Six locality comparisons in the original activity units")
    ax.invert_yaxis()
    ax.grid(axis="x", alpha=0.18)
    fig.text(0.5, 0.015, "Pointwise 95% CIs, not simultaneous intervals. Holm-adjusted p-values are reported separately.\nUncertainty reflects electrode variation within the tested batches.", ha="center", fontsize=9)
    fig.tight_layout(rect=[0, 0.08, 1, 1])
    fig.savefig(out / "locality_effects.png", bbox_inches="tight")
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))
    axes[0].scatter(model.fittedvalues, model.resid, color="#007f83", s=40)
    axes[0].axhline(0, color="#777777", linestyle="--")
    axes[0].set(xlabel="Fitted SA (cell mean)", ylabel="Residual SA", title="Residuals versus fitted values")
    stats.probplot(model.resid, dist="norm", plot=axes[1])
    axes[1].set_title("Normal Q–Q plot (descriptive)")
    fig.text(0.5, 0.01, "With two tests per cell, residuals form opposite pairs by construction. These plots cannot establish independence or normality.", ha="center", fontsize=8.5)
    fig.tight_layout(rect=[0, 0.07, 1, 1])
    fig.savefig(out / "residual_diagnostics.png", bbox_inches="tight")
    plt.close(fig)


def run_analysis(data_path: Path = ROOT / "data/catalyst_activity.csv", output_path: Path = ROOT / "outputs") -> dict:
    df = load_data(data_path)
    model = fit_model(df)
    table = anova_table(model, df.SA.to_numpy())
    manual = manual_decomposition(df)
    np.testing.assert_allclose(table.SS, manual.loc[table.index, "SS_manual"], rtol=1e-10, atol=1e-9)
    np.testing.assert_allclose(table.SS.sum(), np.sum((df.SA - df.SA.mean()) ** 2))
    np.testing.assert_allclose(model.fittedvalues, df.groupby(FACTORS, observed=False).SA.transform("mean"))
    assert model.df_resid == 12 and np.linalg.matrix_rank(model.model.exog) == 12
    cells = cell_summary(df, model)
    effects = simple_locality_effects(model)
    # Verify contrast direction against direct arithmetic, independent of design coding.
    wide = cells.pivot(index=["Carbon", "Loading_level"], columns="Locality", values="mean")
    direct = [wide.loc[(r.Carbon, r.Loading_level), "OUT"] - wide.loc[(r.Carbon, r.Loading_level), "IN"] for r in effects.itertuples()]
    np.testing.assert_allclose(effects.estimate, direct)
    interactions = interaction_contrasts(model)
    average = average_effects(model)
    log_model = ols(FORMULA.replace("SA ~", "np.log(SA) ~"), data=df).fit()
    log_table = anova_table(log_model, np.log(df.SA).to_numpy())
    tables_path = output_path / "tables"
    tables_path.mkdir(parents=True, exist_ok=True)
    for filename, value, index in [
        ("anova.csv", table, True), ("manual_ss_verification.csv", manual, True),
        ("cell_summary.csv", cells, False), ("locality_contrasts.csv", effects, False),
        ("interaction_contrasts.csv", interactions, False), ("average_effects.csv", average, False),
        ("log_sa_sensitivity.csv", log_table, True),
        ("sa_ratio_audit.csv", df[["Sample", "Replicate", "SA", "SA_ratio_audit", "SA_audit_difference"]], False),
        ("measured_loading_by_sample.csv", df[["Sample", "Carbon", "Loading_level", "Locality", "Loading"]].drop_duplicates(), False),
    ]:
        value.to_csv(tables_path / filename, index=index, float_format="%.10g")
    figures(df, model, cells, effects, output_path / "figures")
    info = {"observations": len(df), "cells": len(cells), "electrode_tests_per_cell": 2,
            "scope": SCOPE, "formula": FORMULA, "residual_df": int(model.df_resid),
            "residual_ms": model.mse_resid, "residual_sd": np.sqrt(model.mse_resid),
            "R_squared_descriptive": model.rsquared,
            "max_absolute_sa_ratio_discrepancy": df.SA_audit_difference.abs().max(),
            "design_rank": int(np.linalg.matrix_rank(model.model.exog)),
            "rank_after_adding_actual_loading": int(np.linalg.matrix_rank(np.column_stack([model.model.exog, df.Loading]))),
            "versions": {"numpy": np.__version__, "pandas": pd.__version__, "scipy": scipy.__version__,
                         "statsmodels": statsmodels.__version__, "matplotlib": matplotlib.__version__}}
    (output_path / "run_info.json").write_text(json.dumps(info, indent=2) + "\n", encoding="utf-8")
    return {"data": df, "model": model, "anova": table, "cells": cells, "locality": effects,
            "interaction": interactions, "average": average, "log_anova": log_table, "info": info}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=ROOT / "data/catalyst_activity.csv")
    parser.add_argument("--out", type=Path, default=ROOT / "outputs")
    args = parser.parse_args()
    results = run_analysis(args.data, args.out)
    print(SCOPE)
    print(results["anova"].to_string(float_format=lambda x: f"{x:.6g}"))
    print("\nOUT - IN comparisons:\n", results["locality"].to_string(index=False))
    print("\nIndependent SS and contrast arithmetic checks passed.")
