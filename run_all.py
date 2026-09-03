"""
Run every experiment in the study and write results/ and figures/.

    python run_all.py
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from offtarget_panel.model import (
    PanelModel,
    concentrated_panel,
    PIC50_FLOOR,
    panel_npv_from_hits,
    spread_panel,
    to_pic50,
)

ROOT = Path(__file__).resolve().parent
FIG = ROOT / "figures"
RES = ROOT / "results"
FIG.mkdir(exist_ok=True)
RES.mkdir(exist_ok=True)

T = 1000
PANEL = 50
HIT_RATE = 0.003
NU = 5.0
SEED = 20260903

plt.rcParams.update(
    {
        "figure.dpi": 150,
        "savefig.bbox": "tight",
        "font.size": 9,
        "axes.grid": True,
        "grid.alpha": 0.3,
        "axes.spines.top": False,
        "axes.spines.right": False,
    }
)


def base_model(rho: float, **kw) -> PanelModel:
    return PanelModel(n_targets=T, rho_global=rho, nu=NU, hit_rate=HIT_RATE, **kw)


# ----------------------------------------------------------------------
# Experiment 0: what the simulated affinity landscape looks like
# ----------------------------------------------------------------------


def experiment_0() -> None:
    rng = np.random.default_rng(SEED)
    fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.2))

    for ax, rho in zip(axes[:2], [0.0, 0.6]):
        m = base_model(rho)
        a, s = np.sqrt(m.rho_global), np.sqrt(1 - m.rho_global)
        for i in range(6):
            phi = rng.standard_normal()
            zc = a * phi + s * rng.standard_t(m.nu, T) * np.sqrt((m.nu - 2) / m.nu)
            ax.plot(np.arange(1, T + 1), np.sort(to_pic50(zc))[::-1], lw=1.2, alpha=0.85)
        ax.axhline(to_pic50(m.theta), color="k", lw=1, ls="--")
        ax.text(T * 0.02, to_pic50(m.theta) + 0.15, "liability threshold", fontsize=7)
        ax.axhline(PIC50_FLOOR, color="0.3", lw=0.8)
        ax.text(T * 0.02, PIC50_FLOOR + 0.06, "assay floor: no measurable affinity", fontsize=6.5)
        ax.set_ylim(3.8, 9.2)
        ax.set_title(f"$\\rho$ = {rho:.1f}   (6 compounds)")
        ax.set_xscale("log")
        ax.set_xlabel("off-targets ranked by affinity")
        ax.set_ylabel("pIC50")

    m = base_model(0.6)
    a, s = np.sqrt(m.rho_global), np.sqrt(1 - m.rho_global)
    phi = rng.standard_normal(300_000)
    zz = a * phi + s * rng.standard_t(m.nu, 300_000) * np.sqrt((m.nu - 2) / m.nu)
    axes[2].hist(to_pic50(zz), bins=180, range=(3.8, 10), density=True, color="0.4")
    axes[2].axvline(to_pic50(m.theta), color="k", lw=1, ls="--")
    axes[2].set_yscale("log")
    axes[2].set_xlabel("pIC50")
    axes[2].set_ylabel("density (log)")
    axes[2].set_title("marginal affinity distribution")

    fig.suptitle("Simulated affinity landscape across 1000 physiological off-targets", y=1.03)
    fig.savefig(FIG / "fig0_landscape.png")
    plt.close(fig)


# ----------------------------------------------------------------------
# Experiment 1: how much does a clean 50-target panel tell you?
# ----------------------------------------------------------------------


def experiment_1() -> pd.DataFrame:
    rows = []
    for rho in np.round(np.arange(0.0, 0.95, 0.05), 3):
        m = base_model(float(rho))
        rows.append(
            {
                "rho": rho,
                "theta": m.theta,
                "base_rate_clean": m.p_clean(T),
                "p_clean_panel": m.p_clean(PANEL),
                "npv": m.npv(PANEL),
                "expected_missed": m.expected_missed_hits(PANEL),
                "residual_risk_ratio": m.residual_risk_ratio(PANEL),
                "variance_to_mean": m.variance_to_mean(),
            }
        )
    df = pd.DataFrame(rows)
    df["lift"] = df["npv"] / df["base_rate_clean"]
    df.to_csv(RES / "exp1_rho_sweep.csv", index=False)

    fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.2))
    axes[0].plot(df["rho"], df["npv"], "o-", ms=3, label="NPV of clean 50-panel")
    axes[0].plot(df["rho"], df["base_rate_clean"], "s--", ms=3, color="0.5",
                 label="base rate (no screen)")
    axes[0].set_xlabel("shared variance share $\\rho$")
    axes[0].set_ylabel("P(clean on unscreened 950)")
    axes[0].legend(fontsize=7)
    axes[0].set_title("Negative predictive value")

    axes[1].plot(df["rho"], df["residual_risk_ratio"], "o-", ms=3, color="C1")
    axes[1].axhline(1.0, color="0.5", ls="--", lw=1)
    axes[1].set_ylim(0, 1.1)
    axes[1].set_xlabel("shared variance share $\\rho$")
    axes[1].set_ylabel("residual risk ratio")
    axes[1].set_title("Per-target risk on the unscreened 950")

    axes[2].plot(df["rho"], df["variance_to_mean"], "o-", ms=3)
    axes[2].axhline(1.0, color="0.5", ls="--", lw=1)
    axes[2].set_yscale("log")
    axes[2].set_xlabel("shared variance share $\\rho$")
    axes[2].set_ylabel("variance / mean of hit counts (log)")
    axes[2].set_title("Observable signature of $\\rho$")

    fig.savefig(FIG / "fig1_rho_sweep.png")
    plt.close(fig)
    return df


# ----------------------------------------------------------------------
# Experiment 1b: does the base liability rate change the conclusion?
# ----------------------------------------------------------------------


def experiment_1b() -> pd.DataFrame:
    rows = []
    for p in [0.0003, 0.001, 0.003, 0.01, 0.03]:
        for rho in [0.0, 0.3, 0.6]:
            m = PanelModel(n_targets=T, rho_global=rho, nu=NU, hit_rate=p)
            rows.append(
                {
                    "hit_rate": p,
                    "rho": rho,
                    "expected_liabilities": p * T,
                    "p_clean_panel": m.p_clean(PANEL),
                    "base_rate_clean": m.p_clean(T),
                    "npv": m.npv(PANEL),
                    "residual_risk_ratio": m.residual_risk_ratio(PANEL),
                }
            )
    df = pd.DataFrame(rows)
    df.to_csv(RES / "exp1b_hit_rate.csv", index=False)

    fig, ax = plt.subplots(figsize=(4.6, 3.2))
    for rho, g in df.groupby("rho"):
        ax.plot(g["expected_liabilities"], g["residual_risk_ratio"], "o-", ms=3,
                label=f"$\\rho$={rho}")
    ax.axhline(1.0, color="0.5", ls="--", lw=1)
    ax.set_xscale("log")
    ax.set_ylim(0, 1.1)
    ax.set_xlabel("expected liabilities per compound across 1000 targets")
    ax.set_ylabel("residual risk ratio")
    ax.set_title("A panel informs only when it can fire")
    ax.legend(fontsize=7)
    fig.savefig(FIG / "fig1b_hit_rate.png")
    plt.close(fig)
    return df


# ----------------------------------------------------------------------
# Experiment 2: panel size
# ----------------------------------------------------------------------


def experiment_2() -> pd.DataFrame:
    sizes = [5, 10, 20, 44, 50, 100, 200, 400, 800]
    rows = []
    for rho in [0.0, 0.3, 0.5, 0.7, 0.85]:
        m = base_model(rho)
        for n in sizes:
            rows.append(
                {
                    "rho": rho,
                    "panel_size": n,
                    "npv": m.npv(n),
                    "expected_missed": m.expected_missed_hits(n),
                    "residual_risk_ratio": m.residual_risk_ratio(n),
                }
            )
    df = pd.DataFrame(rows)
    df.to_csv(RES / "exp2_panel_size.csv", index=False)

    fig, axes = plt.subplots(1, 2, figsize=(8, 3.2))
    for rho, g in df.groupby("rho"):
        axes[0].plot(g["panel_size"], g["npv"], "o-", ms=3, label=f"$\\rho$={rho}")
        axes[1].plot(g["panel_size"], g["residual_risk_ratio"], "o-", ms=3, label=f"$\\rho$={rho}")
    axes[0].set_xscale("log")
    axes[0].set_xlabel("targets screened")
    axes[0].set_ylabel("P(rest clean | panel clean)")
    axes[0].legend(fontsize=7)
    axes[0].set_title("Returns to panel size")
    axes[1].set_xscale("log")
    axes[1].set_xlabel("targets screened")
    axes[1].set_ylabel("residual risk ratio")
    axes[1].set_ylim(0, 1.1)
    axes[1].axhline(1.0, color="0.5", ls="--", lw=1)
    axes[1].set_title("Diminishing returns to screening")
    fig.savefig(FIG / "fig2_panel_size.png")
    plt.close(fig)
    return df


# ----------------------------------------------------------------------
# Experiment 3: tail weight of the idiosyncratic term
# ----------------------------------------------------------------------


def experiment_3() -> pd.DataFrame:
    rows = []
    for nu in [3.0, 4.0, 5.0, 8.0, 15.0, 50.0, 200.0]:
        for rho in [0.3, 0.5, 0.7]:
            m = PanelModel(n_targets=T, rho_global=rho, nu=nu, hit_rate=HIT_RATE)
            rows.append(
                {
                    "nu": nu,
                    "rho": rho,
                    "npv": m.npv(PANEL),
                    "base_rate_clean": m.p_clean(T),
                    "residual_risk_ratio": m.residual_risk_ratio(PANEL),
                }
            )
    df = pd.DataFrame(rows)
    df.to_csv(RES / "exp3_tail_weight.csv", index=False)

    fig, ax = plt.subplots(figsize=(4.4, 3.2))
    for rho, g in df.groupby("rho"):
        ax.plot(g["nu"], g["residual_risk_ratio"], "o-", ms=3, label=f"$\\rho$={rho}")
    ax.set_xscale("log")
    ax.set_ylim(0, 1.05)
    ax.set_xlabel("idiosyncratic degrees of freedom $\\nu$  (low = heavier tail)")
    ax.set_ylabel("residual risk ratio")
    ax.set_title("Heavy idiosyncratic tails blunt the panel")
    ax.legend(fontsize=7)
    fig.savefig(FIG / "fig3_tail_weight.png")
    plt.close(fig)
    return df


# ----------------------------------------------------------------------
# Experiment 4: family structure and panel composition (Monte Carlo)
# ----------------------------------------------------------------------


def experiment_4(n_compounds: int = 200_000) -> pd.DataFrame:
    rng = np.random.default_rng(SEED + 1)
    rows = []
    configs = [
        ("global only", 0.6, 0.0),
        ("half global, half family", 0.3, 0.3),
        ("family only", 0.0, 0.6),
    ]
    for name, rg, rf in configs:
        m = PanelModel(
            n_targets=T, rho_global=rg, rho_family=rf, n_families=20, nu=NU, hit_rate=HIT_RATE
        )
        hits = m.simulate_hits(n_compounds, rng)
        for design, panel in [
            ("spread over 20 families", spread_panel(T, 20, PANEL)),
            ("concentrated in 1 family", concentrated_panel(T, 20, PANEL)),
        ]:
            r = panel_npv_from_hits(hits, panel)
            r.update({"structure": name, "rho_global": rg, "rho_family": rf, "design": design})
            rows.append(r)
        del hits
    df = pd.DataFrame(rows)
    df.to_csv(RES / "exp4_family_structure.csv", index=False)

    fig, ax = plt.subplots(figsize=(5.4, 3.2))
    labels = [c[0] for c in configs]
    x = np.arange(len(labels))
    w = 0.36
    for i, design in enumerate(df["design"].unique()):
        sub = df[df["design"] == design].set_index("structure").loc[labels]
        ax.bar(x + (i - 0.5) * w, sub["residual_risk_ratio"], w, label=design)
    ax.axhline(1.0, color="k", ls="--", lw=1, label="no information")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=7)
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("residual risk ratio")
    ax.set_title("Panel composition matters only under family structure")
    ax.legend(fontsize=7)
    fig.savefig(FIG / "fig4_family_structure.png")
    plt.close(fig)
    return df


# ----------------------------------------------------------------------
# Experiment 5: Monte Carlo check of the quadrature, and the 44-target diagnostic
# ----------------------------------------------------------------------


def experiment_5(n_compounds: int = 200_000) -> pd.DataFrame:
    rng = np.random.default_rng(SEED + 2)
    rows = []
    for rho in [0.0, 0.3, 0.6]:
        m = base_model(rho)
        hits = m.simulate_hits(n_compounds, rng)
        mc = panel_npv_from_hits(hits, np.arange(PANEL))
        counts = hits[:, :44].sum(axis=1)
        rows.append(
            {
                "rho": rho,
                "npv_analytic": m.npv(PANEL),
                "npv_mc": mc["npv"],
                "rrr_analytic": m.residual_risk_ratio(PANEL),
                "rrr_mc": mc["residual_risk_ratio"],
                "base_rate_analytic": m.p_clean(T),
                "base_rate_mc": mc["base_rate"],
                "vmr_analytic_1000": m.variance_to_mean(),
                "vmr_mc_44": float(counts.var() / counts.mean()),
                "frac_clean_on_44": float((counts == 0).mean()),
                "n_clean_panel": mc["n_clean_panel"],
            }
        )
        del hits
    df = pd.DataFrame(rows)
    df.to_csv(RES / "exp5_mc_check.csv", index=False)
    return df


# ----------------------------------------------------------------------
# Experiment 6: reading rho off a panel you have already run
# ----------------------------------------------------------------------


def experiment_6() -> pd.DataFrame:
    """
    The quantities above are only useful if rho can be estimated.  It can:
    shared structure over-disperses the per-compound hit count relative to the
    binomial count expected under independence, and the hit counts from any
    panel already run are enough to measure that.
    """
    rows = []
    for n_panel in [20, 44, 100]:
        for rho in np.round(np.arange(0.0, 0.95, 0.05), 3):
            m = base_model(float(rho))
            rows.append(
                {
                    "n_panel": n_panel,
                    "rho": rho,
                    "mean_hits_in_panel": m.hit_count_moments(n_panel)[0],
                    "vmr_in_panel": m.variance_to_mean(n_panel),
                    "residual_risk_ratio": m.residual_risk_ratio(PANEL),
                }
            )
    df = pd.DataFrame(rows)
    df.to_csv(RES / "exp6_dispersion_diagnostic.csv", index=False)

    fig, axes = plt.subplots(1, 2, figsize=(8, 3.2))
    for n_panel, g in df.groupby("n_panel"):
        axes[0].plot(g["rho"], g["vmr_in_panel"], "o-", ms=3, label=f"{n_panel} targets")
    axes[0].axhline(1.0, color="0.5", ls="--", lw=1)
    axes[0].set_xlabel("shared variance share $\\rho$")
    axes[0].set_ylabel("variance / mean of hit counts in panel")
    axes[0].set_title("What over-dispersion you would measure")
    axes[0].legend(fontsize=7)

    g = df[df["n_panel"] == 44]
    axes[1].plot(g["vmr_in_panel"], g["residual_risk_ratio"], "o-", ms=3, color="C2")
    axes[1].set_xlabel("measured variance / mean in a 44-target panel")
    axes[1].set_ylabel("implied residual risk ratio")
    axes[1].set_ylim(0, 1.05)
    axes[1].set_title("From a measurement to a claim")
    fig.savefig(FIG / "fig5_dispersion_diagnostic.png")
    plt.close(fig)
    return df


def main() -> None:
    experiment_0()
    d1 = experiment_1()
    d1b = experiment_1b()
    d2 = experiment_2()
    d3 = experiment_3()
    d4 = experiment_4()
    d5 = experiment_5()
    d6 = experiment_6()

    summary = {
        "panel_size": PANEL,
        "n_targets": T,
        "marginal_hit_rate": HIT_RATE,
        "nu": NU,
        "base_rate_clean_all_targets": float(d1.loc[d1.rho == 0.0, "base_rate_clean"].item()),
    }
    (RES / "summary.json").write_text(json.dumps(summary, indent=2))

    pd.set_option("display.width", 200)
    print("\n=== Experiment 1: rho sweep ===")
    print(d1.round(4).to_string(index=False))
    print("\n=== Experiment 1b: base liability rate ===")
    print(d1b.round(4).to_string(index=False))
    print("\n=== Experiment 2: panel size ===")
    print(d2.round(4).to_string(index=False))
    print("\n=== Experiment 3: tail weight ===")
    print(d3.round(4).to_string(index=False))
    print("\n=== Experiment 4: family structure ===")
    print(d4.round(4).to_string(index=False))
    print("\n=== Experiment 5: Monte Carlo check ===")
    print(d5.round(4).to_string(index=False))
    print("\n=== Experiment 6: dispersion diagnostic (44-target panel) ===")
    print(d6[d6.n_panel == 44].round(4).to_string(index=False))


if __name__ == "__main__":
    main()
