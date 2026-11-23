"""Plots the figures of docs/impurity/impurity.tex from the data written by `make doc impurity compute`
(kicked_top_spin_chain with every namelist in docs/impurity/figures/configs_spin_chain/), and writes
params.tex: the parameters of these namelists and the numbers quoted in the text, as LaTeX macros,
so that impurity.tex has no hard-coded values. One run per namelist, ordered by delta = deltaz (+ deltax)."""

import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

FIG_DIR = Path(__file__).resolve().parent
ROOT = FIG_DIR.parents[2]
CONFIG_DIR = FIG_DIR / "configs_spin_chain"
LABELS = ["A", "B", "C"]  # macro suffixes of the runs (delta = 0, small, large), LaTeX names have no digits
DEFAULTS = {"deltaz": 0.0, "deltax": 0.0, "imp_seed": 12345, "butterfly_effect": False, "s0": -1, "pol": "z"}

plt.rcParams.update({"font.size": 10, "axes.labelsize": 11, "savefig.bbox": "tight"})


def value(s):
    """Python value of one Fortran namelist item."""
    s = s.strip()
    if s[0] in "'\"":
        return s[1:-1]
    if s.lower() in (".true.", "t"):
        return True
    if s.lower() in (".false.", "f"):
        return False
    try:
        return int(s)
    except ValueError:
        return float(s.lower().replace("d", "e"))


def read_namelist(path):
    """All items of all groups of a namelist file as one dict (group names are not needed here)."""
    cfg = dict(DEFAULTS)
    for line in path.read_text().splitlines():
        line = re.sub(r"!(?=(?:[^']*'[^']*')*[^']*$).*", "", line).strip()  # comment outside quotes
        if "=" not in line or line.startswith("&"):
            continue
        key, rhs = line.split("=", 1)
        items = [value(x) for x in rhs.split(",") if x.strip()]
        cfg[key.strip().lower()] = items[0] if len(items) == 1 else items
    return cfg


def runs():
    """Config of every namelist, with its data directory, ordered by impurity width."""
    cfgs = [read_namelist(p) for p in sorted(CONFIG_DIR.glob("*.nml"))]
    for c in cfgs:
        c["data"] = ROOT / c["out_dir"] / "h_1" / "chain"
    if len({c["data"].resolve() for c in cfgs}) != len(cfgs):
        raise SystemExit(f"plot_impurity.py: the namelists in {CONFIG_DIR} must have different out_dir")
    return sorted(cfgs, key=lambda c: (c["deltaz"] + c["deltax"], c["deltaz"]))


def load(c, name):
    return np.loadtxt(c["data"] / name, ndmin=2)


def label(c):
    if c["deltaz"] == c["deltax"]:
        return rf"$\delta = {c['deltaz']:g}$"
    return rf"$\delta_z = {c['deltaz']:g}$, $\delta_x = {c['deltax']:g}$"


def fidelity(cs):
    """Fig. [fig:fidelity]: F(t) of all runs."""
    fig, ax = plt.subplots(figsize=(7, 3.5))
    for c in cs:
        f = load(c, "fidelity.dat")
        ax.plot(f[:, 0], f[:, 1], "o-", ms=2.5, lw=1, label=label(c))
    ax.axvline(cs[0]["t_r"], color="0.5", ls="--", lw=0.8)
    ax.set_ylim(-0.02, 1.02)
    ax.set_xlabel(r"$t$")
    ax.set_ylabel(r"$F(t)$")
    ax.legend(fontsize=8)
    fig.savefig(FIG_DIR / "impurity_fidelity.pdf")
    plt.close(fig)


def spin(cs):
    """Figs. [fig:J_a], [fig:J2_a]: <J_a>/j and <J_a^2>/j^2, a = x, y, z, one figure each."""
    for col, (a, power) in enumerate([(a, p) for p in (1, 2) for a in "xyz"], start=1):
        fig, ax = plt.subplots(figsize=(7, 3))
        for c in cs:
            s = load(c, "total_spin.dat")
            ax.plot(s[:, 0], s[:, col] / (c["n_q"] / 2) ** power, "o-", ms=2.5, lw=1, label=label(c))
        ax.axvline(cs[0]["t_r"], color="0.5", ls="--", lw=0.8)
        ax.set_xlabel(r"$t$")
        ax.set_ylabel(rf"$\langle J_{a}\rangle/j$" if power == 1 else rf"$\langle J_{a}^2\rangle/j^2$")
        ax.legend(fontsize=8)
        fig.savefig(FIG_DIR / f"impurity_J{'2' if power == 2 else ''}{a}.pdf")
        plt.close(fig)


def husimi(cs):
    """Fig. [fig:husimi]: Husimi distribution at t = 0, t_r, 2 t_r, one row per run."""
    t_r = cs[0]["t_r"]
    steps = [0, t_r, 2 * t_r]
    fig, axes = plt.subplots(len(cs), len(steps), figsize=(3.4 * len(steps), 2.9 * len(cs)),
                             sharex=True, sharey=True, squeeze=False)
    for row, c in zip(axes, cs):
        d = load(c, "husimi.dat")
        n_phi = c["n_phi"]
        for ax, t in zip(row, steps):
            b = d[d[:, 0] == t]
            shape = (-1, n_phi)
            ax.pcolormesh(b[:, 1].reshape(shape), b[:, 2].reshape(shape), b[:, 5].reshape(shape),
                          shading="nearest", cmap="viridis", rasterized=True)
            ax.set_xlim(-np.pi, np.pi)
            ax.set_ylim(-1, 1)
            ax.set_xticks([-np.pi, 0, np.pi], [r"$-\pi$", "0", r"$\pi$"])
            ax.set_title(f"{label(c)}, $t = {t}$", fontsize=9)
            ax.label_outer()
        row[0].set_ylabel(r"$Z=\cos\theta$")
    for ax in axes[-1]:
        ax.set_xlabel(r"$\phi$")
    fig.savefig(FIG_DIR / "impurity_husimi.pdf", dpi=200)
    plt.close(fig)


def sci(x):
    """x in LaTeX math, e.g. 3.2\\times10^{-14}."""
    m, e = f"{x:.1e}".split("e")
    return rf"{m}\times10^{{{int(e)}}}"


def params(cs):
    """params.tex: common parameters (must agree in all namelists) and per-run values A, B, C, ..."""
    common = ["n_q", "eps", "h", "t_r", "imp_seed", "s0", "pol", "butterfly_effect", "n_theta", "n_phi"]
    for k in common:
        if any(c[k] != cs[0][k] for c in cs):
            raise SystemExit(f"plot_impurity.py: {k} differs between the namelists in {CONFIG_DIR}")
    c0 = cs[0]
    if any("theta0" in c or "phi0" in c for c in cs):
        raise SystemExit("plot_impurity.py: impurity.tex assumes the initial state is given by s0 and pol, not theta0, phi0")
    if c0["butterfly_effect"]:
        raise SystemExit("plot_impurity.py: impurity.tex assumes butterfly_effect = .false.")
    macros = {
        "nq": c0["n_q"], "jval": f"{c0['n_q'] / 2:g}", "epsval": f"{c0['eps']:g}", "hval": f"{c0['h']:g}",
        "kval": f"{c0['n_q'] * c0['eps']:g}", "tr": c0["t_r"], "trtwo": 2 * c0["t_r"], "impseed": c0["imp_seed"],
        "initpol": c0["pol"], "initsign": "+" if c0["s0"] == -1 else "-",
        "ntheta": c0["n_theta"], "nphi": c0["n_phi"], "nruns": len(cs),
    }
    for L, c in zip(LABELS, cs):
        f = load(c, "fidelity.dat")[:, 1]
        macros[f"deltaz{L}"] = f"{c['deltaz']:g}"
        macros[f"deltax{L}"] = f"{c['deltax']:g}"
        macros[f"Ftr{L}"] = f"{f[c['t_r']]:.3g}"
        macros[f"Fend{L}"] = f"{f[-1]:.3g}"
        macros[f"Fdef{L}"] = sci(abs(1 - f[-1]))
        summary = dict(line.split() for line in (c["data"] / "summary.dat").read_text().splitlines()
                       if not line.startswith("#"))
        macros.setdefault("Qzero", f"{float(summary['husimi_norm(0)']):.4g}")
        macros[f"Qtr{L}"] = f"{float(summary['husimi_norm(t_r)']):.3g}"
        macros[f"Qend{L}"] = f"{float(summary['husimi_norm(2t_r)']):.3g}"
    lines = ["% Written by figures/plot_impurity.py from figures/configs_spin_chain/*.nml; do not edit."]
    lines += [rf"\newcommand{{\{k}}}{{{v}}}" for k, v in macros.items()]
    rows = [rf"  \deltaz{L} & \deltax{L} & \Ftr{L} & \Fend{L} & $\Fdef{L}$ & \Qtr{L} & \Qend{L} \\" for L in LABELS[:len(cs)]]
    lines += [r"\newcommand{\runrows}{%"] + rows + ["}"]
    (FIG_DIR / "params.tex").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    cs = runs()
    if len(cs) != 3 or cs[0]["deltaz"] != 0 or cs[0]["deltax"] != 0:
        raise SystemExit("plot_impurity.py: impurity.tex expects 3 namelists, the first with deltaz = deltax = 0")
    params(cs)
    fidelity(cs)
    spin(cs)
    husimi(cs)
