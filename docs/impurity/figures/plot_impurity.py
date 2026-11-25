"""Plots the figures of docs/impurity/impurity.tex from the data written by `make doc impurity compute`
(kicked_top_spin_chain with every namelist in docs/impurity/figures/configs_spin_chain/), and writes
params.tex: the parameters of these namelists and the numbers quoted in the text, as LaTeX macros,
so that impurity.tex has no hard-coded values. One run per namelist. The runs are grouped by the initial state
(theta0, phi0), in the order of the namelist file names, and every group has the same impurity widths, ordered by deltaz (+ deltax).
The namelists with butterfly_effect = .true. are the butterfly runs (no impurity, one qubit rotated at t_r): one figure
of F(t) per initial state and rotation axis, one line per rotation angle. The classical kicked top (kicked_top with every
namelist in docs/impurity/figures/configs/) gives the phase space at the chain's h and k = n_q*eps: one 'grid' namelist
for the Poincare section and one 'single' namelist per initial state for its orbit and Lyapunov exponent."""

import re
from fractions import Fraction
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

FIG_DIR = Path(__file__).resolve().parent
ROOT = FIG_DIR.parents[2]
CONFIG_DIR = FIG_DIR / "configs_spin_chain"
CLASSICAL_DIR = FIG_DIR / "configs"
LABELS = ["A", "B", "C"]  # macro suffixes of the impurity widths (delta = 0, small, large), LaTeX names have no digits
STATES = ["One", "Two"]  # macro suffixes of the initial states
ALPHAS = ["A", "B", "C", "D", "E", "F", "G", "H"]  # macro suffixes of the rotation angles of the butterfly runs
AXES = "xyz"
DEFAULTS = {"deltaz": 0.0, "deltax": 0.0, "imp_seed": 12345, "butterfly_effect": False, "i_flip": 1, "flip_axis": "x",
            "flip_angle": np.pi}

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
    """Configs of all namelists, with their data directory: the impurity runs as one list per initial state (in the
    order of the first namelist file name of each state), each ordered by impurity width, and the butterfly runs as
    a dict {(state, axis): runs ordered by flip_angle}."""
    cfgs = [read_namelist(p) for p in sorted(CONFIG_DIR.glob("*.nml"))]
    for c in cfgs:
        if "theta0" not in c or "phi0" not in c:
            raise SystemExit(f"plot_impurity.py: every namelist in {CONFIG_DIR} must set theta0 and phi0")
        c["data"] = ROOT / c["out_dir"] / "h_1" / "chain"
    if len({c["data"].resolve() for c in cfgs}) != len(cfgs):
        raise SystemExit(f"plot_impurity.py: the namelists in {CONFIG_DIR} must have different out_dir")
    imp = [c for c in cfgs if not c["butterfly_effect"]]
    bf = [c for c in cfgs if c["butterfly_effect"]]
    states = list(dict.fromkeys((c["theta0"], c["phi0"]) for c in imp))
    groups = [sorted([c for c in imp if (c["theta0"], c["phi0"]) == s], key=lambda c: (c["deltaz"] + c["deltax"], c["deltaz"]))
              for s in states]
    if any((c["theta0"], c["phi0"]) not in states or c["flip_axis"] not in AXES for c in bf):
        raise SystemExit(f"plot_impurity.py: every butterfly namelist in {CONFIG_DIR} must use an initial state of the"
                         " impurity runs and flip_axis 'x', 'y' or 'z'")
    butterfly = {(s, a): sorted([c for c in bf if (c["theta0"], c["phi0"]) == s and c["flip_axis"] == a],
                                key=lambda c: c["flip_angle"]) for s in states for a in AXES}
    return groups, butterfly


def load(c, name):
    return np.loadtxt(c["data"] / name, ndmin=2)


def label(c):
    if c["deltaz"] == c["deltax"]:
        return rf"$\delta = {c['deltaz']:g}$"
    return rf"$\delta_z = {c['deltaz']:g}$, $\delta_x = {c['deltax']:g}$"


def angle(x):
    """x in LaTeX math as a rational multiple of pi if it is one (e.g. \\pi/4), else as a number."""
    f = Fraction(x / np.pi).limit_denominator(64)
    if abs(float(f) * np.pi - x) > 1e-12:
        return f"{x:g}"
    if f == 0:
        return "0"
    num = {1: "", -1: "-"}.get(f.numerator, str(f.numerator))
    return rf"{num}\pi" + (f"/{f.denominator}" if f.denominator != 1 else "")


def state_label(c):
    return rf"$(\theta_0,\phi_0) = ({angle(c['theta0'])}, {angle(c['phi0'])})$"


def classical(groups):
    """Configs of the classical runs: the Poincare section and one orbit per initial state (in the order of groups).
    h and k must be those of the chain, k = n_q*eps."""
    cfgs = [read_namelist(p) for p in sorted(CLASSICAL_DIR.glob("*.nml"))]
    for c in cfgs:
        c["data"] = ROOT / c["out_dir"]
    c0 = groups[0][0]
    if any(abs(c["h"] - c0["h"]) > 1e-12 or abs(c["k"] - c0["n_q"] * c0["eps"]) > 1e-12 for c in cfgs):
        raise SystemExit(f"plot_impurity.py: every namelist in {CLASSICAL_DIR} must have h and k = n_q*eps of the chain")
    section = [c for c in cfgs if c["mode"] == "grid"]
    single = {(c["theta0"], c["phi0"]): c for c in cfgs if c["mode"] == "single"}
    states = [(g[0]["theta0"], g[0]["phi0"]) for g in groups]
    if len(section) != 1 or any(s not in single for s in states):
        raise SystemExit(f"plot_impurity.py: {CLASSICAL_DIR} needs one namelist with mode = 'grid' and one with"
                         " mode = 'single' for every initial state (theta0, phi0) of the chain")
    return section[0], [single[s] for s in states]


def phase_space(section, orbits):
    """Fig. [fig:classical]: Poincare section of the classical kicked top and the orbits of the two initial states."""
    fig, ax = plt.subplots(figsize=(6.5, 4.2))
    s = np.loadtxt(section["data"] / "phase_space.dat", ndmin=2)
    ax.plot(s[:, 2], s[:, 3], ",", color="k", alpha=0.4, rasterized=True)
    for c, color in zip(orbits, ["C0", "C3"]):
        o = np.loadtxt(c["data"] / "phase_space.dat", ndmin=2)[: section["n_kicks"] + 1]
        lam = np.loadtxt(c["data"] / "lyapunov.dat", ndmin=2)[0, 3]
        ax.plot(o[:, 2], o[:, 3], ".", color=color, ms=1.5, rasterized=True)
        ax.plot(o[0, 2], o[0, 3], "*", color=color, mec="k", ms=13, clip_on=False, label=rf"{state_label(c)}, $\Lambda = {num(lam)}$")
    ax.set_xlim(-np.pi, np.pi)
    ax.set_ylim(-1, 1)
    ax.set_xticks([-np.pi, -np.pi / 2, 0, np.pi / 2, np.pi], [r"$-\pi$", r"$-\pi/2$", "0", r"$\pi/2$", r"$\pi$"])
    ax.set_xlabel(r"$\phi$")
    ax.set_ylabel(r"$Z=\cos\theta$")
    ax.legend(loc="upper right", fontsize=8, framealpha=0.9)
    fig.savefig(FIG_DIR / "impurity_classical.pdf", dpi=200)
    plt.close(fig)


def fidelity(cs, S):
    """Fig. [fig:fidelity]: F(t) of all runs of one initial state."""
    fig, ax = plt.subplots(figsize=(4.6, 3))
    for c in cs:
        f = load(c, "fidelity.dat")
        ax.plot(f[:, 0], f[:, 1], "o-", ms=2, lw=1, label=label(c))
    ax.axvline(cs[0]["t_r"], color="0.5", ls="--", lw=0.8)
    ax.set_ylim(-0.02, 1.02)
    ax.set_xlabel(r"$t$")
    ax.set_ylabel(r"$F(t)$")
    ax.set_title(state_label(cs[0]), fontsize=10)
    ax.legend(fontsize=8)
    fig.savefig(FIG_DIR / f"impurity_fidelity_{S}.pdf")
    plt.close(fig)


def butterfly_fidelity(cs, S, a):
    """Fig. [fig:butterfly_a]: F(t) of the butterfly runs of one initial state and rotation axis a."""
    fig, ax = plt.subplots(figsize=(4.6, 3))
    for c in cs:
        f = load(c, "fidelity.dat")
        ax.plot(f[:, 0], f[:, 1], "o-", ms=2, lw=1, label=rf"$\alpha = {angle(c['flip_angle'])}$")
    ax.axvline(cs[0]["t_r"], color="0.5", ls="--", lw=0.8)
    ax.set_ylim(-0.02, 1.02)
    ax.set_xlabel(r"$t$")
    ax.set_ylabel(r"$F(t)$")
    ax.set_title(rf"{state_label(cs[0])}, rotation about ${a}$", fontsize=10)
    ax.legend(fontsize=7)
    fig.savefig(FIG_DIR / f"impurity_butterfly_{a}_{S}.pdf")
    plt.close(fig)


def spin(cs, S):
    """Figs. [fig:S_a], [fig:J2_a]: <S_a>/j and <S_a^2>/j^2, a = x, y, z, one figure each, for one initial state."""
    for col, (a, power) in enumerate([(a, p) for p in (1, 2) for a in "xyz"], start=1):
        fig, ax = plt.subplots(figsize=(4.6, 2.6))
        for c in cs:
            s = load(c, "total_spin.dat")
            ax.plot(s[:, 0], s[:, col] / (c["n_q"] / 2) ** power, "o-", ms=2, lw=1, label=label(c))
        ax.axvline(cs[0]["t_r"], color="0.5", ls="--", lw=0.8)
        ax.set_xlabel(r"$t$")
        ax.set_ylabel(rf"$\langle S_{a}\rangle/j$" if power == 1 else rf"$\langle S_{a}^2\rangle/j^2$")
        ax.set_title(state_label(cs[0]), fontsize=10)
        ax.legend(fontsize=7)
        fig.savefig(FIG_DIR / f"impurity_J{'2' if power == 2 else ''}{a}_{S}.pdf")
        plt.close(fig)


def husimi(cs, S):
    """Fig. [fig:husimi]: Husimi distribution at t = 0, t_r, 2 t_r, one row per run of one initial state."""
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
    fig.suptitle(state_label(cs[0]), fontsize=11)
    fig.savefig(FIG_DIR / f"impurity_husimi_{S}.pdf", dpi=200)
    plt.close(fig)


def sci(x):
    """x in LaTeX math, e.g. 3.2\\times10^{-14}."""
    m, e = f"{x:.1e}".split("e")
    return rf"{m}\times10^{{{int(e)}}}"


def num(x):
    """x in LaTeX math with 3 digits, in scientific notation below 10^-2."""
    return sci(x) if 0 < abs(x) < 1e-2 else f"{x:.3g}"


def params(groups, butterfly, section, orbits):
    """params.tex: common parameters (must agree in all namelists), the impurity widths A, B, C (must agree for
    all initial states), the rotation angles A, B, ... of the butterfly runs (must agree for all initial states and
    axes) and the per-run values, e.g. \\FendBTwo for width B and initial state Two."""
    cs = [c for g in groups for c in g]
    bfs = [c for g in butterfly.values() for c in g]
    for k in ["n_q", "eps", "h", "t_r", "n_theta", "n_phi"]:
        if any(c[k] != cs[0][k] for c in bfs):
            raise SystemExit(f"plot_impurity.py: {k} differs between the namelists in {CONFIG_DIR}")
    if any(c["deltaz"] != 0 or c["deltax"] != 0 or c["i_flip"] != bfs[0]["i_flip"] for c in bfs):
        raise SystemExit("plot_impurity.py: the butterfly runs must have deltaz = deltax = 0 and the same i_flip")
    alphas = [c["flip_angle"] for c in next(iter(butterfly.values()))]
    if not alphas or len(alphas) > len(ALPHAS) or any([c["flip_angle"] for c in g] != alphas for g in butterfly.values()):
        raise SystemExit(f"plot_impurity.py: every initial state and axis needs the same 1 .. {len(ALPHAS)} flip_angle values")
    common = ["n_q", "eps", "h", "t_r", "imp_seed", "butterfly_effect", "n_theta", "n_phi"]
    for k in common:
        if any(c[k] != cs[0][k] for c in cs):
            raise SystemExit(f"plot_impurity.py: {k} differs between the namelists in {CONFIG_DIR}")
    for g in groups:
        if [(c["deltaz"], c["deltax"]) for c in g] != [(c["deltaz"], c["deltax"]) for c in groups[0]]:
            raise SystemExit(f"plot_impurity.py: the initial states in {CONFIG_DIR} must have the same deltaz, deltax")
    c0 = cs[0]
    if c0["butterfly_effect"]:
        raise SystemExit("plot_impurity.py: impurity.tex assumes butterfly_effect = .false.")
    macros = {
        "nq": c0["n_q"], "jval": f"{c0['n_q'] / 2:g}", "epsval": f"{c0['eps']:g}", "hval": f"{c0['h']:g}",
        "kval": f"{c0['n_q'] * c0['eps']:g}", "tr": c0["t_r"], "trtwo": 2 * c0["t_r"], "impseed": c0["imp_seed"],
        "ntheta": c0["n_theta"], "nphi": c0["n_phi"], "nruns": len(groups[0]),
        "iflip": bfs[0]["i_flip"], "nalpha": len(alphas),
        "norbits": section["n_theta"] * section["n_phi"], "nkicksSection": section["n_kicks"],
        "nkicksState": orbits[0]["n_kicks"],
    }
    if any(c["n_kicks"] != orbits[0]["n_kicks"] for c in orbits):
        raise SystemExit(f"plot_impurity.py: the 'single' namelists in {CLASSICAL_DIR} must have the same n_kicks")
    for L, x in zip(ALPHAS, alphas):
        macros[f"alpha{L}"] = angle(x)
    for L, c in zip(LABELS, groups[0]):
        macros[f"deltaz{L}"] = f"{c['deltaz']:g}"
        macros[f"deltax{L}"] = f"{c['deltax']:g}"
    for S, g in zip(STATES, groups):
        macros[f"theta{S}"] = angle(g[0]["theta0"])
        macros[f"phi{S}"] = angle(g[0]["phi0"])
        macros[f"Lambda{S}"] = num(np.loadtxt(orbits[STATES.index(S)]["data"] / "lyapunov.dat", ndmin=2)[0, 3])
        for L, c in zip(LABELS, g):
            f = load(c, "fidelity.dat")[:, 1]
            macros[f"Ftr{L}{S}"] = num(f[c['t_r']])
            macros[f"Fend{L}{S}"] = num(f[-1])
            macros[f"Fdef{L}{S}"] = sci(abs(1 - f[-1]))
            summary = dict(line.split() for line in (c["data"] / "summary.dat").read_text().splitlines()
                           if not line.startswith("#"))
            macros.setdefault("Qzero", f"{float(summary['husimi_norm(0)']):.4g}")
    lines = ["% Written by figures/plot_impurity.py from figures/configs_spin_chain/*.nml and figures/configs/*.nml;"
             " do not edit."]
    lines += [rf"\newcommand{{\{k}}}{{{v}}}" for k, v in macros.items()]
    (FIG_DIR / "params.tex").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    groups, butterfly = runs()
    if len(groups) != 2 or any(len(g) != 3 or g[0]["deltaz"] != 0 or g[0]["deltax"] != 0 for g in groups):
        raise SystemExit("plot_impurity.py: impurity.tex expects 2 initial states (theta0, phi0) with 3 namelists each,"
                         " the first with deltaz = deltax = 0")
    section, orbits = classical(groups)
    params(groups, butterfly, section, orbits)
    phase_space(section, orbits)
    for S, g in zip(STATES, groups):
        for a in AXES:
            butterfly_fidelity(butterfly[(g[0]["theta0"], g[0]["phi0"]), a], S, a)
        fidelity(g, S)
        spin(g, S)
        husimi(g, S)
