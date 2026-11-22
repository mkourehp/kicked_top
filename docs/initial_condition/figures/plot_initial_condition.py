"""Plots the four figures of docs/initial_condition/initial_condition.tex from the data written by
`make doc initial_condition compute`: classical Poincare section, quantum Husimi distribution,
spin-chain fidelity echo and spin-chain Husimi distribution, all at h = pi/2, k = 2.4
(eps = k/(2j) = 0.2, j = 6, n_q = 12) for one initial point in an island and one in the chaotic sea."""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

FIG_DIR = Path(__file__).resolve().parent
DATA_DIR = FIG_DIR / "data"
STATES = [("island", "island", "C0"), ("sea", "chaotic sea", "C3")]

plt.rcParams.update({"font.size": 10, "axes.labelsize": 11, "savefig.bbox": "tight"})


def load(run, name):
    """Numeric data of output file name of run."""
    return np.loadtxt(DATA_DIR / run / name, ndmin=2)


def phase_axes(ax):
    ax.set_xlim(-np.pi, np.pi)
    ax.set_ylim(-1, 1)
    ax.set_xticks([-np.pi, -np.pi / 2, 0, np.pi / 2, np.pi], [r"$-\pi$", r"$-\pi/2$", "0", r"$\pi/2$", r"$\pi$"])
    ax.set_xlabel(r"$\phi$")
    ax.set_ylabel(r"$Z=\cos\theta$")


def husimi_panels(blocks, titles, name):
    """One row per initial state, one panel per time; blocks[row][col] are rows (t, phi, Z, X, Y, Q)."""
    fig, axes = plt.subplots(len(blocks), len(blocks[0]), figsize=(3.2 * len(blocks[0]), 2.9 * len(blocks)),
                             sharex=True, sharey=True, squeeze=False)
    for r, row in enumerate(blocks):
        for c, b in enumerate(row):
            ax = axes[r, c]
            n_phi = len(np.unique(b[:, 1]))
            ax.pcolormesh(b[:, 1].reshape(-1, n_phi), b[:, 2].reshape(-1, n_phi), b[:, 5].reshape(-1, n_phi),
                          shading="nearest", cmap="viridis", rasterized=True)
            phase_axes(ax)
            ax.set_title(f"{STATES[r][1]}, {titles[c]}")
            ax.label_outer()
    fig.savefig(FIG_DIR / name, dpi=150)


def section():
    """Fig. [fig:section]: Poincare section and the two orbits started at the initial points."""
    fig, ax = plt.subplots(figsize=(6, 4))
    s = load("section", "phase_space.dat")
    ax.plot(s[:, 2], s[:, 3], ",", color="k", alpha=0.5, rasterized=True)
    for run, label, color in STATES:
        o = load(f"{run}_classical", "phase_space.dat")
        lam = load(f"{run}_classical", "lyapunov.dat")[0, 3]
        ax.plot(o[:, 2], o[:, 3], ".", color=color, ms=1.5, rasterized=True)
        ax.plot(o[0, 2], o[0, 3], "*", color=color, mec="k", ms=12, label=rf"{label}, $\Lambda = {lam:.2g}$")
    phase_axes(ax)
    ax.legend(loc="upper right", fontsize=8, framealpha=0.9)
    fig.savefig(FIG_DIR / "section.pdf", dpi=200)


def quantum_husimi():
    """Fig. [fig:quantum-husimi]: Husimi distribution of the kicked top, j = 6."""
    blocks, titles = [], None
    for run, _, _ in STATES:
        d = load(f"{run}_quantum", "husimi.dat")
        kicks = np.unique(d[:, 0]).astype(int)
        blocks.append([d[d[:, 0] == n] for n in kicks])
        titles = [f"$n = {n}$" for n in kicks]
    husimi_panels(blocks, titles, "quantum_husimi.pdf")


def chain_fidelity():
    """Fig. [fig:chain-fidelity]: fidelity F(t) = |<psi_0|psi(t)>|^2, forward to t_r, then backward."""
    fig, ax = plt.subplots(figsize=(6, 3.5))
    for run, label, color in STATES:
        f = load(f"{run}_chain/h_1/chain", "fidelity.dat")
        ax.plot(f[:, 0], f[:, 1], "o-", color=color, ms=2.5, lw=1, label=label)
    t_r = int(round(f[-1, 0] / 2))
    ax.axvline(t_r, color="0.5", ls="--", lw=0.8)
    ax.set_xlabel("$t$")
    ax.set_ylabel("$F(t)$")
    ax.set_ylim(0, 1.05)
    ax.legend()
    fig.savefig(FIG_DIR / "chain_fidelity.pdf", dpi=200)


def chain_husimi():
    """Fig. [fig:chain-husimi]: Husimi distribution of the collective spin of the chain, n_q = 12."""
    blocks, titles = [], None
    for run, _, _ in STATES:
        d = load(f"{run}_chain/h_1/chain", "husimi.dat")
        ts = np.unique(d[:, 0]).astype(int)
        blocks.append([d[d[:, 0] == t] for t in ts])
        titles = [f"$t = {t}$" for t in ts]
    husimi_panels(blocks, titles, "chain_husimi.pdf")


if __name__ == "__main__":
    section()
    quantum_husimi()
    chain_fidelity()
    chain_husimi()
