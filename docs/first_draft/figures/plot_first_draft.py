"""Figures and numbers of docs/first_draft/first_draft.tex from the data written by
`make doc first_draft compute` (kicked_top_scan with configs_scan/scan.nml):
initial-condition scan, phase portrait with the two selected states, spin-j scan series,
chain echo fidelity, Husimi snapshots, flip-angle/axis and t_r scans.
Writes figures/*.pdf and figures/numbers.tex (one \\newcommand per number quoted in the draft)."""

import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

FIG_DIR = Path(__file__).resolve().parent
DATA_DIR = FIG_DIR / "data"
CONFIG = FIG_DIR / "configs_scan" / "scan.nml"
STATES = [("island", "island", "C0"), ("sea", "chaotic sea", "C3")]

plt.rcParams.update({"font.size": 10, "axes.labelsize": 11, "savefig.bbox": "tight"})


def namelist(path):
    """Values of the scalar entries of a namelist file as {group: {name: value}} (strings without quotes)."""
    groups, group = {}, None
    for line in path.read_text().splitlines():
        line = line.split("!")[0].strip()
        if line.startswith("&"):
            group = groups.setdefault(line[1:].split()[0].lower(), {})
            line = line[1 + len(line[1:].split()[0]):]
        for name, value in re.findall(r"(\w+)\s*=\s*('[^']*'|[^,\s/]+)", line):
            group[name.lower()] = value.strip("'")
    return groups


# Defaults of src/scan/scan_config_mod.f90 for entries the namelist omits.
DEFAULTS = {
    "system": {"h": str(np.pi / 2), "k": "2.4", "n_q": "12"},
    "scan": {"n_theta": "40", "n_phi": "80", "n_lyap": "2000", "n_renorm": "10", "t_scan": "200",
             "lam_regular": "0.01", "lam_chaotic": "0.1", "sat_tol": "0.05"},
    "classical": {"n_kicks": "1000", "n_theta": "12", "n_phi": "12"},
    "chain": {"t_r": "30", "i_flip": "1", "flip_axis": "y", "flip_angle": str(np.pi),
              "n_angle": "0", "tr_step": "0", "tr_max": "0"},
}
NML = {g: {**d, **namelist(CONFIG).get(g, {})} for g, d in DEFAULTS.items()}
CFG = {**NML["system"], **NML["scan"], **NML["chain"]}
H, K, NQ = float(CFG["h"]), float(CFG["k"]), int(CFG["n_q"])
J = NQ / 2
T_R = int(CFG["t_r"])


def load(path):
    return np.loadtxt(DATA_DIR / path, ndmin=2)


def summary(path):
    """name -> value of a chain summary.dat."""
    out = {}
    for line in (DATA_DIR / path).read_text().splitlines():
        if not line.startswith("#"):
            name, value = line.split()
            out[name] = float(value)
    return out


def selected():
    """name -> row of selected.dat: ic theta0 phi0 Z0 Lambda F(t_scan) PR(t_scan) F_late PR_early PR_late t_E t_scan/t_E F_random."""
    out = {}
    for line in (DATA_DIR / "selected.dat").read_text().splitlines():
        if not line.startswith("#"):
            w = line.split()
            out[w[0]] = np.array([float(x) for x in w[1:]])
    return out


SEL = selected()


def phase_axes(ax, labels=True):
    ax.set_xlim(-np.pi, np.pi)
    ax.set_ylim(-1, 1)
    ax.set_xticks([-np.pi, -np.pi / 2, 0, np.pi / 2, np.pi], [r"$-\pi$", r"$-\pi/2$", "0", r"$\pi/2$", r"$\pi$"])
    if labels:
        ax.set_xlabel(r"$\phi$ (rad)")
        ax.set_ylabel(r"$Z=\cos\theta$")


def mark(ax, size=11):
    for run, label, color in STATES:
        r = SEL[run]
        ax.plot(r[2], r[3], "*", color=color, mec="k", mew=0.8, ms=size, label=label)


def scan_maps():
    """Fig. [fig:scan]: Lambda, F_late, PR_late and PR_late/PR_early over the grid of initial states."""
    d = load("scan.dat")
    n_phi = len(np.unique(np.round(d[:, 2], 12)))
    grid = lambda c: d[:, c].reshape(-1, n_phi)
    phi, z = grid(2), grid(3)
    panels = [(grid(4), r"(a) $\Lambda$", "magma", None),
              (grid(7), r"(b) $\overline{F}_{\rm late}$", "viridis", (0, 1)),
              (grid(9) / (2 * J + 1), r"(c) $\overline{P}_{\rm late}/(2j+1)$", "viridis", (0, 0.6)),
              (grid(9) / grid(8), r"(d) $\overline{P}_{\rm late}/\overline{P}_{\rm early}$", "coolwarm", (0.8, 1.2))]
    fig, axes = plt.subplots(2, 2, figsize=(9, 6.2), sharex=True, sharey=True)
    for ax, (v, title, cmap, lim) in zip(axes.flat, panels):
        im = ax.pcolormesh(phi, z, v, shading="nearest", cmap=cmap, rasterized=True,
                           vmin=None if lim is None else lim[0], vmax=None if lim is None else lim[1])
        fig.colorbar(im, ax=ax, pad=0.02)
        phase_axes(ax)
        mark(ax)
        ax.set_title(title)
        ax.label_outer()
    axes[0, 0].legend(loc="lower left", fontsize=8, framealpha=0.9)
    fig.savefig(FIG_DIR / "scan.pdf", dpi=200)


def section():
    """Fig. [fig:section]: Poincare section with the two selected initial states and their orbits."""
    fig, ax = plt.subplots(figsize=(6, 4))
    s = load("section/phase_space.dat")
    ax.plot(s[:, 2], s[:, 3], ",", color="k", alpha=0.5, rasterized=True)
    for run, label, color in STATES:
        o = load(f"{run}/classical/phase_space.dat")
        ax.plot(o[:, 2], o[:, 3], ".", color=color, ms=1.2, alpha=0.6, rasterized=True)
        ax.plot(o[0, 2], o[0, 3], "*", color=color, mec="k", ms=13, label=label)
    phase_axes(ax)
    ax.legend(loc="upper right", fontsize=8, framealpha=0.9)
    fig.savefig(FIG_DIR / "section.pdf", dpi=200)


def scan_series():
    """Fig. [fig:scan-series]: F(t) and PR(t) of the two selected states under the spin-j Floquet operator."""
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.3))
    t_scan = None
    for run, label, color in STATES:
        d = load(f"{run}/scan_series.dat")
        t_scan = int(d[-1, 0])
        axes[0].plot(d[:, 0], d[:, 1], "-", color=color, lw=0.8, label=label)
        axes[1].plot(d[:, 0], d[:, 2], "-", color=color, lw=0.8, label=label)
    axes[0].axhline(1 / (2 * J + 1), color="k", ls=":", lw=1, label=r"$1/(2j+1)$")
    axes[1].axhline(J + 1, color="k", ls=":", lw=1, label=r"$j+1$")
    t_e = SEL["sea"][10]
    for ax in axes:
        ax.axvspan(t_scan / 2, 3 * t_scan / 4, color="0.9", lw=0)
        ax.axvspan(3 * t_scan / 4, t_scan, color="0.8", lw=0)
        ax.axvline(t_e, color="C3", ls="--", lw=0.8)
        ax.set_xlabel("$t$ (kicks)")
        ax.set_xlim(0, t_scan)
    axes[0].set_ylabel("$F(t)$")
    axes[0].set_ylim(0, 1.02)
    axes[1].set_ylabel("$P(t)$")
    axes[1].set_ylim(0, 2 * J + 1)
    axes[0].legend(fontsize=8, loc="upper right")
    axes[1].legend(fontsize=8, loc="lower right")
    fig.savefig(FIG_DIR / "scan_series.pdf", dpi=200)


def echo_fidelity():
    """Fig. [fig:echo]: chain fidelity F(t), forward to t_r, flip, backward to 2 t_r; with and without flip."""
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.3))
    for run, label, color in STATES:
        f = load(f"{run}/chain/fidelity.dat")
        g = load(f"{run}/chain_no_flip/fidelity.dat")
        for ax in axes:
            ax.plot(f[:, 0], f[:, 1], "o-", color=color, ms=2.5, lw=1, label=f"{label}, flip")
            ax.plot(g[:, 0], g[:, 1], "--", color=color, lw=1, label=f"{label}, no flip")
    for ax in axes:
        ax.axvline(T_R, color="0.5", ls=":", lw=0.8)
        ax.axhline(1 / (2 * J + 1), color="k", ls=":", lw=0.8)
        ax.set_xlabel("$t$ (kicks)")
        ax.set_xlim(0, 2 * T_R)
    axes[0].set_ylim(0, 1.05)
    axes[0].set_ylabel("$F(t)$")
    axes[0].set_title("(a) linear scale")
    axes[1].set_yscale("log")
    axes[1].set_ylim(1e-6, 2)
    axes[1].set_title("(b) logarithmic scale")
    axes[0].legend(fontsize=7, loc="lower left")
    fig.savefig(FIG_DIR / "echo_fidelity.pdf", dpi=200)


def echo_husimi():
    """Fig. [fig:husimi]: Husimi distribution of the symmetric part of the chain state at t = 0, t_r/2, t_r, 3t_r/2, 2t_r."""
    blocks, titles = [], None
    for run, _, _ in STATES:
        d = load(f"{run}/chain/husimi.dat")
        ts = np.unique(d[:, 0]).astype(int)
        blocks.append([d[d[:, 0] == t] for t in ts])
        titles = [f"$t = {t}$" for t in ts]
    fig, axes = plt.subplots(2, len(titles), figsize=(3.0 * len(titles), 5.4), sharex=True, sharey=True)
    for r, row in enumerate(blocks):
        vmax = max(b[:, 5].max() for b in row)
        for c, b in enumerate(row):
            ax = axes[r, c]
            n_phi = len(np.unique(np.round(b[:, 1], 12)))
            im = ax.pcolormesh(b[:, 1].reshape(-1, n_phi), b[:, 2].reshape(-1, n_phi), b[:, 5].reshape(-1, n_phi),
                               shading="nearest", cmap="viridis", vmin=0, vmax=vmax, rasterized=True)
            phase_axes(ax)
            ax.plot(SEL[STATES[r][0]][2], SEL[STATES[r][0]][3], "+", color="w", ms=9, mew=1.2)
            ax.set_title(f"{STATES[r][1]}, {titles[c]}", fontsize=9)
            ax.label_outer()
        fig.colorbar(im, ax=axes[r, :].tolist(), pad=0.01, label="$Q$")
    fig.savefig(FIG_DIR / "echo_husimi.pdf", dpi=150)


def echo_scans():
    """Fig. [fig:echo-scans]: F(2t_r) versus flip angle and axis, and versus t_r."""
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.3))
    markers = {"x": "o", "y": "s", "z": "^"}
    for run, label, color in STATES:
        # Only scans enabled in the config (stale files of disabled scans are ignored).
        path = DATA_DIR / run / "flip_scan.dat"
        if int(CFG["n_angle"]) >= 2 and path.exists():
            rows = [l.split() for l in path.read_text().splitlines() if l.strip() and not l.startswith("#")]
            for a in "xyz":
                v = np.array([[float(x) for x in r[1:]] for r in rows if r[0] == a])
                axes[0].plot(v[:, 0], v[:, 1], markers[a], color=color, mfc="none", ms=5, label=f"{label}, $a={a}$")
                axes[0].plot(v[:, 0], v[:, 2], "-", color=color, lw=0.6)
        path = DATA_DIR / run / "tr_scan.dat"
        v = load(f"{run}/tr_scan.dat") if int(CFG["tr_step"]) > 0 and path.exists() else np.empty((0, 0))
        if v.size:
            axes[1].semilogy(v[:, 0], v[:, 2], "o-", color=color, ms=3, lw=1, label=f"{label}, $F(2t_r)$")
            axes[1].semilogy(v[:, 0], v[:, 4], "s--", color=color, mfc="none", ms=3, lw=0.8,
                             label=f"{label}, symmetric weight")
    axes[0].set_xticks([0, np.pi / 4, np.pi / 2, 3 * np.pi / 4, np.pi],
                       ["0", r"$\pi/4$", r"$\pi/2$", r"$3\pi/4$", r"$\pi$"])
    axes[0].set_xlabel(r"flip angle $\vartheta$ (rad)")
    axes[0].set_ylabel("$F(2t_r)$")
    axes[0].set_ylim(0, 1.05)
    axes[0].set_title(f"(a) $t_r = {T_R}$")
    axes[0].legend(fontsize=6.5, ncol=2, loc="lower left")
    axes[1].set_xlabel("$t_r$ (kicks)")
    axes[1].set_title(r"(b) $\vartheta = \pi$, $a = " + CFG["flip_axis"] + "$")
    axes[1].legend(fontsize=7, loc="lower right")
    fig.savefig(FIG_DIR / "echo_scans.pdf", dpi=200)


def fmt(x, digits=2):
    """x with `digits` significant digits, in LaTeX."""
    if x == 0:
        return "0"
    e = int(np.floor(np.log10(abs(x))))
    if -3 <= e <= 3:
        return f"{x:.{max(digits - 1 - e, 0)}f}"
    m = x / 10**e
    return rf"{m:.{digits - 1}f}\times 10^{{{e}}}"


def numbers():
    """figures/numbers.tex: \\newcommand{\\Num...}{value} for every number quoted in the draft."""
    d = load("scan.dat")
    out = {
        "Hval": fmt(H, 3), "Kval": fmt(K, 3), "NqVal": str(NQ), "Jval": fmt(J, 3) if J % 1 else str(int(J)),
        "EpsVal": fmt(K / NQ, 3), "Ngrid": str(d.shape[0]), "NthetaScan": CFG["n_theta"], "NphiScan": CFG["n_phi"],
        "NthetaSection": NML["classical"]["n_theta"], "NphiSection": NML["classical"]["n_phi"],
        "NkicksSection": NML["classical"]["n_kicks"],
        "Tscan": str(int(load("island/scan_series.dat")[-1, 0])), "Tr": str(T_R),
        "Iflip": CFG["i_flip"], "FlipAxis": CFG["flip_axis"],
        "Frandom": fmt(1 / (2 * J + 1), 2), "PRrandom": fmt(J + 1, 3), "PRcoherent": fmt((4 * J + 1) / (2 * J + 1), 3),
        "ChaoticFraction": fmt(np.mean(d[:, 4] > float(CFG["lam_chaotic"])), 2),
        "LamRegular": fmt(float(CFG["lam_regular"]), 1), "LamChaotic": fmt(float(CFG["lam_chaotic"]), 1),
        "SatTol": fmt(float(CFG["sat_tol"]), 1), "Nlyap": CFG["n_lyap"],
    }
    for run, tag in [("island", "Island"), ("sea", "Sea")]:
        r = SEL[run]
        out.update({f"{tag}ThetaZero": fmt(r[1], 3), f"{tag}PhiZero": fmt(r[2], 3), f"{tag}Zzero": fmt(r[3], 2),
                    f"{tag}Lambda": fmt(r[4], 2), f"{tag}Fscan": fmt(r[5], 2), f"{tag}PRscan": fmt(r[6], 2),
                    f"{tag}Flate": fmt(r[7], 2), f"{tag}PRearly": fmt(r[8], 3), f"{tag}PRlate": fmt(r[9], 3),
                    f"{tag}PRratio": fmt(r[9] / r[8], 3)})
        s, s0 = summary(f"{run}/chain/summary.dat"), summary(f"{run}/chain_no_flip/summary.dat")
        out.update({f"{tag}Ftr": fmt(s["F(t_r)"], 2), f"{tag}Fecho": fmt(s["F(2t_r)"], 2),
                    f"{tag}FechoPred": fmt(s["F(2t_r)_predicted"], 2),
                    f"{tag}FechoNoFlipDev": fmt(abs(1 - s0["F(2t_r)"]), 1),
                    f"{tag}SymWeight": fmt(s["symmetric_weight(2t_r)"], 2),
                    f"{tag}SymWeightPred": fmt(s["sym_weight_predicted"], 2),
                    f"{tag}HusimiNormEnd": fmt(s["husimi_norm(2t_r)"], 2),
                    f"{tag}PRtr": fmt(s["participation(t_r)"], 2)})
    out["SeaTehrenfest"] = fmt(SEL["sea"][10], 2)
    out["SeaTscanOverTe"] = fmt(SEL["sea"][11], 2)
    lines = ["% Generated by plot_first_draft.py from figures/data; do not edit."]
    lines += [rf"\newcommand{{\Num{name}}}{{{value}}}" for name, value in out.items()]
    (FIG_DIR / "numbers.tex").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    scan_maps()
    section()
    scan_series()
    echo_fidelity()
    echo_husimi()
    echo_scans()
    numbers()
