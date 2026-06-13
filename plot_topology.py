"""Render a sample SAGIN topology to docs/topology.png.

Uses the same geometry + greedy best-gain association as the experiment, so the
figure is a faithful, reproducible snapshot of one random scenario (fixed seed).
Two panels: a top-down (x-y) view with association links, and a side elevation
(x-altitude) view that makes the three layers explicit.

    python plot_topology.py [--seed N]
"""

from __future__ import annotations

import argparse
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from channel import path_gain_db, random_positions
from config import Config

LEO_C, HAPS_C, GU_C = "#1f6fbf", "#ff9933", "#2e8b57"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int)
    args = ap.parse_args()

    cfg = Config()
    seed = args.seed if args.seed is not None else cfg.seed
    rng = np.random.default_rng(seed)

    leo = random_positions(rng, cfg.n_leo, cfg.coverage_km, cfg.height_leo_km)
    haps = random_positions(rng, cfg.n_haps, cfg.coverage_km, cfg.height_haps_km)
    gu = random_positions(rng, cfg.n_gu, cfg.coverage_km, cfg.height_gu_km)

    # Greedy best-large-scale-gain association (same rule as the greedy baselines).
    haps_to_leo = path_gain_db(
        leo, haps, cfg.freq_hz, cfg.antenna_gain_leo_haps_db
    ).argmax(axis=0)
    gu_to_haps = path_gain_db(
        haps, gu, cfg.freq_hz, cfg.antenna_gain_haps_gu_db
    ).argmax(axis=0)

    fig, (axt, axs) = plt.subplots(1, 2, figsize=(13, 5.6))

    # ---- top-down (x-y) with association links ----
    for h in range(cfg.n_haps):
        l = haps_to_leo[h]
        axt.plot(
            [haps[0, h], leo[0, l]],
            [haps[1, h], leo[1, l]],
            color=LEO_C,
            lw=0.8,
            alpha=0.45,
            zorder=1,
        )
    for g in range(cfg.n_gu):
        h = gu_to_haps[g]
        axt.plot(
            [gu[0, g], haps[0, h]],
            [gu[1, g], haps[1, h]],
            color=GU_C,
            lw=0.8,
            alpha=0.45,
            zorder=1,
        )
    axt.scatter(
        leo[0],
        leo[1],
        marker="*",
        s=300,
        c=LEO_C,
        edgecolors="k",
        linewidths=0.5,
        label=f"LEO ×{cfg.n_leo}  (300 km)",
        zorder=3,
    )
    axt.scatter(
        haps[0],
        haps[1],
        marker="^",
        s=130,
        c=HAPS_C,
        edgecolors="k",
        linewidths=0.5,
        label=f"HAPS ×{cfg.n_haps}  (20 km)",
        zorder=3,
    )
    axt.scatter(
        gu[0],
        gu[1],
        marker="o",
        s=48,
        c=GU_C,
        edgecolors="k",
        linewidths=0.4,
        label=f"Ground user ×{cfg.n_gu}  (0 km)",
        zorder=3,
    )
    axt.set_xlabel("x (km)")
    axt.set_ylabel("y (km)")
    axt.set_title("Top-down view — greedy best-gain association")
    axt.set_xlim(-5, 105)
    axt.set_ylim(-5, 105)
    axt.set_aspect("equal", "box")
    axt.legend(loc="upper right", fontsize=9, framealpha=0.9)
    axt.grid(True, ls=":", alpha=0.4)

    # ---- side elevation (x vs altitude) ----
    for z, c in ((cfg.height_leo_km, LEO_C), (cfg.height_haps_km, HAPS_C), (0, GU_C)):
        axs.axhline(z, color=c, ls="--", lw=0.7, alpha=0.5)
    axs.scatter(
        leo[0],
        leo[2],
        marker="*",
        s=300,
        c=LEO_C,
        edgecolors="k",
        linewidths=0.5,
        zorder=3,
    )
    axs.scatter(
        haps[0],
        haps[2],
        marker="^",
        s=130,
        c=HAPS_C,
        edgecolors="k",
        linewidths=0.5,
        zorder=3,
    )
    axs.scatter(
        gu[0], gu[2], marker="o", s=48, c=GU_C, edgecolors="k", linewidths=0.4, zorder=3
    )
    axs.set_xlabel("x (km)")
    axs.set_ylabel("altitude (km)")
    axs.set_title("Side elevation — three layers")
    axs.set_xlim(-5, 105)
    axs.set_ylim(-18, 330)
    axs.grid(True, ls=":", alpha=0.3)

    fig.suptitle(
        f"Sample SAGIN topology — {cfg.n_leo} LEO, {cfg.n_haps} HAPS, {cfg.n_gu} GU  (seed {seed})",
        fontsize=13,
        fontweight="bold",
    )
    fig.tight_layout(rect=[0, 0, 1, 0.95])

    out = os.path.join(os.path.dirname(__file__), "docs", "topology.png")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print("wrote", out)


if __name__ == "__main__":
    main()
