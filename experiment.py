"""End-to-end experiment: train Q-learning vs baselines over many topologies.

Usage
-----
    python experiment.py            # full run  (see Config defaults)
    python experiment.py --quick    # fast smoke run
    python experiment.py --seed 7 --episodes 30000 --topologies 20

Outputs (written to ./results/):
    learning_curve.png   averaged training curve with +/- 1 std band
    comparison.png       grouped bars: RL vs baselines (tier1 / tier2 / total)
    results.csv          per-topology numbers for every method
    summary.txt          mean +/- std table printed at the end
"""

from __future__ import annotations

import argparse
import csv
import os
from dataclasses import replace

import numpy as np

import baselines
from channel import make_topology
from config import Config
from trainer import eval_policy, train_one_topology

RESULTS_DIR = os.path.join(os.path.dirname(__file__), "results")


def parse_args() -> Config:
    p = argparse.ArgumentParser(
        description="SAGIN association+power Q-learning experiment"
    )
    p.add_argument("--quick", action="store_true", help="fast low-fidelity run")
    p.add_argument("--seed", type=int)
    p.add_argument("--episodes", type=int)
    p.add_argument("--topologies", type=int)
    p.add_argument("--eval-samples", type=int)
    a = p.parse_args()

    cfg = Config().quick() if a.quick else Config()
    overrides = {}
    if a.seed is not None:
        overrides["seed"] = a.seed
    if a.episodes is not None:
        overrides["episodes"] = a.episodes
    if a.topologies is not None:
        overrides["n_topologies"] = a.topologies
    if a.eval_samples is not None:
        overrides["eval_fading_samples"] = a.eval_samples
    return replace(cfg, **overrides) if overrides else cfg


def run(cfg: Config):
    master = np.random.default_rng(cfg.seed)
    curves = np.empty((cfg.n_topologies, cfg.episodes))
    methods = ["Q-Learning", "Random", "MaxPower+Greedy", "BestUniform+Greedy"]
    totals = {m: [] for m in methods}
    per_tier = {m: [] for m in methods}

    for t in range(cfg.n_topologies):
        rng = np.random.default_rng(master.integers(2**31))
        topo = make_topology(cfg, rng)

        curve, policy = train_one_topology(cfg, topo, rng)
        curves[t] = curve

        rl = eval_policy(cfg, topo, policy, rng, cfg.eval_fading_samples)
        base = baselines.run_all(cfg, topo, rng, cfg.eval_fading_samples)
        results = {"Q-Learning": rl, **base}

        for m in methods:
            r1, r2 = results[m]
            per_tier[m].append((r1, r2))
            totals[m].append(r1 + r2)

        print(
            f"topo {t + 1:>2}/{cfg.n_topologies}  "
            + "  ".join(f"{m}={sum(results[m]):6.2f}" for m in methods)
        )

    return curves, totals, per_tier, methods


def summarise(totals, methods) -> str:
    lines = ["method                 total sum-rate (bits/s/Hz)", "-" * 52]
    for m in methods:
        arr = np.array(totals[m])
        lines.append(f"{m:<22} {arr.mean():7.3f}  +/- {arr.std():5.3f}")
    rl = np.array(totals["Q-Learning"]).mean()
    best_base = max(np.array(totals[m]).mean() for m in methods if m != "Q-Learning")
    lines.append("-" * 52)
    lines.append(
        f"Q-Learning gain over best baseline: {100 * (rl / best_base - 1):+.1f}%"
    )
    return "\n".join(lines)


def make_plots(cfg, curves, per_tier, methods):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # ---- learning curve (smoothed mean +/- std across topologies) ----
    win = max(1, cfg.episodes // 200)
    kernel = np.ones(win) / win
    smoothed = np.array([np.convolve(c, kernel, mode="valid") for c in curves])
    mean = smoothed.mean(axis=0)
    std = smoothed.std(axis=0)
    x = np.arange(mean.size)

    fig, ax = plt.subplots(figsize=(7, 4.2))
    ax.plot(x, mean, color="#1f77b4", label="Q-Learning (mean)")
    ax.fill_between(
        x, mean - std, mean + std, color="#1f77b4", alpha=0.2, label="+/- 1 std"
    )
    ax.set_xlabel("episode")
    ax.set_ylabel("total sum-rate (bits/s/Hz)")
    ax.set_title("Q-learning convergence (averaged over topologies)")
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(RESULTS_DIR, "learning_curve.png"), dpi=130)
    plt.close(fig)

    # ---- grouped bar comparison ----
    tiers = ["tier1 (LEO->HAPS)", "tier2 (HAPS->GU)", "total"]
    vals = {m: np.array(per_tier[m]) for m in methods}  # [topo, 2]
    means = {
        m: [vals[m][:, 0].mean(), vals[m][:, 1].mean(), vals[m].sum(axis=1).mean()]
        for m in methods
    }
    errs = {
        m: [vals[m][:, 0].std(), vals[m][:, 1].std(), vals[m].sum(axis=1).std()]
        for m in methods
    }

    fig, ax = plt.subplots(figsize=(8, 4.5))
    width = 0.2
    xpos = np.arange(len(tiers))
    for i, m in enumerate(methods):
        ax.bar(
            xpos + (i - 1.5) * width, means[m], width, yerr=errs[m], capsize=3, label=m
        )
    ax.set_xticks(xpos)
    ax.set_xticklabels(tiers)
    ax.set_ylabel("sum-rate (bits/s/Hz)")
    ax.set_title("Q-learning vs baselines")
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(RESULTS_DIR, "comparison.png"), dpi=130)
    plt.close(fig)


def save_csv(totals, per_tier, methods):
    path = os.path.join(RESULTS_DIR, "results.csv")
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["method", "topology", "tier1", "tier2", "total"])
        for m in methods:
            for i, (r1, r2) in enumerate(per_tier[m]):
                w.writerow([m, i, f"{r1:.4f}", f"{r2:.4f}", f"{r1 + r2:.4f}"])


def main():
    cfg = parse_args()
    os.makedirs(RESULTS_DIR, exist_ok=True)
    print(
        f"config: n_leo={cfg.n_leo} n_haps={cfg.n_haps} n_gu={cfg.n_gu} "
        f"episodes={cfg.episodes} topologies={cfg.n_topologies} seed={cfg.seed}\n"
    )

    curves, totals, per_tier, methods = run(cfg)
    summary = summarise(totals, methods)
    print("\n" + summary)

    save_csv(totals, per_tier, methods)
    with open(os.path.join(RESULTS_DIR, "summary.txt"), "w") as f:
        f.write(summary + "\n")
    make_plots(cfg, curves, per_tier, methods)
    print(f"\nartifacts written to {RESULTS_DIR}/")


if __name__ == "__main__":
    main()
