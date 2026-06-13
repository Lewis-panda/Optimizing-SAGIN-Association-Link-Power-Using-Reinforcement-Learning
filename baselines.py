"""Non-learning baselines, evaluated on the same topology / fading as the agent.

    random            : random association + random power  (lower bound / sanity)
    max_power_greedy  : best mean-gain association, everyone at P_max
    best_uniform      : best mean-gain association, single best *common* power level

``best_uniform`` is the strong baseline: if RL cannot beat the best single power
applied to everyone, its per-node power control buys nothing. Beating it is the
whole point -- and it is only possible because interference makes P_max suboptimal.
"""

from __future__ import annotations

import numpy as np

from channel import Topology, noise_power_w, sample_gain
from env import sum_rate


def _greedy_assoc(topo: Topology) -> tuple[np.ndarray, np.ndarray]:
    """Each receiver associates with the transmitter of highest large-scale gain."""
    a_haps = topo.pg_leo_haps_db.argmax(axis=0)  # [n_haps] -> LEO index
    a_gu = topo.pg_haps_gu_db.argmax(axis=0)  # [n_gu]   -> HAPS index
    return a_haps, a_gu


def _mean_rate(cfg, topo, a_ha, a_lp, a_ga, a_hp, rng, n_samples):
    noise_w = noise_power_w(cfg.bandwidth_hz, cfg.noise_figure_db)
    r1 = np.empty(n_samples)
    r2 = np.empty(n_samples)
    for k in range(n_samples):
        g1 = sample_gain(topo.pg_leo_haps_db, rng)
        g2 = sample_gain(topo.pg_haps_gu_db, rng)
        r1[k], _ = sum_rate(a_lp, a_ha, g1, noise_w)
        r2[k], _ = sum_rate(a_hp, a_ga, g2, noise_w)
    return float(r1.mean()), float(r2.mean())


def random_policy(cfg, topo, rng, n_samples) -> tuple[float, float]:
    powers = cfg.power_levels_dbm
    noise_w = noise_power_w(cfg.bandwidth_hz, cfg.noise_figure_db)
    r1 = np.empty(n_samples)
    r2 = np.empty(n_samples)
    for k in range(n_samples):
        a_ha = rng.integers(0, cfg.n_leo, cfg.n_haps)
        a_ga = rng.integers(0, cfg.n_haps, cfg.n_gu)
        p_l = powers[rng.integers(0, cfg.n_power_levels, cfg.n_leo)]
        p_h = powers[rng.integers(0, cfg.n_power_levels, cfg.n_haps)]
        g1 = sample_gain(topo.pg_leo_haps_db, rng)
        g2 = sample_gain(topo.pg_haps_gu_db, rng)
        r1[k], _ = sum_rate(p_l, a_ha, g1, noise_w)
        r2[k], _ = sum_rate(p_h, a_ga, g2, noise_w)
    return float(r1.mean()), float(r2.mean())


def max_power_greedy(cfg, topo, rng, n_samples) -> tuple[float, float]:
    a_ha, a_ga = _greedy_assoc(topo)
    a_lp = np.full(cfg.n_leo, cfg.power_max_dbm)
    a_hp = np.full(cfg.n_haps, cfg.power_max_dbm)
    return _mean_rate(cfg, topo, a_ha, a_lp, a_ga, a_hp, rng, n_samples)


def best_uniform_greedy(cfg, topo, rng, n_samples) -> tuple[float, float]:
    """Greedy association; sweep the common power level and keep the best per tier."""
    a_ha, a_ga = _greedy_assoc(topo)
    best_r1 = -np.inf
    best_r2 = -np.inf
    for p in cfg.power_levels_dbm:
        # Independent RNG copies so each candidate power sees the same fading draws.
        r1, _ = _mean_rate(
            cfg,
            topo,
            a_ha,
            np.full(cfg.n_leo, p),
            a_ga,
            np.full(cfg.n_haps, p),
            np.random.default_rng(rng.integers(2**31)),
            n_samples,
        )
        # r1 and r2 are coupled in that call; recompute r2 with its own sweep below.
        best_r1 = max(best_r1, r1)
    for p in cfg.power_levels_dbm:
        _, r2 = _mean_rate(
            cfg,
            topo,
            a_ha,
            np.full(cfg.n_leo, p),
            a_ga,
            np.full(cfg.n_haps, p),
            np.random.default_rng(rng.integers(2**31)),
            n_samples,
        )
        best_r2 = max(best_r2, r2)
    return best_r1, best_r2


def run_all(cfg, topo, rng, n_samples) -> dict[str, tuple[float, float]]:
    return {
        "Random": random_policy(cfg, topo, rng, n_samples),
        "MaxPower+Greedy": max_power_greedy(cfg, topo, rng, n_samples),
        "BestUniform+Greedy": best_uniform_greedy(cfg, topo, rng, n_samples),
    }
