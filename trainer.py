"""Training and greedy-policy evaluation for one topology.

Four agent groups act jointly every episode:

    tier 1 reward (LEO->HAPS sum-rate)  drives  HAPS-association + LEO-power
    tier 2 reward (HAPS->GU  sum-rate)  drives  GU-association   + HAPS-power

Small-scale fading is redrawn every episode, so the agents learn *ergodic*
(expected) rates rather than overfitting one frozen channel realisation.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from agents import AgentGroup, epsilon_at
from channel import Topology, noise_power_w, sample_gain
from env import sum_rate


@dataclass
class Policy:
    """The four trained Q-tables plus the states they converged to."""

    haps_assoc: AgentGroup
    leo_power: AgentGroup
    gu_assoc: AgentGroup
    haps_power: AgentGroup
    state: tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]


def train_one_topology(
    cfg, topo: Topology, rng: np.random.Generator
) -> tuple[np.ndarray, Policy]:
    noise_w = noise_power_w(cfg.bandwidth_hz, cfg.noise_figure_db)
    powers = cfg.power_levels_dbm

    haps_assoc = AgentGroup(cfg.n_haps, cfg.n_leo, cfg.n_leo)
    leo_power = AgentGroup(cfg.n_leo, cfg.n_power_levels, cfg.n_power_levels)
    gu_assoc = AgentGroup(cfg.n_gu, cfg.n_haps, cfg.n_haps)
    haps_power = AgentGroup(cfg.n_haps, cfg.n_power_levels, cfg.n_power_levels)

    # Random initial states (== initial actions).
    s_ha = rng.integers(0, cfg.n_leo, cfg.n_haps)
    s_lp = rng.integers(0, cfg.n_power_levels, cfg.n_leo)
    s_ga = rng.integers(0, cfg.n_haps, cfg.n_gu)
    s_hp = rng.integers(0, cfg.n_power_levels, cfg.n_haps)

    curve = np.empty(cfg.episodes)
    for ep in range(cfg.episodes):
        eps = epsilon_at(ep, cfg)

        a_ha = haps_assoc.select(s_ha, eps, rng)
        a_lp = leo_power.select(s_lp, eps, rng)
        a_ga = gu_assoc.select(s_ga, eps, rng)
        a_hp = haps_power.select(s_hp, eps, rng)

        g1 = sample_gain(topo.pg_leo_haps_db, rng)
        g2 = sample_gain(topo.pg_haps_gu_db, rng)

        r1, rates1 = sum_rate(
            powers[a_lp], a_ha, g1, noise_w
        )  # tier 1: HAPS pick a LEO
        r2, rates2 = sum_rate(powers[a_hp], a_ga, g2, noise_w)  # tier 2: GU pick a HAPS

        # Credit assignment:
        #   association agents learn from their OWN link rate (local reward) -- a
        #     receiver can directly tell whether its chosen transmitter is good;
        #   power agents learn from the GLOBAL tier rate so they still feel the
        #     interference externality that raising power inflicts on others.
        # next state == chosen action
        haps_assoc.update(s_ha, a_ha, rates1, a_ha, cfg.alpha, cfg.gamma)
        leo_power.update(s_lp, a_lp, r1, a_lp, cfg.alpha, cfg.gamma)
        gu_assoc.update(s_ga, a_ga, rates2, a_ga, cfg.alpha, cfg.gamma)
        haps_power.update(s_hp, a_hp, r2, a_hp, cfg.alpha, cfg.gamma)

        s_ha, s_lp, s_ga, s_hp = a_ha, a_lp, a_ga, a_hp
        curve[ep] = r1 + r2

    policy = Policy(
        haps_assoc, leo_power, gu_assoc, haps_power, (s_ha, s_lp, s_ga, s_hp)
    )
    return curve, policy


def eval_policy(
    cfg, topo: Topology, policy: Policy, rng: np.random.Generator, n_samples: int
) -> tuple[float, float]:
    """Mean tier-1 / tier-2 rate of the converged greedy policy over fresh fading."""
    noise_w = noise_power_w(cfg.bandwidth_hz, cfg.noise_figure_db)
    powers = cfg.power_levels_dbm
    s_ha, s_lp, s_ga, s_hp = policy.state

    a_ha = policy.haps_assoc.greedy(s_ha)
    a_lp = policy.leo_power.greedy(s_lp)
    a_ga = policy.gu_assoc.greedy(s_ga)
    a_hp = policy.haps_power.greedy(s_hp)

    r1 = np.empty(n_samples)
    r2 = np.empty(n_samples)
    for k in range(n_samples):
        g1 = sample_gain(topo.pg_leo_haps_db, rng)
        g2 = sample_gain(topo.pg_haps_gu_db, rng)
        r1[k], _ = sum_rate(powers[a_lp], a_ha, g1, noise_w)
        r2[k], _ = sum_rate(powers[a_hp], a_ga, g2, noise_w)
    return float(r1.mean()), float(r2.mean())
