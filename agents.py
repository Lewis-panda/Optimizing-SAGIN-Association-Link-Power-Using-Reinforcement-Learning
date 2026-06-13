"""Tabular Q-learning agents (faithful to the original idea, vectorised & fixed).

Every decision-maker is an independent Q-learner. We group identical agents so a
whole group updates with one set of NumPy ops:

    * association agents : state = current association, action = next association
    * power agents       : state = current power index, action = next power index

Using "current action" as the state turns each agent into a small, well-defined
MDP (the original used a quantised copy of the reward as the state, which made
state == reward -- degenerate). Agents share a global per-tier reward, i.e. this
is *independent Q-learning*, the standard multi-agent baseline.
"""

from __future__ import annotations

import numpy as np


class AgentGroup:
    """A batch of ``n_agents`` independent tabular Q-learners with equal shapes."""

    def __init__(self, n_agents: int, n_states: int, n_actions: int):
        self.n_agents = n_agents
        self.n_states = n_states
        self.n_actions = n_actions
        self.q = np.zeros((n_agents, n_states, n_actions))

    def select(
        self, states: np.ndarray, eps: float, rng: np.random.Generator
    ) -> np.ndarray:
        """Epsilon-greedy action for every agent given its current state."""
        idx = np.arange(self.n_agents)
        greedy = self.q[idx, states].argmax(axis=1)
        rand = rng.integers(0, self.n_actions, size=self.n_agents)
        explore = rng.random(self.n_agents) < eps
        return np.where(explore, rand, greedy)

    def greedy(self, states: np.ndarray) -> np.ndarray:
        """Pure-exploitation action (used at evaluation time)."""
        idx = np.arange(self.n_agents)
        return self.q[idx, states].argmax(axis=1)

    def update(
        self,
        states: np.ndarray,
        actions: np.ndarray,
        reward: float,
        next_states: np.ndarray,
        alpha: float,
        gamma: float,
    ) -> None:
        """Q-learning (Bellman) update.

        ``reward`` may be a scalar (shared/global) or a per-agent vector of shape
        ``[n_agents]`` for local credit assignment; NumPy broadcasting handles both.
        """
        idx = np.arange(self.n_agents)
        best_next = self.q[idx, next_states].max(axis=1)
        td_target = reward + gamma * best_next
        self.q[idx, states, actions] += alpha * (
            td_target - self.q[idx, states, actions]
        )


def epsilon_at(ep: int, cfg) -> float:
    """Linearly decay epsilon from eps_start to eps_end over a fraction of training."""
    span = max(1, int(cfg.eps_decay_frac * cfg.episodes))
    frac = min(1.0, ep / span)
    return cfg.eps_start + (cfg.eps_end - cfg.eps_start) * frac
