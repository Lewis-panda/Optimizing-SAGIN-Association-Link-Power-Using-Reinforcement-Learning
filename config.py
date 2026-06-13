"""Central configuration for the SAGIN RL experiment.

All tunable parameters live here so that experiments are reproducible and
self-documenting. Nothing else in the codebase should hard-code a magic number.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class Config:
    # ------------------------------------------------------------------ geometry
    coverage_km: float = 100.0  # square service area side length
    n_leo: int = 5  # LEO satellites  (tier-1 transmitters)
    n_haps: int = 10  # HAPS platforms  (tier-1 RX / tier-2 TX)
    n_gu: int = 15  # ground users    (tier-2 receivers)
    height_leo_km: float = 300.0
    height_haps_km: float = 20.0
    height_gu_km: float = 0.0

    # -------------------------------------------------------------------- radio
    freq_hz: float = 2.0e9  # S-band, typical for NTN access
    bandwidth_hz: float = 10e6
    noise_figure_db: float = 7.0
    # Combined Tx+Rx antenna gain folded into the link budget per tier (dBi).
    # Satellite links rely on highly directional antennas, hence the larger value.
    antenna_gain_leo_haps_db: float = 40.0
    antenna_gain_haps_gu_db: float = 20.0
    # Discrete transmit-power codebook (dBm). 10 dBm = 10 mW, 40 dBm = 10 W.
    power_min_dbm: float = 10.0
    power_max_dbm: float = 40.0
    n_power_levels: int = 10

    # -------------------------------------------------------------- q-learning
    episodes: int = 20000
    alpha: float = 0.1  # learning rate
    gamma: float = 0.9  # discount factor
    eps_start: float = 1.0
    eps_end: float = 0.01
    eps_decay_frac: float = 0.6  # fraction of episodes over which eps decays

    # -------------------------------------------------------------- experiment
    n_topologies: int = 10  # Monte-Carlo over independent random topologies
    eval_fading_samples: int = 200  # fading realisations used to estimate ergodic rate
    seed: int = 0

    @property
    def power_levels_dbm(self) -> np.ndarray:
        return np.linspace(self.power_min_dbm, self.power_max_dbm, self.n_power_levels)

    def quick(self) -> "Config":
        """Return a fast, low-fidelity variant for smoke tests / iteration."""
        from dataclasses import replace

        return replace(
            self,
            episodes=3000,
            n_topologies=3,
            eval_fading_samples=50,
        )
