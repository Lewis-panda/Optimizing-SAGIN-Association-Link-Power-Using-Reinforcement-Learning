"""Physical-layer model: geometry, path loss, fading and the link budget.

This module is deliberately small and pure (no RL, no state) so the radio model
can be unit-tested in isolation. The single most important fix relative to the
original MATLAB code lives here: a channel *gain* is a small linear quantity
(``10**(-PL/10)``), not a path-*loss* in dB. Larger distance -> smaller gain.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

SPEED_OF_LIGHT = 3.0e8  # m/s
BOLTZMANN = 1.380649e-23  # J/K
TEMPERATURE_K = 290.0  # reference noise temperature


def fspl_db(distance_m: np.ndarray | float, freq_hz: float) -> np.ndarray | float:
    """Free-space path loss in dB.

    FSPL(dB) = 20*log10(d) + 20*log10(f) + 20*log10(4*pi/c),  d in metres.
    The third term equals -147.55 dB; it is *included*, not added on top of an
    already-complete expression (the original code did both and cancelled it).
    """
    return (
        20.0 * np.log10(distance_m)
        + 20.0 * np.log10(freq_hz)
        + 20.0 * np.log10(4.0 * np.pi / SPEED_OF_LIGHT)
    )


def noise_power_w(bandwidth_hz: float, noise_figure_db: float) -> float:
    """Thermal noise power in watts: kT B times the receiver noise figure."""
    n0 = BOLTZMANN * TEMPERATURE_K * bandwidth_hz
    return n0 * 10.0 ** (noise_figure_db / 10.0)


def dbm_to_w(power_dbm: np.ndarray | float) -> np.ndarray | float:
    return 10.0 ** ((power_dbm - 30.0) / 10.0)


def random_positions(
    rng: np.random.Generator, n: int, area_km: float, height_km: float
) -> np.ndarray:
    """Return a 3 x n array of (x, y, z) positions in km, z fixed at ``height_km``."""
    xy = rng.uniform(0.0, area_km, size=(2, n))
    z = np.full((1, n), height_km)
    return np.vstack([xy, z])


def path_gain_db(
    tx_pos: np.ndarray, rx_pos: np.ndarray, freq_hz: float, antenna_gain_db: float
) -> np.ndarray:
    """Large-scale channel gain in dB for every (tx, rx) pair.

    Returns an ``[n_tx, n_rx]`` matrix of ``antenna_gain - FSPL`` (a negative dB
    value). Distances are converted km -> m here; the rest of the code never has
    to worry about units again.
    """
    diff = tx_pos[:, :, None] - rx_pos[:, None, :]  # 3 x n_tx x n_rx
    dist_m = np.sqrt((diff**2).sum(axis=0)) * 1.0e3  # n_tx x n_rx, metres
    return -fspl_db(dist_m, freq_hz) + antenna_gain_db


@dataclass
class Topology:
    """A frozen geometry: large-scale gains are fixed; only fading varies per step."""

    pg_leo_haps_db: np.ndarray  # [n_leo, n_haps]
    pg_haps_gu_db: np.ndarray  # [n_haps, n_gu]


def make_topology(cfg, rng: np.random.Generator) -> Topology:
    leo = random_positions(rng, cfg.n_leo, cfg.coverage_km, cfg.height_leo_km)
    haps = random_positions(rng, cfg.n_haps, cfg.coverage_km, cfg.height_haps_km)
    gu = random_positions(rng, cfg.n_gu, cfg.coverage_km, cfg.height_gu_km)
    return Topology(
        pg_leo_haps_db=path_gain_db(
            leo, haps, cfg.freq_hz, cfg.antenna_gain_leo_haps_db
        ),
        pg_haps_gu_db=path_gain_db(haps, gu, cfg.freq_hz, cfg.antenna_gain_haps_gu_db),
    )


def sample_gain(pg_db: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """Instantaneous linear gain = large-scale gain x Rayleigh fading.

    ``|h|^2 ~ Exp(1)`` models Rayleigh small-scale fading power with unit mean.
    """
    fading = rng.exponential(1.0, size=pg_db.shape)
    return 10.0 ** (pg_db / 10.0) * fading
