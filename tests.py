"""Lightweight sanity tests for the radio model -- run: ``python tests.py``.

These guard the exact bugs that broke the original code, so a regression here
means the physics is wrong again.
"""

from __future__ import annotations

import numpy as np

from channel import dbm_to_w, fspl_db, path_gain_db
from env import sum_rate


def test_fspl_constant():
    # The standard constant 20*log10(4*pi/c) must be ~ -147.55 dB.
    const = 20.0 * np.log10(4.0 * np.pi / 3.0e8)
    assert abs(const - (-147.55)) < 0.1, const


def test_path_gain_decreases_with_distance():
    tx = np.array([[0.0], [0.0], [20.0]])  # one tx at 20 km altitude
    near = np.array([[0.0], [0.0], [0.0]])  # directly below
    far = np.array([[80.0], [80.0], [0.0]])  # far corner
    g_near = path_gain_db(tx, near, 2.0e9, 0.0)[0, 0]
    g_far = path_gain_db(tx, far, 2.0e9, 0.0)[0, 0]
    assert g_near > g_far, (g_near, g_far)  # closer => higher gain


def test_gain_is_small_linear():
    # A channel gain must be << 1 in linear scale (it is a loss, not a +147 dB number).
    tx = np.array([[0.0], [0.0], [300.0]])
    rx = np.array([[10.0], [10.0], [0.0]])
    g_lin = 10.0 ** (path_gain_db(tx, rx, 2.0e9, 40.0)[0, 0] / 10.0)
    assert 0.0 < g_lin < 1e-6, g_lin


def test_interference_tradeoff():
    # Two tx, two rx (rx0<-tx0, rx1<-tx1). Raising tx1 power must lower rx0 SINR.
    gain = np.array(
        [
            [1e-12, 5e-13],  # tx0 -> rx0, rx1
            [4e-13, 1e-12],
        ]
    )  # tx1 -> rx0, rx1
    assoc = np.array([0, 1])
    noise = 1e-13
    _, rates_lo = sum_rate(np.array([30.0, 10.0]), assoc, gain, noise)
    _, rates_hi = sum_rate(np.array([30.0, 40.0]), assoc, gain, noise)
    assert rates_hi[0] < rates_lo[0], (
        rates_lo[0],
        rates_hi[0],
    )  # victim hurt by interferer


def test_sum_rate_matches_hand_calc():
    gain = np.array([[1e-12, 5e-13], [4e-13, 1e-12]])
    assoc = np.array([0, 1])
    noise = 1e-13
    p = np.array([30.0, 30.0])
    total, rates = sum_rate(p, assoc, gain, noise)
    p_w = dbm_to_w(p)
    # rx0: desired tx0, interferer tx1
    sinr0 = (p_w[0] * gain[0, 0]) / (p_w[1] * gain[1, 0] + noise)
    sinr1 = (p_w[1] * gain[1, 1]) / (p_w[0] * gain[0, 1] + noise)
    expected = np.log2(1 + sinr0) + np.log2(1 + sinr1)
    assert abs(total - expected) < 1e-9, (total, expected)
    assert abs(rates[0] - np.log2(1 + sinr0)) < 1e-9


def main():
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for t in tests:
        t()
        print(f"PASS  {t.__name__}")
    print(f"\n{len(tests)} sanity tests passed.")


if __name__ == "__main__":
    main()
