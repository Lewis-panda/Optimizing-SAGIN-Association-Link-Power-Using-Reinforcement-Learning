"""The SAGIN sum-rate environment.

Two independent downlink tiers share the same machinery:

    tier 1:  LEO  -> HAPS   (backhaul)
    tier 2:  HAPS -> GU     (access)

Each receiver is served by exactly one transmitter (the *association*). Every
active transmitter reuses the whole band, so a receiver suffers co-channel
interference from *all other* active transmitters. This is what makes power
control a real problem: turning your power up helps you but hurts everyone else,
so "everybody at max power" is generally **not** optimal -- the bug that made the
original reward monotone in power is fixed here.
"""

from __future__ import annotations

import numpy as np

from channel import dbm_to_w


def sum_rate(
    power_dbm: np.ndarray,
    assoc: np.ndarray,
    gain_lin: np.ndarray,
    noise_w: float,
    share_load: bool = True,
) -> tuple[float, np.ndarray]:
    """Spectral efficiency of one tier under universal frequency reuse.

    Parameters
    ----------
    power_dbm : [n_tx]   transmit power chosen by each transmitter
    assoc     : [n_rx]   index of the transmitter serving each receiver
    gain_lin  : [n_tx, n_rx]   instantaneous linear channel gains
    noise_w   : scalar thermal noise power (W)

    share_load : bool
        If True (default) each transmitter splits its band equally among the
        receivers it serves (TDMA/OFDMA within a cell), so a receiver's share is
        ``(1/load) * log2(1+SINR)``. Without this, the optimiser could dump every
        receiver onto a single interference-free transmitter -- physically bogus.

    Returns
    -------
    (total_rate, per_rx_rate) in bits/s/Hz, with

        SINR_j = P[a_j] g[a_j, j] / ( sum_{m != a_j, active} P[m] g[m, j] + N0 ).
    """
    n_tx, n_rx = gain_lin.shape
    p_w = dbm_to_w(power_dbm)  # [n_tx]

    # A transmitter only emits if it serves at least one receiver.
    counts = np.bincount(assoc, minlength=n_tx)  # receivers served per tx
    active = counts > 0

    rx_power = p_w[:, None] * gain_lin  # power from every tx at every rx
    total_at_rx = (rx_power * active[:, None]).sum(axis=0)  # [n_rx]
    desired = rx_power[assoc, np.arange(n_rx)]  # [n_rx]
    interference = total_at_rx - desired  # other active tx only

    sinr = desired / (interference + noise_w)
    rates = np.log2(1.0 + sinr)
    if share_load:
        rates = rates / counts[assoc]  # equal intra-cell resource sharing
    return float(rates.sum()), rates
