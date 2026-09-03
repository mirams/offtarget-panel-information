"""Consistency checks. Run with: python -m pytest tests -q"""

import numpy as np
import pytest

from offtarget_panel.model import (
    PanelModel,
    concentrated_panel,
    panel_npv_from_hits,
    spread_panel,
)


def test_threshold_calibration_hits_target_rate():
    for rho in [0.0, 0.4, 0.8]:
        m = PanelModel(rho_global=rho, hit_rate=0.003)
        assert m._marginal_hit_rate(m.theta) == pytest.approx(0.003, rel=1e-6)


def test_independence_gives_no_information():
    """The premise of the study: with rho = 0 a clean panel says nothing."""
    m = PanelModel(rho_global=0.0)
    for n in [5, 50, 200, 900]:
        assert m.residual_risk_ratio(n) == pytest.approx(1.0, abs=1e-9)
    # and the NPV is just the base rate over the shrunken remainder
    assert m.npv(50) == pytest.approx(m.p_clean(m.n_targets - 50), rel=1e-9)


def test_more_covariance_is_more_informative():
    ratios = [PanelModel(rho_global=r).residual_risk_ratio(50) for r in [0.0, 0.3, 0.6, 0.9]]
    assert all(a > b for a, b in zip(ratios, ratios[1:]))


def test_quadrature_matches_monte_carlo():
    rng = np.random.default_rng(0)
    for rho in [0.0, 0.5]:
        m = PanelModel(rho_global=rho, hit_rate=0.003)
        hits = m.simulate_hits(60_000, rng)
        mc = panel_npv_from_hits(hits, np.arange(50))
        assert mc["npv"] == pytest.approx(m.npv(50), abs=0.01)
        assert mc["residual_risk_ratio"] == pytest.approx(m.residual_risk_ratio(50), abs=0.03)
        assert mc["base_rate"] == pytest.approx(m.p_clean(m.n_targets), abs=0.01)


def test_panel_designs_have_the_requested_size():
    for size in [20, 44, 50, 97]:
        assert len(spread_panel(1000, 20, size)) == size
        assert len(np.unique(spread_panel(1000, 20, size))) == size
        assert len(concentrated_panel(1000, 20, size)) == size


def test_spread_panel_covers_every_family_when_it_can():
    panel = spread_panel(1000, 20, 50)
    assert len(np.unique(panel // 50)) == 20


def test_family_structure_requires_monte_carlo():
    m = PanelModel(rho_global=0.2, rho_family=0.4)
    with pytest.raises(ValueError):
        m.npv(50)
