"""
Generative model for compound-by-off-target affinities.

Latent affinity of compound c at target t, on a standardised pIC50-like scale:

    z_ct = mu_t + sqrt(rho_g) * phi_c
                + sqrt(rho_f) * psi_{c, f(t)}
                + sqrt(1 - rho_g - rho_f) * eps_ct

    phi_c   ~ N(0, 1)          compound-level promiscuity (one global axis)
    psi_cf  ~ N(0, 1)          compound-by-target-family promiscuity
    eps_ct  ~ standardised t_nu   idiosyncratic, heavy tailed

Var(z_ct) = 1 + Var(mu_t) by construction, so rho_g and rho_f are directly
interpretable as variance shares.  A compound is counted as having a
safety-relevant liability at target t ("a hit") when z_ct > theta.  The long
right tail of the affinity distribution comes from the Student-t idiosyncratic
term; most (compound, target) pairs sit far below theta, a few are extreme.

Two framings the user might reach for are the same object here:

  * "targets covary"                 -> Corr(z_ct, z_ct') = rho_g (+ rho_f
                                        within a family)
  * "some compounds are dirty across
     the board"                      -> phi_c shifts the whole affinity
                                        distribution for compound c

For a single factor these are algebraically identical; the interesting choice
is how many factors, not which framing.

With mu_t constant and rho_f = 0 the targets are exchangeable, so everything
we need is a one-dimensional integral over phi and we do not have to simulate.
The Monte Carlo path is kept for the family-structured case and to check the
quadrature.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy import optimize, stats

# Quadrature grid over the standard normal promiscuity axis.  The integrands
# below contain (1 - q)^T with T ~ 1000, which is sharply varying in phi, so we
# use a fine trapezoidal grid rather than a low-order Gauss rule.
_PHI_GRID = np.linspace(-9.0, 9.0, 40001)
_PHI_WEIGHT = stats.norm.pdf(_PHI_GRID)
_PHI_WEIGHT = _PHI_WEIGHT / np.trapezoid(_PHI_WEIGHT, _PHI_GRID)


def _standardised_t_sf(x: np.ndarray, nu: float) -> np.ndarray:
    """Survival function of a Student-t rescaled to unit variance."""
    return stats.t.sf(x * np.sqrt(nu / (nu - 2.0)), df=nu)


def _standardised_t_rvs(size, nu: float, rng: np.random.Generator) -> np.ndarray:
    """Draws from a Student-t rescaled to unit variance."""
    return rng.standard_t(nu, size=size) * np.sqrt((nu - 2.0) / nu)


@dataclass
class PanelModel:
    """One- or two-factor model for a physiological off-target panel."""

    n_targets: int = 1000
    rho_global: float = 0.5
    rho_family: float = 0.0
    n_families: int = 20
    nu: float = 5.0
    hit_rate: float = 0.003  # marginal P(hit) per compound per target
    theta: float = field(init=False)

    def __post_init__(self) -> None:
        if self.nu <= 2.0:
            raise ValueError("nu must exceed 2 for the idiosyncratic term to have finite variance")
        if self.rho_global + self.rho_family >= 1.0:
            raise ValueError("shared variance shares must leave a positive idiosyncratic share")
        if self.n_targets % self.n_families:
            raise ValueError("n_targets must divide evenly into n_families")
        self.theta = self._calibrate_threshold()

    # ------------------------------------------------------------------
    # threshold calibration
    # ------------------------------------------------------------------

    def _shared_sd(self) -> float:
        """SD of the shared (phi + psi) part, which is Gaussian."""
        return np.sqrt(self.rho_global + self.rho_family)

    def _idio_sd(self) -> float:
        return np.sqrt(1.0 - self.rho_global - self.rho_family)

    def _marginal_hit_rate(self, theta: float) -> float:
        """P(z > theta) averaged over compounds, integrating out all shared factors."""
        s = self._shared_sd()
        if s == 0.0:
            return float(_standardised_t_sf(np.array(theta / self._idio_sd()), self.nu))
        q = _standardised_t_sf((theta - s * _PHI_GRID) / self._idio_sd(), self.nu)
        return float(np.trapezoid(q * _PHI_WEIGHT, _PHI_GRID))

    def _calibrate_threshold(self) -> float:
        """Choose theta so the marginal per-target hit rate matches hit_rate."""
        f = lambda th: self._marginal_hit_rate(th) - self.hit_rate
        return float(optimize.brentq(f, 0.0, 40.0, xtol=1e-12, rtol=1e-14))

    # ------------------------------------------------------------------
    # conditional hit probability
    # ------------------------------------------------------------------

    def hit_prob(self, shared: np.ndarray) -> np.ndarray:
        """P(hit at a given target | shared latent value), i.e. q(phi)."""
        shared = np.asarray(shared, dtype=float)
        return _standardised_t_sf(
            (self.theta - self._shared_sd() * shared) / self._idio_sd(), self.nu
        )

    # ------------------------------------------------------------------
    # exact results for the exchangeable (single global factor) case
    # ------------------------------------------------------------------

    def _require_exchangeable(self) -> None:
        if self.rho_family > 0.0:
            raise ValueError(
                "closed-form results assume exchangeable targets; use the Monte Carlo "
                "routines when rho_family > 0"
            )

    def _expect(self, values: np.ndarray) -> float:
        return float(np.trapezoid(values * _PHI_WEIGHT, _PHI_GRID))

    def p_clean(self, n: int) -> float:
        """P(a compound has no hits across n exchangeable targets)."""
        self._require_exchangeable()
        q = self.hit_prob(_PHI_GRID)
        return self._expect(np.power(1.0 - q, n))

    def npv(self, panel_size: int) -> float:
        """P(no hits among the unscreened targets | panel of this size came back clean)."""
        self._require_exchangeable()
        q = self.hit_prob(_PHI_GRID)
        clean_panel = np.power(1.0 - q, panel_size)
        clean_all = np.power(1.0 - q, self.n_targets)
        return self._expect(clean_all) / self._expect(clean_panel)

    def expected_missed_hits(self, panel_size: int) -> float:
        """E[number of hits among unscreened targets | clean panel]."""
        self._require_exchangeable()
        q = self.hit_prob(_PHI_GRID)
        clean_panel = np.power(1.0 - q, panel_size)
        n_rest = self.n_targets - panel_size
        return self._expect(clean_panel * n_rest * q) / self._expect(clean_panel)

    def residual_risk_ratio(self, panel_size: int) -> float:
        """
        Per-target liability rate among unscreened targets after a clean panel,
        divided by the unconditional per-target rate.

        This is the metric to read first.  It is exactly 1.0 when targets are
        independent, for every panel size, because a clean panel then says
        nothing at all about the targets you did not screen.  A value of 0.4
        means the clean panel has cut the expected liability rate elsewhere in
        the proteome by 60%.  Unlike a negative predictive value it does not
        move around simply because the base rate of clean compounds moves.
        """
        self._require_exchangeable()
        q = self.hit_prob(_PHI_GRID)
        w = np.power(1.0 - q, panel_size)
        return self._expect(w * q) / self._expect(w) / self.hit_rate

    def hit_count_moments(self, n: int | None = None) -> tuple[float, float]:
        """
        Mean and variance of the number of hits per compound across n targets.

        Defaults to the full proteome, but the interesting case is n = the size
        of a panel you actually ran, because that is what you can measure.
        """
        self._require_exchangeable()
        q = self.hit_prob(_PHI_GRID)
        T = self.n_targets if n is None else n
        mean = T * self._expect(q)
        # law of total variance for K | phi ~ Binomial(T, q(phi))
        e_var = T * self._expect(q * (1.0 - q))
        var_e = T**2 * (self._expect(q**2) - self._expect(q) ** 2)
        return mean, e_var + var_e

    def variance_to_mean(self, n: int | None = None) -> float:
        """
        Over-dispersion of hit counts relative to the independent (binomial) case.

        Equals 1 - q ~= 1 when targets are independent.  Anything appreciably
        above 1 in real panel data is evidence for shared structure, and this is
        the handle by which rho can be estimated from data you already have.
        """
        mean, var = self.hit_count_moments(n)
        return var / mean

    # ------------------------------------------------------------------
    # Monte Carlo
    # ------------------------------------------------------------------

    def family_index(self) -> np.ndarray:
        per_family = self.n_targets // self.n_families
        return np.repeat(np.arange(self.n_families), per_family)

    def simulate_hits(
        self, n_compounds: int, rng: np.random.Generator, chunk: int = 5000
    ) -> np.ndarray:
        """
        Boolean hit matrix of shape (n_compounds, n_targets).

        Generated in chunks so the full latent matrix is never held at once.
        """
        fam = self.family_index()
        out = np.empty((n_compounds, self.n_targets), dtype=bool)
        a = np.sqrt(self.rho_global)
        b = np.sqrt(self.rho_family)
        s = self._idio_sd()
        for start in range(0, n_compounds, chunk):
            stop = min(start + chunk, n_compounds)
            m = stop - start
            z = a * rng.standard_normal((m, 1))
            if b > 0.0:
                psi = rng.standard_normal((m, self.n_families))
                z = z + b * psi[:, fam]
            z = z + s * _standardised_t_rvs((m, self.n_targets), self.nu, rng)
            out[start:stop] = z > self.theta
        return out


PIC50_CENTRE = 4.4
PIC50_SPREAD = 0.9
PIC50_FLOOR = 4.0


def to_pic50(z: np.ndarray) -> np.ndarray:
    """
    Map the standardised latent scale onto something that looks like assay output.

    Nothing in the analysis depends on this: the model is scale-free and every
    result is driven by where the liability threshold sits in the tail.  The
    transform exists so that figures show what a profiling dataset actually
    looks like, with most (compound, target) pairs left-censored at the assay
    detection limit rather than sitting at a symmetric negative value.
    """
    return np.maximum(PIC50_CENTRE + PIC50_SPREAD * np.asarray(z), PIC50_FLOOR)


def panel_npv_from_hits(hits: np.ndarray, panel: np.ndarray) -> dict[str, float]:
    """
    Empirical negative predictive value of a panel.

    hits  : (n_compounds, n_targets) boolean
    panel : integer indices of the screened targets

    Returns the base rate (unconditional probability a compound is clean
    everywhere), the NPV, and the expected number of missed liabilities.
    """
    n_targets = hits.shape[1]
    mask = np.zeros(n_targets, dtype=bool)
    mask[panel] = True
    panel_hits = hits[:, mask].sum(axis=1)
    rest_hits = hits[:, ~mask].sum(axis=1)

    n_rest = int((~mask).sum())
    overall_rate = float(hits.mean())

    clean_panel = panel_hits == 0
    n_clean_panel = int(clean_panel.sum())
    if n_clean_panel == 0:
        return {
            "base_rate": float((hits.sum(axis=1) == 0).mean()),
            "p_clean_panel": 0.0,
            "npv": np.nan,
            "expected_missed": np.nan,
            "residual_risk_ratio": np.nan,
            "n_clean_panel": 0,
        }
    expected_missed = float(rest_hits[clean_panel].mean())
    return {
        "base_rate": float((hits.sum(axis=1) == 0).mean()),
        "p_clean_panel": float(clean_panel.mean()),
        "npv": float((rest_hits[clean_panel] == 0).mean()),
        "expected_missed": expected_missed,
        "residual_risk_ratio": expected_missed / n_rest / overall_rate,
        "n_clean_panel": n_clean_panel,
    }


def spread_panel(n_targets: int, n_families: int, size: int) -> np.ndarray:
    """
    Panel spread as evenly as possible across families.

    The remainder is distributed one target at a time so the panel really does
    contain `size` targets; rounding the allocation down would silently shrink
    the panel and make design comparisons unfair.
    """
    per_family = n_targets // n_families
    take = np.full(n_families, size // n_families, dtype=int)
    take[: size % n_families] += 1
    if take.max() > per_family:
        raise ValueError("panel larger than the available targets per family")
    idx = []
    for f, k in enumerate(take):
        idx.extend(range(f * per_family, f * per_family + int(k)))
    assert len(idx) == size
    return np.array(idx)


def concentrated_panel(n_targets: int, n_families: int, size: int) -> np.ndarray:
    """Panel drawn from as few families as possible."""
    per_family = n_targets // n_families
    idx = []
    f = 0
    while len(idx) < size:
        idx.extend(range(f * per_family, (f + 1) * per_family))
        f += 1
    return np.array(idx[:size])
