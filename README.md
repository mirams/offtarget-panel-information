# offtarget-panel-information

An illustrative simulation study: **what can a 50-target safety pharmacology panel tell you
about the ~950 physiological off-targets you did not screen?**

The short answer is that under independent target affinities it can tell you nothing at all,
by construction, so any proteome-wide reading of a clean panel is implicitly a claim about
shared compound-level structure. This repo quantifies how much of that structure is needed.

**Read [`report.md`](report.md) first.**

## Contents

| Path | What it is |
|---|---|
| `report.md` | The write-up: model, results, figures, modelling recommendations, limitations |
| `offtarget_panel/model.py` | Generative model, threshold calibration, quadrature results, Monte Carlo sampler |
| `run_all.py` | Every experiment; writes `results/*.csv` and `figures/*.png` |
| `tests/test_model.py` | Consistency checks, including quadrature against Monte Carlo |
| `results/` | Numerical output as CSV |
| `figures/` | Generated figures |

## Model in one line

Latent affinity of compound *c* at target *t* is a shared compound-level promiscuity term, an
optional target-family term, and a heavy-tailed idiosyncratic term; a liability is a threshold
exceedance. The shared variance share `rho` is the parameter of interest.

## Headline result

The **residual risk ratio** is the per-target liability rate among unscreened targets given a
clean panel, divided by the unconditional rate. It equals 1 exactly when targets are
independent, for any panel size.

| `rho` | Residual risk ratio after a clean 50-target panel |
|---|---|
| 0.0 | 1.00 |
| 0.3 | 0.88 |
| 0.5 | 0.70 |
| 0.7 | 0.42 |
| 0.85 | 0.19 |

Roughly half of all affinity variance must be shared before a clean 50-target panel cuts
expected liabilities elsewhere in the proteome by even 30%.

`rho` is estimable from panel data already held: shared structure over-disperses per-compound
hit counts relative to the binomial expectation. See section 9 of the report.

## Running it

```bash
pip install -r requirements.txt
python run_all.py
python -m pytest tests -q
```

All randomness is seeded; the report's numbers reproduce exactly.
