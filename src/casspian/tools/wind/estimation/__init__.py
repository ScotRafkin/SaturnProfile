"""The wind estimation tool, `casspian-wind-fit` (SPEC_19, design note v0.5).

Three parts, of which only the first two grow as datasets are added:

* **Observation types** (`radio.py`, `gradient.py`): each a reader of its standard data, an operator
  with its linearization, and an error covariance.
* **The prior** (`prior.py`): a mean and a covariance.
* **The solver** (`solver.py`): weighted least squares, the posterior, the consistency check and
  leave-one-out. It sees only matrices and group labels.

`state.py` is the basis the shear is carried on, `driver.py` the command and the relinearization.
The package is not `tools/wind/fit.py` because that module is SPEC_01's point-based latitude fit.
"""
