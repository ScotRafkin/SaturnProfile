"""SPEC_19 v0.5 section 5 item 6: the absolute profiles, plotting only, from the existing runs. For each target,
the observed temperature against pressure with its sigma_T band, and the delivered profiles of hypotheses 1, 2
and 3 and the untransferred anchor, 1.3 bar to 0.2 mbar, log pressure.

Run from the repository root, after the runs are built. Writes `reports/figures/step19_1_absolute_profiles.png`.
"""

import sys
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

sys.path.insert(0, str(Path(__file__).parent))
import accept_step19_1 as accept  # noqa: E402
from casspian.tools.lindal import radio_profiles as rp  # noqa: E402
from casspian.tools.wind.estimation import driver  # noqa: E402

OUTPUT = Path("reports/figures/step19_1_absolute_profiles.png")
#: The range of the figure, mbar.
BOTTOM_MBAR, TOP_MBAR = 1300.0, 0.2
#: Target, then its planetocentric latitude for the panel's title.
PLANETOCENTRIC = {"31.2 S": "26.1 S", "74.2 S": "71.1 S", "3.0 S": "2.4 S"}


def main():
    warnings.simplefilter("ignore")
    problem = driver.read_control(accept.CONTROL)
    fig, axes = plt.subplots(1, 3, figsize=(15, 7.5), sharey=True, sharex=True)
    for ax, target in zip(axes, problem.radio.targets):
        table = rp.read_csv(target.observed)
        p_obs, T_obs = table["pressure_mbar"], table["temperature_K"]
        good = np.isfinite(T_obs)
        p_obs, T_obs = p_obs[good], T_obs[good]
        sigma = problem.radio.sigma(p_obs * 100.0)
        ax.fill_betweenx(p_obs, T_obs - sigma, T_obs + sigma, color="0.82", label="sigma_T")
        ax.plot(T_obs, p_obs, color="k", linewidth=1.6, label="observed")
        anchor = None
        for name, label, color in zip(accept.RUNS[target.name], ("hypothesis 1", "hypothesis 2", "hypothesis 3"), ("C2", "C4", "C3")):
            product = driver.tree(Path("forward") / name / "output" / f"{name}_profile.nc")
            root = product.to_dataset(inherit=False)
            ax.plot(root["temperature_K"].values, root["pressure_Pa"].values / 100.0, color=color, linewidth=1.2, label=label)
            anchor = anchor or product["anchors/lindal/inputs/thermo"].to_dataset(inherit=False)
        ax.plot(anchor["temperature_K"].values, anchor["pressure_Pa"].values / 100.0, color="C0", linestyle="--", linewidth=1.2,
                label="anchor, 36.5 N (30.8 N), not transferred")
        ax.set_yscale("log")
        ax.set_title(f"{target.name} ({PLANETOCENTRIC[target.name]} planetocentric)")
        ax.set_xlabel("temperature (K)")
        ax.grid(True, color="0.92", linewidth=0.6)
    axes[0].set_ylim(BOTTOM_MBAR, TOP_MBAR)
    axes[0].set_ylabel("pressure (mbar)")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=len(labels), fontsize=9, frameon=False)
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    OUTPUT.parent.mkdir(exist_ok=True)
    fig.savefig(OUTPUT, dpi=110)
    print(f"wrote {OUTPUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
