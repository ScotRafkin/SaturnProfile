# REPORT 11, Step 2. The Lindal anchor under the Lindal-only wind

Coding agent, 8 October 2026. Specification: `docs/specs/SPEC_11_Lindal_Wind.md` v0.4 §1a. Working
tree on `main` at `7c06fb4` (Step 1 accepted and pushed), uncommitted. Acceptance:
`tests/step11_2/accept_step11_2.py`, **2 of 2** (both reported items); output in
`reports/step11_2/output.txt`. Diagnosis of the round trip:
`tests/step11_2/diagnose_round_trip.py`, output in `reports/step11_2/diagnose_round_trip.txt`. No
regression, and no existing suite rerun, as the author directed.

**In short.**
- **The run.** `lindal_iris_ii_110` converges at the 1e-8 tolerance in 5 passes, and its pressure
  identity is at most 0.074 K.
- **The round trip** misses the IRIS fit's change by -1.54, +0.77 and -1.13 K at 110, 290 and
  730 mbar, and each miss is explained (§3). The transfer reproduces the change its own wind implies
  to within 0.06 to 0.23 K. The larger parts are the comparison itself and the equatorial bridge:
  - the fit's slope, which the wind is built from, does not integrate to the fit's change in value
    (1.26 K at 110 mbar);
  - inside the equatorial band the shear is the bridge, not the data (up to 0.85 K at 730 mbar).

## 1. The changes

| File | Change |
|---|---|
| `tools/wind/shear.py` | The case `lindal_iris`, whose function calls `lindal_wind.construct` and returns its whole dataset (`Case.whole`). Its inputs `temperatures`, `composition`, `gravity` and `rotation` are path keys; `build` reads them into memory and records them in `input_hashes`, so `construct` still touches no file. The tool's refusals are unchanged. |
| `tools/run/run_inputs.py` | The sections run gravity, rotation, wind, composition, shear. The driver's refusal names the sections in that order. |
| `tools/lindal/lindal_wind.py` | `build`, `main` and the `[lindal_wind]` keys removed. `title` is now among the globals replaced, since the cloud wind's would be false of this wind. |
| `forward/lindal_iris_ii_110/` | Made by `casspian-new-run` from `lindal_transfer`. In the build file, `[shear]` is the `lindal_iris` case at `top_level_Pa = 11000`, `above_top = "relaxed"`. In the namelist, `[target]` is -26.4 deg and the description and header say so. Nothing else differs. |

**Suites edited, not rerun** (the author's direction; named here so none is edited silently):

- **`step05_3` check 1** expects the driver's refusal to name the sections. The expected string goes
  from `['gravity', 'rotation', 'wind', 'shear', 'composition']` to
  `['gravity', 'rotation', 'wind', 'composition', 'shear']`. Its regular expression edits the build
  file's text, whose layout is unchanged, so nothing else in the check moves.
- **`step11_1`** built its winds through `lindal_wind.build`, which is removed. It now writes a
  `[shear]` control with `case = "lindal_iris"` and calls `shear.build`. The winds are the same
  constructions, and the title and `created_by` change.

Both compile. Neither has been run since the edit.

## 2. Acceptance

The run: 1149 latitudes by 82 geopotentials. The outer loop's largest change in `ln p` falls from
9.0e-2 to 8.3e-9 in 5 passes (tolerance 1e-8).

1. **The round trip** (reported). The anchor is at 30.806 deg planetocentric, 36.30 planetographic.
   The target is at -26.40 planetocentric, -31.49 planetographic by the wind tool's rule (the
   specification states -31.35). The IRIS fit is taken at its defaults at the two planetographic
   latitudes.

   | Level | Delivered `T(target) - T(anchor)` | IRIS fit's change | Difference |
   |---|---|---|---|
   | 110 mbar | +3.72 K | +5.27 K (81.42 to 86.69) | -1.54 K |
   | 290 mbar | +1.29 K | +0.52 K (96.00 to 96.51) | +0.77 K |
   | 730 mbar | -5.35 K | -4.22 K (119.30 to 115.08) | -1.13 K |

   All three are over 0.5 K and are explained in §3.
2. **For the author** (`target_profile.png`):
   - **The delivered profile** at -26.4 deg with the anchor's. The IRIS fit's values at the two
     latitudes are drawn as dots; the IRIS and Lindal temperatures differ, and only changes are
     compared.
   - **Its lapse rate** with the dry adiabat, `g / c_p`. Here `g = dPhi/dz` along the target's
     vertical, and `c_p` is taken with frozen molecular rotation: 7/2 R per mole of H2, 5/2 R of He
     and 4 R of NH3. The adiabat is 0.73 K/km.
     - The delivered lapse rate is above it at 15 levels between 601 and 1258 mbar, reaching about
       1.0 K/km.
     - The anchor's own lapse rate there is 0.75 to 0.89 K/km, above the same adiabat (0.78 K/km
       at the anchor). So most of the excess is the adiabat's heat capacity, not the transfer. The
       transfer adds up to about 0.15 K/km.
     - The lapse rate is a difference between adjacent levels and is noisy at depth.
   - **The pressure identity in kelvins.** At fixed refractivity a residual `d` is `T d`. It is at
     most 0.074 K (residual -6.6e-4) at 660 mbar, with an RMS of 0.035 K. The largest values sit
     between 300 and 1000 mbar, where the shear changes.

## 3. The round trip explained

`diagnose_round_trip.py` reads the run's own wind and splits the path from the anchor to the target
into the equatorial band (`|phi_c| < 5 deg`) and the rest. On each piece it integrates two
gradients of `ln T` along the level:
- the IRIS fit's slope;
- the gradient the wind implies by the balance the tool inverted, on the tool's geometry.

Each is expressed in kelvins as `T_anchor x` the integral. The difference splits into four parts
that add to the reported one:

| Part | 110 mbar | 290 mbar | 730 mbar |
|---|---|---|---|
| (A) the fit's integrated slope against its change in value | -1.26 | +0.34 | -0.23 |
| (B) inside the band: the bridge against the data | -0.14 | +0.43 | -0.85 |
| (C) outside the band: the wind read back against the data | -0.09 | -0.21 | +0.18 |
| (D) the transfer against the balance, linearized | -0.06 | +0.21 | -0.23 |
| **Sum, against the reported difference** | **-1.55 (-1.54)** | **+0.77 (+0.77)** | **-1.13 (-1.13)** |

- **(A) is the comparison, not the model.**
  - The wind is built from `fit`'s slope. The round trip compares with the change in `fit`'s value.
    For a local linear fit the slope at each latitude is that window's regression slope, not the
    derivative of the fitted curve, and where the temperature curves the two differ.
  - Along the path the difference builds up in two places. One is the equatorial structure
    (+0.73 K by the equator at 110 mbar). The other is the last few degrees before the target, where
    the 110 mbar temperatures climb steeply toward the data gap at -31.5 deg (-0.6 K over the last
    5 deg).
  - At 110 mbar the slope integrates to +4.08 K where the value changes by +5.26 K.
- **(B) is ruling 2's bridge.** Inside the band the shear is PCHIP between the edges, not the
  inversion of the data. The data's gradient there implies +0.81 K at 730 mbar, and the bridge
  -0.04 K.
- **(C) is the wind as the model reads it.** `WindField` reads `u` and `du/dln p` linearly between
  the file's nodes, against the tool's continuous construction: 0.09 to 0.21 K.
- **(D) is the transfer.** Its isobar geometry is the model's own mesh against the tool's
  approximation, and the balance is linearized as `T_anchor x` the integral: 0.06 to 0.23 K.

## 4. For the author

1. **The round trip's benchmark.** As written, the round trip mixes (A) into the model's error.
   There are two ways to make the comparison like for like:
   - compare with the integral of the fit's slope, which is what the wind encodes; or
   - build the wind from the derivative of the fit's value, so that slope and value agree by
     construction.

   The first changes only the acceptance. The second changes SPEC_10's `fit` or the tool's use of
   it.
2. **The band.** (B) is the only part that the data could not set. It is up to 0.85 K at 730 mbar,
   on a path that crosses the equator.
3. **The adiabat.** The frozen-rotation `c_p` puts the anchor's own deep profile above the adiabat.
   A heat capacity for equilibrium or normal hydrogen at 100 to 150 K would draw the adiabat
   elsewhere. Which one the comparison should use is the author's.

## 5. Notes

- **Not committed.** The run's inputs and product are under `forward/lindal_iris_ii_110/inputs/`
  and `output/`, which git ignores for every run but the two registered. Only its two control files
  would be committed.
- **The relaxation.** The run was built from an uncommitted tree, so its inputs carry `-dirty`. The
  acceptance replaces `lib.control._refuse_dirty_commit` for its own process, as `step05_3` does.
  After the acceptance commit the run rebuilds clean.
