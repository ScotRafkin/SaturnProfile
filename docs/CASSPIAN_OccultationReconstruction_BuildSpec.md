# CASSPIAN Occultation Reconstruction Toolkit — Build Specification

**Document version:** 0.1
**Target implementer:** Claude Code agent in VS Code
**Purpose:** Build a modular Python toolkit to reconstruct atmospheric profiles from radio occultation data, starting with Voyager 2 (Lindal et al. 1985), in a form that supports arbitrary composition assumptions and downstream uncertainty propagation.

---

## 0. Guiding Principles (read first)

1. **Deterministic modules only.** Every module in this spec is a pure, deterministic function or class. No random number generation, no Monte Carlo, no loops over ensembles. Orchestration (including Monte Carlo) is a separate future layer that will call these modules. Do not build it here.

2. **Single source of truth for constants.** All physical constants and species properties live in ONE module (`constants.py`). No constant is ever hard-coded anywhere else.

3. **Composition flows as mole fractions.** Composition is always represented as a set of species mole fractions as a function of height. Derived quantities (mean molecular weight, mean refractivity) are computed from mole fractions by utility functions, never stored or passed independently, to guarantee consistency.

4. **Provenance travels with data.** Every quantity carries a flag indicating whether it is measured, derived, or assumed.

5. **Standard interchange format is netCDF.** All profile data reads and writes go through netCDF with a defined schema (Section 2).

6. **Each module is independently testable.** Provide a unit test for each module. Tests must run without network access.

7. **No prose postambles in code.** Docstrings and comments only. No attribution to any AI tool in code or commits.

8. **Formatting standard:** Do not use em dashes or en dashes in comments, docstrings, or documentation prose.

---

## 1. Repository Structure

```
occultation_recon/
├── occultation_recon/
│   ├── __init__.py
│   ├── constants.py            # Module 0: all physical constants and species data
│   ├── composition_utils.py    # Module U: derived quantities from mole fractions
│   ├── io_standard.py          # Module 2: read/write standard netCDF format
│   ├── number_density.py       # Module 3: ensure n(h) is available
│   ├── composition.py          # Module 4: composition profiles as mole fractions
│   ├── refractivity.py         # Module 5: refractivity from n and composition
│   └── build_table1_netcdf.py  # Module 1: one-time Table I digitization
├── data/
│   ├── lindal1985_table1.csv   # raw digitized Table I (input to Module 1)
│   └── lindal1985_voyager2.nc  # standard-format output of Module 1
├── tests/
│   ├── test_constants.py
│   ├── test_composition_utils.py
│   ├── test_io_standard.py
│   ├── test_number_density.py
│   ├── test_composition.py
│   └── test_refractivity.py
├── docs/
│   └── this_spec.md
├── pyproject.toml
└── README.md
```

Use Python 3.11+. Dependencies: numpy, scipy, netCDF4 (or xarray, implementer's choice but be consistent), matplotlib for the graphical outputs. Use pytest for tests.

---

## 2. The Standard netCDF Format (schema)

This schema is the interchange contract between all modules. Every profile dataset conforms to it.

### Dimensions
- `level`: the number of vertical levels in the profile

### Coordinate variable
- `height` (float64, dimension `level`): geometric altitude in meters, referenced to a stated datum. Monotonic. This is typically the measured independent variable.

### Data variables (all dimension `level`, all float64, all optional except height)
- `pressure`: pressure in pascals
- `temperature`: temperature in kelvin
- `number_density`: total number density in m^-3
- `refractivity`: refractivity N = (n_refractive - 1) * 1e6, dimensionless (scaled)
- `nh3_mixing_ratio`: ammonia mole fraction (or ppm, see attribute), optional
- One variable per species mole fraction if composition is stored: `x_H2`, `x_He`, `x_CH4`, etc.

### Required per-variable attributes
- `units`: string (e.g., "Pa", "K", "m", "m-3")
- `provenance`: string, one of `"measured"`, `"derived"`, `"assumed"`
- `long_name`: human-readable description

### Required global attributes
- `title`: dataset title
- `source`: origin of the data (e.g., "Lindal et al. 1985, AJ 90, 1136, Table I, Voyager 2 ingress")
- `height_datum`: description of the altitude reference (e.g., "1-bar pressure level per source")
- `latitude_deg`: planetographic or planetocentric latitude of the profile (state which in the attribute string)
- `longitude_deg`: if applicable
- `history`: creation provenance string
- `composition_note`: free text describing the composition assumption under which any derived quantities were produced

### Notes
- A dataset need not contain all variables. A raw Table I dataset contains height, pressure, temperature, nh3_mixing_ratio. A synthetic dataset might contain height and number_density only. Downstream modules inspect which variables are present and act accordingly.
- Provide a small documented Python dataclass or dict schema in `io_standard.py` that defines this contract in code.

---

## 3. Module 0: constants.py

Single source of truth. No computation, only data and simple named constants.

### Physical constants (SI, CODATA values)
```python
BOLTZMANN = 1.380649e-23        # J/K, exact
ATOMIC_MASS_UNIT = 1.66053906660e-27  # kg
AVOGADRO = 6.02214076e23        # 1/mol, exact
LOSCHMIDT = 2.6867811e25        # m^-3, number density at 273.15 K, 101325 Pa
GAS_CONSTANT = 8.314462618      # J/(mol K)
```

### Species data table
Provide a dictionary keyed by species name. Each entry contains:
- `molar_mass_kg_per_mol`: molar mass in kg/mol
- `refractivity_stp`: refractivity N_STP = (n_refractive - 1) * 1e6 for the pure gas at STP (273.15 K, 101325 Pa), at radio/microwave frequency
- `refractivity_reference`: string citing the source of the refractivity value
- `is_polar`: bool, True for species with significant orientational polarizability (e.g., NH3), which have frequency-dependent refractivity

Populate with these species. Use these approximate values but mark each with its reference string, and leave clear TODO markers where a more authoritative value should be substituted:

```python
SPECIES = {
    "H2":  {"molar_mass_kg_per_mol": 2.01588e-3,  "refractivity_stp": 136.0, "is_polar": False, ...},
    "He":  {"molar_mass_kg_per_mol": 4.002602e-3, "refractivity_stp": 35.0,  "is_polar": False, ...},
    "CH4": {"molar_mass_kg_per_mol": 16.0425e-3,  "refractivity_stp": 430.0, "is_polar": False, ...},
    "NH3": {"molar_mass_kg_per_mol": 17.0305e-3,  "refractivity_stp": 375.0, "is_polar": True,  ...},
}
```

IMPORTANT: The refractivity_stp values above are approximate placeholders. Add a clearly marked comment block instructing that these should be replaced with authoritative radio-frequency refractivity values from a spectroscopic reference before scientific use. The CH4 and NH3 values in particular carry uncertainty and NH3 is frequency-dependent.

### Derived helper (convenience, not a substitute for the utility module)
Provide a function to convert `refractivity_stp` to per-particle refractivity volume:
```python
def per_particle_refractivity(species_name):
    # returns refractivity per unit number density such that
    # N = number_density * per_particle_refractivity summed over species
    # per_particle = refractivity_stp / LOSCHMIDT
```

### Test (test_constants.py)
- Assert all species have all required keys.
- Assert molar masses are positive and in plausible range.
- Assert per_particle_refractivity(H2) recovers refractivity_stp when multiplied by LOSCHMIDT.

---

## 4. Module U: composition_utils.py

Derived quantities from mole fractions. Pure functions. This is where mean molecular weight and mean refractivity are computed, guaranteeing they always come from the same mole fractions.

### Functions

```python
def mean_molar_mass(mole_fractions: dict) -> float:
    """
    mole_fractions: dict mapping species name -> mole fraction (scalar).
    Returns mean molar mass in kg/mol = sum_i x_i * M_i.
    Must validate that mole fractions sum to 1 within a tolerance (e.g. 1e-6);
    raise or warn if not, implementer's choice, but do not silently renormalize
    without flagging.
    """

def mean_molecular_weight_amu(mole_fractions: dict) -> float:
    """Mean molecular weight in atomic mass units (dimensionless mu)."""

def mean_refractivity_per_particle(mole_fractions: dict) -> float:
    """
    Returns the composition-weighted per-particle refractivity R_mean such that
    refractivity N = number_density * R_mean.
    R_mean = sum_i x_i * per_particle_refractivity(species_i)
    """
```

### Vectorized variants
Provide variants that accept mole fractions as a dict of species -> array over height, returning arrays. These are what the height-dependent composition case uses.

```python
def mean_molar_mass_profile(mole_fraction_profiles: dict) -> np.ndarray
def mean_refractivity_per_particle_profile(mole_fraction_profiles: dict) -> np.ndarray
```

### Test (test_composition_utils.py)
- Pure H2 returns H2 molar mass and H2 per-particle refractivity.
- 50/50 H2/He returns the arithmetic mean of the two.
- Mole fractions not summing to 1 triggers the validation path.
- Profile (array) variants agree with scalar variants applied level by level.

---

## 5. Module 1: build_table1_netcdf.py

One-time digitization of Lindal et al. 1985 Table I into the standard netCDF format.

### Input
A CSV file `data/lindal1985_table1.csv` containing the Table I data. The implementer should create this CSV from the following data (Voyager 2 ingress, from Table I of the paper). The table has TWO halves printed side by side in the paper; concatenate them into one monotonic profile ordered by increasing pressure.

Columns: `pressure_mbar`, `temperature_K`, `nh3_ppm`, `altitude_km`

Data (pressure in mbar, T in K, NH3 in ppm where given else blank, altitude in km relative to 1 bar):
```
0.20,138.7,,376.7
0.25,143.2,,363.9
0.32,147.3,,350.7
0.40,146.2,,337.4
0.50,143.8,,324.3
0.63,142.2,,311.4
0.79,141.2,,298.6
1.00,141.6,,285.8
1.26,142.5,,273.0
1.58,140.3,,260.2
2.00,137.7,,247.7
2.51,136.9,,235.3
3.16,133.4,,223.1
3.98,129.0,,211.3
5.01,123.5,,199.9
6.31,117.7,,189.0
7.94,112.0,,178.7
10.00,104.8,,168.9
12.59,99.5,,159.8
15.85,95.0,,151.0
19.95,90.5,,142.8
23.99,87.6,,136.4
28.84,85.3,,130.2
34.67,83.9,,124.1
41.69,83.1,,118.1
50.12,82.3,,112.2
60.26,82.0,,106.3
72.44,82.4,,100.4
87.10,82.7,,94.4
100.00,83.4,,90.0
114.82,84.3,,85.5
131.83,85.2,,80.9
151.36,86.3,,76.3
173.78,87.9,,71.6
199.53,89.8,,66.9
229.09,91.4,,62.0
251.19,92.4,,58.7
275.42,93.5,,55.4
302.00,94.8,,52.0
331.13,96.5,,48.6
363.08,98.4,,45.1
398.11,100.7,,41.6
436.52,103.2,,37.9
478.63,106.0,,34.2
501.19,107.5,,32.3
524.81,109.1,,30.3
549.54,110.7,,28.4
575.44,112.4,,26.4
602.56,114.1,,24.4
630.96,115.8,,22.3
660.69,117.6,,20.2
691.83,119.4,,18.1
724.44,121.3,,16.0
758.58,123.2,,13.8
794.33,125.0,,11.6
831.76,127.0,2.6,9.3
870.96,128.9,5.1,7.0
912.01,130.8,7.1,4.7
954.99,132.8,9.5,2.4
1000.00,134.8,10.9,0.0
1047.13,136.8,,-2.4
1096.48,138.9,21.2,-4.9
1148.15,140.9,41.9,-7.4
1202.26,142.9,49.2,-9.9
1258.93,145.0,66.9,-12.5
1298.48,146.2,,-14.1
```

### Processing
1. Read the CSV.
2. Convert units: pressure mbar to Pa (multiply by 100), altitude km to m (multiply by 1000), NH3 ppm to mole fraction (multiply by 1e-6) where present.
3. Order by increasing pressure (equivalently decreasing altitude).
4. Write to `data/lindal1985_voyager2.nc` conforming to the Section 2 schema.

### Provenance flags for this dataset
- `height`: provenance = "measured" (from spacecraft geometry and occultation ray tracing)
- `pressure`: provenance = "derived" (Lindal hydrostatic integration)
- `temperature`: provenance = "derived" (Lindal, from number density and pressure via ideal gas law)
- `nh3_mixing_ratio`: provenance = "derived" (from microwave absorption)

### Global attributes for this dataset
- title: "Voyager 2 Radio Occultation Ingress Profile, Saturn"
- source: "Lindal, Sweetnam, Eshleman 1985, AJ 90, 1136, Table I"
- height_datum: "1-bar pressure level as defined in source"
- latitude_deg: "36.5 N planetographic (ingress, evening terminator); see source Section on Voyager 2"
- composition_note: "Source p and T derived assuming 94% H2 / 6% He by number, mean molecular mass 2.135 amu, per Lindal et al. 1985"

### Test (test_io_standard.py covers reading; add a check here)
- After writing, read back and assert pressure is monotonic increasing.
- Assert the 1000 mbar level has altitude 0.0 within rounding (it is the 1-bar reference).
- Assert temperature at 1000 mbar is 134.8 K.

---

## 6. Module 2: io_standard.py

Read and write the standard netCDF format.

### Functions
```python
def write_profile(path: str, data: dict, attributes: dict) -> None:
    """
    data: dict mapping variable name -> (array, units, provenance, long_name)
          must include 'height'
    attributes: dict of global attributes
    Writes a netCDF file conforming to the Section 2 schema.
    """

def read_profile(path: str) -> ProfileData:
    """
    Returns a ProfileData object (dataclass) exposing:
      .height (np.ndarray)
      .variables: dict of name -> np.ndarray
      .provenance: dict of name -> str
      .units: dict of name -> str
      .attributes: dict of global attributes
      .has(name) -> bool   convenience to check variable presence
    """
```

### Validation on read (keep reasonable, do not over-engineer)
- Require `height` present.
- Require height monotonic (increasing or decreasing); record direction.
- Warn (do not fail) if a variable lacks a provenance attribute.
- Check units strings are present for all variables.

### Test (test_io_standard.py)
- Round trip: write a small synthetic profile, read it back, assert arrays and attributes match.
- Reading a file without height raises a clear error.
- has() correctly reports presence and absence.

---

## 7. Module 3: number_density.py

Ensure total number density n(h) is available, computing it from whatever the source provides.

### Primary function
```python
def ensure_number_density(profile: ProfileData,
                          composition=None) -> np.ndarray:
    """
    Return total number density in m^-3 at each level of the profile.

    Logic (dispatch on what the profile contains):
      Case A: profile has 'number_density' already -> return it directly.
      Case B: profile has 'pressure' and 'temperature' -> compute
              n = pressure / (BOLTZMANN * temperature).  This is the Lindal case.
      Case C: profile has 'refractivity' and composition is provided ->
              compute n = refractivity / mean_refractivity_per_particle(composition).
      Otherwise: raise a clear error stating what is missing.

    The 'composition' argument is only needed for Case C. For the initial
    build, implement Case A and Case B. Stub Case C with a clear
    NotImplementedError message so it is obvious where to add it later.
    """
```

### Notes
- Case B is the exact algebraic identity n = p / (k_B T). It carries no integration and no composition assumption. This is the fundamental recovery step for Lindal data.
- Do NOT integrate anything in this module. This module only applies local algebraic relations.

### Test (test_number_density.py)
- Case B: feed a profile with known p and T, assert n = p/(kB T) exactly.
- Case A: feed a profile that already has number_density, assert pass-through.
- Feed a profile with neither, assert the error is raised.
- Sanity: for Lindal 1000 mbar / 134.8 K, assert n is approximately 5.37e25 m^-3 (compute the expected value in the test from the identity, do not hard-code blindly).

---

## 8. Module 4: composition.py

Produce composition as mole fractions versus height. Multiple constructors for different sources. All return the SAME representation: a dict mapping species name -> np.ndarray over the profile levels (or a callable/interpolator evaluated onto the profile's height grid).

### Constructors
```python
def constant_composition(species_fractions: dict, n_levels: int) -> dict:
    """
    species_fractions: dict species -> scalar mole fraction, must sum to 1.
    Returns dict species -> array of length n_levels, each constant.
    Validate the sum to 1 within tolerance.
    """

def lindal_composition(n_levels: int) -> dict:
    """
    The composition assumed by Lindal et al. 1985 for their retrieval:
    94% H2, 6% He by number (their adopted refractivity mixture).
    Returns constant profiles. This is the composition that reproduces
    the measured refractivity when combined with Lindal's recovered n(h).
    """

def analytical_composition(height: np.ndarray, spec) -> dict:
    """
    Composition from an analytical prescription (e.g. a smooth transition
    in mole fraction with height). 'spec' defines the functional form.
    Stub this with a simple example (e.g. constant below a knee, linear
    ramp above) and a clear extension point. Not critical for first build.
    """

def composition_from_file(path: str, height: np.ndarray) -> dict:
    """
    Read a composition profile from an external file (e.g. Julie Moses
    photochemical model output) and interpolate onto the given height grid.
    Stub with a documented expected file format and a NotImplementedError,
    to be completed when a real Moses file is in hand.
    """
```

### Requirements
- Every constructor returns mole fractions that sum to 1 at every level (validate).
- The representation must be directly consumable by composition_utils functions.
- Do not compute mean molecular weight or refractivity here. That is composition_utils' job. This module only produces mole fractions.

### Test (test_composition.py)
- constant_composition with 0.94 H2 / 0.06 He returns arrays of the right length that sum to 1.
- lindal_composition returns 0.94 / 0.06.
- A composition not summing to 1 triggers validation.
- Output feeds mean_molar_mass_profile and mean_refractivity_per_particle_profile without error and gives sensible numbers (mu approx 2.135 amu for the Lindal case).

---

## 9. Module 5: refractivity.py

The pure physical relation between number density, composition, and refractivity. This is the fundamental data product generator.

### Core function (belongs conceptually to the utility layer, place it here or in composition_utils, implementer's choice, but it must be a single pure function)
```python
def refractivity_from_number_density(number_density: np.ndarray,
                                     composition: dict) -> np.ndarray:
    """
    N(h) = number_density(h) * mean_refractivity_per_particle(composition, h)

    Pure physical relation. Agnostic about the source of its inputs.
      - Feed Lindal's n(h) and Lindal's composition -> recover measured N(h) (the anchor).
      - Feed any synthetic n(h) and any composition -> that atmosphere's N(h).
    Returns refractivity N = (n_refractive - 1) * 1e6, dimensionless.
    """
```

Also provide the inverse, which downstream reconstruction will use:
```python
def number_density_from_refractivity(refractivity: np.ndarray,
                                     composition: dict) -> np.ndarray:
    """
    n(h) = refractivity(h) / mean_refractivity_per_particle(composition, h)
    The inverse relation. Used to reinterpret a fixed measured refractivity
    under an alternative composition.
    """
```

### Product generation function
```python
def generate_refractivity_product(profile_path: str,
                                  composition,
                                  output_path: str,
                                  plot_path: str = None) -> None:
    """
    High level driver for the 'fundamental data product':
      1. Read the profile (io_standard.read_profile).
      2. Ensure number density (number_density.ensure_number_density).
      3. Compute refractivity (refractivity_from_number_density).
      4. Write an output netCDF containing height, number_density, refractivity,
         and the composition mole fractions used, with full provenance and a
         composition_note global attribute documenting the assumption.
      5. If plot_path given, produce a figure: refractivity vs height (log x
         if it spans decades) and number density vs height, with a title
         stating the composition assumption. Save to plot_path.
    This function documents its input assumptions in both the netCDF metadata
    and the figure, per the requirement that the product be self-documenting.
    """
```

### Test (test_refractivity.py)
- refractivity_from_number_density then number_density_from_refractivity round trips to the original n within floating point tolerance for a fixed composition.
- For the Lindal profile with Lindal composition, the recovered refractivity is positive, monotonic-ish with height, and of plausible magnitude (spot check one level by hand in the test).
- Feeding a different composition changes N by the ratio of mean refractivities, and the test asserts that ratio explicitly (this is the key composition-sensitivity behavior).

---

## 10. Build and Validation Sequence

Implement and test in this order. Each step must pass its test before proceeding.

1. `constants.py` + `test_constants.py`
2. `composition_utils.py` + `test_composition_utils.py`
3. `io_standard.py` + `test_io_standard.py`
4. `build_table1_netcdf.py` (create the CSV, run it, produce `lindal1985_voyager2.nc`)
5. `number_density.py` + `test_number_density.py`
6. `composition.py` + `test_composition.py`
7. `refractivity.py` + `test_refractivity.py`
8. End to end: run `generate_refractivity_product` on the Lindal netCDF with `lindal_composition`, produce the output netCDF and the plot. Visually inspect that number density decreases smoothly with height and refractivity spans the expected range.

---

## 11. Explicit Non-Goals (do NOT build these yet)

- No hydrostatic integration. That is the next phase, building on the refractivity product.
- No Monte Carlo or ensemble machinery. That is orchestration, layered on later.
- No top-boundary-pressure handling. Belongs to the hydrostatic phase.
- No gravity model. Belongs to the hydrostatic phase.
- No radiative transfer. Separate effort entirely.

Keep this build strictly to the deterministic modules that produce the fundamental refractivity data product from occultation input.

---

## 12. Coding Standards

- Type hints on all public functions.
- Docstrings on all modules, classes, and public functions, stating inputs, outputs, units, and assumptions.
- Units stated explicitly in every docstring and every netCDF variable attribute.
- No hard-coded physical constants outside constants.py.
- Each module runnable and testable in isolation.
- pytest for all tests; tests must not require network access.
- Prose in comments and docstrings must not use em dashes or en dashes.
- Do not include AI attribution in code comments or commit messages.

