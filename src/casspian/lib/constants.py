"""Physical constants of nature, CODATA 2018.

SPEC_00 section 4: constants of nature live here and nowhere else. Nothing Saturn specific
belongs in this module. Every product records `codata_release` in its global attributes so that
the values a file was computed with are traceable without being editable.

Manuscript equations implemented: none. This module holds values only.

Naming follows SPEC_00 section 3.1: nothing in the package is named `R`, so the molar gas
constant is spelled out. The mean specific gas constant `R_bar` and the mean molecular
refractivity `script_R` are formed in `lib.thermo` from composition, never here.
"""

CODATA_RELEASE = "2018"

# CODATA name: Boltzmann constant. Exact by the 2019 SI definition. J K^-1.
BOLTZMANN_CONSTANT = 1.380649e-23

# CODATA name: Avogadro constant. Exact by the 2019 SI definition. mol^-1.
AVOGADRO_CONSTANT = 6.02214076e23

# CODATA name: atomic mass constant. kg. Relative standard uncertainty 3.0e-10.
ATOMIC_MASS_CONSTANT = 1.66053906660e-27

# CODATA name: molar gas constant. Exact, the product of the two constants above.
# J mol^-1 K^-1.
MOLAR_GAS_CONSTANT = 8.314462618

# CODATA name: Newtonian constant of gravitation. m^3 kg^-1 s^-2.
# Relative standard uncertainty 2.2e-5, the least well known constant used here.
NEWTONIAN_CONSTANT_OF_GRAVITATION = 6.67430e-11

# CODATA name: Loschmidt constant (273.15 K, 101 325 Pa). Exact, being
# STANDARD_PRESSURE / (BOLTZMANN_CONSTANT * STANDARD_TEMPERATURE). m^-3.
# Used to convert published per molecule refractivities, which are quoted as
# refractivity of a gas at standard temperature and pressure, into m^3 per molecule
# (SPEC_00 section 5).
LOSCHMIDT_CONSTANT = 2.686780111e25

# CODATA name: standard-state temperature. Exact, definitional. K.
STANDARD_TEMPERATURE = 273.15

# CODATA name: standard-state pressure. Exact, definitional. Pa.
STANDARD_PRESSURE = 101325.0
