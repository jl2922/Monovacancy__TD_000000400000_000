# Configuration File For:
# Vacancy Migration Energy (VME) and Formation Energy (VFE) Test Driver
# Author: Junhao Li <streaver91@gmail.com>

DEBUG = 1

# Accuracy related parameters
FMAX_TOL = 1e-3 # absolute
FIRE_MAX_STEPS = 200
MDMIN_MAX_STEPS = 200
MIN_ATOMS = 100
NEB_POINTS = 20
UNCERT_STEPS = 20
EPS = 1e-10
STRESS_DX = 1e-3
# Parameters for Debugging
if DEBUG == 1:
    FMAX_TOL = 1e-2 # absolute

# Logs output
FIRE_LOGFILE = 'fire.log'
MDMIN_LOGFILE = 'mdmin.log'
OUTPUT_BASIS = True

# Extrapolation Parameters
FITS_CNT = [2, 3] # Number of data points used for each fitting
FITS_ORDERS = [
    [0, -3],
    [0, -3],
] # Number of orders included in each fitting
# Fit Results Used (Corresponding to the above)
FITS_VFE_VALUE = 0 # Vacancy Formation Energy
FITS_VFE_UNCERT = [1]
FITS_VME_VALUE = 0 # Vacancy Migration Energy
FITS_VME_UNCERT = [1]

# Conversion constants
e = 1.602176565e-19
A = 1.0e-10
eVoverA32GPa = e / A**3 / 1.0e9 # Convert eV to GPa

# Output Configuration
UNIT_ENERGY = 'eV'
UNIT_LENGTH = 'angstrom'
UNIT_ANGLE = 'degree'
UNIT_PRESSURE = 'GPa'
UNIT_VOLUME = UNIT_LENGTH + '^3'
SPACE_GROUPS = {
    'fcc': 'Fm-3m',
    'bcc': 'Im-3m',
    'sc': 'Pm-3m',
    'diamond': 'Fd-3m',
    'hcp': 'P63/mmc',
}
WYCKOFF_CODES = {
    'fcc': ['4a'],
    'bcc': ['2a'],
    'sc': ['1a'],
    'diamond': ['8a'],
    'hcp': ['2d'],
}
WYCKOFF_SITES = {
    'fcc': [[0.0, 0.0, 0.0]],
    'bcc': [[0.0, 0.0, 0.0]],
    'sc': [[0.0, 0.0, 0.0]],
    'diamond': [[0.0, 0.0, 0.0]],
    'hcp': [[2.0 / 3.0, 1.0 / 3.0, 0.25]],
}
VFE_PROP_ID = 'tag:staff@noreply.openkim.org,2016-04-06:property/monovacancy-neutral-formation-free-energy-crystal-npt'
VME_PROP_ID = 'tag:staff@noreply.openkim.org,2016-04-06:property/monovacancy-neutral-migration-energy-crystal-npt'
