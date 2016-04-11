# Configuration File For:
# Vacancy Migration Energy (VME) and Formation Energy (VFE) Test Driver
# Author: Junhao Li <streaver91@gmail.com>

DEBUG = 1

# Accuracy related parameters
FMAX_TOL = 1.0e-3 # absolute
FIRE_MAX_STEPS = 500
MDMIN_MAX_STEPS = 500
MIN_ATOMS = 100
NEB_POINTS = 15
UNCERT_STEPS = 20
EPS = 1.0e-6
STRESS_DX = 1.0e-3
NUM_SIZES = 3

# Parameters for Debugging
if DEBUG == 1:
    FMAX_TOL = 50.0e-3 # absolute
    NEB_POINTS = 15

# Logs output
FIRE_LOGFILE = 'fire.log'
MDMIN_LOGFILE = 'mdmin.log'
SAVE_BASIS = True

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
PROPERTY_DEFINITIONS = [
    'monovacancy-neutral-formation-free-energy-crystal-npt',
    'monovacancy-neutral-migration-energy-crystal-npt',
    'monovacancy-neutral-elastic-dipole-tensor-crystal-npt',
    'monovacancy-neutral-migration-saddle-point-elastic-dipole-tensor-crystal-npt',
    'monovacancy-neutral-defect-strain-tensor-crystal-npt',
    'monovacancy-neutral-relaxation-volume-crystal-npt'
]
