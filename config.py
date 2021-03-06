# Configuration File for Vacancy Test Driver
# Author: Junhao Li <streaver91@gmail.com>

DEBUG = 1

# Accuracy related parameters.
FMAX_TOL = 10.0e-3 # absolute
FIRE_MAX_STEPS = 1000
MDMIN_MAX_STEPS = 1000
MIN_ATOMS = 100
NEB_POINTS = 12
EPS = 1.0e-6
STRESS_DX = 1.0e-3
NUM_SIZES = 3
MIGRATION = (1.0, 0.0, 0.0) # Migration along x by one conventional lattice.
STRAIN = (0.0, 0.0, 0.0)

# Logs output
FIRE_LOGFILE = 'output/fire.log'
MDMIN_LOGFILE = 'output/mdmin.log'
SAVE_BASIS = True
SAVE_JSON = True
REMOVE_CIFS = True
OUTPUT_INSTANCES = False

# Conversion constants.
eV = 1.602176565e-19
angstrom = 1.0e-10
toGPa = eV / angstrom**3 / 1.0e9 # Convert eV to GPa.
voigtEncode = ((0, 1, 2, 1, 2), (0, 1, 2, 0, 0))

# Output information.
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
ALIASES = {
    'vacancy-formation-energy': 'vacancy-formation-free-energy',
    'basis-atom-species': 'species',
    'basis-short-name': 'short-name',
}
