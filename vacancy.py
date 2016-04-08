# Vacancy Class Definition File For:
# Vacancy Migration Energy (VME) and Formation Energy (VFE) Test Driver
# Author: Junhao Li <streaver91@gmail.com>

# Python Modules
import sys
import re
import json
import math
from collections import OrderedDict
from scipy.optimize import fmin
from scipy.interpolate import interp1d
import numpy as np

# ASE Modules
try:
    from ase.lattice import bulk
    print 'Imported bulk from ase.lattice' # For ASE version 3.9
except ImportError:
    from ase.structure import bulk
    print 'Imported bulk from ase.structure' # For ASE version <= 3.8
from ase.optimize import FIRE, QuasiNewton, MDMin
from ase.data import chemical_symbols
from ase.data import reference_states
from ase import Atoms, Atom
from ase.io import write
from ase.constraints import FixAtoms
from ase.neb import NEB

# KIM Modules
from kimcalculator import *
from kimservice import KIM_API_get_data_double

# Load Configuration
import config as C

class Vacancy(object):
    # Class for calculating vacancy formation energy and relaxation volume
    def __init__(self, elem, model, lattice, latticeConsts, migration):
        # Process Inputs
        self.elem = elem
        self.model = model
        self.lattice = lattice
        self.latticeConsts = np.array(latticeConsts)
        self.migration = np.array(migration) # Vacancy migration vector
        self._printInputs()

        # Create basis for constructing supercell
        self.basis = self._createBasis()

        # Initialize minimization
        # self.FIREUncert = 0

        # Determine size of the supercell
        # nAtoms = atoms.get_number_of_atoms()
        # factor = math.pow(8 / nAtoms, 0.333)
        # self.cellSizeMin = int(math.ceil(factor * C.CELL_SIZE_MIN))
        # self.cellSizeMax = self.cellSizeMin + 2
        # print 'Cell Size Min:', self.cellSizeMin
        # print 'Cell Size Max:', self.cellSizeMax
        # print 'Smallest System Size:', nAtoms * self.cellSizeMin**3
        # print 'Largest System Size:', nAtoms * self.cellSizeMax**3
        # print 'Model Cutoff:', KIM_API_get_data_double(atoms.calc.pkim, 'cutoff')[0]

    def _printInputs(self):
        # Print out the inputs
        print 'Inputs:'
        print 'Element: ', self.elem
        print 'Model: ', self.model
        print 'Lattice: ', self.lattice
        print 'Lattice Constants: ', self.latticeConsts

    def _createBasis(self):
        # Create basic atoms
        # Cache to local scope
        elem = self.elem
        lattice = self.lattice
        latticeConsts = self.latticeConsts
        # Check whether lattice is valid
        supportedLattices = ['sc', 'fcc', 'bcc', 'diamond', 'hcp']
        assert lattice in supportedLattices, 'lattice not supported'
        # Create basis according to structure
        if lattice == 'hcp':
            # Create cubic lattice to accommodate MI_OPBC
            a = latticeConsts[0]
            c = latticeConsts[1]
            rt3 = np.sqrt(3.0)
            basis = Atoms(
                elem + '4',
                positions = np.array([
                    [0.0, 0.0, 0.0],
                    [0.5, 0.5 * rt3, 0.0],
                    [0.5, 0.5 / 3.0 * rt3, 0.5],
                    [0.0, (0.5 + 0.5 / 3.0) * rt3, 0.5]
                ]),
                cell = np.array([1.0, rt3, 1.0]),
                pbc = [1, 1, 1]
            )
            basis.set_cell(np.array([a, a * rt3, c]), scale_atoms = True)
        else:
            basis = bulk(
                elem,
                a = latticeConsts[0],
                crystalstructure = lattice,
                cubic = True,
            )
        if C.OUTPUT_BASIS == True:
            write('output/basis.cif', basis, format = 'cif')
        return basis

    def _createSupercell(self, size):
        # Create supercell w/ basis repeated 'size' times in each direction
        supercell = self.basis.copy()
        supercell.set_calculator(KIMCalculator(self.model))
        supercell *= (size, size, size)
        return supercell


    def _relaxAtoms(self, atoms, tol = C.FIRE_TOL, steps = C.FIRE_MAX_STEPS):
        fire = FIRE(atoms, logfile = C.FIRE_LOGFILE)
        fire.run(fmax = tol, steps = steps)
        fireSteps = fire.get_number_of_steps()
        return fireSteps < steps ? True : False

    def _relaxPath(self, images):
        neb = NEB(images)
        neb.interpolate()
        mdmin = MDMin(neb, logfile = C.MDMIN_LOGFILE)
        mdmin.run(fmax = C.MDMIN_TOL, steps = C.MDMIN_MAX_STEPS)
        mdminSetps = mdmin.get_number_of_steps()
        return mdminSetps < C.MDMIN_MAX_STEPS ? True : False

    def _findAtomId(self, positions, position):
        return 0

    def _getSizeResult(self, size):
        print 'Calculating Supercell of Size', size, '...'

        # Object for output
        res = {}

        # Setup supercells
        supercell = self._createSupercell(size)
        en0 = supercell.get_potential_energy()
        nAtoms = supercell.get_number_of_atoms()
        del supercell[0]

        # Initial and final unrelaxed state for migration
        initial = supercell.copy()
        initial.set_calculator(KIMCalculator(self.model))
        final = supercell.copy()
        finalPositions = final.get_positions()
        # Assuming migration is from (0,0,0)
        movedAtomId = self._findAtomId(finalPositions, self.migration)
        finalPositions[movedAtomId] = np.array([0.0, 0.0, 0.0])
        final.set_positions(finalPositions)
        final.set_calculator(KIMCalculator(self.model))

        # Relax initial and final state
        tmp = self._relaxAtoms(initial)
        assert tmp, 'Maximum FIRE Steps Reached.'
        tmp = self._relaxAtoms(final)
        assert tmp, 'Maximum FIRE Steps Reached.'

        # Interpolate between initial and final
        images = []
        for i in range(C.NEB_POINTS):
            image = initial.copy()
            image.set_calculator(KIMCalculator(self.model))
            images.append(image)
        images.append(final.copy())
        tmp = self._relaxPath(images)
        assert tmp, 'Maximum MDMin Steps Reached.'

        # Obtaining vacancy formation energy
        relaxedInitial = initial.copy()
        relaxedInitial.set_calculator(KIMCalculator(self.model))
        en1 = relaxedInitial.get_potential_energy()
        enCoh = en0 / nAtoms
        res['cohesive-energy'] = enCoh
        res['vacancy-formation-energy'] = en1 - enCoh * (nAtoms - 1)

        # Obtain fmax induced error (md error)
        self._relaxAtoms(
            relaxedInitial,
            tol = C.FIRE_TOL * C.EPS,
            steps = C.UNCERT_STEPS
        )
        en1Refined = relaxedInitial.get_potential_energy()
        res['md-uncert'] = np.abs(en1Refined - en1)

        # Obtain elastic dipole tensor


        # Obtain vacancy migration energy


        # Obtain vacancy
        sys.exit(0)
        # removedAtomPosition = atoms.get_positions()[0]


        # Get Initial State


        # Unrelaxed Energy
        # enUnrelaxed = initial.get_potential_energy()
        # print 'Assert n initial:', initial.get_number_of_atoms()
        # unrelaxedFormationEnergy = enUnrelaxed - enAtoms * (nAtoms - 1) / nAtoms

        # Output Before Calculation

        # print 'Unrelaxed Energy w/ Vacancy:', enUnrelaxed
        # print 'Unrelaxed Formation Energy:', unrelaxedFormationEnergy

        # Test Moving Epsilon
        # for stepLenFactor in range(10):
        #     stepLen = 0.01 * 2**stepLenFactor # A
        #     tmpPositions = initial.get_positions()
        #     tmpPositions[0][0] = tmpPositions[0][0] + stepLen
        #     initial.set_positions(tmpPositions)
        #     enMoved = initial.get_potential_energy()
        #     print 'Moving Energy (', str(stepLen), 'A)', enMoved - enUnrelaxed
        #     tmpPositions[0][0] = tmpPositions[0][0] - stepLen
        #     initial.set_positions(tmpPositions)
        # print tmpPositions

        # Relaxation
        # dyn = FIRE(initial, logfile = C.FIRE_LOG)
        # dyn.run(fmax = C.FIRE_TOL, steps = C.FIRE_MAX_STEPS)
        # dynSteps = dyn.get_number_of_steps()
        # assert dynSteps < C.FIRE_MAX_STEPS, 'FIRE Steps Maximum Reached.'
        # en1 = initial.get_potential_energy()

        # Perform more steps and obtain FIRE Uncertainty
        # self.en1Uncert =

        # Get Final State
        tmpPositions = atoms.get_positions()
        tmpPositions[0] = removedAtomPosition
        atoms.set_positions(tmpPositions)
        final = atoms.copy()
        final.set_calculator(KIMCalculator(self.model))
        dyn = FIRE(final, logfile = C.FIRE_LOG)
        dyn.run(fmax = C.FIRE_TOL, steps = C.FIRE_MAX_STEPS)
        if dyn.get_number_of_steps() >= C.FIRE_MAX_STEPS:
            print '[ERR1047]FIRE Failed: Steps Exceeds Maximum.'
            print 'Vacancy Formation Energy May Not Be Accurate.'

        # Calculate VFE
        enInitial = initial.get_potential_energy()
        formationEnergy = enInitial - enAtoms * (nAtoms - 1) / nAtoms

        # Spline Interpolation
        x = np.arange(0, C.NEB_POINTS + 2)
        y = np.array([image.get_potential_energy() for image in images])
        f = interp1d(x, -y, kind = 'cubic') # f(x) = -y(x) in order to use fmin
        xmax = fmin(f, x[C.NEB_POINTS / 2 + 1], ftol = C.FMIN_FTOL)
        ymax = -f(xmax)
        enSaddle = ymax[0]
        migrationEnergy = enSaddle - enInitial
        
        # Output Results
        print 'Formation Energy:', formationEnergy
        print 'Migration Energy:', migrationEnergy
        return migrationEnergy, formationEnergy

    def _getAngle(self, vec1, vec2):
        # Get angle between two vectors in degrees (always between 0 - 180)
        vec1Unit = vec1 / np.linalg.norm(vec1)
        vec2Unit = vec2 / np.linalg.norm(vec2)
        angle = np.arccos(np.dot(vec1Unit, vec2Unit))
        if np.isnan(angle):
            return 0.0
        angle = angle * 180.0 / np.pi
        return angle

    def _getFit(self, xdata, ydata, orders):
        # Polynomial Fitting with Specific Orders
        A = []
        print '[Fitting]'
        print 'xdata:', xdata
        print 'ydata:', ydata
        print 'Orders:', orders
        for order in orders:
            A.append(np.power(xdata * 1.0, order))
        A = np.vstack(A).T
        res = np.linalg.lstsq(A, ydata)
        print 'Results:', res
        return res[0]

    def _extrapolate(self, sizes, valueBySize, valueFitId, uncertFitsId, systemUncert):
        # Extrapolate to Obtain VFE and VME of Infinite Size
        naSizes = np.array(sizes)
        naValueBySize = np.array(valueBySize)
        valueFitsBySize = []
        for i in range(len(C.FITS_CNT)):
            cnt = C.FITS_CNT[i]
            orders = C.FITS_ORDERS[i]
            print 'Fitting w/', cnt, 'points, including orders', orders
            valueFits = []
            for j in range(0, len(sizes) - cnt + 1):
                xdata = naSizes[j:(j + cnt)]
                valueFits.append(self._getFit(
                    xdata,
                    valueBySize[j:(j + cnt)],
                    orders
                )[0])
            valueFitsBySize.append(valueFits)

        # Get Source Value
        valueFitCnt = C.FITS_CNT[valueFitId]
        maxSizeId = len(sizes) - 1
        sourceValue = valueFitsBySize[valueFitId][maxSizeId - valueFitCnt + 1]

        # Get Source Uncertainty (Statistical)
        sourceUncert = 0
        for uncertFitId in uncertFitsId:
            uncertFitCnt = C.FITS_CNT[uncertFitId]
            uncertValue = valueFitsBySize[uncertFitId][maxSizeId - uncertFitCnt + 1]
            sourceUncert = max([abs(uncertValue - sourceValue), sourceUncert])

        # Include systematic error, assuming it is independent of statistical errors
        sourceUncert = math.sqrt(sourceUncert**2 + systemUncert**2)
        return sourceValue, sourceUncert

    def _getCrystalInfo(self):
        unitBulk = self.atoms
        unitCell = unitBulk.get_cell()
        hostInfo = OrderedDict([
            ('host-cauchy-stress', V([0, 0, 0, 0, 0, 0], C.UNIT_PRESSURE)),
            ('host-short-name', V([self.lattice])),
            ('host-a', V(np.linalg.norm(unitCell[0]), C.UNIT_LENGTH)),
            ('host-b', V(np.linalg.norm(unitCell[1]), C.UNIT_LENGTH)),
            ('host-c', V(np.linalg.norm(unitCell[2]), C.UNIT_LENGTH)),
            ('host-alpha', V(self._getAngle(unitCell[1], unitCell[2]), C.UNIT_ANGLE)),
            ('host-beta', V(self._getAngle(unitCell[2], unitCell[0]), C.UNIT_ANGLE)),
            ('host-gamma', V(self._getAngle(unitCell[0], unitCell[1]), C.UNIT_ANGLE)),
            ('host-space-group', V(C.SPACE_GROUPS[self.lattice])),
            ('host-wyckoff-multiplicity-and-letter', V(C.WYCKOFF_CODES[self.lattice])),
            ('host-wyckoff-coordinates', V(C.WYCKOFF_SITES[self.lattice])),
            ('host-wyckoff-species', V([self.elem] * len(C.WYCKOFF_CODES[self.lattice]))),
        ])
        reservoirInfo = OrderedDict([
            ('reservoir-cohesive-potential-energy', V(unitBulk.get_potential_energy(), C.UNIT_ENERGY)),
            ('reservoir-short-name', V([self.lattice])),
            ('reservoir-cauchy-stress', V([0.0, 0.0, 0.0, 0.0, 0.0, 0.0], C.UNIT_PRESSURE)),
            ('reservoir-a', V(np.linalg.norm(unitCell[0]), C.UNIT_LENGTH)),
            ('reservoir-b', V(np.linalg.norm(unitCell[1]), C.UNIT_LENGTH)),
            ('reservoir-c', V(np.linalg.norm(unitCell[2]), C.UNIT_LENGTH)),
            ('reservoir-alpha', V(self._getAngle(unitCell[1], unitCell[2]), C.UNIT_ANGLE)),
            ('reservoir-beta', V(self._getAngle(unitCell[2], unitCell[0]), C.UNIT_ANGLE)),
            ('reservoir-gamma', V(self._getAngle(unitCell[0], unitCell[1]), C.UNIT_ANGLE)),
            ('reservoir-space-group', V(C.SPACE_GROUPS[self.lattice])),
            ('reservoir-wyckoff-multiplicity-and-letter', V(C.WYCKOFF_CODES[self.lattice])),
            ('reservoir-wyckoff-coordinates', V(C.WYCKOFF_SITES[self.lattice])),
            ('reservoir-wyckoff-species', V([self.elem] * len(C.WYCKOFF_CODES[self.lattice]))),
        ])
        return hostInfo, reservoirInfo

    def run(self):
        # Determine sizes of the supercell
        nBasisAtoms = self.basis.get_number_of_atoms()
        minSize = np.ceil(np.power(C.MIN_ATOMS * 1.0 / nBasisAtoms, 0.333))
        sizes = np.arange(minSize, minSize + 2.0, 1.0)

        # Obtain results for each size
        for i in range(sizes.shape[0]):
            size = sizes[i].astype(int)
            sizeResult = self._getSizeResult(size)

        sys.exit(0)



    def getResult(self):
        # Calculate VME and VFE for Each Size
        sizes = []
        migrationEnergyBySize = []
        formationEnergyBySize = []
        print '\n[Calculation]'
        for size in range(self.cellSizeMin, self.cellSizeMax + 1):
            migrationEnergy, formationEnergy = self._getResultsForSize(size)
            sizes.append(size)
            migrationEnergyBySize.append(migrationEnergy)
            formationEnergyBySize.append(formationEnergy)



        # For Debugging Output
        # sizes = [4, 5, 6]
        # migrationEnergyBySize = [
            # 0.64576855097573116,
            # 0.64448594514897195,
            # 0.64412217147128104,
        # ]
        # formationEnergyBySize = [
            # 0.67914728564926463,
            # 0.67877480646416188,
            # 0.67873178748595819,
        # ]

        print '\n[Calculation Results Summary]'
        print 'Size\tMigrationEnergy\tFormationEnergy'
        for i in range(len(sizes)):
            print [sizes[i], migrationEnergyBySize[i], formationEnergyBySize[i]]

        # Extrapolate for VFE and VME of Infinite Size
        print '\n[Extrapolation]'
        VMEValue, VMEUncert = self._extrapolate(
            sizes,
            migrationEnergyBySize,
            C.FITS_VME_VALUE,
            C.FITS_VME_UNCERT,
            self.FIREUncert * 1.414, # Assuming MDMin Uncertainty Same as FIRE
        )
        VFEValue, VFEUncert = self._extrapolate(
            sizes,
            formationEnergyBySize,
            C.FITS_VFE_VALUE,
            C.FITS_VFE_UNCERT,
            self.FIREUncert,
        )

        # Print Main Results
        print 'Vacancy Migration Energy:', [VMEValue, VMEUncert]
        print 'Vacancy Formation Energy:', [VFEValue, VFEUncert]
        print 'FIRE Uncertainty:', self.FIREUncert

        # Prepare Output
        migrationEnergyResult = OrderedDict([
            ('property-id', C.VME_PROP_ID),
            ('instance-id', 1),
            ('vacancy-migration-energy', V(VMEValue, C.UNIT_ENERGY, VMEUncert)),
            ('host-missing-atom-start', V(1)),
            ('host-missing-atom-end', V(1)),
        ])
        formationEnergyResult = OrderedDict([
            ('property-id', C.VFE_PROP_ID),
            ('instance-id', 2),
            ('relaxed-formation-potential-energy', V(VFEValue, C.UNIT_ENERGY, VFEUncert)),
            ('host-removed-atom', V(1)),
        ])
        hostInfo, reservoirInfo = self._getCrystalInfo()
        migrationEnergyResult.update(hostInfo)
        formationEnergyResult.update(hostInfo)
        formationEnergyResult.update(reservoirInfo)
        return [migrationEnergyResult, formationEnergyResult]
