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
    def __init__(self, elem, model, lattice, latticeConsts):
        # Process Inputs
        self._elem = elem
        self._model = model
        self._lattice = lattice
        self._latticeConsts = np.array(latticeConsts)
        _printInputs()

        self._basis = _createBasis()
        self.FIREUncert = 0
        sys.exit(0)
        
        # Determine size of the supercell
        nAtoms = atoms.get_number_of_atoms()
        factor = math.pow(8 / nAtoms, 0.333)
        self.cellSizeMin = int(math.ceil(factor * C.CELL_SIZE_MIN))
        self.cellSizeMax = self.cellSizeMin + 2
        print 'Cell Size Min:', self.cellSizeMin
        print 'Cell Size Max:', self.cellSizeMax
        print 'Smallest System Size:', nAtoms * self.cellSizeMin**3
        print 'Largest System Size:', nAtoms * self.cellSizeMax**3
        print 'Model Cutoff:', KIM_API_get_data_double(atoms.calc.pkim, 'cutoff')[0]


    def _createBasis():
        # Create basic atoms
        # Cache to local scope
        lattice = self._lattice
        latticeConsts = self._latticeConsts
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
                    [0.0, 0.5 + 0.5 / 3.0 * rt3, 0.5]
                ]),
                cell = np.array([1, rt3, 1]),
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
        atoms.set_calculator(KIMCalculator(self._model))
        if C.OUTPUT_BASIS == True:
            write('basis.exyz', basis, format = 'extxyz')
        return basis

    def _printInputs():
        # Print out the inputs
        print 'Inputs:'
        print 'Element: ', self.elem
        print 'Model: ', self.model
        print 'Lattice: ', self.lattice
        print 'Lattice Constants: ', self.latticeConsts

    def _createSupercell(self, size):
        # Create supercell w/ basic atoms repeated [size] times in each direction
        superAtoms = self.atoms.copy()
        superAtoms.set_calculator(KIMCalculator(self.model))
        superAtoms *= (size, size, size)
        return superAtoms

    def _getResultsForSize(self, size):
        print '========';
        print 'Calculating Size', size, '...'
        atoms = self._createSupercell(size)
        enAtoms = atoms.get_potential_energy()
        nAtoms = atoms.get_number_of_atoms()
        removedAtomPosition = atoms.get_positions()[0]
        del atoms[0]

        # Get Initial State
        initial = atoms.copy()
        initial.set_calculator(KIMCalculator(self.model))

        # Unrelaxed Energy
        enUnrelaxed = initial.get_potential_energy()
        print 'Assert n initial:', initial.get_number_of_atoms()
        unrelaxedFormationEnergy = enUnrelaxed - enAtoms * (nAtoms - 1) / nAtoms

        # Output Before Calculation
        print 'Cohesive Energy:', enAtoms
        print 'Number of Atoms:', nAtoms
        print 'Cohesive Energy Per Atom:', enAtoms / nAtoms
        print 'Unrelaxed Energy w/ Vacancy:', enUnrelaxed
        print 'Unrelaxed Formation Energy:', unrelaxedFormationEnergy

        # Test Moving Epsilon
        for stepLenFactor in range(10):
            stepLen = 0.01 * 2**stepLenFactor # A
            tmpPositions = initial.get_positions()
            tmpPositions[0][0] = tmpPositions[0][0] + stepLen
            initial.set_positions(tmpPositions)
            enMoved = initial.get_potential_energy()
            print 'Moving Energy (', str(stepLen), 'A)', enMoved - enUnrelaxed
            tmpPositions[0][0] = tmpPositions[0][0] - stepLen
            initial.set_positions(tmpPositions)
        print tmpPositions

        # Relaxation
        dyn = FIRE(initial, logfile = C.FIRE_LOG)
        dyn.run(fmax = C.FIRE_TOL, steps = C.FIRE_MAX_STEPS)
        if dyn.get_number_of_steps() >= C.FIRE_MAX_STEPS:
            print '[ERR1047]FIRE Failed: Steps Exceeds Maximum.'
            print 'Vacancy Formation Energy May Not Be Accurate.'

        # Get FIRE Uncertainty
        enTmp = initial.get_potential_energy()
        dyn.run(fmax = C.FIRE_TOL * C.EPS, steps = C.UNCERT_STEPS)
        enUncert = abs(enTmp - initial.get_potential_energy())
        self.FIREUncert = max([enUncert, self.FIREUncert])

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

        # Calculate VME
        # Apply NEB
        images = [initial]
        for i in range(C.NEB_POINTS):
            images.append(initial.copy())
        images.append(final)
        for image in images:
            image.set_calculator(KIMCalculator(self.model))
        neb = NEB(images)
        neb.interpolate()
        minimizer = MDMin(neb)
        minimizer.run(fmax = C.MDMIN_TOL, steps = C.MDMIN_MAX_STEPS)
        if minimizer.get_number_of_steps() >= C.MDMIN_MAX_STEPS:
            print '[ERR1048]NEB Failed: Steps Exceeds Maximum.'
            print 'Vacancy Migration Energy May Not Be Accurate.'

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
        print '\nFit with Size:', xdata
        print 'Orders:', orders
        for order in orders:
            A.append(np.power(xdata * 1.0, -order))
        A = np.vstack(A).T
        print 'Matrix A (Ax = y):\n', A
        print 'Data for Fitting:', ydata
        res = np.linalg.lstsq(A, ydata)
        print 'Fitting Results:', res
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

    def getResults(self):
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
