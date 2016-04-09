# Vacancy Class Definition File For:
# Vacancy Migration Energy (VME) and Formation Energy (VFE) Test Driver
# Author: Junhao Li <streaver91@gmail.com>

# Python Modules
import sys
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
import functions as F

class Vacancy(object):
    # Class for calculating vacancy formation energy and relaxation volume
    def __init__(self, elem, model, lattice, latticeConsts, migration):
        # Inputs description:
        # latticeConsts specified in angstrom
        # migration vector specified in fractional coordinates

        # Process Inputs
        self.elem = elem
        self.model = model
        self.lattice = lattice
        self.latticeConsts = np.array(latticeConsts)
        self.migration = np.array(migration) # Vacancy migration vector
        self._printInputs()

        # Create basis for constructing supercell
        self.basis = self._createBasis()

    def _printInputs(self):
        # Print out the inputs
        print '[Inputs]'
        print 'Element: ', self.elem
        print 'Model: ', self.model
        print 'Lattice: ', self.lattice
        print 'Lattice Constants: ', self.latticeConsts
        print 'Migration Vector:', self.migration

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


    def _relaxAtoms(self, atoms, tol = C.FMAX_TOL, steps = C.FIRE_MAX_STEPS):
        # Relax atoms with FIRE algorithm
        fire = FIRE(atoms, logfile = C.FIRE_LOGFILE)
        fire.run(fmax = tol, steps = steps)
        fireSteps = fire.get_number_of_steps()
        return fireSteps < steps

    def _relaxPath(self, images, tol = C.FMAX_TOL, steps = C.MDMIN_MAX_STEPS):
        # Relax migration path with NEB
        neb = NEB(images)
        neb.interpolate()
        mdmin = MDMin(neb, logfile = C.MDMIN_LOGFILE)
        mdmin.run(fmax = tol, steps = steps)
        mdminSetps = mdmin.get_number_of_steps()
        return mdminSetps < C.MDMIN_MAX_STEPS

    def _findAtomId(self, positions, position):
        # Find the id of positions closest to position
        positionsCopy = positions.copy()
        positionsCopy -= position
        dist = np.sum(positionsCopy**2, axis = 0)
        return np.argmin(dist)

    def _getInitialFinal(self, supercell):
        # Setup image before migration
        initial = supercell.copy()
        del initial[0]

        # Setup image after migration
        final = initial.copy()
        finalPositions = final.get_positions()
        cellSize = np.diagonal(self.basis.get_cell())
        migration = cellSize * np.array(self.migration)
        movedAtomId = self._findAtomId(finalPositions, migration)
        finalPositions[movedAtomId] = np.array([0.0, 0.0, 0.0])
        final.set_positions(finalPositions)

        # Relax both images
        initial.set_calculator(KIMCalculator(self.model))
        tmp = self._relaxAtoms(initial)
        assert tmp, 'Maximum FIRE Steps Reached (initial).'
        final.set_calculator(KIMCalculator(self.model))
        tmp = self._relaxAtoms(final)
        assert tmp, 'Maximum FIRE Steps Reached. (final)'

        return initial, final

    def _getImages(self, initial, final):
        # Create Interpolated images between initial and final
        images = []
        for i in range(C.NEB_POINTS):
            image = initial.copy()
            image.set_calculator(KIMCalculator(self.model))
            images.append(image)
        finalCopy = final.copy()
        finalCopy.set_calculator(KIMCalculator(self.model))
        images.append(finalCopy)
        tmp = self._relaxPath(images)
        assert tmp, 'Maximum MDMin Steps Reached. (migration path)'
        return images

    def _getFmaxUncert(self, image):
        imageCopy = image.copy()
        imageCopy.set_calculator(KIMCalculator(self.model))
        en1 = imageCopy.get_potential_energy()
        self._relaxAtoms(
            imageCopy,
            tol = C.FMAX_TOL * C.EPS,
            steps = C.UNCERT_STEPS
        )
        en1Refined = imageCopy.get_potential_energy()
        return np.abs(en1Refined - en1)

    def _getSaddleImage(self, images):
        # Obtain saddle point energy of vacancy migration
        nImages = len(images)
        xmid = nImages / 2
        x = np.arange(0, nImages, 1.0)
        y = np.array([image.get_potential_energy() for image in images])
        pos = np.array([image.get_positions() for image in images])
        fy = interp1d(x, y, kind = 'cubic') # fy = -y
        fpos = interp1d(x, pos, kind = 'cubic', axis = 0)
        xmax = F.fmax_jh(fy, x[xmid])[0]
        saddleImage = images[0].copy()
        saddleImage.set_positions(fpos(xmax))
        saddleImage.set_calculator(KIMCalculator(self.model))
        return saddleImage

    def _getStressTensor(self, atoms):
        # Obtain stress tensor for orthorhombic atoms
        stress = np.zeros((3, 3))
        atomsCopy = atoms.copy()
        atomsCopy.set_calculator(KIMCalculator(self.model))
        cellCopy = atoms.get_cell().copy()
        dx = C.STRESS_DX
        dcell = np.eye(3) * dx
        areas = atoms.get_volume() / np.diagonal(cellCopy)
        for i in range(3):
            atomsCopy.set_cell(cellCopy + dcell[i], scale_atoms = True)
            eninc = atomsCopy.get_potential_energy()
            atomsCopy.set_cell(cellCopy - dcell[i], scale_atoms = True)
            endec = atomsCopy.get_potential_energy()
            stress[i][i] = (eninc - endec) / (2.0 * dx) / areas[i]
        return stress

    def _getSizeResult(self, size):
        print '[Calculating Supercell of Size', size, ']'

        # Setup supercell
        supercell = self._createSupercell(size)
        en0 = supercell.get_potential_energy()
        nAtoms = supercell.get_number_of_atoms()
        print 'Number of atoms:', nAtoms

        # Initial and final unrelaxed state for migration
        initial, final = self._getInitialFinal(supercell)

        # Interpolate between initial and final
        images = self._getImages(initial, final)

        # Object for output
        res = {}

        # Obtain fmax induced error (md error)
        res['fmax-uncert'] = self._getFmaxUncert(initial)

        # Obtaining vacancy formation energy
        en1 = initial.get_potential_energy()
        res['vacancy-formation-energy'] = en1 - en0 * (nAtoms - 1) / nAtoms

        # Obtain elastic dipole tensor
        nd = 1.0 / supercell.get_volume()  # defect concentration
        stress0 = self._getStressTensor(supercell)
        stress1 = self._getStressTensor(initial)
        res['elastic-dipole-tensor'] = (stress1 - stress0) / nd

        # Obtain vacancy migration energy
        saddleImage = self._getSaddleImage(images)
        en2 = saddleImage.get_potential_energy()
        res['vacancy-migration-energy'] = en2 - en1

        # Obtain saddle point elastic dipole tensor
        stress2 = self._getStressTensor(saddleImage)
        res['saddle-point-elastic-dipole-tensor'] = (stress2 - stress0) / nd

        print res
        return res

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
        sizes = np.arange(minSize, minSize + 3.0, 1.0)

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
        print 'Size    MigrationEnergy    FormationEnergy'
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
