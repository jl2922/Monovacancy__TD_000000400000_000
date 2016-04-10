# Vacancy Class Definition File For:
# Vacancy Migration Energy (VME) and Formation Energy (VFE) Test Driver
# Author: Junhao Li <streaver91@gmail.com>

# Python Modules
import sys
from scipy.optimize import fmin
from scipy.interpolate import interp1d
import numpy as np
import json

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
import ase.units as units

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
        # print 'cell:', cellCopy
        # print 'volume:', atoms.get_volume()
        areas = atoms.get_volume() / np.diagonal(cellCopy)
        print atoms
        # print areas
        for i in range(3):
            dcell = np.zeros((3, 3))
            dcell[i][i] = dx
            atomsCopy.set_cell(cellCopy - 2 * dcell, scale_atoms = True)
            enM2 = atomsCopy.get_potential_energy()
            atomsCopy.set_cell(cellCopy - dcell, scale_atoms = True)
            enM1 = atomsCopy.get_potential_energy()
            atomsCopy.set_cell(cellCopy + 2 * dcell, scale_atoms = True)
            enP2 = atomsCopy.get_potential_energy()
            atomsCopy.set_cell(cellCopy + dcell, scale_atoms = True)
            enP1 = atomsCopy.get_potential_energy()
            derivative = (-enP2 + 8 * enP1 - 8 * enM1 + enM2) / (12 * dx)
            stress[i][i] = derivative / areas[i]
        return stress

    def _getElasticStiffness(self, atoms):
        # Obtain the top left 3*3 elastic stiffness tensor
        atomsCopy = atoms.copy()
        # atomsCopy.set_cell(atomsCopy.get_cell() * 1.0, scale_atoms = True)
        cellCopy = atomsCopy.get_cell().copy()
        atomsCopy.set_calculator(KIMCalculator(self.model))

        # Obtain elastic stiffness tensor
        EST = np.zeros((3, 3))
        dx = C.STRESS_DX
        dStrain = dx / np.diagonal(cellCopy)
        print 'atoms:', atomsCopy
        for i in range(3):
            dcell = np.zeros((3, 3))
            dcell[i][i] = dx
            atomsCopy.set_cell(cellCopy + dcell, scale_atoms = True)
            stressPlus = self._getStressTensor(atomsCopy)
            atomsCopy.set_cell(cellCopy - dcell, scale_atoms = True)
            stressMinus = self._getStressTensor(atomsCopy)
            print 'stress:', stressPlus
            print stressMinus
            dStress = np.diagonal(stressPlus - stressMinus) / 2.0
            print dStress
            print dStrain
            EST[:, i] = dStress / dStrain
        print 'est:', EST
        # print EST / units.GPa
        return EST

    def _getElasticCompliance(self, atoms):
        # Obtain the top left 3*3 elastic compliance tensor
        EST = self._getElasticStiffness(atoms)
        # print EST
        ECT = np.linalg.inv(EST)
        # print ECT
        return ECT

    def _getDefectStrainTensor(self, EDT, ECT):
        # Compute defect strain tensor from:
        # EDT: elastic dipole tensor
        # ECT: elastic compliance tensor
        DST = np.zeros((3, 3))
        EDTDiag = np.diagonal(EDT)
        for i in range(3):
            DST[i][i] = np.sum(np.diagonal(EDT) * ECT[i])
        return DST

    def _getSizeResult(self, size):
        print '[Calculating Supercell of Size', size, ']'
        F.clock('start')

        # Setup supercell
        supercell = self._createSupercell(size)
        en0 = supercell.get_potential_energy()
        nAtoms = supercell.get_number_of_atoms()
        print 'Number of atoms:', nAtoms

        # Initial and final unrelaxed state for migration
        initial, final = self._getInitialFinal(supercell)
        F.clock('initial and final obtained')

        # Interpolate between initial and final
        images = self._getImages(initial, final)
        F.clock('migration path obtained')

        # Object for output
        res = {}

        # Obtain fmax induced error (md error)
        res['fmax-uncert'] = self._getFmaxUncert(initial)
        F.clock('fmax uncert')

        # sys.exit(0)

        # Obtaining vacancy formation energy
        en1 = initial.get_potential_energy()
        VFE = en1 - en0 * (nAtoms - 1) / nAtoms
        res['vacancy-formation-energy'] = VFE

        # Obtain dipole strain / stress tensor and relxation volume
        nd = 1.0 / supercell.get_volume()  # defect concentration
        stress0 = self._getStressTensor(supercell)
        stress1 = self._getStressTensor(initial)
        EDT = (stress1 - stress0) / nd
        DST = self._getDefectStrainTensor(EDT, self.ECT)
        RV = np.trace(DST)
        res['elastic-dipole-tensor'] = EDT
        res['defect-strain-tensor'] = DST
        res['vacancy-relaxation-volume'] = RV

        # Obtain vacancy migration energy
        saddleImage = self._getSaddleImage(images)
        en2 = saddleImage.get_potential_energy()
        VME = en2 - en1
        res['vacancy-migration-energy'] = VME

        # Obtain saddle point elastic dipole tensor
        stress2 = self._getStressTensor(saddleImage)
        SPEDT = (stress2 - stress0) / nd
        res['saddle-point-elastic-dipole-tensor'] = SPEDT

        F.clock('finish')
        print res
        return res

    def _getFit(self, xdata, ydata, orders):
        # Polynomial Fitting with Specific Orders
        # Deal with higher order ydate
        if len(ydata.shape) > 2:
            ydataShape = ydata.shape
            ydataNew = ydata.reshape((ydataShape[0], -1))
            res = self._getFit(xdata, ydataNew, orders)
            print ydataShape[1:]
            res = res.reshape((len(orders), ) + ydataShape[1:])
            return res

        # Return fitted results Corresponding to the orders as np array
        A = []
        print '[Fitting]'
        print 'xdata:', xdata
        print 'ydata:', ydata
        print 'Orders:', orders
        for order in orders:
            A.append(np.power(xdata * 1.0, order))
        A = np.vstack(A).T
        res = np.linalg.lstsq(A, ydata)
        res = res[0]
        print 'Results:', res
        return res

    def _extrapolate(self, sizes, sizeResults):
        properties = sizeResults[0].keys()
        propValues = {}
        nSizes = sizes.shape[0]
        for prop in properties:
            tmp = np.array([sizeResults[i][prop] for i in range(nSizes)])
            propValues[prop] = tmp

        print propValues

        def _checkProp(prop):
            assert prop in properties, prop + ' not found'

        def _extProp(prop, orders, nPoints):
            # Always use the last few points
            assert nPoints <= nSizes, 'not enough points for extrapolation'
            sizesNew = sizes[-nPoints:]
            valuesNew = propValues[prop][-nPoints:]
            fitRes = self._getFit(sizesNew, valuesNew, orders)
            return fitRes[0]

        def _format(value, unit = '', uncert = ''):
            if type(value) == np.ndarray:
                value = value.tolist()
            res = {
                'source-value': value
            }
            if unit != '':
                res['source-unit'] = unit
            if uncert != '':
                if type(uncert) == np.ndarray:
                    uncert = uncert.tolist()
                res['source-std-uncert-value'] = uncert
            return res

        def _getValueUncert(prop, orders, nPts, nPtsUncert, sysUncert):
            print 'extrapolating', prop
            _checkProp(prop)
            val = _extProp(prop, [0, -3], nPts)
            val2 = _extProp(prop, [0, -3], nPtsUncert)
            uncert = np.sqrt(np.abs(val - val2)**2 + sysUncert**2)
            return val, uncert

        _checkProp('fmax-uncert')
        fmaxUncert = np.max(propValues['fmax-uncert'])

        res = {} # Object for return

        # Extrapolate vacancy formation energy
        VFE, VFEUncert = _getValueUncert(
            'vacancy-formation-energy',
            [0, -3], 3, 2, fmaxUncert
        )
        res['vacancy-formation-energy'] = _format(VFE, 'eV', VFEUncert)

        # Extrapolate vacancy migration energy
        VME, VMEUncert = _getValueUncert(
            'vacancy-migration-energy',
            [0, -3], 3, 2, fmaxUncert
        )
        res['vacancy-migration-energy'] = _format(VME, 'eV', VMEUncert)

        # Extrapolate elastic dipole tensors
        EDT, EDTUncert = _getValueUncert(
            'elastic-dipole-tensor',
            [0, -3], 3, 2, fmaxUncert
        )
        res['elastic-dipole-tensor'] = _format(EDT, 'eV', EDTUncert)
        SPEDT, SPEDTUncert = _getValueUncert(
            'saddle-point-elastic-dipole-tensor',
            [0, -3], 3, 2, fmaxUncert
        )
        res['saddle-point-elastic-dipole-tensor'] = _format(SPEDT, 'eV', SPEDTUncert)

        # Extrapolate defect strain tensor
        DST, DSTUncert = _getValueUncert(
            'defect-strain-tensor',
            [0, -3], 3, 2, 0.0
        )
        DSTUncert = np.sqrt(DSTUncert**2 + (fmaxUncert / VFE * DST)**2)
        res['defect-strain-tensor'] = _format(DST, 'angstrom^3', DSTUncert)

        # Extrapolate vacancy relxation volume
        RV = np.array(DST).trace()
        RVUncert = np.sqrt(np.sum(np.diagonal(np.array(DSTUncert))**2))
        res['vacancy-relxation-volume'] = _format(RV, 'angstrom^3', RVUncert)

        resStr = json.dumps(res, separators = (' ',' '), indent = 2)
        print resStr

    def run(self):
        # Determine sizes of the supercell
        nBasisAtoms = self.basis.get_number_of_atoms()
        minSize = np.ceil(np.power(C.MIN_ATOMS * 1.0 / nBasisAtoms, 0.333))
        sizes = np.arange(minSize, minSize + C.NUM_SIZES, 1.0)

        # Obtain common variables
        supercell = self._createSupercell(minSize.astype(int))
        self.ECT = self._getElasticCompliance(self.basis)
        print self.ECT
        
        # Obtain results for each size
        sizeResults = []
        for i in range(sizes.shape[0]):
            size = sizes[i].astype(int)
            sizeResults.append(self._getSizeResult(size))

        self.res = self._extrapolate(sizes, sizeResults)


    def getResult(self):
        return self.res
