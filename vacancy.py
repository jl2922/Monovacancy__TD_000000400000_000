# Vacancy Class Definition File For:
# Vacancy Migration Energy (VME) and Formation Energy (VFE) Test Driver
# Author: Junhao Li <streaver91@gmail.com>

# Python Modules
import sys
from scipy.optimize import fmin
from scipy.interpolate import interp1d
import numpy as np
import json
import tarfile

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

        # For storing cif file paths to be compressed
        self._cifs = []

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

        res = {}
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
            res['a'] = _format(a, 'angstrom')
            res['b'] = _format(a, 'angstrom')
            res['c'] = _format(c, 'angstrom')
            res['alpha'] = _format(90.0, 'degree')
            res['beta'] = _format(90.0, 'degree')
            res['gamma'] = _format(120.0, 'degree')
            basisCoords = [[0.0, 0.0, 0.0], [2.0 / 3, 1.0 / 3, 0.5]]
            res['basis-atom-coordinates'] = _format(basisCoords)
            res['basis-atom-species'] = _format([elem] * 2)
        else:
            a = latticeConsts[0]
            basis = bulk(
                elem,
                a = a,
                crystalstructure = lattice,
                cubic = True,
            )
            res['a'] = _format(a, 'angstrom')
            res['b'] = _format(a, 'angstrom')
            res['c'] = _format(a, 'angstrom')
            res['alpha'] = _format(90.0, 'degree')
            res['beta'] = _format(90.0, 'degree')
            res['gamma'] = _format(90.0, 'degree')
            basisCoords = basis.get_positions() / a
            basisNAtoms = basis.get_number_of_atoms()
            res['basis-atom-coordinates'] = _format(basisCoords)
            res['basis-atom-species'] = _format([elem] * basisNAtoms)

        res['basis-short-name'] = _format([lattice])
        res['space-group'] = _format(C.SPACE_GROUPS[lattice])
        res['wyckoff-multiplicity-and-letter'] = _format(C.WYCKOFF_CODES[lattice])
        res['wyckoff-coordinates'] = _format(C.WYCKOFF_SITES[lattice])
        res['wyckoff-species'] = _format([elem])

        self._res = res

        if C.SAVE_BASIS == True:
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
        # Relax migration path with nudged elastic band
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
        nAtoms = supercell.get_number_of_atoms()
        # Setup image before migration
        initial = supercell.copy()
        del initial[0]
        self._saveAtoms('unrelaxedInitial%d' % nAtoms, initial)

        # Setup image after migration
        final = initial.copy()
        finalPositions = final.get_positions()
        cellSize = np.diagonal(self.basis.get_cell())
        migration = cellSize * np.array(self.migration)
        movedAtomId = self._findAtomId(finalPositions, migration)
        finalPositions[movedAtomId] = np.array([0.0, 0.0, 0.0])
        final.set_positions(finalPositions)
        self._saveAtoms('unrelaxedFinal%d' % nAtoms, initial)

        # Relax both images
        initial.set_calculator(KIMCalculator(self.model))
        tmp = self._relaxAtoms(initial)
        assert tmp, 'Maximum FIRE Steps Reached (initial).'
        self._saveAtoms('relaxedInitial%d' % nAtoms, initial)
        final.set_calculator(KIMCalculator(self.model))
        tmp = self._relaxAtoms(final)
        assert tmp, 'Maximum FIRE Steps Reached. (final)'
        self._saveAtoms('relaxedFinal%d' % nAtoms, initial)

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
        # Obtain fmax induced uncertainty
        # Relax images for a few more steps
        imageCopy = image.copy()
        imageCopy.set_calculator(KIMCalculator(self.model))
        enOrigin = imageCopy.get_potential_energy()
        stressOrigin = self._getStressTensor(imageCopy)
        self._relaxAtoms(
            imageCopy,
            tol = C.FMAX_TOL * C.EPS,
            steps = C.UNCERT_STEPS
        )

        # Obtain energy uncertainty
        en1Refined = imageCopy.get_potential_energy()
        enUncert = np.abs(en1Refined - enOrigin)

        # Obtain stress uncertainty
        stressRefined = self._getStressTensor(imageCopy)
        stressUncert = np.abs(stressRefined - stressOrigin)
        edtUncert = stressUncert * imageCopy.get_volume()

        return enUncert, edtUncert


    def _getSaddleImage(self, images):
        # Obtain saddle point energy of vacancy migration
        nAtoms = images[0].get_number_of_atoms() + 1
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
        self._saveAtoms('relaxedSaddle%d' % nAtoms, saddleImage)
        return saddleImage

    def _getStressTensor(self, atoms):
        # Obtain stress tensor for orthorhombic atoms
        cellCopy = atoms.get_cell().copy()
        areas = atoms.get_volume() / np.diagonal(cellCopy)
        stress = np.zeros((3, 3))

        def energyToCell(cell):
            # Potential energy partial cell size
            # Warning: will cause atoms to change
            atoms.set_cell(cell, scale_atoms = True)
            return atoms.get_potential_energy()

        # Loop through each direction
        for i in range(3):
            dEnergy = F.partialD(energyToCell, cellCopy, (i, i))
            stress[i][i] = dEnergy / areas[i]

        # Restore to original cell
        atoms.set_cell(cellCopy, scale_atoms = True)

        return stress

    def _getElasticStiffness(self, atoms):
        # Obtain the top left 3*3 elastic stiffness tensor
        cellCopy = atoms.get_cell().copy()
        cellDiag = np.diagonal(cellCopy)

        # Obtain elastic stiffness tensor
        EST = np.zeros((3, 3))

        def stressToCell(cell):
            # Partial stress partial cell size
            # Warning: will cause atoms to change
            atoms.set_cell(cell, scale_atoms = True)
            return self._getStressTensor(atoms)

        # Loop through strain in each direction
        for i in range(3):
            dStress = F.partialD(stressToCell, cellCopy, (i, i))
            dStress = np.diagonal(dStress)
            EST[:, i] = dStress * cellDiag

        # Restore to original cell
        atoms.set_cell(cellCopy, scale_atoms = True)

        return EST

    def _getElasticCompliance(self, atoms):
        # Obtain the top left 3*3 elastic compliance tensor
        EST = self._getElasticStiffness(atoms)
        ECT = np.linalg.inv(EST)
        return ECT

    def _getDefectStrainTensor(self, EDT, ECT):
        # Compute defect strain tensor from:
        # EDT: elastic dipole tensor
        # ECT: elastic compliance tensor
        DST = np.zeros((3, 3))
        EDTDiag = np.diagonal(EDT)
        for i in range(3):
            DST[i][i] = np.sum(EDTDiag * ECT[i, :])
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
        fmaxEnergyUncert, fmaxEDTUncert = self._getFmaxUncert(initial)
        res['fmax-energy-uncert'] = fmaxEnergyUncert
        res['fmax-edt-uncert'] = fmaxEDTUncert
        F.clock('fmax uncert obtained')

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

        # Print results for size
        F.printDict(res)

        F.clock('size finished')
        return res

    def _getFit(self, xdata, ydata, orders):
        # Polynomial Fitting with Specific Orders
        # Deal with higher order ydate
        if len(ydata.shape) > 2:
            ydataShape = ydata.shape
            ydataNew = ydata.reshape((ydataShape[0], -1))
            res = self._getFit(xdata, ydataNew, orders)
            res = res.reshape((len(orders), ) + ydataShape[1:])
            return res

        # Return fitted results Corresponding to the orders as np array
        A = []
        # print '[Fitting]'
        # print 'xdata:', xdata
        # print 'ydata:', ydata
        # print 'Orders:', orders
        for order in orders:
            A.append(np.power(xdata * 1.0, order))
        A = np.vstack(A).T
        res = np.linalg.lstsq(A, ydata)
        res = res[0]
        # print 'Results:', res
        return res

    def _extrapolate(self, sizes, sizeResults):
        properties = sizeResults[0].keys()
        propValues = {}
        nSizes = sizes.shape[0]
        for prop in properties:
            tmp = np.array([sizeResults[i][prop] for i in range(nSizes)])
            propValues[prop] = tmp

        F.printDict(propValues)

        def _checkProp(prop):
            assert prop in properties, prop + ' not found'

        def _extProp(prop, orders, nPoints):
            # Always use the last few points
            assert nPoints <= nSizes, 'not enough points for extrapolation'
            sizesNew = sizes[-nPoints:]
            valuesNew = propValues[prop][-nPoints:]
            fitRes = self._getFit(sizesNew, valuesNew, orders)
            return fitRes[0]

        def _getValueUncert(prop, orders, nPts, nPtsUncert, sysUncert):
            F.clock('extrapolating ' + prop)
            _checkProp(prop)
            val = _extProp(prop, orders, nPts)
            val2 = _extProp(prop, orders, nPtsUncert)
            uncert = np.sqrt(np.abs(val - val2)**2 + sysUncert**2)
            return val, uncert

        _checkProp('fmax-energy-uncert')
        fmaxEnergyUncert = np.max(propValues['fmax-energy-uncert'][-3:])
        fmaxEDTUncert = np.max(propValues['fmax-edt-uncert'][-3:], axis = 0)

        res = self._res # Object for return

        # Extrapolate vacancy formation energy
        VFE, VFEUncert = _getValueUncert(
            'vacancy-formation-energy',
            [0, -3], 3, 2, fmaxEnergyUncert
        )
        res['vacancy-formation-energy'] = _format(VFE, 'eV', VFEUncert)

        # Extrapolate vacancy migration energy
        VME, VMEUncert = _getValueUncert(
            'vacancy-migration-energy',
            [0, -3], 3, 2, fmaxEnergyUncert
        )
        res['vacancy-migration-energy'] = _format(VME, 'eV', VMEUncert)

        # Extrapolate elastic dipole tensors
        EDT, EDTUncert = _getValueUncert(
            'elastic-dipole-tensor',
            [0, -3], 3, 2, fmaxEDTUncert
        )
        res['elastic-dipole-tensor'] = _format(EDT, 'eV', EDTUncert)
        SPEDT, SPEDTUncert = _getValueUncert(
            'saddle-point-elastic-dipole-tensor',
            [0, -3], 3, 2, fmaxEDTUncert
        )
        res['saddle-point-elastic-dipole-tensor'] = _format(SPEDT, 'eV', SPEDTUncert)

        # Extrapolate defect strain tensor
        DST, DSTUncert = _getValueUncert(
            'defect-strain-tensor',
            [0, -3], 3, 2, 0.0
        )
        DSTSysUncert = fmaxEDTUncert / (EDT + C.EPS) * DST
        DSTUncert = np.sqrt(DSTUncert**2 + DSTSysUncert**2)
        res['defect-strain-tensor'] = _format(DST, 'angstrom^3', DSTUncert)

        # Extrapolate vacancy relxation volume
        RV = np.array(DST).trace()
        RVUncert = np.sqrt(np.sum(np.diagonal(np.array(DSTUncert))**2))
        res['vacancy-relxation-volume'] = _format(RV, 'angstrom^3', RVUncert)

        F.printDict(res)

    def run(self):
        # Determine sizes of the supercell
        nBasisAtoms = self.basis.get_number_of_atoms()
        minSize = np.ceil(np.power(C.MIN_ATOMS * 1.0 / nBasisAtoms, 0.333))
        sizes = np.arange(minSize, minSize + C.NUM_SIZES, 1.0)

        # Obtain common variables
        supercell = self._createSupercell(minSize.astype(int))
        self.ECT = self._getElasticCompliance(supercell)
        print self.ECT

        # Obtain results for each size
        sizeResults = []
        for i in range(sizes.shape[0]):
            size = sizes[i].astype(int)
            sizeResults.append(self._getSizeResult(size))

        self.res = self._extrapolate(sizes, sizeResults)

    def _saveAtoms(self, filename, atoms):
        path = 'output/%s.cif' % filename
        self._cifs.append(path)
        write(path, atoms, format = 'cif')

    def _addInfo(self):
        res = self._res
        res['missing-atom-wyckoff-site-start'] = _format(0)
        res['missing-atom-wyckoff-site-end'] = _format(0)
        res['vacancy-short-name-start'] = _format('tetrahedral')
        res['vacancy-short-name-end'] = _format('tetrahedral')
        res['vacancy-migration-direction'] = _format(self.migration)
        res['missing-atom-species'] = _format(self.elem)
        res['vacancy-position-start'] = _format(np.zeros(3))
        res['vacancy-position-end'] = _format(self.migration)
        res['vacancy-short-name-start'] = _format('')
        res['cauchy-stress'] = _format(np.zeros(6), 'GPa')
        res['temperature'] = _format(0.0, 'K')

        # Supercell structure
        tar = tarfile.open("output/structure.tgz", "w:gz")
        for cif in self._cifs:
            tar.add(cif)
        tar.close()
        structuralProps = [
            'vacancy-unrelaxed-structure-start',
            'vacancy-unrelaxed-structure-end',
            'vacancy-relxed-structure-start',
            'vacancy-relaxed-structure-end',
            'vacancy-relaxed-structure-saddle-point',
        ]
        for prop in structuralProps:
            res[prop] = _format('structure.tgz')

    def getResult(self):
        self._addInfo()
        return self._res
