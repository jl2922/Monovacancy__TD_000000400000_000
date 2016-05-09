# Vacancy Class Definition File for Vacancy Test Driver
# Author: Junhao Li <streaver91@gmail.com>

# Python Modules
import os
import sys
from scipy import interpolate
import numpy as np
import tarfile

# ASE Modules
try:
    from ase.lattice import bulk
    print 'Imported bulk from ase.lattice' # For ASE version 3.9
except ImportError:
    from ase.structure import bulk
    print 'Imported bulk from ase.structure' # For ASE version <= 3.8
from ase.optimize import FIRE, MDMin
from ase import Atoms
from ase.io import write
from ase.neb import NEB

# KIM Modules
from kimcalculator import KIMCalculator
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
        self.strain = np.array(C.STRAIN)

        # For storing results object
        self._res = {}

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

        res = self._res
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
            basis.set_cell(basis.get_cell() * (1 + self.strain), scale_atoms = True)
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
            print basis.get_cell()
            basis.set_cell(basis.get_cell() * (1 + self.strain), scale_atoms = True)
            print basis.get_cell()
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
        print 'fire finished in %d steps' % fireSteps
        return fireSteps < steps

    def _relaxPath(self, images, tol = C.FMAX_TOL, steps = C.MDMIN_MAX_STEPS):
        # Relax migration path with nudged elastic band
        neb = NEB(images)
        neb.interpolate()
        mdmin = MDMin(neb, logfile = C.MDMIN_LOGFILE)
        mdmin.run(fmax = tol, steps = steps)
        mdminSetps = mdmin.get_number_of_steps()
        print 'mdmin finished in %d steps' % mdminSetps
        return mdminSetps < C.MDMIN_MAX_STEPS

    def _findAtomId(self, positions, position):
        # Find the id of positions closest to position
        positionsCopy = positions.copy()
        positionsCopy -= position
        dist = np.sum(positionsCopy**2, axis = 1)
        return np.argmin(dist)

    def _getInitialFinal(self, supercell):
        nAtoms = supercell.get_number_of_atoms()
        # Setup unrelaxed image before migration
        initial = supercell.copy()
        del initial[0]
        self._saveAtoms('unrelaxedInitial%d' % nAtoms, initial)

        # Setup unrelaxed image after migration
        final = initial.copy()
        finalPositions = final.get_positions()
        cellSize = np.diagonal(self.basis.get_cell())
        migration = cellSize * np.array(self.migration)
        movedAtomId = self._findAtomId(finalPositions, migration)
        print 'atom #%d moved from ' % movedAtomId, finalPositions[movedAtomId]
        finalPositions[movedAtomId] = np.array([0.0, 0.0, 0.0])
        final.set_positions(finalPositions)
        self._saveAtoms('unrelaxedFinal%d' % nAtoms, final)

        # Relax both images
        initial.set_calculator(KIMCalculator(self.model))
        tmp = self._relaxAtoms(initial)
        assert tmp, 'Maximum FIRE Steps Reached (initial).'
        self._saveAtoms('relaxedInitial%d' % nAtoms, initial)
        final.set_calculator(KIMCalculator(self.model))
        tmp = self._relaxAtoms(final)
        assert tmp, 'Maximum FIRE Steps Reached. (final)'
        self._saveAtoms('relaxedFinal%d' % nAtoms, final)

        return initial, final

    def _getImages(self, initial, final):
        # Create and relax Interpolated images between initial and final
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

    def _getSaddleImage(self, images, saveImage = True):
        # Obtain saddle point energy of vacancy migration
        nAtoms = images[0].get_number_of_atoms() + 1
        nImages = len(images)
        xmid = nImages / 2
        x = np.arange(0, nImages, 1.0)
        y = np.array([image.get_potential_energy() for image in images])
        pos = np.array([image.get_positions() for image in images])
        fy = interpolate.interp1d(x, y, kind = 'cubic')
        fpos = interpolate.interp1d(x, pos, kind = 'cubic', axis = 0)
        xmax = F.fmax_jh(fy, x[xmid], ftol = 1.0e-10)[0]
        saddleImage = images[0].copy()
        saddleImage.set_positions(fpos(xmax))
        saddleImage.set_calculator(KIMCalculator(self.model))
        if saveImage:
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

    def _getDefectStrainTensor(self, EDT):
        # Compute defect strain tensor from:
        # EDT: elastic dipole tensor
        # ECT: elastic compliance tensor
        DST = np.zeros((3, 3))
        EDTDiag = np.diagonal(EDT)
        for i in range(3):
            DST[i][i] = np.sum(EDTDiag * self.ECT[i, :])
        return DST

    def _getSizeResult(self, size):
        # Procedure:
        # 1. Create supercell
        # 2. Obtain initial and final configuration
        # 3. Obtain path with NEB
        # 4. Record values needed, v
        # 5. Run 2 and 3 with more precision for a few more steps
        # 6. Record refined values, vRefined, for res
        # 6. Use abs(v - vRefined) for res uncertainty

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

        # Record values needed
        en1 = initial.get_potential_energy()
        VFE = en1 - en0 * (nAtoms - 1) / nAtoms  # vacancy formation energy
        nd = 1.0 / supercell.get_volume()  # defect concentration
        stress0 = self._getStressTensor(supercell)
        print 'STRESS:'
        print stress0
        stress1 = self._getStressTensor(initial)
        print stress1
        EDT = (stress1 - stress0) / nd  # elastic dipole tensor
        DST = self._getDefectStrainTensor(EDT)  # defect strain tensor
        VRV = np.trace(DST)  # vacancy relaxation volume
        saddleImage = self._getSaddleImage(images)
        en2 = saddleImage.get_potential_energy()
        VME = en2 - en1  # vacancy migration energy
        stress2 = self._getStressTensor(saddleImage)
        SPEDT = (stress2 - stress0) / nd  # saddle point elastic dipole tensor
        SPDST = self._getDefectStrainTensor(SPEDT)

        # Refine results with high precision for a few more steps
        initialRefined = initial.copy()
        initialRefined.set_calculator(KIMCalculator(self.model))
        self._relaxAtoms(initialRefined, tol = C.FMAX_TOL * C.EPS, steps = 5)
        self._relaxPath(images, tol = C.FMAX_TOL * C.EPS, steps = 5)
        F.clock('results refined')

        # Record refined values
        en1Refined = initialRefined.get_potential_energy()
        VFERefined = en1Refined - en0 * (nAtoms - 1) / nAtoms
        stress1Refined = self._getStressTensor(initialRefined)
        EDTRefined = (stress1Refined - stress0) / nd
        DSTRefined = self._getDefectStrainTensor(EDTRefined)
        VRVRefined = np.trace(DSTRefined)
        saddleImageRefined = self._getSaddleImage(images, False)
        en2Refined = saddleImageRefined.get_potential_energy()
        VMERefined = en2Refined - en1Refined
        stress2Refined = self._getStressTensor(saddleImageRefined)
        SPEDTRefined = (stress2Refined - stress0) / nd
        SPDSTRefined = self._getDefectStrainTensor(SPEDTRefined)

        # Object for output
        sizeRes = {
            'number-of-atoms': nAtoms
        }

        # Save refined results as final results
        # Save difference as uncertainty
        def saveSizeRes(prop, valueInit, valueRefined, sysUncert = 0.0):
            sizeRes[prop] = valueInit
            propUncert = np.abs(valueRefined - valueInit)
            propUncert = np.sqrt(propUncert**2 + sysUncert**2)
            sizeRes[prop + '-uncert'] = propUncert

        saveSizeRes('vacancy-formation-energy', VFE, VFERefined)
        saveSizeRes('elastic-dipole-tensor', EDT, EDTRefined)
        saveSizeRes('defect-strain-tensor', DST, DSTRefined)
        saveSizeRes('vacancy-relaxation-volume', VRV, VRVRefined)
        VMESysUncert = np.abs(VFERefined - VFE)
        SPEDTSysUncert = np.abs(EDTRefined - EDT)
        SPDSTSysUncert = np.abs(DSTRefined - DST)
        saveSizeRes('vacancy-migration-energy', VME, VMERefined, VMESysUncert)
        saveSizeRes('saddle-point-elastic-dipole-tensor', SPEDT, SPEDTRefined, SPEDTSysUncert)
        saveSizeRes('saddle-point-defect-strain-tensor', SPDST, SPDSTRefined, SPDSTSysUncert)

        # Print results for size
        F.printDict(sizeRes)

        F.clock('size finished')
        return sizeRes

    def _getFit(self, xdata, ydata, orders):
        # Polynomial Fitting with Specific Orders
        # Deal with higher order ydate
        if len(ydata.shape) > 2:
            ydataShape = ydata.shape
            ydataNew = ydata.reshape((ydataShape[0], -1))
            fitRes = self._getFit(xdata, ydataNew, orders)
            fitRes = fitRes.reshape((len(orders), ) + ydataShape[1:])
            return fitRes

        # Return fitted results Corresponding to the orders as np array
        A = []
        for order in orders:
            A.append(np.power(xdata * 1.0, order))
        A = np.vstack(A).T
        fitRes = np.linalg.lstsq(A, ydata)
        fitRes = fitRes[0]
        return fitRes

    def _extrapolate(self, sizes, sizeResults):
        # Extrapolate results for each size and obtain dilute limit
        # First pack sizeResults into propValues
        properties = sizeResults[0].keys()
        propValues = {}
        nSizes = sizes.shape[0]
        for prop in properties:
            tmp = np.array([sizeResults[i][prop] for i in range(nSizes)])
            propValues[prop] = tmp
        print 'size calculation summary:'
        F.printDict(propValues)

        def extProp(prop, orders, nPoints):
            # Always use the last few points
            assert nPoints <= nSizes, 'not enough points for extrapolation'
            sizesNew = sizes[-nPoints:]
            valuesNew = propValues[prop][-nPoints:]
            fitRes = self._getFit(sizesNew, valuesNew, orders)
            return fitRes[0]

        def getValueUncert(prop, orders, sysUncert = 0.0):
            # Obtain value and uncertainty
            # Uncertainty consist of three parts:
            # 1. statistical uncertainty from extrapolation
            # 2. uncertainty of the original value
            # 2. additional systematic uncertainty
            print 'extrapolating ' + prop
            assert prop in properties, prop + ' not found'

            # Obtain number of points for extrapolation
            nPts1 = len(orders)
            nPts2 = nPts1 + 1

            # Obtain statistical uncertainty
            val1 = extProp(prop, orders, nPts1)
            val2 = extProp(prop, orders, nPts2)
            statUncert = np.abs(val1 - val2)

            # Obtain value uncertainty
            valUncert = 0.0
            valUncertKey = prop + '-uncert'
            if valUncertKey in propValues:
                valUncert = propValues[valUncertKey]
                valUncert = np.max(valUncert[-nPts2:], axis = 0)
            # Final results
            val = (val1 + val2) / 2
            uncert = np.sqrt(statUncert**2 + valUncert**2 + sysUncert**2)
            return val, uncert

        res = self._res # Object for storing results

        # Extrapolate vacancy formation energy
        VFE, VFEUncert = getValueUncert('vacancy-formation-energy', [0, -3])
        res['vacancy-formation-energy'] = _format(VFE, 'eV', VFEUncert)

        # Extrapolate vacancy migration energy
        VME, VMEUncert = getValueUncert('vacancy-migration-energy', [0, -3])
        res['vacancy-migration-energy'] = _format(VME, 'eV', VMEUncert)

        # Extrapolate elastic dipole tensors
        EDT, EDTUncert = getValueUncert('elastic-dipole-tensor', [0, -3])
        res['elastic-dipole-tensor'] = _format(EDT, 'eV', EDTUncert)
        SPEDT, SPEDTUncert = getValueUncert('saddle-point-elastic-dipole-tensor', [0, -3])
        res['saddle-point-elastic-dipole-tensor'] = _format(SPEDT, 'eV', SPEDTUncert)

        # Extrapolate defect strain tensor
        DST, DSTUncert = getValueUncert('defect-strain-tensor', [0, -3])
        res['defect-strain-tensor'] = _format(DST, 'angstrom^3', DSTUncert)
        SPDST, SPDSTUncert = getValueUncert('saddle-point-defect-strain-tensor', [0, -3])
        res['saddle-point-defect-strain-tensor'] = _format(SPDST, 'angstrom^3', SPDSTUncert)

        # Extrapolate vacancy relxation volume
        VRV, VRVUncert = getValueUncert('vacancy-relaxation-volume', [0, -3])
        res['vacancy-relaxation-volume'] = _format(VRV, 'angstrom^3', VRVUncert)

        print 'extrapolated results:'
        F.printDict(res)
        F.clock('extrapolation finished')

    def run(self):
        # Determine sizes of the supercell
        nBasisAtoms = self.basis.get_number_of_atoms()
        minSize = np.ceil(np.power(C.MIN_ATOMS * 1.0 / nBasisAtoms, 0.333))
        sizes = np.arange(minSize, minSize + C.NUM_SIZES, 1.0)

        # Obtain common variables
        supercell = self._createSupercell(minSize.astype(int))
        EST = self._getElasticStiffness(supercell)
        ECT = np.linalg.inv(EST)
        self.ECT = ECT

        # Add EST and ECT to res
        res = self._res
        res['partial-elastic-stiffness-tensor'] = _format(EST)
        res['partial-elastic-compliance-tensor'] = _format(ECT)
        enCoh = supercell.get_potential_energy() / supercell.get_number_of_atoms()
        res['free-energy-per-atom'] = _format(enCoh, 'eV')

        # Obtain results for each size
        sizeResults = []
        for i in range(sizes.shape[0]):
            size = sizes[i].astype(int)
            sizeResults.append(self._getSizeResult(size))

        self._extrapolate(sizes, sizeResults)

    def _saveAtoms(self, filename, atoms):
        path = 'output/%s.cif' % filename
        self._cifs.append(path)
        write(path, atoms, format = 'cif')

    def _packStructure(self):
        tar = tarfile.open("output/structure.tgz", "w:gz")
        for cif in self._cifs:
            tar.add(cif)
            os.remove(cif)
        tar.close()
        structuralProps = [
            'vacancy-unrelaxed-structure-start',
            'vacancy-unrelaxed-structure-end',
            # 'vacancy-unrelaxed-structure-saddle-point', # maynot be accurate
            'vacancy-relaxed-structure-start',
            'vacancy-relaxed-structure-end',
            'vacancy-relaxed-structure-saddle-point',
        ]
        res = self._res
        for prop in structuralProps:
            res[prop] = _format('structure.tgz')

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
        res['vacancy-position-saddle-point'] = _format(self.migration / 2)
        res['vacancy-short-name-start'] = _format('')
        res['cauchy-stress'] = _format(np.zeros(6), 'GPa')
        res['temperature'] = _format(0.0, 'K')

        self._packStructure()

        aliases = C.ALIASES
        for key in res:
            # Static vacancy structure same as migration start
            if key[-6:] == '-start':
                aliases[key] = key[:-6]
        for originalName in aliases:
            aliasName = aliases[originalName]
            res[aliasName] = res[originalName]

    def getResult(self):
        self._addInfo()
        return self._res
