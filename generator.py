"""
Generate test_generator.json for tests generation

Date: 2015/08/24
Author: Junhao Li <streaver91@gmail.com>
"""

import json
from ase.data import chemical_symbols
from ase.data import reference_states
from random import randint

lattices = ['fcc', 'bcc', 'sc', 'diamond', 'hcp']

# Parameters for Debugging
lattices = ['fcc']
# chemical_symbols = ['Al', 'Ni', 'Si']

with open('test_generator.json', 'w') as f:
    for pk, elem in enumerate(chemical_symbols):
        for lattice in lattices:
            kimnum = '{:012d}'.format(randint(0, 10**12 - 1))

            # hcp lattice will use a different lattice constant driver
            if lattice == 'hcp':
                latticeProperty = 'structure-hexagonal-crystal-npt'
            else:
                latticeProperty = 'structure-cubic-crystal-npt'

            f.write(json.dumps({
                'elem': elem,
                'lattice': lattice,
                'kimnum': kimnum,
                'latticeProperty': latticeProperty,
            }) + '\n')
