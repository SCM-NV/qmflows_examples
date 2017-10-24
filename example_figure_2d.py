from scm.plams import Molecule
from qmflows import (adf, run, templates)

molecule = Molecule('uranyl.xyz')
inp = templates.geometry

# Generic Settings
inp.basis = 'TZ2P'
inp.functional = 'bp86'

# Specific
inp.specific.adf.relativistic = "scalar Zora"

result = run(inp, molecule)
