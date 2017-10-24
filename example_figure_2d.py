from scm.plams import Molecule
from qmflows import (adf, run, templates)

molecule = Molecule('files/uranyl.xyz')
inp = templates.geometry

# Generic Settings
inp.basis = 'TZ2P'
inp.functional = 'bp86'

# Specific
inp.specific.adf.charge = "2  0"
inp.specific.adf.relativistic = "scalar Zora"
inp.specific.adf.geometry.convergence = 'grad=1e-4'

result = run(adf(inp, molecule))
