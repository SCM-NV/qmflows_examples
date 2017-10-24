from scm.plams import Molecule
from qmflows import (adf, run, templates)

molecule = Molecule('files/uranyl.xyz')
inp = templates.geometry

# Generic Settings
# inp.basis = 'TZ2P'
inp.basis = 'DZP'
inp.functional = 'bp86'

# Specific
inp.specific.adf.relativistic = "scalar Zora"
inp.specific.adf.scf.iterations = 500

result = run(adf(inp, molecule))
