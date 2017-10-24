from scm.plams import Molecule
from qmflows import (adf, run, templates)

molecule = Molecule('files/uranyl.xyz')

inp = templates.geometry
result = run(adf(inp, molecule))

print(result.molecule)
