from scm.plams import Molecule
from qmflows import (adf, dftb, run, templates)

molecule = Molecule('acetonitrile.xyz')

inp = templates.geometry

# Geometry optimization with DFTB (ADF)
dftb_opt = dftb(inp, molecule)

# Geometry optimization with ADF
inp.functional = "PBE"
inp.basis = "DZP"

adf_opt = adf(inp, dftb_opt.molecule)

# Running the results
result = run(adf_opt.molecule)

print(result)
