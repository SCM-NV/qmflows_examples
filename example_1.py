from scm.plams import Molecule
from qmflows import (
    adf, dftb, orca, run, templates, Settings)

molecule = Molecule('files/acetonitrile.xyz')

inp = templates.geometry

# Geometry optimization with DFTB (ADF)
dftb_opt = dftb(inp, molecule)

# Geometry optimization with ADF
inp.functional = "pbe"
inp.basis = "DZP"
adf_opt = adf(inp, dftb_opt.molecule)

# Excited state geometry optimization
cis = Settings()
cis.basis = 'def2TZVP'
cis.specific.orca.cis.Nroots = 1
cis.specific.orca.cis.iroot = 1
cis.specific.orca.main = 'RHF tightscf'

orca_cis = orca(cis, adf_opt.molecule, job_name='cis')

# Running the results
result = run(orca_cis.molecule, folder='cis')

print(result)
