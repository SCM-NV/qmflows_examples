from scm.plams import Molecule
from qmflows import (
    adf, dftb, orca, run, templates, Settings)

molecule = Molecule('files/chlorophyll_a.xyz')

inp = templates.geometry

# Geometry optimization with DFTB (ADF)
dftb_opt = dftb(inp, molecule)

# Geometry optimization with ADF
inp.functional = "pbe"
inp.basis = "DZP"
adf_opt = adf(inp, dftb_opt.molecule)

# single point calculation using EOMS-CCSD
ccsd = Settings()
ccsd.specific.orca.mdci.nroots = 10
ccsd.specific.orca.main = 'Pal4 RHF EOM-CCSD def2-TZVP tightscf'

orca_ccsd = orca(ccsd, adf_opt.molecule, job_name='ccsd')

# Running the results
run(orca_ccsd, folder='ccsd')

