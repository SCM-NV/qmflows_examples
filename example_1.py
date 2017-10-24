from scm.plams import Molecule
from qmflows import (
    adf, dftb, gamess, run, templates, Settings)

molecule = Molecule('files/acetonitrile.xyz')

inp = templates.geometry

# Geometry optimization with DFTB (ADF)
dftb_opt = dftb(inp, molecule)

# Geometry optimization with ADF
inp.functional = "pbe"
inp.basis = "DZP"
adf_opt = adf(inp, dftb_opt.molecule)

# Single and double excitations CISD using Gamess
cisd = Settings()
cisd.specific.gamess.control.cityp = 'GUGA'
cisd.specific.gamess.control.runtype = 'energy'
cisd.specific.gamess.cidrt.group = 'c3v'
cisd.specific.gamess.cidrt.iexcit = 2
cisd.specific.gamess.basis.gbasis = 'sto'
cisd.specific.gamess.basis.ngauss = 3

gamess_cisd = gamess(cisd, adf_opt.molecule)

# Running the results
result = run(gamess_cisd.energy)

print(result)
