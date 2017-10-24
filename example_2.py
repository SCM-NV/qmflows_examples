from noodles import gather
from qmflows import (Settings, adf, dftb, run, templates)
from qmflows.molkit import from_smiles

# List of Molecules to simulate
smiles = [
        'C1=CC2=C(C=C1)C=CC=C2',
        'C1CC2(CCCCC2)C=C1',
        'C1=CC2=CC=CC=CC2=C1',
        'C1=CC=CC=C1',
        'C1CCC=CC1',
        'C1=CC=C(C=C1)OC1=CC=CC=C1',
        'C1C=CC=CC=C1'
]

# transform the string into a format understandable by Qmflows
molecules = map(from_smiles, smiles)

# Used DFTB to optimize the geometries
dftb_jobs = [dftb(templates.geometry, mol) for mol in molecules]
optimized_mols = [job.molecule for job in dftb_jobs]

# Settings for ADF
s = Settings()
s.basis = 'DZP'
s.functional = 'b3lyp'
s.specific.adf.scf.converge = 1e-6
s.specific.adf.symmetry = 'nosym'

# Compute the single point calculation
singlepoints = [adf(s, mol) for mol in optimized_mols]

# Extract the HOMO-LUMO values
homos_lumos = [gather(sp.homo, sp.lumo) for sp in singlepoints]

# Run the computation
results = run(gather(*homos_lumos))

# Compute the Gap
gaps = [(lumo - homo) * 27.21 for homo, lumo in results]


print("HOMO-LUMO gap (eV)", gaps)
