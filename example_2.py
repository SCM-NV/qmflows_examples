from noodles import gather
from qmflows import (Settings, adf, dftb, run, templates)
from qmflows.molkit import from_smiles

# List of Molecules to simulate
smiles = {'naphthalene': 'C1=CC2=C(C=C1)C=CC=C2',
          'azulene': 'C1=CC2=CC=CC=CC2=C1',
          'cycloheptatriene': 'C1=CC=C(C=C1)OC1=CC=CC=C1',
          'furaldehyde': 'O=CC1=COC=C1',
          'methylthiophene': 'CC1=C(C=CS1)C#C'}


# transform the string into a format understandable by Qmflows
molecules = {name: from_smiles(smile) for name, smile in smiles.items()}

# Used DFTB to optimize the geometries
dftb_jobs = {name: dftb(templates.geometry, mol, name='dftb_{}'.format(name))
             for name, mol in molecules.items()}
optimized_mols = {name: job.molecule for name, job in dftb_jobs.items()}

# Settings for ADF
s = Settings()
s.basis = 'DZP'
s.functional = 'pbe'
s.specific.adf.scf.converge = 1e-6
s.specific.adf.symmetry = 'nosym'

# Compute the single point calculation
singlepoints = {name: adf(s, mol, name='adf_{}'.format(name))
                for name, mol in optimized_mols.items()}

# Extract the HOMO-LUMO values
homos_lumos = {name: gather(sp.homo, sp.lumo)
               for name, sp in singlepoints.items()}

# Run the computation
results = run(gather(**homos_lumos))

# Compute the Gap
print(results)
# gaps = [(lumo - homo) * 27.21 for homo, lumo in results]


# print("HOMO-LUMO gap (eV)", gaps)
