from noodles import (gather, schedule)
from qmflows import (Settings, adf, dftb, run, templates)
from qmflows.molkit import from_smiles
from scipy.constants import physical_constants

eV = physical_constants['electron volt-hartree relationship'][0]


def filter_homo_lumo_lower_than(jobs, x):
    """
    Filter the `jobs` that fulfill that the HOMO-LUMO gap
    is lower than x
    """
    interesting = [(j.job_name, (j.lumo - j.homo) / eV) for
                   j in jobs if (j.lumo - j.homo) / eV < x]

    return interesting


# List of Molecules to simulate
smiles = {'naphthalene': 'C1=CC2=C(C=C1)C=CC=C2',
          'azulene': 'C1=CC2=CC=CC=CC2=C1',
          'cycloheptatriene': 'C1=CC=C(C=C1)OC1=CC=CC=C1',
          'furaldehyde': 'O=CC1=COC=C1',
          'methylthiophene': 'CC1=C(C=CS1)C#C'}


# transform the string into a format understandable by Qmflows
molecules = {name: from_smiles(smile) for name, smile in smiles.items()}

# Used DFTB to optimize the geometries
dftb_jobs = {name: dftb(templates.geometry, mol, job_name='dftb_{}'.format(name))
             for name, mol in molecules.items()}
optimized_mols = {name: job.molecule for name, job in dftb_jobs.items()}

# Settings for ADF
s = Settings()
s.basis = 'DZP'
s.functional = 'pbe'
s.specific.adf.scf.converge = 1e-6
s.specific.adf.symmetry = 'nosym'

# Compute the single point calculation
singlepoints = [adf(s, mol, job_name='adf_{}'.format(name))
                for name, mol in optimized_mols.items()]

# Filter results with HOMO-LUMO gap lower than 3 eV
filter_schedule = schedule(filter_homo_lumo_lower_than)
interesting = filter_schedule(gather(*singlepoints), 3)

# Run the computation
results = run(interesting)
print(results)
