from noodles import (gather, schedule)
from plams import Molecule
from qmworks import (adf, dftb, orca, run, templates)
import os

# skip DFTB preoptimization in case DFTB fails
def skip_dftb(job, mol):
    """
    return optimized molecule in case of successful otherwise
    used unoptimized molecule.
    """
    if job.status not in ['failed', 'crashed']:
        mol = job.molecule
    return mol

# Create Molecule object from xyz coordinates
path  = 'files/byphenyls/'
content = os.listdir(path)
files = [os.path.join(path, f) for f in content]
names = [os.path.splitext(f)[0] for f in content]
molecules = {n: Molecule(f) for f, n in zip(files, names)}

# DFTB optimization
dftb_jobs = {
    n: dftb(templates.geometry, mol, job_name='dftb_{}'.format(n))
    for n, mol in molecules.items()}

# Check preoptimized molecules
schedule_skip = schedule(skip_dftb)
opt_mols = {
    n: schedule_skip(job, molecules[n]) for n, job in dftb_jobs.items()}

# DFT optimization with Orca
s = templates.geometry
s.functional = 'BP86'
s.basis = 'TZV(P)'
dft_jobs = {n: orca(s, mol, job_name='orca_{}'.format(n)) for n, mol in opt_mols.items()}

# TD-DFT ADF with COSMO
s = templates.singlepoint
s.functional = 'camy-b3lyp'
s.basis = 'TZ2P'
s.specific.adf.Solvation.solvent = 'name=Dichloromethane emp=0.0 neql=2.028'

td_dft_jobs = {n: adf(s, mol, job_name='td_dft_{}'.format(n)) for n, mol in dft_jobs.items()}

energies = [(n, j.energy) for n, j in td_dft_jobs.items()]

result = run(gather(*energies), folder='biphenyl')

print(result)
