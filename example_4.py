from noodles import (gather, schedule)
from qmflows import (adf, dftb, orca, run, templates)
from scm.plams import Molecule
import os

@schedule
def skip_dftb(job, mol):
    """
    return optimized molecule in case of successful otherwise
    used default geometry.
    """
    if job.status not in ['failed', 'crashed']:
        mol = job.molecule
    return mol

# Create Molecule object from xyz coordinates
path  = 'files/biphenyls/'
content = os.listdir(path)
files = [os.path.join(path, f) for f in content]
names = [os.path.splitext(f)[0] for f in content]
molecules = {n: Molecule(f) for f, n in zip(files, names)}

# DFTB optimization
dftb_jobs = {
    n: dftb(templates.geometry, mol, job_name='dftb_{}'.format(n))
    for n, mol in molecules.items()}

# Check preoptimized molecules
opt_mols = {
    n: skip_dftb(job, molecules[n]) for n, job in dftb_jobs.items()}

# DFT optimization with Orca
s = templates.geometry
s.functional = 'BP86'
s.basis = 'def2TZVP'
dft_jobs = {n: orca(s, mol, job_name='orca_{}'.format(n)) for n, mol in opt_mols.items()}

# TD-DFT ADF with COSMO
s = templates.singlepoint
s.functional = 'camy-b3lyp'
s.basis = 'TZ2P'
s.specific.adf.Solvation.solvent = 'name=Dichloromethane emp=0.0 neql=2.028'

td_dft_jobs = {
    n: adf(s, job.molecule, job_name='td_dft_{}'.format(n))
    for n, job in dft_jobs.items()}

energies = [gather(n, j.energy) for n, j in td_dft_jobs.items()]

result = run(gather(*energies), n_processes=4, folder='biphenyl')

print(result)
