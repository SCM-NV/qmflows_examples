from noodles import (gather, lift, schedule)
from qmflows import (Settings, adf, dftb, run, templates)
from qmflows.molkit import from_smiles
from scipy.constants import physical_constants
from scm.plams import KFReader
import numpy as np

eV = physical_constants['electron volt-hartree relationship'][0]

@schedule
def filter_homo_lumo_lower_than(jobs, x):
    """
    Filter the `jobs` that fulfill that the HOMO-LUMO gap
    is lower than x
    """
    return [j for j in jobs if (j.lumo - j.homo) / eV < x]


@schedule
def iterate_over_jobs(promises, fun, inp, prop, prefix=None):
    """
    Iterate over a list of `promised` jobs, calling function `fun`
    with input `inp` and property `prop` from previous jobs.
    """
    return gather(*[fun(inp, getattr(job, prop),
                job_name='{}_{}'.format(prefix, get_job_name(job)))
            for job in promises])


def extract_tdfft_excitations(job):
    """
    Extract the excitation Energies for the `job` using the KF binary
    file from ADF.
    """
    kf = KFReader(job.kf.path)
    energies = np.array(kf.read('All excitations', 'All Sing-Sing excitations'))
    
    return energies / eV


def get_job_name(job):
    """get the base job name"""
    name = job.job_name
    return name.split('_')[1]


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

# Settings for ADF SAOP single point calculation
s = Settings()
s.basis = 'DZP'
s.specific.adf.basis.core = None
s.specific.adf.xc.model = 'saop'
s.specific.adf.scf.converge = 1e-6
s.specific.adf.symmetry = 'nosym'

# Compute the single point calculation
singlepoints = [adf(s, mol, job_name='adf_{}'.format(name))
                for name, mol in optimized_mols.items()]

# Filter results with HOMO-LUMO gap lower than 3 eV
candidates = filter_homo_lumo_lower_than(gather(*singlepoints), 3)

# Optimize the selected candidates
inp = templates.geometry
inp.functional = 'pbe'
inp.basis = 'DZP'
inp.specific.adf.scf.converge = 1e-6
inp.specific.adf.symmetry = 'nosym'

opt_jobs = iterate_over_jobs(candidates, adf, inp, 'molecule', prefix='pbe')

# Single point TD-DFT calculations
s.specific.adf.excitations.allowed = ''
s.specific.adf.excitations.lowest = 10

td_dft_jobs = iterate_over_jobs(opt_jobs, adf, s, 'molecule', prefix='tddft')

# Filter The optimize molecules based on TD-DFT
candidates_td_dft = filter_homo_lumo_lower_than(td_dft_jobs, 3)

# Run the computation
results = run(candidates_td_dft, folder='screening')
for r in results:
    excitations  = extract_tdfft_excitations(r)
    print('Molecule: ', get_job_name(r))
    print('Excitations (eV):\n', np.array2string(excitations, precision=2))



