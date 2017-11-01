"""
This example generates an approximate TS for rotation in 2-methyl-biphenyl
using DFTB, and performs a full TS optimization in Orca.
It illustrates the use of a hessian from one package, DFTB in this case,
to initialize a TS optimization in another, i.e. Orca in this case
"""


from qmflows import (
    orca, dftb, templates, molkit, Dihedral, run, Settings)

# Generate 2-Methyl-biphenyl molecule
mol = molkit.from_smiles('c1ccccc1c2ccccc2C')

# Define dihedral angle
dihe = Dihedral(1, 6, 7, 12)

# Constrained geometry optimization with DFTB
# The dihedral is set to 0.0 to obtain an approximate TS
s1 = Settings()
s1.constraint.update(dihe.get_settings(0.0))
dftb_opt = dftb(templates.geometry.overlay(s1), mol)

# Calculate the DFTB hessian
dftb_freq = dftb(templates.freq, dftb_opt.molecule)

# Transition state calculation using the DFTB hessian as starting point
s2 = Settings()
s2.inithess = dftb_freq.hessian
orca_ts = orca(templates.ts.overlay(s2), dftb_opt.molecule)
orca_freq = orca(templates.freq, orca_ts.molecule)

# Execute the workflow
result = run(orca_freq)

# Analyse the result
ts_dihe = round(dihe.get_current_value(result.molecule))
frequencies = [f for f in result.frequencies if abs(f) > 1e-3]

print('Dihedral angle (degrees): {:.0f}'.format(ts_dihe))
print('Three lowest frequencies: ', frequencies[:3])
