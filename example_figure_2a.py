from qmflows.molkit import from_smiles

# List of Molecules to simulate
smiles = ['C1=CC2=C(C=C1)C=CC=C2', 'C1CC2(CCCCC2)C=C1']

# Transform the smiles into a format understandable by QMFlows
molecules = [from_smiles(s) for s in smiles]
