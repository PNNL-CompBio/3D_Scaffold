import os
# File path to your FASTA file
path_to_file = '500_random_SMILES.txt'
save_dir = 'murkowski_scaffolds'
# Open file with "with" statement to avoid problems with access
# to original file (in case computer hangs
# or there will be any other problem)
with open(path_to_file, mode='r') as handle:
    for line in handle:
        funcional_group_smiles = line.strip()
        print(funcional_group_smiles)
