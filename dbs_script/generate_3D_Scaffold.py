from ase.db import connect
import os
import re
from ase.io.extxyz import read_xyz
import shutil
import tempfile
import argparse


parser = argparse.ArgumentParser()
parser.add_argument('--xyz_path', type = str, help='Path to directory containing xyz_files for creating database',
                    default='./xyz_files/')

args = parser.parse_args()

tmpdir = tempfile.mkdtemp('temporary')

with connect(os.path.join('scaffold3D.db')) as con:
    ordered_files = sorted(os.listdir(args.xyz_path),
                           key=lambda x: (int(re.sub('\D', '', x)), x))
    for i, xyzfile in enumerate(ordered_files):
        xyzfile = os.path.join(args.xyz_path, xyzfile)

        if (i + 1) % 100 == 0:
            print('Parsed: {:6d}'.format(i + 1))
        properties = {}
        tmp = os.path.join(tmpdir, 'tmp.xyz')

        with open(xyzfile, 'r') as f:
            lines = f.readlines()
            l = lines[1].split()[2:]
            smiles_string_molecule = lines[-2].strip().split("\t")[0]
            smiles_string_func_grp = lines[-1].strip().split("\t")[0]
            properties['Smiles_String'] = smiles_string_molecule
            properties['Functional_Group'] = smiles_string_func_grp
            with open(tmp, "wt") as fout:
                for line in lines:
                    fout.write(line.replace('*^', 'e'))

        with open(tmp, 'r') as f:
            ats = list(read_xyz(f, 0))[0]
        print(properties)
        con.write(ats, data=properties)
print('Done.')
shutil.rmtree(tmpdir)


os.system('python {} {}'.format('preprocess_dataset.py', './scaffold3D.db'))
