import os
import ase
#from ase.io import write
from ase.db import connect
import sys

generated_molecules = [i.rstrip() for i in os.listdir(sys.argv[1]) if i.endswith('.db')][0] 
if not os.path.exists(sys.argv[1]+"/generated_xyz_files"):
    os.mkdir(sys.argv[1]+"/generated_xyz_files")
with connect(sys.argv[1] + '/' + generated_molecules) as conn:
    for i in range(conn.count()):
        molecule = conn.get(i+1).toatoms()
        ase.io.write(sys.argv[1]+"/generated_xyz_files/"+str(i)+'.xyz', molecule, format='xyz', parallel=True, append=False)
