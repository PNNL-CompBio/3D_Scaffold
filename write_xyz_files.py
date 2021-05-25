import os
import ase
from ase import io
from ase.io import read
from ase.db import connect

generated_molecules = [i.rstrip() for i in os.listdir('.') if i.endswith('.db')][0] 

with connect(generated_molecules) as conn:
    for i in range(conn.count()):
        molecule = conn.get(i+1).toatoms()
        ase.io.write(str(i)+'.xyz', molecule, format='xyz', parallel=True, append=False)
