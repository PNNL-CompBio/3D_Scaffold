****3D-Scaffold: Deep Generative Framework to Generate Molecule's 3D-Cordinates with Desired Scaffolds.****
![molecule](Visual.gif)




3D-Scaffold generates entire molecules with desired scaffolds in an autoregressive fashion, placing one atom after another in 3d Euclidean space of scaffolds. It only uses the 3D positions and types of atoms in a molecule along with the SMILES string of scaffold as an input and generates novel therapeutic candidates as output. 3D-Scaffold framework is built on top of G-SchNet [1].

This work is published as [3D-Scaffold: A Deep Learning Framework to Generate 3D Coordinates of Drug-like Molecules with Desired Scaffolds](https://pubs.acs.org/doi/10.1021/acs.jpcb.1c06437). Citation information at the bottom.

**Requirements**

gschnet, 
schnetpack 0.3, 
pytorch >= 1.2, 
python >= 3.7,
ASE >= 3.17.0, 
Open Babel 2.41, 
rdkit >= 2019.03.4.0 

Please install dependencies and requirements for running the model from Reference 1.

![molecule2](Architecture.png)

# Generate a database:

To generate a database from set of XYZ files:

```bash
python generate_3D_Scaffold.py --xyz_path ./xyz_files
```

Move the generated database file 'Scaffold3D.db' to `./data/` directory

If there is no `Scaffold3D.db` in the directory, the code will automatically download the QM9 dataset for use in training.

# Training a model:

To train a model with same hyperparameters as described in paper;

```bash
python ./scaffold3D.py train 3D_Scaffold ./data/ ./model --split 2000 500 --cuda --batch_size 5 --draw_random_samples 5 --features 64 --interactions 6 --max_epochs 1000
```

# Generate molecules from trained model:

Once a model has been trained, you can use the trained model to generate molecules based on the desired scaffold (the scaffold C=CC(=O)N in the case of the example).

```bash
python ./scaffold3D.py generate 3D_Scaffold  ./model/ 100 --functional_group 'C=CC(=O)N' --chunk_size 100 --max_length 65 --file_name scaffold
```
The output of this script will be a dictionary of molecules saved as `{--file_name}.dict`.

## Filter the generated molecules:

Once the molecules have been generated and stored as a dictionary, we can filter the duplicates/unique molecules.

```bash
python filter_generated.py ./model/generated/
```

## Write generated molecules in to xyz file

To view the molecules that we have generated as xyz files, we can run the `write_xyz.py` script to convert the `scaffold3D.db` into XYZ files. This will create a new directory `/generated_xyz_files` in the directory which you run the code.

```bash
python write_xyz.py ./models/generated/
```

`./models/generated/generated_xyz_files` will be made in this example

**References;**
1. Gebauer, N.; Gastegger, M.; Sch ̈utt, K. Symmetry-adapted generation of 3d point setsfor  the  targeted  discovery  of  molecules.  Advances  in  Neural  Information  ProcessingSystems. 2019; pp 7566–7578.

# How to Cite:

```bibtex
@article{Joshi2021,
author = {Joshi, Rajendra P. and Gebauer, Niklas W. A. and Bontha, Mridula and Khazaieli, Mercedeh and James, Rhema M. and Brown, James B. and Kumar, Neeraj},
doi = {10.1021/acs.jpcb.1c06437},
issn = {1520-6106},
journal = {The Journal of Physical Chemistry B},
month = {nov},
number = {44},
pages = {12166--12176},
title = {{3D-Scaffold: A Deep Learning Framework to Generate 3D Coordinates of Drug-like Molecules with Desired Scaffolds}},
url = {https://pubs.acs.org/doi/10.1021/acs.jpcb.1c06437},
volume = {125},
year = {2021}
}
```