import numpy as np 
import os
import scm.plams as plams


## ================================================================= ##
# MOLECULE_FUNCS
#module containing functions for manipulating molecules

structures_folder = os.getcwd() + r'\structures\\'


def find_mol(mol):
		'''
		Method that loads a given molecule. Must be either a path to
		a file, if it is not found it will search in structures_folder,
		otherwise it will search pubchem for the structure

		mol - string specifying desired molecule
		'''

		

		#check if the file exists on disk
		if os.path.isfile(mol):
			molecule = plams.Molecule(mol)
			name = mol.split('\\')[-1][:-4]

		else:
			#else look in strucutres_folder
			if os.path.isfile(structures_folder + mol + '.xyz'):
				molecule = plams.Molecule(structures_folder + mol + '.xyz')
				name = mol

			#if still not found search pubchem
			else:
				import pubchempy as pcp

				mols = pcp.get_compounds(mol, 'name', record_type='3d')

				#check if there was a match
				if len(mols) == 0:
					raise Exception(f'Molecule {mol} not found on PubChem or on disk.')

				#save xyz file to disk
				else:
					coords = np.asarray([[a.x, a.y, a.z] for a in mols[0].atoms])
					coords = np.where(coords == None, 0, coords).astype(float)
					elements = np.asarray([a.element for a in mols[0].atoms])

					name = mol.lower()
					mol = structures_folder + f'{name}.xyz'

					with open(mol, 'w+') as f:
						f.write(f'{len(elements)}\n')
						f.write('Downloaded from PubChem using PucbChemPy\n')
						for i, e in enumerate(elements):
							f.write(f'{e: <2} \t {coords[i][0]: >8.5f} \t {coords[i][1]: >8.5f} \t {coords[i][2]: >8.5f}\n')

					molecule = plams.Molecule(mol)

		return molecule, name



def get_coords(mol):
	'''
	Function that gets the coordinates of atoms in molecule

	mol - molecules object or string
	'''

	if type(mol) is str: molecule, _ = find_mol(mol)
	coords = np.asarray([atom.coords for atom in molecule.atoms])

	return coords



def mirror(mol):
	'''
	Function that mirrors the molecule

	mol - molecules object or string
	'''

	if type(mol) is str: molecule, _ = find_mol(mol)
	for atom in molecule.atoms:
		atom.x *= -1


	return molecule



def save_to_xyz(mol, path, comment=''):
	'''
	Function that writes a molecule to xyz format

	mol - molecule object
	path - path to write xyz file to. If not a valid path is given it will write to a default directory
	'''

	if not os.path.exists(path): path = structures_folder + f'{path}.xyz'

	elements = [a.symbol for a in mol.atoms]
	coords = [a.coords for a in mol.atoms]

	with open(path, 'w+') as f:
		f.write(f'{len(elements)}\n')
		f.write('\n')
		for i, e in enumerate(elements):
			f.write(f'{e: <2} \t {coords[i][0]: >8.5f} \t {coords[i][1]: >8.5f} \t {coords[i][2]: >8.5f}\n')

