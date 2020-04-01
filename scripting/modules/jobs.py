import numpy as np
import scm.plams as plams
import os, time


## ================================================================= ##
# JOBS
#classes and functions used to run jobs for ADF

#default settings:
structures_folder = os.getcwd() + r'\structures\\'


class JobQueue(list):
	def __init__(self, jobs=[]):
		self.jobs = jobs


	def append(self, job):
		self.jobs.append(job)


	def run(self):
		path = os.getcwd()+r'\RUNS'
		folder = time.strftime("%d-%m-%Y", time.localtime())

		plams.init(path=path, folder=folder)

		results = []
		for job in self.jobs:
			res = job.run(False)
			results.append(res)

		plams.finish()

		return results



class Job:
	def __init__(self, mol, settings=None, job_name=None):
		self._find_mol(mol)
		self._job_name = job_name

		if settings is None:
			self.settings = plams.Settings()
			self._set_std_settings()

		else: self.settings = settings


	def _find_mol(self, mol):
		'''
		Method that loads a given molecule. Must be either a path to
		a file, if it is not found it will search in structures_folder,
		otherwise it will search pubchem for the structure

		mol - string specifying desired molecule
		'''

		#check if the file exists on disk
		if os.path.isfile(mol):
			self.mol = plams.Molecule(mol)
			self.name = mol.split('\\')[-1][:-4]

		else:
			#else look in strucutres_folder
			if os.path.isfile(structures_folder + mol + '.xyz'):
				self.mol = plams.Molecule(structures_folder + mol + '.xyz')
				self.name = mol

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

					self.name = mol.lower()
					mol = structures_folder + f'{self.name}.xyz'

					with open(mol, 'w+') as f:
						f.write(f'{len(elements)}\n')
						f.write('Downloaded from PubChem using PucbChemPy\n')
						for i, e in enumerate(elements):
							f.write(f'{e: <2} \t {coords[i][0]: >8.5f} \t {coords[i][1]: >8.5f} \t {coords[i][2]: >8.5f}\n')

					self.mol = plams.Molecule(mol)


	def run(self, init=True):
		'''
		Method that runs this job
		'''
		if init: plams.init(path=os.getcwd()+r'\RUNS', folder=time.strftime("%d-%m-%Y", time.localtime()))

		s = self.settings
		job = plams.ADFJob(molecule=self.mol, name=self.name, settings=s)
		results = job.run()

		if init: plams.finish()

		return results




class DFTJob(Job):
	'''
	Class used for geometry optimization + frequency jobs using DFT
	'''

	def _set_std_settings(self):
		'''
		Method that specifies standard settings for a DFT geometry optimization + freqs job
		'''
		self.settings.input.Basis.type = 'DZP'
		self.settings.input.Basis.core = 'None'
		self.settings.input.XC.GGA = 'PBE'
		self.settings.input['Relativistic Scalar'] = 'ZORA'
		self.settings.input.Geometry
		self.settings.input.AnalyticalFreq 
		# self.settings.input.NumericalQuality = 'Excellent'
		self.settings.input.SYMMETRY = 'NOSYM'




