import numpy as np
import scm.plams as plams
# import modules.molecule_funcs as mf
import os, time


## ================================================================= ##
# JOBS
#classes and functions used to run jobs for ADF

#default settings:
structures_folder = os.getcwd() + r'\structures\\'


class JobQueue(list):
	def __init__(self, jobs=[], run_path=os.getcwd()+r'\RUNS'):
		self.jobs = jobs
		self.run_path = run_path


	def append(self, job):
		self.jobs.append(job)


	def run(self):
		folder = time.strftime("[%H:%M](%d-%m-%Y)", time.localtime())

		plams.init(path=self.run_path, folder=folder)

		kffiles = []
		for job in self.jobs:
			res = job.run(False)
			kffiles.append(res.KFPATH)

		plams.finish()

		return kffiles


	def clear(self):
		self.jobs = []



class Job:
	def __init__(self, mol_file, job_name, settings=None, geo_opt=False):
		self.job_name = job_name
		self.mol = plams.Molecule(mol_file)
		#sort atoms to prevent errors in ADF
		self.mol.atoms.sort(key=lambda x: x.symbol)

		self.settings = plams.Settings()
		self._set_std_settings(geo_opt)

		if settings is not None:
			self.settings.update(settings)



class DFTJob(Job):
	'''
	Class used for geometry optimization + frequency jobs using DFT
	'''

	def _set_std_settings(self, geo_opt=False):
		'''
		Method that specifies standard settings for a DFT geometry optimization + freqs job
		'''

		if geo_opt: self.settings.input.Geometry

		self.settings.input.Basis.type = 'DZP'
		self.settings.input.Basis.core = 'None'
		# self.settings.input.XC.GGA = 'BP86'
		# self.settings.input['Relativistic Scalar'] = 'ZORA'
		self.settings.input.AnalyticalFreq 
		# self.settings.input.NumericalQuality = 'Excellent'
		self.settings.input.SYMMETRY = 'NOSYM'
		self.settings.input.VCD = 'Yes'


	def run(self, init=True, path=None):
		'''
		Method that runs this job
		'''
		if init: 
			if path is None:
				plams.init(path=os.getcwd()+r'\RUNS', folder=time.strftime("[%H:%M](%d-%m-%Y)", time.localtime()))
			else:
				plams.init(path=path)

		s = self.settings
		job = plams.ADFJob(molecule=self.mol, name=self.job_name, settings=s)
		results = job.run()
		results.dir = '\\'.join(results._kfpath().split('\\')[:-2])
		results.KFPATH = results._kfpath()

		if init: plams.finish()

		return results._kfpath()



class DFTBJob(Job):
	'''
	Class used for geometry optimization + frequency jobs using DF tight binding methods
	'''

	def _set_std_settings(self, geo_opt=False):
		'''
		Method that specifies standard settings for a DFTB geometry optimization + freqs job
		'''

		if geo_opt: self.settings.input.ams.Task = 'GeometryOptimization'
		else: self.settings.input.ams.Task = 'SinglePoint'

		self.settings.input.ams.Properties.NormalModes = 'Yes'
		self.settings.input.DFTB
		self.settings.input.DFTB.Model = 'GFN1-xTB'
		self.settings.input.DFTB.ResourcesDir = 'GFN1-xTB'
		self.settings.input.DFTB.Properties.VCD = 'Yes'	


	def run(self, init=True, path=None):
		'''
		Method that runs this job
		'''
		if init: 
			if path is None:
				plams.init(path=os.getcwd()+r'\RUNS', folder=time.strftime("[%H:%M](%d-%m-%Y)", time.localtime()))
			else:
				plams.init(path=path)

		s = self.settings
		job = plams.AMSJob(molecule=self.mol, name=self.job_name, settings=s)
		results = job.run()
		results.dir = '\\'.join(results.rkfpath('dftb').split('\\')[:-2])
		results.KFPATH = results.rkfpath('dftb')

		if init: plams.finish()

		self.results = results

		return results.rkfpath('dftb')

