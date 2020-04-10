import numpy as np
import scm.plams as plams
import modules.molecule_funcs as mf
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

	def clear(self):
		self.jobs = []



class Job:
	def __init__(self, mol, settings=None, job_name=None):
		self.mol, self.name = mf.find_mol(mol)
		if job_name is None:
			self._job_name = self.name
		else:
			self._job_name = job_name

		if settings is None:
			self.settings = plams.Settings()
			self._set_std_settings()

		else: self.settings = settings



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
		self.settings.input.VCD = 'Yes'


	def run(self, init=True):
		'''
		Method that runs this job
		'''
		if init: plams.init(path=os.getcwd()+r'\RUNS', folder=time.strftime("%d-%m-%Y", time.localtime()))

		s = self.settings
		job = plams.ADFJob(molecule=self.mol, name=self._job_name, settings=s)
		results = job.run()
		results.dir = '\\'.join(results._kfpath().split('\\')[:-2])
		results.KFPATH = results._kfpath()

		if init: plams.finish()

		self.results = results

		return results



class DFTBJob(Job):
	'''
	Class used for geometry optimization + frequency jobs using DF tight binding methods
	'''

	def _set_std_settings(self):
		'''
		Method that specifies standard settings for a DFTB geometry optimization + freqs job
		'''

		self.settings.input.ams.Task = 'GeometryOptimization'
		self.settings.input.ams.Properties.NormalModes = 'Yes'
		self.settings.input.DFTB
		self.settings.input.DFTB.Model = 'GFN1-xTB'
		# self.settings.input.DFTB.ResourcesDir = 'DFTB.org/3ob-freq-1-2'
		self.settings.input.DFTB.ResourcesDir = 'GFN1-xTB'
		self.settings.input.DFTB.Properties.VCD = 'Yes'	


	def run(self, init=True):
		'''
		Method that runs this job
		'''
		if init: plams.init(path=os.getcwd()+r'\RUNS', folder=time.strftime("%d-%m-%Y", time.localtime()))

		s = self.settings
		job = plams.AMSJob(molecule=self.mol, name=self._job_name, settings=s)
		results = job.run()
		results.dir = '\\'.join(results.rkfpath('dftb').split('\\')[:-2])
		results.KFPATH = results.rkfpath('dftb')

		if init: plams.finish()

		self.results = results

		return results

