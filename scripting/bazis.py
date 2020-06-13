import modules.jobs as jobs
import scm.plams as plams
import os, time, shutil


# #DFT + DFTB
# plams.init(path=os.getcwd()+r'/RUNS', folder=time.strftime("[%H-%M](%d-%m-%Y)", time.localtime()))
# kf_dir = os.getcwd() + f'/RUNS/#KFFiles/{time.strftime("[%H-%M](%d-%m-%Y)", time.localtime())}'
# try: os.mkdir(kf_dir)
# except: pass


# dataset_name = 'second_set'
# dataset_dir = os.getcwd() + f'/structures/{dataset_name}'
# mol_dirs = os.listdir(dataset_dir)
# mol_names = [d + '_' + p.strip('.xyz') for d in mol_dirs for p in os.listdir(dataset_dir + '/' + d)]
# mol_paths = [dataset_dir + '/' + d + '/' + p for d in mol_dirs for p in os.listdir(dataset_dir + '/' + d)]





# print(f'Starting singlepoint DFT and DFTB vibration calculations for {len(mol_names)} molecules from the {dataset_name} dataset:')
# print('Molecules:')
# [print('\t'+i+1,n) for i,n in enumerate(mol_names)]



# for n, p in zip(mol_names, mol_paths):
# 	kf_DFT = jobs.DFTJob(p, job_name=f'{n}_DFT', geo_opt=True).run(init=False)
# 	shutil.copy2(kf_DFT, f'{kf_dir}/{n}_DFT.t21')

# 	kf_DFTB = jobs.DFTBJob(p, job_name=f'{n}_DFTB', geo_opt=True).run(init=False)
# 	shutil.copy2(kf_DFTB, f'{kf_dir}/{n}_DFTB.rkf')






plams.init(path=os.getcwd()+r'/RUNS', folder=time.strftime("[%H-%M](%d-%m-%Y)", time.localtime()))
kf_dir = os.getcwd() + f'/RUNS/#KFFiles/{time.strftime("[%H-%M](%d-%m-%Y)", time.localtime())}'
try: os.mkdir(kf_dir)
except: pass


dataset_name = 'second_set'
dataset_dir = os.getcwd() + f'/structures/{dataset_name}'
mol_dirs = os.listdir(dataset_dir)
mol_names = [d + '_' + p.strip('.xyz') for d in mol_dirs for p in os.listdir(dataset_dir + '/' + d)]
mol_paths = [dataset_dir + '/' + d + '/' + p for d in mol_dirs for p in os.listdir(dataset_dir + '/' + d)]


functionals = ['PBE']


print(f'Starting singlepoint DFT and DFTB vibration calculations for {len(mol_names)} molecules from the {dataset_name} dataset:')
print('Molecules:')
[print('\t'+str(i+1),n) for i,n in enumerate(mol_names)]


for f in functionals:
	print(f'Starting geo-opt DFT vibration calculations for {len(mol_names)} molecules from the {dataset_name} dataset with functional {f}')

	for n, p in zip(mol_names, mol_paths):
		settings = plams.Settings()


		if f == 'PBE':
			settings.input.XC.GGA = 'PBE'
			
			kf_DFT = jobs.DFTJob(p, job_name=f'{n}_DFT', geo_opt=True).run(init=False)
			shutil.copy2(kf_DFT, f'{kf_dir}/{n}_DFT.t21')

		if f == 'DFTB3':
			settings.input.DFTB.Model = 'DFTB3'
			settings.input.DFTB.ResourcesDir = 'DFTB.org/3ob-3-1'

			kf_DFTB = jobs.DFTBJob(p, job_name=f'{n}_{f}_DFTB', geo_opt=True).run(init=False)
			shutil.copy2(kf_DFTB, f'{kf_dir}/{n}_{f}_DFTB.rkf')

		if f == 'DFTB3_freq':
			settings.input.DFTB.Model = 'DFTB3'
			settings.input.DFTB.ResourcesDir = 'DFTB.org/3ob-freq-1-2'

			kf_DFTB = jobs.DFTBJob(p, job_name=f'{n}_{f}_DFTB', geo_opt=True).run(init=False)
			shutil.copy2(kf_DFTB, f'{kf_dir}/{n}_{f}_DFTB.rkf')
