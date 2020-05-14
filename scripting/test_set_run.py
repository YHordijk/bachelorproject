import modules.molecule_funcs as mf
import modules.jobs as jb
import modules.ir as ir
import modules.plot as plot
import modules.sinkhorn_algorithm as sink
import os


### ============== setup jobs ============== ###
# root = os.getcwd() + r'\structures\hydrocarbons'
# spectra_path = r"C:\Users\Yuman\Desktop\Programmeren\bachelorproject\scripting\RUNS\#KFFiles\hydrocarbons"

# mol_names = ['methane', 'ethane', 'propane', 'butane', 'pentane', 'hexane', 'heptane', 'octane']
# mol_names = ['butane', 'pentane', 'hexane', 'heptane', 'octane']

# files = [mf.get_mol_paths(mol_name, root)[0] for mol_name in mol_names]


mol_names = ['AI4','AI5','AI6']
files = [r"C:\Users\Yuman\Desktop\Programmeren\bachelorproject\scripting\structures\testSet\Aminoindan_CONF\4.xyz",
		 r"C:\Users\Yuman\Desktop\Programmeren\bachelorproject\scripting\structures\testSet\Aminoindan_CONF\5.xyz",
		 r"C:\Users\Yuman\Desktop\Programmeren\bachelorproject\scripting\structures\testSet\Aminoindan_CONF\6.xyz"]


dft_q = jb.JobQueue()
for i in range(len(files)):
	dft_q.append(jb.DFTJob(files[i], mol_names[i]+'_DFT'))

mol_names = ['AI1','AI2','AI3','AI4','AI5','AI6']
files = [r"C:\Users\Yuman\Desktop\Programmeren\bachelorproject\scripting\structures\testSet\Aminoindan_CONF\1.xyz",
		 r"C:\Users\Yuman\Desktop\Programmeren\bachelorproject\scripting\structures\testSet\Aminoindan_CONF\2.xyz",
		 r"C:\Users\Yuman\Desktop\Programmeren\bachelorproject\scripting\structures\testSet\Aminoindan_CONF\3.xyz",
		 r"C:\Users\Yuman\Desktop\Programmeren\bachelorproject\scripting\structures\testSet\Aminoindan_CONF\4.xyz",
		 r"C:\Users\Yuman\Desktop\Programmeren\bachelorproject\scripting\structures\testSet\Aminoindan_CONF\5.xyz",
		 r"C:\Users\Yuman\Desktop\Programmeren\bachelorproject\scripting\structures\testSet\Aminoindan_CONF\6.xyz"]

dftb_q = jb.JobQueue()
for i in range(len(files)):
	dftb_q.append(jb.DFTBJob(files[i], mol_names[i]+'_DFTB'))



### ============== run jobs ============== ###
kffiles_dftb = dftb_q.run()
kffiles_dft = dft_q.run()

# dft_kf_path = r"C:\Users\Yuman\Desktop\Programmeren\bachelorproject\scripting\RUNS\#KFFiles\hydrocarbons\DFT\\"
# dftb_kf_path = r"C:\Users\Yuman\Desktop\Programmeren\bachelorproject\scripting\RUNS\#KFFiles\hydrocarbons\DFTB\\"
# kffiles_dft = [dft_kf_path + f for f in os.listdir(dft_kf_path)]
# kffiles_dftb = [dftb_kf_path + f for f in os.listdir(dftb_kf_path)]


### ============== generate IR spectra ============== ###
# for name, dft_kf, dftb_kf in zip(sorted(mol_names), kffiles_dft, kffiles_dftb):
# 	irs = (ir.get_spectrum_from_kf(dft_kf), ir.get_spectrum_from_kf(dftb_kf))
# 	res = sink.sinkhorn(*irs, 0.0004)
# 	plot.plot_sink_results(res, save_to=spectra_path+ rf'\{name}.png', show_plot=False)


print(kffiles_dft)
print(kffiles_dftb)