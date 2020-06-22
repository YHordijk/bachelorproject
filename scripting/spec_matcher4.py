import modules.ir as ir
import modules.plot as plot
import modules.sinkhorn_algorithm as sink
import modules.colour_maps as cmap
import os
import numpy as np
import modules.comp_funcs as comp_funcs
import matplotlib.pyplot as plt 
import moviepy.editor as mvp



def setup():
	global animation_folder
	animation_folder = os.getcwd() + r'\animation\\'

	#get files named the same as name
	animation_files = os.listdir(animation_folder)
	animation_files = [f for f in animation_files if f.startswith(name)]
	global index

	if len(animation_files) > 0:
		#get index of folders:
		animation_files_indices = [int(f.split('_')[-1].strip('.gif')) for f in animation_files]
		#index for new folder is then the max index plus one
		index = max(animation_files_indices) + 1

	else:
		index = 1

	global frames_folder
	frames_folder = os.getcwd() + r'\animation_frames\\'
	#specify target frames_folder
	frames_folder = frames_folder + name + rf'_{index}\\'
	#create folder if it does not exist yet
	if not os.path.exists(frames_folder): os.mkdir(frames_folder)



def make_gif(ds, titles=None):
	frames = []
	if titles is None: titles = ['' for _ in range(len(ds))]

	vmin = min([np.min(d) for d in ds])
	vmax = max([np.max(d) for d in ds])

	for i, d, title in zip(range(len(ds)), ds, titles):
		plt.figure(1)
		plt.imshow(d)
		plt.gca().invert_yaxis()
		plt.title(title)
		plt.colorbar()
		

		# plt.figure(2)
		# b = np.zeros_like(d)
		# b[np.arange(len(d)), d.argmin(1)] = 1
		# plt.imshow(b)
		# plt.gca().invert_yaxis()
		# plt.show()

		plt.savefig(frames_folder + f'{i}.png')
		frames.append(frames_folder + f'{i}.png')

		plt.close()

	clip = mvp.ImageSequenceClip(frames, fps=14)
	clip.write_gif(animation_folder + f'{name}_{index}.gif', fps=15)




def main(functionals, func, category='Aminoindan', bins=600, reg_m=10**2, **kwargs):

			#####################
			#### SETUP START ####
			#####################


	############################################ Directory with KF files


	kf_dir = r"D:\Users\Yuman\Desktop\Programmeren\bachelorproject\scripting\RUNS\#KFFiles\functionals"



			#####################
			####  SETUP END  ####
			#####################


	func_best = {comp_funcs.wasserstein_distance: 				min, 
				 comp_funcs.wasserstein_distance_unbalanced: 	min, 
				 comp_funcs.freq_int_wasserstein:				min,
				 comp_funcs.l2:									min, 
				 comp_funcs.diagonality: 						max, 
				 comp_funcs.bhattacharyya: 						min, 
				 comp_funcs.correlation:						min, 
				 comp_funcs.chi_square: 						min, 
				 comp_funcs.kl_divergence: 						min}


	#get kffiles
	kff = [f for f in os.listdir(kf_dir)]
	kff = list(filter(lambda x: functionals[0] in x or functionals[1] in x, kff))


	kff_of_cat = [kf_dir + '\\' + f for f in kff if f.startswith(category)]


	#get the spectra
	ir1 = [ir.get_spectrum_from_kf(f,width=50,n=bins) for f in kff_of_cat if functionals[0] in f]
	ir2 = [ir.get_spectrum_from_kf(f,width=50,n=bins) for f in kff_of_cat if functionals[1] in f]

	peaksa = [ir.get_freqs_intens(f) for f in kff_of_cat if functionals[0] in f]
	freqsa, intensa = [list(c) for c in zip(*peaksa)]
	peaksa = [list(zip(freqsa[i], intensa[i])) for i in range(len(freqsa))]

	peaksb = [ir.get_freqs_intens(f) for f in kff_of_cat if functionals[1] in f]
	freqsb, intensb = [list(c) for c in zip(*peaksb)]
	peaksb = [list(zip(freqsb[i], intensb[i])) for i in range(len(freqsb))]


	d = func(ir1, ir2, reg_m=reg_m, peaksa=peaksa, peaksb=peaksb, **kwargs)
	return d.astype(float)





category = ''
# category = 'Aminoindan'
# category = 'ISO34_E'
# category = 'ISO34_P'
# category = 'SCONF'


# func = comp_funcs.wasserstein_distance
# func = comp_funcs.wasserstein_distance_unbalanced
func = comp_funcs.freq_int_wasserstein
# func = comp_funcs.l2
# func = comp_funcs.diagonality
# func = comp_funcs.bhattacharyya
# func = comp_funcs.correlation
# func = comp_funcs.chi_square
# func = comp_funcs.kl_divergence




title = category + '_' + str(func).split()[1] 






global name
name = 'spec_match_reg_m'
setup()
ran = (0, 4)
titles = []
ds = []
for param in np.linspace(*ran, 75):
	print(param)
	d = main(['LDA_DFT','DFTB3_DFTB'], func=func, category=category, bins=601, int_weight=0, freq_exp=param)
	# d = 

	titles.append(f'{title}\nparam={param:.2f}')
	ds.append(d)


make_gif(ds, titles)

# plt.figure(1)
# plt.imshow(ds[0])
# plt.gca().invert_yaxis()
# plt.colorbar()
# plt.show()
