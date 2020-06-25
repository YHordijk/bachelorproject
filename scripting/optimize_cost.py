

import modules.ir as ir
import modules.plot as plot
import modules.sinkhorn_algorithm as sink
import modules.colour_maps as cmap
import os, time
import numpy as np
import modules.comp_funcs as comp_funcs
import matplotlib.pyplot as plt 
import moviepy.editor as mvp




category = ''
# category = 'Aminoindan'
# category = 'ISO34_E'
# category = 'ISO34_P'
# category = 'SCONF'
# category = 'ISO34'


# func = comp_funcs.wasserstein_distance
# func = comp_funcs.wasserstein_distance_unbalanced
func = comp_funcs.freq_int_wasserstein
# func = comp_funcs.l2
# func = comp_funcs.diagonality
# func = comp_funcs.bhattacharyya
# func = comp_funcs.correlation
# func = comp_funcs.chi_square
# func = comp_funcs.kl_divergence


functionals = ['LDA_DFT', 'DFTB3_DFTB']


bins = 20


title = category + '_' + str(func).split()[1] 


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



def make_gif(frames):
	clip = mvp.ImageSequenceClip(frames, fps=14)
	clip.write_gif(animation_folder + f'{name}_{index}.gif', fps=14)


def normalize(y): 
	y = np.asarray(y)
	return (y-y.min())/(y-y.min()).max()


def dist(P):
	j = np.ones(P.shape[0])
	r = np.arange(P.shape[0])
	r2 = r**2

	n = j@P@j.T
	sum_x = r@P@j.T
	sum_y = j@P@r.T
	sum_x2 = r2@P@j.T
	sum_y2 = j@P@r2.T 
	sum_xy = r@P@r.T

	return (n*sum_xy - sum_x*sum_y)/(np.sqrt(n*sum_x2 - sum_x**2) * np.sqrt(n*sum_y2 - sum_y**2))

i = 0

def objective(x, plot=True):
	global i 
	i += 1
	d, C = func(ir1, ir2, peaksa=peaksa, peaksb=peaksb, freq_weight=x[0], freq_exp=x[1], int_weight=x[2], int_exp=x[3])
	d_new = (1 - d)/np.sum(1 - d)

	b = np.zeros_like(d)
	b[np.arange(len(d)), d.argmin(1)] = 1

	if plot:
		plt.figure()
		plt.suptitle(rf'a={x[0]:.6f}, b={x[2]:.6f}' + '\n' + rf'$\alpha$={x[1]:.6f}, $\beta$={x[3]:.6f}')
		plt.subplot(2,2,1)
		plt.imshow(d)
		plt.gca().set_title('Distance matrix')
		# plt.colorbar()

		plt.subplot(2,2,2)
		plt.gca().set_title('Compressed distance matrix')
		plt.imshow(b)

		plt.subplot(2,2,3)
		plt.gca().set_title('Example Cost matrix')
		plt.imshow(C)
		plt.show()

		plt.savefig(frames_folder + f'{i}.png')

		global frames
		frames.append(frames_folder + f'{i}.png')
	# fun = 1-dist(b)

	fun = np.sum(np.eye(d.shape[0]) * (d_new))

	with open('progress.txt', 'a') as file:
		print(f'x: {x}, fun: {fun}', file=file)
		print(f'x: {x}, fun: {fun}')

	return fun

name = 'optimize_params'

setup()
frames = []

import scipy.optimize as opt

bound = opt.Bounds(0, np.inf)





params = [0.000040, 1.955098, 0.000282, 1.990915]


# objective([10000000,4,0,2])
objective(params)


# opt.minimize(objective, [0, 2, 1, 2], bounds=bound)




# opt.minimize(objective, [-2.56491310e-06,  1e+00,  2.04119914e-03,  0.5e+00,])

# objective([ 2.02016817e-09,  2.17189317e+00, -5.29376073e-09,  2.69305257e+00])


# make_gif(frames)