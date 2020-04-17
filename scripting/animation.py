import numpy as np 
import matplotlib.pyplot as plt
import modules.histograms as hist
import modules.sinkhorn_algorithm as sink
import modules.barycenter as bc
import modules.plot as plot
import modules.jobs as jobs
import modules.ir as ir
import math, os
import moviepy.editor as mvp
import modules.ir as ir


## ================================================================= ##
# ANIMATION
#script used for generating frames of pyplot figures
#for example, we can animate changing epsilon or the exponent in the 
#ground cost matrix function



## ======== SETUP ======== ##
def setup(name):
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



## ======== FRAME COMBINATION ======== ##
def make_gif(name):
	frames = [frames_folder + f for f in os.listdir(frames_folder)]
	clip = mvp.ImageSequenceClip(frames, fps=14)
	clip.write_gif(animation_folder + f'{name}_{index}.gif', fps=14)



## ======== HIST SELECTION ======== ##
#some histogram functions
def func(x):
	return np.where(x < 0.5, x, 1-x)

def func2(x):
	return 0.5 - func(x)

#choose two histograms a and b:

a = hist.gaussian(400, 0.4, 0.045, 0)
# a = hist.gaussian(400, 0.2, 0.06) + hist.gaussian(400, 0.5, 0.06)*2 + hist.gaussian(400, 0.8, 0.06)
# a = hist.slater(400, 0.5, 30, 0)
# a = hist.from_func(400, lambda x: 1-x**2)
# a = hist.from_func(400, func)
# a = hist.from_func(400, lambda x: 1-x,0)
# a = hist.from_func(400, lambda x: np.cos(x*6*3.14)+1)
# a = hist.from_func(400, lambda x: x)
# a = hist.dirac_delta(400, 0.5)
# a = hist.from_func(400, lambda x: ((x-0.5)*10)**4) + 3*hist.gaussian(400, 0.7, 0.1)
# r = jobs.DFTJob('l-alanine', job_name='l-alanine_DFT').run(); a = ir.get_spectrum(r, xlim=(0,2000), n=400)


b = hist.gaussian(400, 0.7, 0.08, 0)
# b = hist.gaussian(400, 0.8, 0.05, 0) + hist.gaussian(400, 0.2, 0.05, 0)
# b = hist.gaussian(400, 0.5, 0.05, 0)
# b = hist.gaussian(400, 0.2, 0.06)*2 + hist.gaussian(400, 0.5, 0.06) + hist.gaussian(400, 0.8, 0.06)*2
# b = hist.from_func(400, lambda x: np.sin(x*10*3.14)+1)
# b = hist.from_func(400, func2)
# b = hist.from_func(400, lambda x: x,0)
# b = hist.from_func(400, lambda x: x**0)
# b = hist.from_func(400, lambda x: np.cos(x*5*3.14)+1)
# r = jobs.DFTBJob('l-alanine', job_name='l-alanine_DFT').run(); b = ir.get_spectrum(r, xlim=(0,2000), n=400)





## ======== FRAME GENERATION ======== ##

# for i in np.linspace(0,6, 75):
# 	# e = (math.exp(i*0.001)-1)*4
# 	e = 0.3
# 	converge_thresh = 10**-15
	
# 	cost_fn = lambda x1, x2: abs(x1-x2)**i #function used for calculating ground cost matrix
# 	error_fn = lambda x1, x2: np.sqrt(np.sum((x1-x2)**2)) #function used for calculating errors

# 	plot_subtitle = fr'ground cost matrix function exponent factor $\alpha = {i:.3f}$'
# 	save_to = frames_folder + f'{i:.3f}.png'

# 	res = sink.sinkhorn(a, b, e, converge_thresh=converge_thresh, cost_fn=cost_fn, error_fn=error_fn)
# 	plot.plot_results(res, show_plot=False, save_to=save_to, plot_subtitle=plot_subtitle)



# for i in np.linspace(0, 1, 50):
# 	# e = (math.exp(i*0.001)-1)*4
# 	e = 0.3
# 	converge_thresh = 10**-15
	
# 	cost_fn = lambda x1, x2: abs(x1-x2)**2 #function used for calculating ground cost matrix
# 	error_fn = lambda x1, x2: np.sqrt(np.sum((x1-x2)**2)) #function used for calculating errors

# 	plot_subtitle = fr'ground cost matrix function exponent factor $\alpha = {i:.3f}$'
# 	save_to = frames_folder + f'{i:.3f}.png'

# 	res = sink.sinkhorn(a, b, e, converge_thresh=converge_thresh, cost_fn=cost_fn, error_fn=error_fn)
# 	bc_map = res.bc_map

# 	plot.plot_transport(a, b, bc_map, weight=i, show_plot=False, save_to=save_to)

# setup('interpolations')
# for i in np.linspace(0, 1, 50):
# 	save_to = frames_folder + f'{i:.3f}.png'

# 	bc_map = bc.barycenter(np.vstack((a,b)), (i, 1-i))

# 	plot.plot_hists((a, b, bc_map), labels=('a', 'b', 'barycenter interp'), title=rf'Interpolations with $\alpha={i:.3f}$', show_plot=False, save_to=save_to)
# make_gif()


setup('analytical_interp')
for i in np.linspace(0, 1, 50):
	save_to = frames_folder + f'{i:.3f}.png'

	bc_map = bc.barycenter(np.vstack((a,b)), (i, 1-i))
	analytical = hist.gaussian(400, 0.7-0.3*i, 0.08-0.035*i)

	plot.plot_hists((a, b, analytical, bc_map), labels=('a', 'b', 'analytical interp', 'barycenter interp'), title=rf'Interpolations with $\alpha={i:.3f}$', show_plot=False, save_to=save_to)
make_gif('analytical_interp')


