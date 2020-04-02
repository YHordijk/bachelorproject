import numpy as np 
import matplotlib.pyplot as plt
import modules.histograms as hist
import modules.sinkhorn_algorithm as sink
import modules.plot as plot
import math, os
import moviepy.editor as mvp


## ================================================================= ##
# ANIMATION
#script used for generating frames of pyplot figures
#for example, we can animate changing epsilon or the exponent in the 
#ground cost matrix function
#then use https://ezgif.com/maker/ to animate the frames




#some histogram functions
def func(x):
	return np.where(x < 0.5, x, 1-x)

def func2(x):
	return 0.5 - func(x)

### SETUP
#choose two histograms a and b:
# a = hist.gaussian(600, 0.2, 0.06)
# a = hist.gaussian(600, 0.2, 0.06) + hist.gaussian(600, 0.5, 0.06)*2 + hist.gaussian(600, 0.8, 0.06)
# a = hist.slater(600, 0.5, 30, 0)
# a = hist.from_func(400, lambda x: 1-x**2)
# a = hist.from_func(400, func)
a = hist.from_func(500, lambda x: 1-x,0)
# a = hist.from_func(600, lambda x: np.cos(x*10*3.14)+1)
# a = hist.from_func(600, lambda x: x)
# a = hist.dirac_delta(600, 0.5)
a = hist.from_func(600, lambda x: ((x-0.5)*10)**4) + 3*hist.gaussian(600, 0.7, 0.1)


# b = hist.gaussian(400, 0.5, 0.05, 0)
b = hist.gaussian(600, 0.2, 0.06)*2 + hist.gaussian(600, 0.5, 0.06) + hist.gaussian(600, 0.8, 0.06)*2
# b = hist.from_func(600, lambda x: np.sin(x*10*3.14)+1)
# b = hist.from_func(400, func2)
# b = hist.from_func(500, lambda x: x,0)
# b = hist.from_func(600, lambda x: x**0)
b = hist.from_func(600, lambda x: np.cos(x*5*3.14)+1)



animation_name = 'bc_transport'

animation_folder = os.getcwd() + r'\animation\\'

#get files named the same as animation_name
animation_files = os.listdir(animation_folder)
animation_files = [f for f in animation_files if f.startswith(animation_name)]
if len(animation_files) > 0:
	#get index of folders:
	animation_files_indices = [int(f.split('_')[-1].strip('.gif')) for f in animation_files]
	#index for new folder is then the max index plus one
	index = max(animation_files_indices) + 1

else:
	index = 1

frames_folder = os.getcwd() + r'\animation_frames\\'
#specify target frames_folder
frames_folder = frames_folder + animation_name + rf'_{index}\\'
#create folder if it does not exist yet
if not os.path.exists(frames_folder): os.mkdir(frames_folder)




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



for i in np.linspace(0, 1, 50):
	# e = (math.exp(i*0.001)-1)*4
	e = 0.3
	converge_thresh = 10**-15
	
	cost_fn = lambda x1, x2: abs(x1-x2)**2 #function used for calculating ground cost matrix
	error_fn = lambda x1, x2: np.sqrt(np.sum((x1-x2)**2)) #function used for calculating errors

	plot_subtitle = fr'ground cost matrix function exponent factor $\alpha = {i:.3f}$'
	save_to = frames_folder + f'{i:.3f}.png'

	res = sink.sinkhorn(a, b, e, converge_thresh=converge_thresh, cost_fn=cost_fn, error_fn=error_fn)
	bc_map = res.bc_map

	plot.plot_transport(a, b, bc_map, weight=i, show_plot=False, save_to=save_to)


# make gif
frames = [frames_folder + f for f in os.listdir(frames_folder)]
clip = mvp.ImageSequenceClip(frames, fps=14)
clip.write_gif(animation_folder + f'{animation_name}_{index}.gif', fps=14)


