import numpy as np 
import matplotlib.pyplot as plt
import matplotlib.colors as clr
from itertools import zip_longest


## ================================================================= ##
# plot



def plot_hists(hists, labels=None, axislabels=('Bins', 'Mass'), xlim=None,
			   title='Histogram h', invert_xaxis=False, invert_yaxis=False,
			   save_to=None, show_plot=True, shade=False):
	'''
	Function that plots multiple histograms

	hists - histograms
	'''

	for h in hists:
		h /= np.sum(np.absolute(h))


	plt.figure(figsize=(16,9))
	plt.title(title)
	plt.xlabel(axislabels[0])
	plt.ylabel(axislabels[1])
	if xlim is None: plt.xticks([]); xlim=(0,1)
	plt.yticks([])

	if invert_xaxis: plt.gca().invert_xaxis()
	if invert_yaxis: plt.gca().invert_yaxis()

	# color_array = clr.makeMappingArray(len(h), np.linspace(0,1,len(h)))

	for h, l in zip_longest(hists, labels):
		plt.plot(np.linspace(xlim[0],xlim[1],len(h)), h, linewidth=2, label=l)
		if shade: plt.fill_between(np.linspace(xlim[0],xlim[1],len(h)), h)

	plt.legend()

	#save plot if a file path is given
	if save_to is not None:
		plt.savefig(save_to)

	#show plot if needed
	if show_plot:
		plt.show()
	#otherwise clear and close plot to save memory
	else:
		plt.clf()
		plt.cla()
		plt.close()


def plot_hist(h, axislabels=('Bins', 'Mass'), xlim=None, title='Histogram h',
			  invert_xaxis=False, invert_yaxis=False, save_to=None, show_plot=True):
	'''
	Function that plots a single histogram

	h - histogram
	'''

	h /= abs(np.sum(h))


	plt.figure(figsize=(16,9))
	plt.title(title)
	plt.xlabel(axislabels[0])
	plt.ylabel(axislabels[1])
	if xlim is None: plt.xticks([]); xlim=(0,1)
	plt.yticks([])

	if invert_xaxis: plt.gca().invert_xaxis()
	if invert_yaxis: plt.gca().invert_yaxis()

	plt.plot(np.linspace(xlim[0],xlim[1],len(h)), h, 'blue', linewidth=1)
	plt.fill_between(np.linspace(xlim[0],xlim[1],len(h)), h, facecolor='lightblue')

	#save plot if a file path is given
	if save_to is not None:
		plt.savefig(save_to)

	#show plot if needed
	if show_plot:
		plt.show()
	#otherwise clear and close plot to save memory
	else:
		plt.clf()
		plt.cla()
		plt.close()



def plot_transport(a, b, bc_map, weight=0.5, labels=('Bins', 'Mass'), xlim=None, 
				   save_to=None, show_plot=True, plot_subtitle=None):
	'''
	Function that plots transport between a and b using some barycentric map

	a, b - histograms
	bc_map - barycentric mapping (n-sized ndarray)
	weight - weight given to barycentric mapping
	'''

	a /= np.sum(a)
	b /= np.sum(b)

	plt.figure(figsize=(16,9))
	plt.title(f'Transport of histogram $a$ to $b$ using barycentric mapping with weighting $w={weight:.2f}$')
	plt.xlabel(labels[0])
	plt.ylabel(labels[1])
	if xlim is None: plt.xticks([]); xlim=(0,1)
	plt.yticks([])

	plt.plot(np.linspace(0,1,len(a)), a, 'blue', linewidth=1, label='$a$')
	plt.plot(np.linspace(0,1,len(b)), b, 'green', linewidth=1, label='$b$')

	#calculate transport histogram:
	t = weight*a*bc_map + (1-weight)*a
	# t = weight*a + (1-weight)*b
	t /= np.sum(t)

	plt.plot(np.linspace(0,1,len(t)), t, 'black', linewidth=3, label='Transported histogram')

	plt.legend()

	#save plot if a file path is given
	if save_to is not None:
		plt.savefig(save_to)

	#show plot if needed
	if show_plot:
		plt.show()
	#otherwise clear and close plot to save memory
	else:
		plt.clf()
		plt.cla()
		plt.close()




def plot_results(res, labels=('Bins', 'Mass'), xlim=None, save_to=None, show_plot=True, plot_subtitle=None):
	#get results
	a, b = res.a, res.b
	P = res.P
	v, u = res.v, res.u
	error_a, error_b = res.error_a, res.error_b
	W = res.W
	bc_map = res.bc_map
	epsilon = res.epsilon
	converge_thresh = res.converge_thresh
	In = np.ones((len(a),1))
	Im = np.ones((len(b),1))

	### PLOTTING
	fig = plt.figure(figsize=(16,9))
	title = f'Sinkhorn algorithm errors and coupling matrix for two histograms with $\epsilon={epsilon:.3f}$\nWasserstein distance $W_\epsilon={W:.2f}$; Converged after {len(error_a)} iterations with threshold$={converge_thresh}$'
	if plot_subtitle is None:
		plt.suptitle(title)
	else:
		plt.suptitle(f'{title}\n{plot_subtitle}')

	
	#plot the histograms:
	#histogram a
	plt.subplot(2,2,1)
	plt.title(f'Histogram $a$ and $P_\epsilon a$, $n={len(a)}$')
	plt.xlabel(labels[0])
	plt.ylabel(labels[1])
	plt.plot(range(len(a)), a, 'blue', linewidth=2, label='$a$')
	plt.plot(range(len(a)), P@b/np.sum(P@b), 'red', linewidth=1, label='$P_\epsilon I_n$')
	plt.fill_between(range(len(a)), a, facecolor='lightblue')
	plt.yticks([])
	plt.xlim(0,len(a))
	plt.ylim(0,max(max(P@b/np.sum(P@b)), max(a))*1.1)
	plt.legend()

	#histogram b
	plt.subplot(2,2,2)
	plt.title(f'Histogram $b$ and $P_\epsilon^⊤ b$, $n={len(b)}$')
	plt.xlabel('Bin')
	plt.ylabel('Mass')
	plt.plot(range(len(b)), b, 'blue', linewidth=2, label='$b$')
	plt.plot(range(len(b)), P.T@a/np.sum(P.T@a), 'red', linewidth=1, label='$P_\epsilon^⊤ I_m$')
	plt.fill_between(range(len(b)), b, facecolor='lightblue')
	plt.yticks([])
	plt.xlim(0,len(b))
	plt.ylim(0,max(max(P.T@a/np.sum(P.T@a)), max(b))*1.1)
	plt.legend()

	#Plot error
	plt.subplot(2,2,3)
	plt.title('Log errors during algorithm')
	plt.xlabel('Iteration')
	plt.ylabel('$log_{10}(Err)$')
	plt.plot(range(len(error_a)), np.log(error_a)/np.log(10), label='$log_{10}(Err(P_\epsilon I  , a$))')
	plt.plot(range(len(error_b)), np.log(error_b)/np.log(10), label='$log_{10}(Err(P_\epsilon ^⊤I, b$))')
	plt.legend()

	#Plot coupling matrix
	plt.subplot(2,2,4)
	plt.title('Coupling matrix with barycentric map ($P_\epsilon$)')
	plt.xlabel('Bin in histogram b')
	plt.ylabel('Bin in histogram a')
	plt.imshow(np.log(P+1e-5), aspect='auto')

	#Plot baryocentric map:
	# plt.plot(bc_map, range(len(a)), 'r', linewidth=3)
	plt.gca().invert_yaxis()

	#save plot if a file path is given
	if save_to is not None:
		plt.savefig(save_to)

	#show plot if needed
	if show_plot:
		plt.show()
	#otherwise clear and close plot to save memory
	else:
		plt.clf()
		plt.cla()
		plt.close()