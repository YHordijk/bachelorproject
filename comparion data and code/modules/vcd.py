import numpy as np
import scm.plams as plams
import modules.histograms as hist



## ================================================================= ##
# VCD
#classes and functions used to calculate VCD spectra based on results
#from ADF runs




def get_freqs_intens(result):
	'''
	Function to get frequencies and intensities from result objects

	results - plams.Results object
	'''

	#load kffile
	kff = plams.KFFile(result.KFPATH)
	#read freqs and intens from kffile
	freqs = kff.read_section('Vibrations')['Frequencies[cm-1]']
	if result.KFPATH.endswith('.rkf'):
		rotstr = kff.read_section('Vibrations')['RotationalStrength']
	else:
		rotstr = kff.read_section('Freq Symmetry')['VCD rotational strength_A']

	return freqs, rotstr



def get_spectrum(result, n=600, width=8, xlim=(0,4000)):
	'''
	Function to generate a spectrum histogram from results using
	lorentzian line shape.

	results - plams.Results object
	n - size of spectrum
	width - width of lorentzian peaks, defaults to 8 cm-1 (same as in ADF)

	I will try to add support for gaussian line shape as well
	'''

	freqs, rotstr = get_freqs_intens(result)
	spectrum = np.zeros(n)

	#project width onto new range
	width = (width-xlim[0])/(xlim[1]-xlim[0])

	for f, r in zip(freqs, rotstr):
		#we must project f from [xlim[0], xlim[1]] onto the range [0,1] since the histograms 
		#are calculated on that range
		spectrum += hist.lorentzian(n, (f-xlim[0])/(xlim[1]-xlim[0]), width) * r

	return spectrum