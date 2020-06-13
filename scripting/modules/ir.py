import numpy as np
import scm.plams as plams
try:
	import modules.histograms as hist
except:
	import histograms as hist



## ================================================================= ##
# IR
#classes and functions used to calculate IR spectra based on results
#from ADF runs




def get_freqs_intens(kf):
	'''
	Function to get frequencies and intensities from result objects

	kf - path to KF file
	'''

	#load kffile
	kff = plams.KFFile(kf)
	#read freqs and intens from kffile

	try:
		freqs = kff.read_section('Freq Symmetry')['Frequencies_A']
		# intens = kff.read_section('Freq Symmetry')['IR intensities_A']
		intens = [1 for _ in freqs]
	except:
		freqs = kff.read_section('Vibrations')['Frequencies[cm-1]']
		# intens = kff.read_section('Vibrations')['Intensities[km/mol]']
		intens = [1 for _ in freqs]

	return freqs, intens


def get_spectrum(result, n=600, xlim=(0,4000), width=50):
	'''
	Function to generate a spectrum histogram from results using
	lorentzian line shape.

	results - plams.Results object
	n - size of spectrum
	width - width of lorentzian peaks, defaults to 50cm-1 (same as in ADF)

	I will try to add support for gaussian line shape as well


	##DEPRECATED##
	'''

	freqs, intens = get_freqs_intens(result.KFPATH)
	spectrum = np.zeros(n)

	#project width onto new range
	width = width/(xlim[1]-xlim[0])

	for f, i in zip(freqs, intens):
		#we must project f from [xlim[0], xlim[1]] onto the range [0,1] since the histograms 
		#are calculated on that range
		spectrum += hist.lorentzian(n, (f-xlim[0])/(xlim[1]-xlim[0]), width) * i

	return spectrum


def get_spectrum_from_kf(kf, n=600, xlim=(0,4000), width=50):
	'''
	Same function as above but for paths to kf files.
	'''
	
	freqs, intens = get_freqs_intens(kf)
	spectrum = np.zeros(n)

	#project width onto new range
	width = width/(xlim[1]-xlim[0])

	for f, i in zip(freqs, intens):
		#we must project f from [xlim[0], xlim[1]] onto the range [0,1] since the histograms 
		#are calculated on that range
		spectrum += hist.lorentzian(n, (f-xlim[0])/(xlim[1]-xlim[0]), width) * i
		# spectrum += hist.dirac_delta(n, (f-xlim[0])/(xlim[1]-xlim[0])) * i

	return spectrum