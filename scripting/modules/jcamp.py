import numpy as np
import matplotlib.pyplot as plt
import ir



## ================================================================= ##
# jcamp
#functions to read .jdx spectroscopic data files


def read_jdx(file):
	'''
	Function that reads a jdx file and returns the data and the range
	of data and number of points.
	'''

	start_data = 100000
	with open(file, 'r') as f:
		for i, line in enumerate(f.readlines()):
			line = line.strip('\n')

			#gather meta info
			if line.startswith('##FIRSTX'):
				firstx = float(line.split('=')[1])
			if line.startswith('##LASTX'):
				lastx = float(line.split('=')[1])
			if line.startswith('##NPOINTS'):
				npoints = int(line.split('=')[1])
			if line.startswith('##YFACTOR'):
				yfactor = float(line.split('=')[1])

			#gather data
			if line.startswith('##XYDATA'):
				start_data = i+1
				data_array = np.empty(npoints)

			if i >= start_data:
				if line.startswith('##END'): break

				delta_i = i - start_data
				
				data = [float(d) for d in line.split()]
				for j, d in enumerate(data[1:]):
					data_array[delta_i*10+j] = d*yfactor

	return data_array, (firstx, lastx)




plt.title('Lactic acid IR spectrum')
data, ran = read_jdx(r"C:\Users\Yuman\Downloads\50-21-5-IR.jdx")
ir_dft = ir.get_spectrum_from_kf(r"C:\Users\Yuman\Desktop\Programmeren\bachelorproject\scripting\RUNS\#KFFiles\DFT\lactic acid.t21", n=len(data), xlim=ran)
ir_dftb = ir.get_spectrum_from_kf(r"C:\Users\Yuman\Desktop\Programmeren\bachelorproject\scripting\RUNS\#KFFiles\DFTB\lactic acid.rkf", n=len(data), xlim=ran)

plt.title('Butane IR spectrum')
data, ran = read_jdx(r"C:\Users\Yuman\Downloads\106-97-8-IR.jdx")
ir_dft = ir.get_spectrum_from_kf(r"C:\Users\Yuman\Desktop\Programmeren\bachelorproject\scripting\RUNS\#KFFiles\hydrocarbons\DFT\butane_DFT.t21", n=len(data), xlim=ran)
ir_dftb = ir.get_spectrum_from_kf(r"C:\Users\Yuman\Desktop\Programmeren\bachelorproject\scripting\RUNS\#KFFiles\hydrocarbons\DFTB\butane_DFTB.rkf", n=len(data), xlim=ran)

data /= np.sum(data)
ir_dft /= np.sum(ir_dft)
ir_dftb /= np.sum(ir_dftb)

plt.plot(np.linspace(*ran, len(data)), data, label='Experimental', linewidth=2)
plt.fill_between(np.linspace(*ran,len(data)), data, facecolor='lightblue')
plt.plot(np.linspace(*ran, len(data)), ir_dft, label='DFT (accurate)', linewidth=2)
plt.plot(np.linspace(*ran, len(data)), ir_dftb, label='DFTB', linewidth=2)
plt.yticks([])
plt.xlabel('$cm^-1$')

plt.legend()

plt.show()