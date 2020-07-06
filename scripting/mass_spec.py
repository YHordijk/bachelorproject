import modules.massbank_api as mb
import modules.comp_funcs as comp_funcs
import matplotlib.pyplot as plt 
import os
import numpy as np


# peak_list = []
# for i in range(1,256):
# 	ID = f'TY{"0"*(6-len(str(i)))}{i}'
# 	try:
# 		peaks = mb.get_peaks(ID=ID)
# 		save_file = os.getcwd() + rf'\mass_data\{ID}.txt'
# 		mb.save_peaks(save_file, peaks)
# 	except:
# 		pass
	

# print(mb.load_peaks(r"D:\Users\Yuman\Desktop\Programmeren\bachelorproject\scripting\mass_data\TY000001.txt"))

mass_dir = rf"D:\Users\Yuman\Desktop\Programmeren\bachelorproject\scripting\mass_data\\"
peaks = []
IDs = []
for file in os.listdir(mass_dir):
	IDs.append(file.strip('.txt'))
	peaks.append(mb.load_peaks(mass_dir+file))

input = [mb.get_peaks('PR307425    ')]

d = comp_funcs.freq_int_wasserstein(peaksb=peaks, peaksa=input, int_weight=0, freq_weight=10000, int_exp=2, freq_exp=2)
plt.subplot(1,2,1)
plt.imshow(d)
b = np.zeros_like(d)
b[np.arange(len(d)), d.argmin(1)] = 1
plt.subplot(1,2,2)
plt.imshow(b)
# plt.show()


mol_index = b.argmax()
mol = IDs[mol_index]
print(mol, mb.get_mol(mol))
