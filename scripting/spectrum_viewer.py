import numpy as np
import modules.histograms as hist
import modules.sinkhorn_algorithm as sink
import modules.plot as plot
import modules.barycenter as bc
import modules.ir as ir
import os

kf_path = r"C:\Users\Yuman\Desktop\Programmeren\bachelorproject\scripting\RUNS\#KFFiles\[21-05](13-05-2020)\\"
kfs = [kf_path + d for d in os.listdir(kf_path)]

for kf in kfs:
	spec = ir.get_spectrum_from_kf(kf, xlim=(0,4000), n=1000)
	plot.plot_hist(spec, title=kf.split('\\')[-1][:-4], axislabels=('Wavenumber ($cm^{-1}$)', 'Intensity'), xlim=(0,4000), invert_xaxis=True)

