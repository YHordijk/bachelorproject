import modules.ir as ir
import modules.plot as plot
import modules.sinkhorn_algorithm as sink
import modules.colour_maps as cmap
import os
import numpy as np
import modules.comp_funcs as comp_funcs
import matplotlib.pyplot as plt 
import moviepy.editor as mvp




a = r"D:\Users\Yuman\Desktop\Programmeren\bachelorproject\scripting\RUNS\#KFFiles\functionals\Aminoindan_CONF_4_LDA_DFT.t21"
b = r"D:\Users\Yuman\Desktop\Programmeren\bachelorproject\scripting\RUNS\#KFFiles\functionals\ISO34_E4_DFTB3_DFTB.rkf"

fa, ia = ir.get_freqs_intens(a)
fb, ib = ir.get_freqs_intens(b)

facta = fa.sum()
factb = fb.sum()

fa = fa/fa.sum()
fb = fb/fb.sum()



res = sink.sinkhorn(fa, fb)
plot.plot_hists([fa, fb, fa@res.P, fb@res.P.T], ['fa', 'fb', 'fa@res.P', 'fb@res.P.T'], axislabels=('Bin', 'Frequency'), scatter=10, line=True)

#transport fb
fb = fb@res.P.T * factb







# res = sink.sinkhorn(ia, ib)
# plot.plot_hists([ia, ib, ia@res.P, ib@res.P.T], ['ia', 'ib', 'ia@res.P', 'ib@res.P.T'], axislabels=('Bin', 'Intensity'), scatter=10, line=True)

#transport fb
# ib = ib@res.P.T




sa = ir.get_spectrum_from_freqs_intens(fb,ib)
sb = ir.get_spectrum_from_kf(b)

plot.plot_hists([sa, sb], ['sa', 'sb'])
