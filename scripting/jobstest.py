import modules.jobs as jobs
import modules.ir as ir
import modules.vcd as vcd
import modules.plot as plot
import modules.molecule_funcs as mf
import modules.sinkhorn_algorithm as sink
import scm.plams as plams
import os, time



# mol = mf.mirror('lactic acid')
# mf.save_to_xyz(mol, 'mirror lactic acid')


#prepare jobs
q = jobs.JobQueue()
mols = ['lactic acid', 'l-alanine']
for mol in mols:
	q.append(jobs.DFTJob(mol))
	q.append(jobs.DFTBJob(mol))

#run jobs and get results
results = q.run()

xlim = (0,2000)

#generate spectra
s = []
for r in results:
	s.append(ir.get_spectrum(r, xlim=(0,4000)))

# a = vcd.get_spectrum(results[0], xlim=xlim)
# b = vcd.get_spectrum(results[1], xlim=xlim)
# a = ir.get_spectrum(results[0], xlim=(0,4000))
# b = ir.get_spectrum(results[1], xlim=(0,4000))

plot.plot_hists((s[0], s[1]), show_plot=False, axislabels=('Frequency [$cm^{-1}$]', 'Rotational strength'), labels=('DFT', 'DFTB'), xlim=xlim, invert_xaxis=True, title=f'IR spectra of {mols[0]} using DFT and DFTB', save_to=fr'{results[0].dir}\IR_{mols[0]}.png')
plot.plot_hists((s[2], s[3]), show_plot=False, axislabels=('Frequency [$cm^{-1}$]', 'Rotational strength'), labels=('DFT', 'DFTB'), xlim=xlim, invert_xaxis=True, title=f'IR spectra of {mols[1]} using DFT and DFTB', save_to=fr'{results[0].dir}\IR_{mols[1]}.png')

# plot.plot_hist(a, show_plot=False, axislabels=('Frequency [$cm^{-1}$]', 'Rotational strength'), xlim=xlim, invert_xaxis=True, title=f'VCD spectrum of {mols[0]}', save_to=fr'{results[0].dir}\{mols[0]}.png')
# plot.plot_hist(b, show_plot=False, axislabels=('Frequency [$cm^{-1}$]', 'Rotational strength'), xlim=xlim, invert_xaxis=True, title=f'VCD spectrum of {mols[1]}', save_to=fr'{results[0].dir}\{mols[1]}.png')

# perform wasserstein distance and plot the spectra
res1 = sink.sinkhorn(s[0], s[1], epsilon=0.3)
res1 = sink.sinkhorn(s[2], s[3], epsilon=0.3)

plot.plot_results(res1, save_to=fr'{results[0].dir}\results_{mols[0]}.png')

plot.plot_results(res2, save_to=fr'{results[0].dir}\results_{mols[1]}.png')

