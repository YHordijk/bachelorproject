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
mols = ['l-alanine', 'water', 'dichloromethane', 'chloroethane']
for mol in mols:
	print(mol)
	q.clear()
	q.append(jobs.DFTJob(mol, job_name=f'{mol}_DFT'))
	q.append(jobs.DFTBJob(mol, job_name=f'{mol}_DFTB'))

	#run jobs and get results
	results = q.run()

	xlim = (0,4000)

	#generate spectra
	irs = []
	vcds = []
	for r in results:
		irs.append(ir.get_spectrum(r, xlim=xlim))
		vcds.append(vcd.get_spectrum(r, xlim=xlim))


	# a = vcd.get_spectrum(results[0], xlim=xlim)
	# b = vcd.get_spectrum(results[1], xlim=xlim)
	# a = ir.get_spectrum(results[0], xlim=(0,4000))
	# b = ir.get_spectrum(results[1], xlim=(0,4000))

	plot.plot_hists((irs[0], irs[1]), show_plot=False, axislabels=('Frequency [$cm^{-1}$]', 'Rotational strength'), labels=('DFT', 'DFTB'), xlim=xlim, invert_xaxis=True, title=f'IR spectra of {mol} using DFT and DFTB', save_to=fr'{results[0].dir}\IR_{mol}.png')
	# plot.plot_hists((irs[2], irs[3]), show_plot=False, axislabels=('Frequency [$cm^{-1}$]', 'Rotational strength'), labels=('DFT', 'DFTB'), xlim=xlim, invert_xaxis=True, title=f'IR spectra of {mols[1]} using DFT and DFTB', save_to=fr'{results[0].dir}\IR_{mols[1]}.png')

	plot.plot_hists((vcds[0], vcds[1]), show_plot=False, axislabels=('Frequency [$cm^{-1}$]', 'Intensity'), labels=('DFT', 'DFTB'), xlim=xlim, invert_xaxis=True, title=f'VCD spectra of {mol} using DFT and DFTB', save_to=fr'{results[0].dir}\VCD_{mol}.png')
	# plot.plot_hists((vcds[2], vcds[3]), show_plot=False, axislabels=('Frequency [$cm^{-1}$]', 'Intensity'), labels=('DFT', 'DFTB'), xlim=xlim, invert_xaxis=True, title=f'VCD spectra of {mols[1]} using DFT and DFTB', save_to=fr'{results[0].dir}\VCD_{mols[1]}.png')

	# plot.plot_hist(a, show_plot=False, axislabels=('Frequency [$cm^{-1}$]', 'Rotational strength'), xlim=xlim, invert_xaxis=True, title=f'VCD spectrum of {mols[0]}', save_to=fr'{results[0].dir}\{mols[0]}.png')
	# plot.plot_hist(b, show_plot=False, axislabels=('Frequency [$cm^{-1}$]', 'Rotational strength'), xlim=xlim, invert_xaxis=True, title=f'VCD spectrum of {mols[1]}', save_to=fr'{results[0].dir}\{mols[1]}.png')

	# perform wasserstein distance and plot the spectra
	res1 = sink.sinkhorn(irs[0], irs[1], epsilon=0.3)
	# res2 = sink.sinkhorn(irs[2], irs[3], epsilon=0.3)
	plot.plot_results(res1, show_plot=False, save_to=fr'{results[0].dir}\results_ir_{mol}.png')
	# plot.plot_results(res2, save_to=fr'{results[0].dir}\results_ir_{mols[1]}.png')

	# res1 = sink.sinkhorn(vcds[0], vcds[1], epsilon=0.3)
	# # res2 = sink.sinkhorn(vcds[2], vcds[3], epsilon=0.3)
	# plot.plot_results(res1, save_to=fr'{results[0].dir}\results_vcd_{mol}.png')
	# # plot.plot_results(res2, save_to=fr'{results[0].dir}\results_vcd_{mols[1]}.png')