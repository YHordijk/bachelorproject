import modules.jobs as jobs
import modules.ir as ir
import modules.plot as plot
import modules.sinkhorn_algorithm as sink
import scm.plams as plams
import os, time

#prepare jobs
q = jobs.JobQueue()

for mol in ['water', 'l-alanine']:
	q.append(jobs.DFTJob(mol))

#run jobs and get results
results = q.run()


#generate ir spectra
a = ir.get_spectrum(results[0])
b = ir.get_spectrum(results[1])


#perform wasserstein distance and plot the spectra
res = sink.sinkhorn(a, b, epsilon=0.4)
plot.plot(res)

