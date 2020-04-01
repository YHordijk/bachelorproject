import modules.vg.jobs as jobs
import scm.plams as plams
import os, time


q = jobs.JobQueue()

for mol in ['l-alanine', 'water', 'benzene', 'ammonia']:
	q.append(jobs.DFTJob(mol))

q.run()
