import modules.vg.jobs as jobs

j = jobs.DFTJob('l-alanine')
res = j.run()
print(res)