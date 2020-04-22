# Import required packages
import glob,os,os.path,sys
import numpy as np
import importlib
import yaml

from optparse import OptionParser
from scm.plams import *

# Import libraries
import method_library as ml
import folder_library as fl

# YAML : PARSER OPTIONS
parser = OptionParser()
parser.add_option('--NFinal', type=int, default=150,
                  help="Number systems to select")
parser.add_option('--Mode', type='string', default='Analyse',
                  help="Number systems to select")
(Opts, args) = parser.parse_args()

# constants 
UnitconversionFactor = 627.509 # Ha -> kcal/mol

allefile   = 'energySummary.txt'

#
#
# Change these to yout own path 
#
#

# Path to reference values
combopath = '/home/boelrijk/benchmark_project/benchmarker/' 
# Path to DIET test set
testpath = '/home/boelrijk/benchmark_project/benchmarker/GNTKM150_noCP/'   
resultpath = '/home/boelrijk/benchmark_project/benchmarker/resultfiles/'
##########################################################################
# ALL METHODS FOR WHICH TO BENCHMARK
##########################################################################

# These are the parameters over which the benchmarker will loop.
# Fill up the lists below the way you like, have a look in method_library.py to see which ones 
# you can pick at the moment.

#modes = ['MOPAC', 'DFT', 'DFTB'] # Think about adding UFF in future
mode = 'DFTB'      # currently xTB is hardcoded in the routine

#folders = ['small']
folders = ['all']

basis =['']. # options = ['SZ', 'DZ']#, 'DZP'] #, 'TZP', 'TZ2P', 'QZ4P']
fitsets = ['']
fitsets = ['Normal']
frozencores = ['']
#frozencores = ['Large']

#frozencores = ['None', 'Small', 'Large']

methods = ml.DFTB_methods

#displist = ['Grimme4']


init(folder=mode) # Fill in folder here yourself depending on what you are benchmarking
#load_all('plams_workdir.004/') # Load from previous directory if you have already done some calculations

##################################################################################################
# PART 0 -- UTILS FOR EVALUATION OF DATA
##################################################################################################

def ParseLine(line,functional,folder):
    """
    Parses a line from the inputfile and returns the ID
    """ 
    ID = ''

    # determine functional
    #functional = functional.strip('  e_dh_nodisp')
                      
    x = line.split('/')[-1].replace('.t21','')#.replace('cp','')
    ID = x.split(functional)[0]
    return ID

def ExtractEnergies(FileName,f,folder):
    """
    Extracts Energies from a .csv file, with the format 
    # header
    name of reaction, res1, res2, res3, ...
    """

    # get headers and number of columns

    allE = {}
    F=open(os.path.join(resultpath, FileName))
    firstline = F.readline().strip()
    F.close() 
    nColumns  = firstline.count(',') + 1
    columnNames = list(firstline.split(','))

    # parse the file and get ID with energies. 
    for i in range(1,nColumns):
        E = {}
        F=open(os.path.join(resultpath, FileName))
        for L in F:
            if L[0]=="#":
                continue 
            T=L.split(',')
            if len(T)>1: 
                ID=ParseLine(T[0], f, folder)
                En=float(T[i].strip('\n'))
                E[ID]=En*UnitconversionFactor
            else:
                print ("Warning! %s not defined"%(T[0]))
        allE[columnNames[i]] = E

    return allE,columnNames

###########################################################################################################

def ReadSystems(folder,functional, basisset, fitset, fc, badvalue,
                NFinal=150, InputFile=None, ReferenceFile=None, resultfile=None, csv_resultfile=None, mode='ref'):
    """
    compares every reaction energy to the tested energy, calculates the deviation
    and calculates the weighted mean absolute deviation (WTMAD2) in the end for all methods.
    The reference data is stored in a yaml file.
    """
        
    WTMAD2 = {}
    Y=yaml.load(open(ReferenceFile))#, Loader=yaml.FullLoader)
    file = open(os.path.join(resultpath, resultfile), 'w')
    csv_file = open(os.path.join(resultpath, csv_resultfile), 'w')
    csv_file.write('#runinfo  type of calculation  WTMAD2  Max error \n')
    EnergiesAll,nEnergies = ExtractEnergies(InputFile,functional,folder)
    for n in nEnergies:    
       if (n == '# filename'): continue            
       Energies = EnergiesAll[n]
       WarningList = []
       MAE, MAEDen = 0., 0.
       Errors = {}
       file.write('type of calculation = ' + str(n) + '\n')
       file.write(functional + ' ' + basisset + ' ' + fitset + ' ' + fc +  '\nSystem                                      ' + 'Weight  ' +
                  'my En.  ' + 'ref. En.' + ' signed diff' + '  unsigned diff')
       energydiff_old = 0
       for System in Y:
          if ( mode == 'ref'): 
             Energy = Y[System]['Energy']
          else:
             Energy = Y[System]['EnergyMP2']
          Weight = Y[System]['Weight']  
          EnergyError = False
          EnergyApprox = 0. 
          for S in Y[System]['Species']: 
             if not S['ID'] in Energies:
                EnergyError = True
             else:
                EnergyApprox +=Energies[S['ID']] *float(S['Count'])
            
          # print out warning if not all single point energies for a given system are available
          # also print warning if some energy difference is way to high (controlled by badvalue)
            
          # also compute maximum (signed) error. save old value here

          EnergyDiff = EnergyApprox - Energy
          if (abs(EnergyDiff) > abs(energydiff_old)): energydiff_old = EnergyDiff
          if EnergyError:
             file.write("\n Skipping system due to absence: %s [Weight: %.2f]"\
                  %(System, Weight))
             WarningList += [S['ID'] for S in Y[System]['Species']]
          elif abs( EnergyDiff * Weight)> badvalue:
             ErrList = [S['ID'] or S in Y[System]['Species']]
             file.write( "\n Major error: %s - %.2f vs %.2f [Weight: %.2f]:"\
                            %(System, EnergyApprox, Energy, Weight))
             file.write("\n Species: "+", ".join(ErrList))
             file.write("\n Energies:"+", ".join(
                         ["%.2f"%(Energies[x]) for x in ErrList] ))
                     
             WarningList += ErrList
          elif not(EnergyError):
             file.write("\n %-50s %8.2f %8.2f %8.2f %8.2f %8.2f"%(System, Weight,EnergyApprox, Energy,abs( EnergyDiff), EnergyDiff))
             Errors[System] = EnergyDiff
             MAE += abs( EnergyDiff ) * Weight
             MAEDen += 1. 
            
       file.write("\n-----------------------------------------------------------------------------------------------------------")
       if len(WarningList)>0: 
          file.write("\n The following systems probably have errors:")
          file.write("\n - [ %s ]"%( " ".join("%s"%W for W in WarningList)))       
       if MAEDen != 0 : 
          file.write("\n NActual = %3d, WTMAD2 = %.3f"\
                        %(len(list(Errors)), MAE/MAEDen))
          csv_file.write(functional + str(basisset) + str(fitset)  + str(fc) + ',' +  str(n) + ',' + str(MAE/MAEDen) + ',' +  str(abs(energydiff_old)) + '\n')     
       else:
          file.write('total fail ...') 
       file.write("\n  MAX = %.3f"%(energydiff_old))
       file.write("\n SMAX = %.3f"%(abs(energydiff_old)))
       file.write('\n\n\n')
            
       # summarize all results 
            
       if (MAEDen != 0): WTMAD2[functional + ' ' + n]=MAE/MAEDen

    return WTMAD2     
                    
    # and print out summary in the end
   
def add_charge_and_multiplicity(mol,name,folder,lastfolder):
   """ Looks up the charge and multiplicity from second line in the .xyz files.
       Can also read from a CHARGE_MULTIPLICITY file if needs be.
   """
   n1,n2 = 0,0
   if (not os.path.isfile(os.path.join(folder,'CHARGE_MULTIPLICITY_' + lastfolder + '.txt'))):
      with open(os.path.join(folder,name +'.xyz') ) as fp:
         lines = fp.readlines()
         secondLine = lines[1]
         
         if (',' in secondLine): 
            char = ','
         else:
            char = ' '
         n1 = int(secondLine.split(char)[0].strip())
         n2 = int(secondLine.split(char)[1].strip()) - 1

   else:
      with open(os.path.join(folder,'CHARGE_MULTIPLICITY_' + lastfolder + '.txt')) as ff:
         lines = ff.readlines()    
         for line in lines:
             if (name == line.split(' ')[0].strip()):
                n1 = int(line.split(' ')[1].strip())
                n2 = int(line.split(' ')[2].strip()) - 1
                break
   if mode == 'DFTB':
       return n1
   return n1,n2


def add_charge(mol,name,folder,lastfolder):
   "Essentially identical to add_charge_and_multiplicity but only returns charge, for DFTB."
   n1 = 0
   if (not os.path.isfile(os.path.join(folder,'CHARGE_MULTIPLICITY_' + lastfolder + '.txt'))):
      with open(os.path.join(folder,name +'.xyz') ) as fp:
         lines = fp.readlines()
         secondLine = lines[1]

         if (',' in secondLine):
            char = ','
         else:
            char = ' '
         n1 = int(secondLine.split(char)[0].strip())
   else:
      with open(os.path.join(folder,'CHARGE_MULTIPLICITY_' + lastfolder + '.txt')) as ff:
         lines = ff.readlines()
         for line in lines:
             if (name == line.split(' ')[0].strip()):
                n1 = int(line.split(' ')[1].strip())
                break
   return n1

#############################################################################################################
# For if you don't have it in your PLAMS version yet.
#def get_timings(self):
#    """get_timings()
#    This functions is for AMSJob results only. ADFJob already has a built in function for this.
#    Return a dictionary with timing statistics of the job execution. Returned dictionary contains keys ``cpu``, ``system`` and ``elapsed``. The values are corresponding timings, expressed in seconds.
#    """
#    ret = {}
#    cpu = self.grep_output('Total cpu time:')
#    system = self.grep_output('Total system time:')
#    elapsed = self.grep_output('Total elapsed time:')
#    ret['elapsed'] = float(elapsed[0].split()[-1])
#    ret['system'] = float(system[0].split()[-1])
#    ret['cpu'] = float(cpu[0].split()[-1])
#    return ret

def runCalculationsMOPAC(folder, method, outputfile, timingfile):
   """
   performs a PLAMS run for MOPAC calculations, by creating all the jobs, adding all the jobs to a MultiJob and running the MultiJob.
   """
   # all settings 

   s = Settings()
   s.runscript.nproc = 1
   s.input.ams.task = 'SinglePoint'
   s.input.MOPAC.Model = method

   # create all jobs
   methodstring = method+folder

   molecules = read_molecules(os.path.join(testpath,folder),formats=['xyz'])
   names = []
   jobs = []
   for name in molecules:
      mol = molecules[name]
      s.input.ams.system.Charge = '0'
      s.input.MOPAC.UnpairedElectrons = '0'
      #if (folder not in neutralfolders):
      n1,n2 = add_charge_and_multiplicity(mol,name,os.path.join(testpath,folder),folder)
      if (n1 != 0 or n2 != 0):
         s.input.MOPAC.UnpairedElectrons = str(n2)
         s.input.ams.system.Charge = str(n1)

      # create a single job and eppend it to the joblist
      job1 = AMSJob(molecule=mol, name=name+methodstring, settings=s)
      names.append(name)
      jobs.append(job1)

   mj = MultiJob(name = method, children = jobs)

   writefile = open(os.path.join(resultpath, outputfile),'w')
   writefile2 = open(os.path.join(resultpath, timingfile), 'w')

   j = JobRunner(parallel=True, maxjobs=24)
   mj.run(jobrunner=j)
   mj.results.wait()
   writefile.write('# filename, e_dh_nodisp\n')
   writefile2.write('# filename, runtime\n')

   for job in jobs:
      bond_energy = job.results.get_energy()
      #job.results.wait()
      runtime = job.results.get_timings()
      #if ((method in dispmethodlist or method in d0methodlist) and not method == 'SDSCAN69'):
      #disp_energy = job.results.readkf('Energy', 'Dispersion Energy')
      #diff = bond_energy
      #if ((method in dispmethodlist or method in d0methodlist) and not method == 'SDSCAN69'):
      writefile.write(job.name + ',' + str(bond_energy)+'\n')  #+ ',' + str(diff) + '\n')
      writefile2.write(job.name + ',' + str(runtime) + '\n')
   writefile.close()
   writefile2.close()

def runCalculationsXTB(folder, method, outputfile, timingfile): 
   """
   performs a PLAMS run for DFTB
   """
   # all settings 

   s = Settings()
   s.runscript.nproc = 1
   s.input.ams.task = 'SinglePoint'
   s.input.DFTB.Model = 'GFN1-xTB'

   methodstring = method+folder

   molecules = read_molecules(os.path.join(testpath,folder),formats=['xyz'])
   names = []
   jobs = []
   for name in molecules:
      mol = molecules[name]
      s.input.ams.system.Charge = '0'
      #if (folder not in neutralfolders):
      n1 = add_charge(mol,name,os.path.join(testpath,folder),folder)
      if (n1 != 0 ):
         s.input.ams.system.Charge = str(n1)

      # create a single job and eppend it to the joblist
      job1 = AMSJob(molecule=mol, name=name+methodstring, settings=s)
      names.append(name)
      jobs.append(job1)

   mj = MultiJob(name = method, children = jobs)

   writefile = open(os.path.join(resultpath, outputfile),'w')
   writefile2 = open(os.path.join(resultpath, timingfile), 'w')

   j = JobRunner(parallel=True, maxjobs=24)
   mj.run(jobrunner=j)
   mj.results.wait()
   writefile.write('# filename, e_dh_nodisp\n')
   writefile2.write('# filename, runtime\n')

   for job in jobs:
      bond_energy = job.results.get_energy()
      runtime = job.results.get_timings()
      writefile.write(job.name + ',' + str(bond_energy)+'\n')
      writefile2.write(job.name + ',' + str(runtime) + '\n')
   writefile.close()
   writefile2.close()


def runCalculationsDFT(folder, method, bas, fit, fc, outputfile, timingfile):
   """
   performs a PLAMS run
   """
   # all settings 
   s = Settings()
   s.runscript.nproc = 1
   s.input.Basis.Core = fc
   s.input.NumericalQuality = fit
   if (method in ml.meta_gga_methods): 
      s.input.XC.metaGGA = method
   elif (method in ml.gga_methods):
      s.input.XC.GGA = method
   elif (method in ml.meta_hybrid_methods):
      s.input.XC.metahybrid = method
   elif (method in ml.hybrid_methods):
      s.input.XC.HYBRID = method
   elif (method in ml.double_hybrid_methods):
      s.input.XC.doublehybrid = method
   elif (method in ml.libxclist):
      s.input.XC.libXC = method
   elif (method in ml.MP2_methods):
      s.input.XC.MP2 = 'True'
      s.input.MP2.Formalism = 'Auto' 
   elif (method in ml.SOS_MP2_methods):
      s.input.XC.MP2 = 'True'
      s.input.XC.MP2.Formalism = 'Auto'
      s.input.XC.EmpiricalScaling = 'SOS'
   elif (method in ml.HF_methods): 
      s.input.XC.HartreeFock = 'True'

   s.input.Basis.Type = bas
   #s.input.XC.Dispersion = disp
   s.input.Symmetry = 'nosym'
   s.input.occupations = 'IntegerAufbau'
   
   #This is legacy folder stuff from Arno, might be useful for full set.
   #if (folder == 'CYCONF' or folder=='PCONF10'):
   #   s.input.UNITS.length = 'BOHR'
   #if (folder in integfolders):
   #   s.input.unrestricted = 'True'
   #if (folder in zorafolders):
   #   s.input.relativistic = 'ZORA MAPA'
   #   s.input.Basis.Type    = 'ZORA/' + bas  
   #if (folder in sosfolders and method in sosmethods):
   #   s.input.MP2.LT = True
   #   s.input.MP2.NLAPLACE = '10'  
 
   # create all jobs
   methodstring = method+bas+fit+fc+folder

   molecules = read_molecules(os.path.join(testpath,folder),formats=['xyz']) 
   names = []
   jobs = []
   for name in molecules:
      mol = molecules[name]  
      s.input.unrestricted = 'False'
      s.input.Charge = '0 0'
      #if (folder not in neutralfolders):
      n1,n2 = add_charge_and_multiplicity(mol,name,os.path.join(testpath,folder),folder)  
      if (n1 != 0 or n2 != 0):
         s.input.unrestricted = 'True'
         s.input.Charge = str(n1) + ' ' + str(n2) 

      # create a single job and eppend it to the joblist
      job1 = ADFJob(molecule=mol, name=name+methodstring, settings=s)
      names.append(name)
      jobs.append(job1)

   mj = MultiJob(name = method, children = jobs)

   writefile = open(os.path.join(resultpath, outputfile),'w')
   writefile2 = open(os.path.join(resultpath, timingfile), 'w') 

   j = JobRunner(parallel=True, maxjobs=24)
   mj.run(jobrunner=j)
   mj.results.wait()
   writefile.write('# filename, e_dh_nodisp\n')
   writefile2.write('# filename, runtime\n')
   
   for job in jobs:
      bond_energy = job.results.get_energy()
      runtime = job.results.get_timings()
      #if ((method in dispmethodlist or method in d0methodlist) and not method == 'SDSCAN69'):
      # Could be added in if you want to do dispersion.
      #disp_energy = job.results.readkf('Energy', 'Dispersion Energy')
      #diff = bond_energy - disp_energy
      #if ((method in dispmethodlist or method in d0methodlist) and not method == 'SDSCAN69'):
      writefile.write(job.name + ',' + str(bond_energy)+'\n')
      writefile2.write(job.name + ',' + str(runtime) + '\n')
   writefile.close()
   writefile2.close()

###########################################################################################################
### Probably best to make Calculations function return the MultiJob and then run it from here ####
### This also allows for recreation of the MultiJob for checking Results ###
### Do this tommorow ###


for folder in folders:
   f = open(resultpath + '/'  + allefile.replace('.txt','') + '_' + folder.replace('HG/','') + mode + '.txt','w')
   f.write(' WTMADs for all tested functionals:\n\n' )
   f.write('\n Method                          WTMADs\n -------------------------------------------')
   for m in methods:
      for basisset in basis:
          for fitset in fitsets:
              for frozencore in frozencores:
                  method = m
                  f.write('\n %s %s %s %s %s' % (folder, method, basisset, fitset, frozencore))
                  print('current method:   ', method)
                  print(folder,method,basisset,fitset, frozencore)

                  outputfile = folder + '_' + method + '_' + basisset + '_' + fitset + '_' + frozencore + '_' + 'file.txt'
                  timingfile = folder + '_' + method + '_' + basisset + '_' + fitset + '_' + frozencore + '_' + 'timing.txt'
                  resultfile = folder + '_' + method + '_' + basisset + '_' + fitset + '_' + frozencore + '_' + 'results.txt'
                  csv_resultfile = folder + '_' + method + '_' + basisset + '_' + fitset + '_' + frozencore + '_' + 'csv_results.txt'

##########################################################################################################Y
# PART 1 -- RUN THE ACTUAL CALCULATION 
###########################################################################################################
                  if mode == 'MOPAC':
                     runCalculationsMOPAC(folder, method, outputfile, timingfile)
                  if mode == 'DFTB':
                     runCalculationsXTB(folder,method,outputfile, timingfile)
                  if mode == 'DFT':
                     runCalculationsDFT(folder, method, basisset, fitset, frozencore,  outputfile, timingfile)
                  else:
                     print('the mode you have selected cannot be run by this script') 

############################################################################################################
# PART 2 -- EVALUATE DATA --> GET MAD (RESULTS ARE WRITTEN TO OUTPUTFILE
############################################################################################################

                  combofile = os.path.join(combopath,'ComboList_'+folder+'.txt')
                  badvalue = 100000000 #set unreasonably big for now to obtain all energy diffs                            
                  WTMADarray = ReadSystems(folder, method, basisset, fitset, frozencore,badvalue, 34, 
                                             InputFile = outputfile, ReferenceFile=combofile, 
                                             resultfile = resultfile, csv_resultfile = csv_resultfile, mode='ref')
                       
                  if outputfile == None:
                      continue
                  for method, energy in WTMADarray.items():
                      f.write(' \n %-38s %8.2f,'%(method, energy))
 
#############################################################################################################
   f.close()
