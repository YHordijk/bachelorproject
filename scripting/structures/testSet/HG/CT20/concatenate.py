import os 


for file in os.listdir('.'):
    if file.endswith('_monA_CT20.xyz'):
       with open(file,'r') as ffA:
          linesA = ffA.readlines()
       with open(file.replace('monA','monB'),'r') as ffB:
          linesB =ffB.readlines()

       N = int(linesA[0].replace('\n','').strip()) +  int(linesB[0].replace('\n','').strip()) 

       with open(file.replace('monA','dim'),'a+') as ffD:
          ffD.write(str(N) + '\n\n')
          for line in linesA[2:]:
              ffD.write(line) 
          for line in linesB[2:]:
              ffD.write(' ' + line)
                   
