import os 


for file in os.listdir('.'):

    # open only complex files 

    if (not 'mon' in file):
       name1 = file.split('_')[1]
       for file1 in os.listdir('.'):
           
            

       N = int(linesA[0].replace('\n','').strip()) +  int(linesB[0].replace('\n','').strip()) 

       with open(file.replace('A',''),'a+') as ffD:
          ffD.write(str(N) + '\n\n')
          for line in linesA[2:]:
              ffD.write(line) 
          for line in linesB[2:]:
              ffD.write(' ' + line)
                   
