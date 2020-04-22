###########################################
#This file contains all the folders in the full Grimme set with some information
#This is still Legacy stuff from Arno's script but did not want to delete it yet.
#Might come in useful later.
###########################################

GMTKN55folders = [
'W4-11', 'G21EA', 'G21IP', 'DIPCS10','PA26','SIE4x4','ALKBDE10','YBDE18','AL2X6','HEAVYSB11','NBPRC','ALK8','RC21','G2RC','BH76RC','FH51','TAUT15','DC13','MB16-43','DARC','RSE43','BSR36','CDIE20','ISO34','ISOL24','C60ISO','PArel','BH76','BHPERI','BHDIV10','INV24','BHROT27','PX13','WCPT18','RG18','ADIM6','S22','S66','HEAVY28','WATER27','CARBHB12','PNICO23','HAL59','AHB21','CHB6','IL16','IDISP','ICONF','ACONF','Amino20x4i','PCONF21','MCONF','SCONF', 'UPU23', 'BUT14DIOL']


# al folders from Hobza test set
HobzaFolders = ['S66','S66x8','X40','X40x10','IonicHydrogenBonds','P26','L7']

# folders for which counterpoise correction should be used
CPfolders = ['L7','S66','S66_01','S66_02','S66_15172223','X40','HAL59','PNICO23','CARBHB12','S22','AHB21','CHB6','IL16C']

# folders with conformatioal energies --> easy to calculate, no BSSE
conformerfolders = [
'PX13','BHROT27','INV24','BHDIV10','PArel','C60ISO','ISOL24','ISO34','CDIE20','TAUT15','PA26','DIPCS10','G21IP','G21EA']

# folders where everything is more difficult
otherfolders = ['WATER27','HEAVY28','ADIM6','RG18','WCPT18','BHPERI10','BH76','RSE43','DARC','MB16-43','DC13','FH51','BH76RC','G2RC','RC21','ALK8','NBPRC','HEAVYSB11','AL2X6','YBDE18','SIE4x4','W4-11']

# folders involving monoatomic species --> special care needs to be taken for all calcultions
#
# --> switch off spherical symmetry  
# --> unrestricted calcultion
# --> enforce integer aufbau

integfolders = ['HG/HW6Cl','HG/HW6F','W4-11','G21EA','G21IP','DIPCS10','SIE4x4','RG18','ALKBDE10']

# heavy elements --> ZORA should be used

zorafolders = ['HEAVY28','HEAVYSB11','HAL59','X40','X40x10']


# legacy --> should not be needed any more. However, when something goes wrong, it cant hurt to 
#            add yout folder to this list (provided your molecules are nuetral)

neutralfolders = ['L7','PX13','BHROT27','S66','S66_01', 'S66_02','S66_15172223','X40','IL16','IL16C','X40x10','S66x8'] # folder where the second line is empty

# folders with large molecules --> when sos-mp2 is used, switch to AO-based algorithm

sosfolders = ['L7','C60ISO']   # folders with very laerge molecules


