######################
# This file contains (nearly) all the methods currently available in the ADF-GUI that can be used to test the Grimme benchmark set. 
# I have left out all dftb and MOPAC methods that do not support all atom types within the benchmark set.
# 
######################
test_dft = ['PBE', 'BP86']

gga_methods = ['BP86', 'PW91', 'mPW', 'PBE', 'RPBE', 'revPBE', 'mPBE', 'PBEsol', 'HTBS', 'BLYP', 'OLYP', 'OPBE', 'KT1', 'KT2', 'N12', 'GAM', 'HTBS']

meta_gga_methods = ['M06L', 'MN15-L', 'MVS', 'TPSS', 'TASKxc', 'SSB-D', 'MS0', 'MS1', 'MS2', 'revTPSS','SCAN']

meta_hybrid_methods = ['M06-2X','M06-HF','TPSSH','M06']

hybrid_methods = ['B3LYP', 'B3LYP*', 'B1LYP', 'KMLYP', 'O3LYP', 'X3LYP', 'BHandH', 'BHandHLYP', 'B1PW91', 'mPW1PW', 'mPW1K', 'PBE0', 'OPBE0', 'S12H']

MOPAC_methods = ['PM7', 'PM6', 'PM3', 'MNDO']

DFTB_methods = ['GFN1-xTB']

MP2_methods = ['MP2']

SOS_MP2_methods = ['SOS-MP2']

double_hybrid_methods = ['rev-DOD-PBEP86', 'rev-DOD-BLYP', 'rev-DOD-PBE', 'B2PLYP', 'B2GPPLYP'] #NOTE: There seem to be a lot more of them on the website, but not available in the GUI, decided to only add the ones in the GUI for now.

HF_methods = ['HartreeFock']

# LIBXC folder is legacy of Arno

libxclist = ['LDA','PW92','TETER93','AM05','BGCP','B97-GGA1','B97-K','BLYP','BP86','EDF1','GAM','HCTH-93','HCTH-120','HCTH-147','HCTH-407','HCTH-407P','HCTH-P14','PBEINT','HTBS','KT2','MOHLYP','MOHLYP2','MPBE','MPW','N12','OLYP','PBE','PBEINT','PBESOL','PW91','Q2D','SOGGA','SOGGA11','TH-FL','TH-FC','TH-FCFO','TH-FCO',
'TH1','TH2','TH3','TH4','XLYP','XPBE','HLE16','M06-L','M11-L','MN12-L','MS0','MS1','MS2','MVS','PKZB','TPSS','HLE17','B1LYP','B1PW91','B1WC','B3LYP','B3LYP*','B3LYP5','B3LYP5','B3P86','B3PW91','B97','B97-1',
'B97-2','B97-3','BHANDH','BHANDHLYP','EDF2','MB3LYP-RC04','MPW1K','MPW1PW','MPW3LYP','MPW3PW','MPWLYP1M','O3LYP','OPBE','PBE0','PBE0-13','REVB3LYP','REVPBE','RPBE','SB98-1A','SB98-1B','SB98-1C','SB98-2A','SB98-2B','SB98-2C','SOGGA11-X','SSB','SSB-D','X3LYP','B86B95','B88B95','BB1K','M05','M05-2X','M06','M06-2X','M06-HF','M08-HX','M08-SO','MPW1B95','MPWB1K','MS2H','MVSH','PW6B95','PW86B95','PWB6K','REVTPSSH','TPSSH','X1B95','XB1K','CAM-B3LYP','CAMY-B3LYP','HJS-PBE','HJS-PBESOL','HJS-B97X','HSE03','HSE06','LRC_WPBE','LRC_WPBEH','LCY-BLYP','LCY-PBE','M11','MN12-SX','N12-SX','TUNED-CAM-B3LYP','WB97','WB97X']
