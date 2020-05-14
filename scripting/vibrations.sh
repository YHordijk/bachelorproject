#!/bin/bash
#SBATCH --nodes 1
#SBATCH --ntasks 12
#SBATCH --partition defq
#SBATCH --time 24:00:00
#SBATCH --mem 60000

# ==============================================================  Comment Start  =========================================================================
# Step 1: SBATCH settings: run the job on 1 node with 12 cores on the 'defq' partition; limit the runtime to 24 hours and use up to 60 GB of memory.
# Step 2: Load ADF and create your own personal scratch directory, located on the node.
# Step 3: Define the jobname (i.e. the input filename minus its extension) and the cartesian coordinates of the input molecule.
# Step 4: Run the ADF job; the provided example is a single point of the H2 molecule at the BP86/TZ2P level.
# Step 5: Copy the resulting TAPE21 file back to its submission directory; copy the TAPE13 + logfile if TAPE21 does not exist (i.e. the job crashed).
# Step 6: Empty your own personal scratch directory, located on the node.
# ===============================================================  Comment End  ==========================================================================



module load adf/current


$ADFBIN/startpython $PWD/bazis.py > python_log.txt
