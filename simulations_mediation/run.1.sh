# BASH SCRIPT FOR SIMULATION STUDY OF MEDIATION ANALYSIS USING MGREML
# AUTHOR: R. DE VLAMING
# DATE: October 11, 2021
# REQUIREMENTS:
#  - clone of MGREML repository in working directory
#  - PLINK 1.90 in working directory
#  - Python 3.x installed with relevant packages (numpy, pandas, tqdm, networkx, scipy)
# RUN: 1

# set input args (run number; sample size; no. of SNPs)
iRun=1
iN=20000
iM=20000
# 1. simulate data
python ./simulatedata.py ${iRun} ${iN} ${iM}
# 2. gzip simulated GRM
gzip ./run.${iRun}.grm
# 3. convert GRM to binary format using PLINK 1.90
# rel-cutoff should not drop any individuals, which is what we want for present purposes
# rel-cutoff only meant here to force PLINK to convert to binary GRM in the first place
./plink --grm-gz run.${iRun} --rel-cutoff 0.99 --make-grm-bin --out run.${iRun}.bin
# 4. remove redundant files
rm ./run.${iRun}.grm.gz
rm ./run.${iRun}.grm.id
# 5. perform MGREML analyses
python ./mgreml/mgreml.py --grm run.${iRun}.bin --pheno no_mediation.run.${iRun}.pheno.txt --covar run.${iRun}.covar.txt --mediation --out no_mediation.run.${iRun}
python ./mgreml/mgreml.py --grm run.${iRun}.bin --pheno partial_mediation.run.${iRun}.pheno.txt --covar run.${iRun}.covar.txt --mediation --out partial_mediation.run.${iRun}
python ./mgreml/mgreml.py --grm run.${iRun}.bin --pheno full_mediation.run.${iRun}.pheno.txt --covar run.${iRun}.covar.txt --mediation --out full_mediation.run.${iRun}

