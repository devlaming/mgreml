###########################################################
# BASH SCRIPT FOR SIMULATION STUDY IN MANUSCRIPT ENTITLED #
# "Multivariate Estimation of Factor Structures of        #
#  Complex Traits Using SNP-based Genomic Relationships"  #
#                                                         #
# AUTHOR: R. DE VLAMING                                   #
# REQUIREMENTS:                                           #
#  - clone of MGREML repository in working directory      #
#  - PLINK 1.90 in working directory                      #
#  - Python 3.x installed with relevant packages          #
#    (numpy, pandas, tqdm, networkx, scipy)               #
# RUN: 1                                                  #
###########################################################

###########################################################
# SIMULATION SETTINGS 1 AND 2                             #
###########################################################

# set input args (run number; sample size; no. of SNPs; no. of traits)
iRun=1
iN=20000
iM=20000
iT=10
# 1. simulate data
python ./simulatedata.py ${iRun} ${iN} ${iM} ${iT}
# 2. gzip simulated GRM
gzip ./run${iRun}.grm
# 3. convert GRM to binary format using PLINK 1.90
# rel-cutoff should not drop any individuals, which is what we want for present purposes
# rel-cutoff only meant here to force PLINK to convert to binary GRM in the first place
./plink --grm-gz run${iRun} --rel-cutoff 0.99 --make-grm-bin --out run${iRun}bin
# 4. remove redundant files
rm ./run${iRun}.grm.gz
rm ./run${iRun}.grm.id
# 5. perform MGREML analyses
# notice: --no-intercept as there is nothing to control for in our model
python ./mgreml/mgreml.py --grm run${iRun}bin --pheno rhoGneg.run${iRun}.pheno.txt --no-intercept --out rhoGneg.run${iRun}.results
python ./mgreml/mgreml.py --grm run${iRun}bin --pheno rhoGzero.run${iRun}.pheno.txt --no-intercept --out rhoGzero.run${iRun}.results
python ./mgreml/mgreml.py --grm run${iRun}bin --pheno rhoGhalf.run${iRun}.pheno.txt --no-intercept --out rhoGhalf.run${iRun}.results
python ./mgreml/mgreml.py --grm run${iRun}bin --pheno rhoGone.run${iRun}.pheno.txt --no-intercept --out rhoGone.run${iRun}.results
python ./mgreml/mgreml.py --grm run${iRun}bin --pheno blocks.run${iRun}.pheno.txt --no-intercept --genetic-model factors.gen.unres.txt --environment-model factors.env.unres.txt --restricted-genetic-model factors.gen.res.txt --restricted-environment-model factors.env.res.txt --out blocks.run${iRun}.results

###########################################################
# SIMULATION SETTING 3                                    #
###########################################################

# set input args (run number; sample size; no. of SNPs; no. of traits)
iRun=1
iN=20000
iM=20000
iT=50
# 1. simulate data
python ./simulatedata.py ${iRun} ${iN} ${iM} ${iT}
# 2. gzip simulated GRM
gzip ./run${iRun}.grm
# 3. convert GRM to binary format using PLINK 1.90
# rel-cutoff should not drop any individuals, which is what we want for present purposes
# rel-cutoff only meant here to force PLINK to convert to binary GRM in the first place
./plink --grm-gz run${iRun} --rel-cutoff 0.99 --make-grm-bin --out run${iRun}bin
# 4. remove redundant files
rm ./run${iRun}.grm.gz
rm ./run${iRun}.grm.id
# 5. perform MGREML analysis
python ./mgreml/mgreml.py --grm run${iRun}bin --pheno rhoGzero.run${iRun}.pheno.txt --restricted-rho-genetic 0 --out rhoGzero.run${iRun}.runtime.comparison.results
