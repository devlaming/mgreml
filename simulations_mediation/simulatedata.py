''' 
simulatedata.py: simulates data for simulation study
of mediation analysis using MGREML
last edit: October 11, 2021
'''
import numpy as np
import pandas as pd
import sys
from numpy.matlib import repmat

def SimulateGenotypes(iN,iM,rng):
    # settings for distr. for drawing SNP allele frequencies
    dBetaParam1 = 0.35
    dBetaParam2 = 0.35
    dMinMAF = 0.05
    # set ploidy of human genome
    iP = 2   
    # initialise allele frequencies
    vAF = rng.beta(dBetaParam1,dBetaParam2,size=iM)
    # keep drawing allele frequencies until all between dMinMAF and 1-dMinMAF
    while (min(vAF) < dMinMAF) or (max(vAF) > (1-dMinMAF)):
        (vIndAFlow,) = np.where(vAF < dMinMAF)
        (vIndAFhigh,) = np.where(vAF > (1-dMinMAF))
        vAF[vIndAFlow] = rng.beta(dBetaParam1,dBetaParam2,size=vIndAFlow.shape)
        vAF[vIndAFhigh] = rng.beta(dBetaParam1,dBetaParam2,size=vIndAFhigh.shape)
    # initialise matrix of genotypes
    mG = np.zeros((iN,iM))
    # for each SNP
    for m in range(0,iM):
        # draw SNP data as binomially distributed according to AF of that SNP
        mG[:,m] = rng.binomial(iP,vAF[m],size=iN)
    # standardise the SNP data by semi-empirical allele frequency
    vEAF = mG.mean(axis=0)/2
    mG = (mG - repmat(2*vEAF,iN,1))*repmat((2*vEAF*(1-vEAF))**(-0.5),iN,1)
    # construct GRM
    mA = (mG@mG.T)/iM
    # return GRM and genotype matrix
    return mG, mA
    
def WriteGRM(mA,iM,lFID,lIID,sPrefix):
    # write GRM to txt-based format
    sFileAID = sPrefix + 'grm.id'
    # get sample size
    iN = mA.shape[0]
    with open(sFileAID, 'w') as oFile:
        for i in range(0,iN):
            oFile.write(lFID[i] + '\t' + lIID[i] + '\n')
    sFileA = sPrefix + 'grm'
    with open(sFileA, 'w') as oFile:
        for i in range(0,iN):
            if i%500 == 0:
                print('now writing grm: at observation ' + str(i) + ' out of ' + str(iN))
            for j in range(0,i+1):
                oFile.write(str(i+1) + '\t' + str(j+1) + '\t' + str(iM) + '\t' + str(mA[i,j]) + '\n')

def SimulatePhenotypes(mG,miID,iRun,rng):
    # set total no. of confounders (including intercept)
    iK = 10
    # set model parameters
    dGM = 3
    dGSM = 4
    dESM = 5
    dB = 1
    dGY = 4
    dEY = 4
    # get no. of SNPs and sample size
    iN = mG.shape[0]
    iM = mG.shape[1]
    # construct factors
    vG = (mG@rng.normal(0,1,(iM,1)))*(iM**(-0.5))
    vGS = (mG@rng.normal(0,1,(iM,1)))*(iM**(-0.5))
    vE = rng.normal(0,1,(iN,1))
    vES = rng.normal(0,1,(iN,1))
    # standardise factors
    vG = (vG - vG.mean())/(vG.std())
    vGS = (vGS - vGS.mean())/(vGS.std())
    vE = (vE - vE.mean())/(vE.std())
    vES = (vES - vES.mean())/(vES.std())
    # construct matrix of confounders
    mX = rng.normal(0,1,(iN,iK))
    # change first column to intercept
    mX[:,0] = 1
    # draw effects of confounders
    vBetaM = rng.normal(0,1,(iK,1))
    vBetaY = rng.normal(0,1,(iK,1))
    # generate mediator
    vM = mX@vBetaM + vG*dGM + vGS*dGSM + vES*dESM
    # generate outcome
    vY_no = mX@vBetaY + vM*0 + vG*dGY + vE*dEY
    vY_partial = mX@vBetaY + vM*dB + vG*dGY + vE*dEY
    vY_full = mX@vBetaY + vM*dB + vG*0 + vE*dEY
    # set phenotype matrices
    mY_no = np.hstack((vM,vY_no))
    mY_partial = np.hstack((vM,vY_partial))
    mY_full = np.hstack((vM,vY_full))
    # set matrix of covariates including mediator
    mXM = np.hstack((mX,vM))
    # generate names of phenotypes and covariates
    lPheno = ['Mediator','Outcome']
    lCovar = ['Covariate ' + str(i) for i in range(1,iK)]
    lCovarM = ['Covariate ' + str(i) for i in range(1,iK)]
    lCovarM.append('Mediator')
    # construct DataFrames phenotypes
    dfY_no = pd.DataFrame(mY_no, index = miID, columns = lPheno)
    dfY_partial = pd.DataFrame(mY_partial, index = miID, columns = lPheno)
    dfY_full = pd.DataFrame(mY_full, index = miID, columns = lPheno)
    # construct DataFrame mediator only
    dfM = pd.DataFrame(vM, index = miID, columns = [lPheno[0]])
    # construct DataFrames outcome only
    dfO_no = pd.DataFrame(vY_no, index = miID, columns = [lPheno[1]])
    dfO_partial = pd.DataFrame(vY_partial, index = miID, columns = [lPheno[1]])
    dfO_full = pd.DataFrame(vY_full, index = miID, columns = [lPheno[1]])
    # construct DataFrame covariates, ignoring the intercept
    # as MGREML takes care of that by itself
    dfX = pd.DataFrame(mX[:,1:], index = miID, columns = lCovar)
    # construct DataFrame covariates + mediator, ignoring the intercept
    # as MGREML takes care of that by itself
    dfXM = pd.DataFrame(mXM[:,1:], index = miID, columns = lCovarM)
    # create filenames
    sFile_no = 'no_mediation.run.' + str(iRun) + '.pheno.txt'
    sFile_partial = 'partial_mediation.run.' + str(iRun) + '.pheno.txt'
    sFile_full = 'full_mediation.run.' + str(iRun) + '.pheno.txt'
    sFile_covar = 'run.' + str(iRun) + '.covar.txt'
    # write phenotypes and covariates to tab separated files
    dfY_no.to_csv(sFile_no, sep='\t')
    dfY_partial.to_csv(sFile_partial, sep='\t')
    dfY_full.to_csv(sFile_full, sep='\t')
    dfX.to_csv(sFile_covar, sep='\t')
    # create filenames for mediator only and covariates with mediator
    sFile = 'run.' + str(iRun) + '.mediator.txt'
    sFile_covar = 'run.' + str(iRun) + '.covar_mediator.txt'
    # write mediator only and covariates with mediator as tab separated files
    dfM.to_csv(sFile, sep='\t')
    dfXM.to_csv(sFile_covar, sep='\t')
    # create filenames for outcome only
    sFile_no = 'no_mediation.run.' + str(iRun) + '.outcome.txt'
    sFile_partial = 'partial_mediation.run.' + str(iRun) + '.outcome.txt'
    sFile_full = 'full_mediation.run.' + str(iRun) + '.outcome.txt'
    # write outcome only tab separated files
    dfO_no.to_csv(sFile_no, sep='\t')
    dfO_partial.to_csv(sFile_partial, sep='\t')
    dfO_full.to_csv(sFile_full, sep='\t')

def SimulateData(iRun,iN,iM):
    # set main seed for np.random
    iMainSeed = 14159703
    # set maximum allowed no. of runs
    iMAXRUNS = 10000
    # terminate if this run beyond max no. allowed
    if iRun >= iMAXRUNS:
        return
    # use main seed to set main random number generator
    rngMain = np.random.default_rng(iMainSeed)
    # set seeds for all runs
    vSeed = rngMain.integers(0,iMainSeed,iMAXRUNS)
    # get seed for this run
    iSeed = vSeed[iRun]
    # set random number generator this run
    rng = np.random.default_rng(iSeed)
    # generate FIDs and IIDs
    lFID = ['FID ' + str(i) for i in range(1,iN+1)]
    lIID = ['IID ' + str(iN+i) for i in range(1,iN+1)]
    lID = [lFID,lIID]
    # make multiindex
    lLabelsFID_IID = ['FID', 'IID']
    miID = pd.MultiIndex.from_arrays(lID,names=lLabelsFID_IID)
    # simulate genetic data and get GRM
    (mG, mA) = SimulateGenotypes(iN,iM,rng)
    # simulate phenotypes and covariates
    SimulatePhenotypes(mG,miID,iRun,rng)
    # write GRM data
    sPrefix = 'run.' + str(iRun) + '.'
    WriteGRM(mA,iM,lFID,lIID,sPrefix)

def main(argv):
    iRun = int(argv[1])
    iN = int(argv[2])
    iM = int(argv[3])
    print('SIMULATING DATA FOR RUN ' + str(iRun))
    print('WITH N = ' + str(iN) + ' AND M = ' + str(iM))
    SimulateData(iRun,iN,iM)

if __name__ == "__main__":
    main(sys.argv)

