''' 
simulatedata.py: simulates data for simulation study of MGREML
last edit: June 4, 2021
'''
import numpy as np
import pandas as pd
import sys, getopt
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

def SimulatePhenotypes(mG,miID,mWeightG,mWeightE,dHSq,sPrefix,rng):
    # get no. of SNPs and sample size
    iN = mG.shape[0]
    iM = mG.shape[1]
    # get no. of factors and traits
    iT = mWeightG.shape[0]
    iFG = mWeightG.shape[1]
    iFE = mWeightE.shape[1]
    # construct factors
    mFG = (mG@rng.normal(0,1,(iM,iFG)))*(iM**(-0.5))
    mFE = rng.normal(0,1,(iN,iFE))
    # standardise factors
    mFG = (mFG - repmat(mFG.mean(axis=0),iN,1))*repmat(mFG.var(axis=0)**(-0.5),iN,1)
    mFE = (mFE - repmat(mFE.mean(axis=0),iN,1))*repmat(mFE.var(axis=0)**(-0.5),iN,1)
    # standardise factor coefficients, so VarG = VarE = 1 for each trait
    mCG = (dHSq**0.5)*(np.diag(((mWeightG**2).sum(axis=1))**(-0.5))@mWeightG)
    mCE = ((1-dHSq)**0.5)*(np.diag(((mWeightE**2).sum(axis=1))**(-0.5))@mWeightE)
    # generate liabilities and phenotypes
    mLiabG = mFG@mCG.T
    mLiabE = mFE@mCE.T
    mY = mLiabG + mLiabE
    # generate names of phenotypes
    lPheno = ['Trait ' + str(t) for t in range(1,iT+1)]
    # construct DataFrame phenotype
    dfY = pd.DataFrame(mY, index = miID, columns = lPheno)
    # write phenotype and covariate data to tab separated files
    sFileY = sPrefix + 'pheno.txt'
    dfY.to_csv(sFileY, sep='\t')
    # compute true heritabilities and store
    mVG = mCG@mCG.T
    mVE = mCE@mCE.T
    vHSq = np.diag(mVG) / np.diag(mVG+mVE)
    dfHSq = pd.DataFrame(vHSq,index=lPheno,columns=['heritability'])
    sHSq = sPrefix + 'true.HSq.txt'
    dfHSq.to_csv(sHSq, sep='\t')
    # compute true correlations and store
    mRhoG = np.outer(np.diag(mVG)**(-0.5),np.diag(mVG)**(-0.5))*(mVG)
    mRhoE = np.outer(np.diag(mVE)**(-0.5),np.diag(mVE)**(-0.5))*(mVE)
    dfRhoG = pd.DataFrame(mRhoG,index=lPheno,columns=lPheno)
    dfRhoE = pd.DataFrame(mRhoE,index=lPheno,columns=lPheno)
    sRhoG = sPrefix + 'true.RhoG.txt'
    sRhoE = sPrefix + 'true.RhoE.txt'
    dfRhoG.to_csv(sRhoG, sep='\t')
    dfRhoE.to_csv(sRhoE, sep='\t')

def SimulateData(iRun=0,iN=20000,iM=50000):
    # set no. of traits, and HSq
    iT = 10
    dHSq = 0.5
    # set main seed for np.random
    iMainSeed = 81239616
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
    ## SIMULATION DESIGN 1: RhoG as low as possible ##
    dRhoG = -1/(iT-1)
    # set rhoG matrix and find EVD
    mRhoG = ((1-dRhoG)*np.eye(iT)) + dRhoG
    (vDG,mPG) = np.linalg.eigh(mRhoG)
    # set negative EVs as result of lacking precision to zero
    vDG[vDG<0] = 0
    # construct genetic and environment factor weights
    mWeightG = mPG@np.diag(np.sqrt(vDG))
    mWeightE = np.eye(iT)
    # set prefix
    sPrefix = 'rhoGneg.run' + str(iRun) + '.'
    # simulate phenotypes
    SimulatePhenotypes(mG,miID,mWeightG,mWeightE,dHSq,sPrefix,rng)
    ## SIMULATION DESIGN 2: RhoG = 0 ##
    dRhoG = 0
    # set rhoG matrix and find EVD
    mRhoG = ((1-dRhoG)*np.eye(iT)) + dRhoG
    (vDG,mPG) = np.linalg.eigh(mRhoG)
    # set negative EVs as result of lacking precision to zero
    vDG[vDG<0] = 0
    # construct genetic and environment factor weights
    mWeightG = mPG@np.diag(np.sqrt(vDG))
    mWeightE = np.eye(iT)
    # set prefix
    sPrefix = 'rhoGzero.run' + str(iRun) + '.'
    # simulate phenotypes
    SimulatePhenotypes(mG,miID,mWeightG,mWeightE,dHSq,sPrefix,rng)
    ## SIMULATION DESIGN 3: RhoG = 0.5 ##
    dRhoG = 0.5
    # set rhoG matrix and find EVD
    mRhoG = ((1-dRhoG)*np.eye(iT)) + dRhoG
    (vDG,mPG) = np.linalg.eigh(mRhoG)
    # set negative EVs as result of lacking precision to zero
    vDG[vDG<0] = 0
    # construct genetic and environment factor weights
    mWeightG = mPG@np.diag(np.sqrt(vDG))
    mWeightE = np.eye(iT)
    # set prefix
    sPrefix = 'rhoGhalf.run' + str(iRun) + '.'
    # simulate phenotypes
    SimulatePhenotypes(mG,miID,mWeightG,mWeightE,dHSq,sPrefix,rng)
    ## SIMULATION DESIGN 4: RhoG = 1 ##
    # construct genetic and environment factor weights
    mWeightG = np.ones((iT,1))
    mWeightE = np.eye(iT)
    # set prefix
    sPrefix = 'rhoGone.run' + str(iRun) + '.'
    # simulate phenotypes
    SimulatePhenotypes(mG,miID,mWeightG,mWeightE,dHSq,sPrefix,rng)
    ## SIMULATION DESIGN 5: RhoG zero between blocks, and free within blocks ##
    # construct genetic and environment factor weights
    mWeightG = rng.normal(0,1,(iT,iT))
    mWeightE = np.eye(iT)
    # make sure 1st half of pheno's is affected only by 1st half of genetic factors
    # and 2nd half of pheno's is affected only by 2nd half of genetic factors
    mWeightG[0:int(iT/2),int(iT/2):] = 0
    mWeightG[int(iT/2):,0:int(iT/2)] = 0
    # set prefix
    sPrefix = 'blocks.run' + str(iRun) + '.'
    # simulate phenotypes
    SimulatePhenotypes(mG,miID,mWeightG,mWeightE,dHSq,sPrefix,rng)
    # write GRM data
    sPrefix = 'run' + str(iRun) + '.'
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

