import numpy as np
import pandas as pd
from numpy.matlib import repmat

def SimulateData():
    # set seed for np.random
    np.random.seed(392993604)
    # set sample size, no. of SNPs, no. of traits and covariates
    iN = 5000
    iM = 10000
    iT = 10
    iK = 10
    # set ploidy of human genomen
    iP = 2    
    # set the number of environment and genetic factors
    iFG = 2
    iFE = iT
    # settings for distr. for drawing SNP allele frequencies
    dBetaParam1 = 0.35
    dBetaParam2 = 0.35
    dMinMAF = 1E-2    
    # initialise allele frequencies
    vAF = np.random.beta(dBetaParam1,dBetaParam2,size=iM)
    # keep updating until all between dMinMAF and 1-dMinMAF
    while (min(vAF) < dMinMAF) or (max(vAF) > (1-dMinMAF)):
        (vIndAFlow,) = np.where(vAF < dMinMAF)
        (vIndAFhigh,) = np.where(vAF > (1-dMinMAF))
        vAF[vIndAFlow] = np.random.beta(dBetaParam1,dBetaParam2,size=vIndAFlow.shape)
        vAF[vIndAFhigh] = np.random.beta(dBetaParam1,dBetaParam2,size=vIndAFhigh.shape)
    # initialise matrix of genotypes
    mG = np.zeros((iN,iM))
    # for each SNP
    for m in range(0,iM):
        # draw SNP data as binomially distributed according to AF of that SNP
        mG[:,m] = np.random.binomial(iP,vAF[m],size=iN)
    # standardise the SNP data by semi-empirical allele frequency
    vEAF = mG.mean(axis=0)/2
    mX = (mG - repmat(2*vEAF,iN,1))*repmat((2*vEAF*(1-vEAF))**(-0.5),iN,1)
    # construct GRM
    mA = (mX@mX.T)/iM
    # construct factors
    mFG = (mX@np.random.normal(0,1,(iM,iFG)))*(iM**(-0.5))
    mFE = np.random.normal(0,1,(iN,iFE))
    # construct a bunch of covariates and their random effects
    mX = np.hstack((np.ones((iN,1)),np.random.normal(0,1,(iN,iK-1))))
    mBeta = np.random.normal(0,1,(iK,iT))
    # construct factor coefficients and phenotypes
    mCG = np.random.normal(0,1,(iT,iFG))
    mCE = np.random.normal(0,1,(iT,iFE))
    mY = mFG@mCG.T + mFE@mCE.T + mX@mBeta
    # generate FIDs and IIDs
    lFID = ['FID ' + str(i) for i in range(1,iN+1)]
    lIID = ['IID ' + str(iN+i) for i in range(1,iN+1)]
    lID = [lFID,lIID]
    # make multiindex
    lLabelsFID_IID = ['FID', 'IID']
    miID = pd.MultiIndex.from_arrays(lID,names=lLabelsFID_IID)
    # generate names of phenotypes
    lPheno = ['Some pheno ' + str(100+t) for t in range(1,iT+1)]
    # generate names of covariates
    lCov = ['my covar ' + str(300+k) for k in range(0,iK)]
    lCov[0] = 'intercept'
    # construct DataFrames
    dfX = pd.DataFrame(mX, index = miID, columns = lCov)
    dfA = pd.DataFrame(mA, index = miID, columns = lID)
    dfY = pd.DataFrame(mY, index = miID, columns = lPheno)
    # write phenotype and covariate data to tab separated files
    sFileX = 'covar.txt'
    sFileY = 'pheno.txt'
    dfX.to_csv(sFileX, sep='\t')
    dfY.to_csv(sFileY, sep='\t')
    # write GRM to txt-based format
    sFileAID = 'data.grm.id'
    with open(sFileAID, 'w') as oFile:
        for i in range(0,iN):
            oFile.write(lFID[i] + '\t' + lIID[i] + '\n')
    sFileA = 'data.grm'
    with open(sFileA, 'w') as oFile:
        for i in range(0,iN):
            if i%500 == 0:
                print(i)
            for j in range(0,i+1):
                oFile.write(str(i+1) + '\t' + str(j+1) + '\t' + str(iM) + '\t' + str(dfA.iloc[i,j]) + '\n')
    # gzip data.grm and convert to binary GRM e.g. using PLINK

# run SimulateData
SimulateData()
