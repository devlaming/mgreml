import numpy as np
import pandas as pd
import pickle
import math
from mgreml.data import aligner

def TestMgremlAligner():
    # set number of traits, factors, observations, SNPs, and covariates
    iT  = 10
    iFG = iT
    iFE = iT
    iN  = 1000
    iS  = 5*iN
    iK  = 2
    # simulate SNPs effects and SNPs (everything normally distributed)
    mBeta = np.random.normal(size=(iS,iFG))
    mG    = np.random.normal(size=(iN,iS))
    # get genetic factors
    mFG   = (mG@mBeta)/math.sqrt(iS)
    # simulate environment factors
    mFE   = np.random.normal(size=(iN,iFE))
    # simulate factor coefficients
    mCG   = np.random.normal(size=(iT,iFG))
    mCE   = np.random.normal(size=(iT,iFE))
    # get genetic and environment covariance matrices
    mVG = mCG@mCG.T
    mVE = mCE@mCE.T
    # simulate covariates and effects
    mX     = np.hstack((np.ones((iN,1)),np.random.normal(size=(iN,iK-1))))
    mAlpha = np.random.normal(size=(iK,iT))
    # compute GRM
    mA = (mG@mG.T)/iS
    # get phenotypes
    mY    = mFG@mCG.T + mFE@mCE.T + mX@mAlpha
    # draw IDs
    lFID = ['FAM_' + str(x) for x in np.random.permutation(iN)]
    lIID = ['IND_' + str(x) for x in np.random.permutation(iN)]
    # set 2d array of IDs
    lIDs = [lFID, lIID]
    # draw variable names
    lPhenoNames = ['trait ' + str(x) for x in np.random.permutation(iT)]
    lCovarNames = ['covariate ' + str(x) for x in np.random.permutation(iK)]
    # put in some missing values in phenotypes and covariates
    mY[np.where(np.random.normal(size=(iN,iT))>3.3)] = None
    mX[np.where(np.random.normal(size=(iN,iK))>3.3)] = None
    # randomly select which covariates affect which traits
    mXbinary = (np.random.normal(size=(iT,iK))>1.5).astype(int)
    dfBinXY = pd.DataFrame(mXbinary, index=lPhenoNames, columns=lCovarNames)
    # set dataframes for GRM, phenotype data and covariates
    dfY  = pd.DataFrame(mY, columns=lPhenoNames, index=lIDs)
    dfX  = pd.DataFrame(mX, columns=lCovarNames, index=lIDs)
    dfA  = pd.DataFrame(mA, columns=lIDs, index=lIDs)    
    # shuffle dataframes for phenotypes and covariates, drop 100 observations each
    vRandIndY = np.random.permutation(len(dfY))
    vRandIndX = np.random.permutation(len(dfX))
    vRandIndY = vRandIndY[0:(iN-100)]
    vRandIndX = vRandIndX[0:(iN-100)]
    dfY = dfY.iloc[vRandIndY]
    dfX = dfX.iloc[vRandIndX]
    # create mgreml data object
    MyMgremlData = aligner.MgremlData(dfY, dfA, dfX)
    # set filenames
    sFilepath = './temp/'
    sFileName = sFilepath + 'MyMgremlData.pkl'
    # write GremlData instance to pickled files
    with open(sFileName, 'wb') as handle:
        pickle.dump(MyMgremlData, handle)
    sFileNameTrueParams = sFilepath + 'MyMgremlData.TrueParams.pkl'
    dictOut = {'mAlpha': mAlpha, 'mVG': mVG, 'mVE': mVE, 'lPhenoNames': lPhenoNames, 'lCovarNames': lCovarNames}
    with open(sFileNameTrueParams, 'wb') as handle:
        pickle.dump(dictOut, handle)
    # done testing
    print('Done testing the MGREML data aligner.')
