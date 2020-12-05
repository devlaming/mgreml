import numpy as np
import pandas as pd
import math
from numpy.matlib import repmat
from tqdm import tqdm
pd.options.mode.chained_assignment = None

class StructuralModel:
    
    # regularisation constant initialisation variance matrix
    dLambdaInit = 1E-6
    # set minimum required squared sum of coefficients for each factor
    # when diagnosing issues
    dSSTOL  = 1E-9
    # set maximum Rsq of coefficients for each factor w.r.t. all other factors
    # when diagnosing issues
    dRSqTOL = 1-(1E-9)
    
    def __init__(self, mdData, dfBinFY = None, sType=''):
        # check if model has been specified
        bModelSpecified = isinstance(dfBinFY, pd.DataFrame)
        # if model specified, set structural model accordingly
        if bModelSpecified:
            # count number of traits and factors
            self.iT = dfBinFY.shape[0]
            self.iF = dfBinFY.shape[1]
            # check specification
            self.CheckModelSpecification(mdData, dfBinFY)
            # initialise structural model
            self.InitialiseStructuralModel(mdData, dfBinFY)
        else: # else
            # count number of traits and factors
            self.iT = mdData.mY.shape[1]
            self.iF = self.iT
            print('Assuming saturated ' + sType + ' model')
            # initialise saturated model
            self.InitialiseSaturatedModel(mdData)
        
    def CheckModelSpecification(self, mdData, dfBinFY):
        # raise error if there are more factors than traits
        if self.iF > self.iT:
            raise ValueError('You have specified more factors than traits. This model is not identified.') 
        # raise error if number of traits does not match no. of trais in data
        if self.iT != mdData.mY.shape[1]:
            raise ValueError('The number of traits in your structural model does not match the number of traits in your data.')
        # if duplicates in index of dfBinFY
        if dfBinFY.index.duplicated().sum() > 0:
            raise ValueError('You have phenotypes with duplicate labels in your structural model.')
        # if duplicates in columns of dfBinFY
        if dfBinFY.columns.duplicated().sum() > 0:
            raise ValueError('You have factors with duplicate labels in your structural model.')
        # get indices for phentoypes
        indY_Y  = pd.Index(mdData.lPhenos)
        indFY_Y = dfBinFY.index
        # all phenotypes in data should refer to phenotypes in dfBinFY; abort if problematic
        if not(indY_Y.isin(indFY_Y).all()):
            raise ValueError('There is at least one phenotype which has not been specified in your structural model.') 
        # all phenotypes in dfBinFY should refer to phenotypes in data; abort if problematic
        if not(indFY_Y.isin(indY_Y).all()):
            raise ValueError('There is at least one phenotype in your structural model for which you have not provided data.') 
        # if dfBinFY is not binary: abort
        if ((dfBinFY==0).sum().sum() + (dfBinFY==1).sum().sum()) != (dfBinFY.shape[0]*dfBinFY.shape[1]):
            raise ValueError('The specification of your structural model does not only comprise zeros and ones, signifying whether a given factor is permitted to affect a given phenotype.')
            
    def InitialiseStructuralModel(self, mdData, dfBinFY):
        # get indices for phentoypes from data
        indY  = pd.Index(mdData.lPhenos)
        # sort phenotypes in factor model specification by order in data
        dfBinFY = dfBinFY.loc[indY]
        # double-check if everything is now lined up
        if not(all(dfBinFY.index == indY)):
            raise ValueError('The phenotype labels in your data and structural model cannot be lined up properly.')
        # read out factor labels
        self.lFactors = dfBinFY.columns.tolist()
        # read out structural model as numpy as array
        mB = np.array(dfBinFY)
        # if rank falls below self.iF: abort
        if np.linalg.matrix_rank(mB) < self.iF:
            raise ValueError('Your structural model is rank deficient. Please change your specification.')
        # get trait and factor indices based on mB
        (self.vIndT,self.vIndF) = np.where(mB==1)
        # get regularised phenotypic covariance matrix
        mV = self.InitialiseV(mdData)
        # initialise parameters
        self.vParam = np.zeros(self.vIndT.shape)
        # for each trait
        for t in range(0,self.iT):
            # find indices of all factors affecting current trait
            vParamIndThisTrait = np.array(np.where(self.vIndT==t)).ravel().astype(int)
            # count number of factors affecting this trait
            iThisF = len(vParamIndThisTrait)
            if iThisF > 0: # if this trait affected by at least 1 factor
                # set coefficient weights such that each phenotype t
                # has implied variance equal to mV[t,t]
                self.vParam[vParamIndThisTrait] = math.sqrt(mV[t,t]/iThisF)
        
    def InitialiseSaturatedModel(self, mdData):
        # get regularised phenotypic covariance matrix
        mV = self.InitialiseV(mdData)
        # get cholesky decomposition (lower triangular)
        mC = np.linalg.cholesky(mV)
        # set trait and factor indices of saturated model
        self.InitialiseIndicesSaturated()
        # set parameters
        self.vParam = mC[self.vIndT,self.vIndF]
        # set labels for the factors
        self.lFactors = ['factor ' + str(x) for x in np.arange(self.iT)]

    def InitialiseV(self, mdData):
        # get regularised phenotypic covariance matrix
        mV = (1-StructuralModel.dLambdaInit)*np.cov(mdData.mY.T) + StructuralModel.dLambdaInit*np.eye(self.iT)
        return mV
    
    def InitialiseIndicesSaturated(self):
        # compute no. of parameters for satured model
        iParam = int((self.iT+1)*self.iT/2)
        # set conformable vectors for trait and factor indices
        self.vIndT = np.zeros(iParam)
        self.vIndF = np.zeros(iParam)
        # set counter for parameters
        iCount = 0
        # for each factor
        for f in range(0,self.iT):
            # for each trait
            for t in range(f,self.iT):
                # update vectors with trait and factor indices
                self.vIndT[iCount] = t
                self.vIndF[iCount] = f
                # update parameter count
                iCount = iCount + 1
        # make sure trait and factor indices are integer vectors
        self.vIndT = self.vIndT.astype(int)
        self.vIndF = self.vIndF.astype(int)
        
    def DiagnoseProblem(self, sType=''):
        # get coefficient matrix
        mC = self.GetC()
        # for each factor
        for f in range(0,self.iF):
            # compute squared sum of coefficients
            dSS = np.power(mC[:,f],2).sum()
            if (dSS < StructuralModel.dSSTOL): # if squared sum too low
                # print warning
                sWarning = 'COEFFICIENTS FOR ' + sType + ' FACTOR ' + str(f) + ' ARE ALL CLOSE TO ZERO. I.E. SQUARED SUM OF COEFFICIENTS LESS THAN ' + str(StructuralModel.dSSTOL)
                raise ValueError(sWarning)
            # if there are multiple factors
            if (self.iF > 1):
                # take submatrix of all coefficients, except for current factor
                mCS = np.hstack((mC[:,0:f],mC[:,f+1:]))
                # get coefficients for current factor
                vC = mC[:,f]
                # compute OLS residual for current factor w.r.t. other factors
                vR = vC - (mCS@(np.linalg.inv(mCS.T@mCS)@(mCS.T@vC)))
                # compute RSq of regression of current factor coefficients on 
                # all other factor coefficients
                dRSq = 1 - ((np.power(vR,2).sum())/(np.power(vC,2).sum()))
                if (dRSq > StructuralModel.dRSqTOL): # if RSq too high
                    # print warning
                    sWarning = 'COEFFICIENTS FOR ' + sType + ' FACTOR ' + str(f) + ' ARE MULTICOLLINEAR W.R.T. COEFFICIENTS FOR OTHER FACTORS. I.E. R-SQUARED EXCEEDS ' + str(100*StructuralModel.dRSqTOL) + '%'
                    raise ValueError(sWarning)
        # for each phenotype
        for t in range(0,self.iT):
            # compute squared sum of coefficients
            dSS = np.power(mC[t,:],2).sum()
            if (dSS < StructuralModel.dSSTOL): # if squared sum too low
                # print warning
                sWarning = sType + ' COEFFICIENTS FOR PHENOTYPE ' + str(t) + ' ARE ALL CLOSE TO ZERO. I.E. SQUARED SUM OF COEFFICIENTS LESS THAN ' + str(StructuralModel.dSSTOL)
                raise ValueError(sWarning)
            # if there are multiple phenotypes
            if (self.iT > 1):
                # take submatrix of all phenotypes, except for current phenotype
                mCS = np.hstack((mC[0:t,:].T,mC[t+1:,:].T))
                # get coefficients for current factor
                vC = mC[t,:].T
                # compute OLS residual for current factor w.r.t. other factors
                vR = vC - (mCS@(np.linalg.inv(mCS.T@mCS)@(mCS.T@vC)))
                # compute RSq of regression of current factor coefficients on 
                # all other factor coefficients
                dRSq = 1 - ((np.power(vR,2).sum())/(np.power(vC,2).sum()))
                if (dRSq > StructuralModel.dRSqTOL): # if RSq too high
                    # print warning
                    sWarning = sType + ' COEFFICIENTS FOR PHENOTYPE ' + str(t) + ' ARE MULTICOLLINEAR W.R.T. COEFFICIENTS FOR OTHER PHENOTYPES. I.E. R-SQUARED EXCEEDS ' + str(100*StructuralModel.dRSqTOL) + '%'
                    raise ValueError(sWarning)        
    
    def GetVandC(self, vNew = None):
        # get matrix of coefficients
        mC = self.GetC(vNew)
        # return V=CC' as variance matrix
        return np.matmul(mC,mC.T), mC
    
    def GetV(self, vNew = None):
        # get matrix of coefficients
        mC = self.GetC(vNew)
        # return V=CC' as variance matrix
        return np.matmul(mC,mC.T)
    
    def GetC(self, vNew = None):
        # set matrix of coefficients TxF,
        # where T=no. of traits and F=no. of factors
        mC = np.zeros((self.iT,self.iF))
        if isinstance(vNew, np.ndarray):
            mC[self.vIndT,self.vIndF] = vNew
        else:
            # set values in accordance with parameters and their indices
            mC[self.vIndT,self.vIndF] = self.vParam
        return mC
    
    def GetFreeCoeffs(self):
        mB = np.zeros((self.iT,self.iF))
        mB[self.vIndT,self.vIndF] = 1
        return mB
    
    def UpdateParams(self, vNew):
        if isinstance(vNew, np.ndarray):
            self.vParam = vNew
        else:
            raise TypeError('The new parameter estimates of the model are invaled.')

class GeneticModel(StructuralModel):
    
    # set weight of genetic factors to square root of 20%
    # as to emulate initial heritabilities of 20%
    dWeight = math.sqrt(0.2)
    # string describing type of model
    sType = 'GENETIC'
    
    def __init__(self, mdData, dfBinFY = None):
        super().__init__(mdData, dfBinFY, GeneticModel.sType)
        self.vParam = GeneticModel.dWeight*self.vParam
        self.lFactors = ['genetic ' + str(x) for x in self.lFactors]
        
    def DiagnoseProblem(self):
        super().DiagnoseProblem(GeneticModel.sType)

class EnvironmentModel(StructuralModel):
    
    # set weight of environment factors to square root of 80%
    # as to emulate initial heritabilities of 80%
    dWeight = math.sqrt(0.8)
    # string describing type of model
    sType = 'ENVIRONMENT'
    
    def __init__(self, mdData, dfBinFY = None):
        super().__init__(mdData, dfBinFY, EnvironmentModel.sType)
        self.vParam = EnvironmentModel.dWeight*self.vParam
        self.lFactors = ['environment ' + str(x) for x in self.lFactors]
        # if less environment factors than traits: crash, as this is incompatible with MGREML
        if self.iF < self.iT:
            raise ValueError('You have specified less environmental factors than traits. This is not permitted in MGREML.')
            
    def DiagnoseProblem(self):
        super().DiagnoseProblem(EnvironmentModel.sType)

class CombinedModel:
    
    def __init__(self, mdData, dfGenBinFY = None, dfEnvBinFY = None):
        self.genmod = GeneticModel(mdData, dfGenBinFY)
        self.envmod = EnvironmentModel(mdData, dfEnvBinFY)
        # count the number of parameters and parameter combinations
        self.iParamsG = self.genmod.vParam.shape[0]
        self.iParamsE = self.envmod.vParam.shape[0]
        self.iParams  = self.iParamsG + self.iParamsE
        self.iParamCombos = int(((self.iParams+1)*self.iParams)/2)
    
    def SplitNewParameters(self, vNew):
        if isinstance(vNew, np.ndarray):
            vNewG = vNew[0:self.iParamsG]
            vNewE = vNew[self.iParamsG:]
            return vNewG, vNewE
        else:
            raise TypeError('The new parameter estimates of the model are invaled.')
    
    def GetVandC(self, vNew = None):
        if isinstance(vNew, np.ndarray):
            (vNewG, vNewE) = self.SplitNewParameters(vNew)
            (mVG, mCG) = self.genmod.GetVandC(vNewG)
            (mVE, mCE) = self.envmod.GetVandC(vNewE)
        else:
            (mVG, mCG) = self.genmod.GetVandC()
            (mVE, mCE) = self.envmod.GetVandC()
        return mVG, mCG, mVE, mCE
        
    def GetParams(self):
        vParam = np.zeros(self.iParams)
        vParam[0:self.iParamsG] = self.genmod.vParam
        vParam[self.iParamsG:]  = self.envmod.vParam
        return vParam
    
    def UpdateParams(self, vNew):
        if isinstance(vNew, np.ndarray):
            (vNewG, vNewE) = self.SplitNewParameters(vNew)
            self.genmod.UpdateParams(vNewG)
            self.envmod.UpdateParams(vNewE)
        else:
            raise TypeError('The new parameter estimates of the model are invaled.')
            
    def GetFreeCoeffs(self):
        mBG = self.genmod.GetFreeCoeffs()
        mBE = self.envmod.GetFreeCoeffs()
        return mBG, mBE

class MgremlModel:
    
    # set lowest eigenvalue permitted in VE matrix without aborting
    dMinEigVal = 1E-12
    
    def __init__(self, mdData, dfGenBinFY = None, dfEnvBinFY = None):
        self.model = CombinedModel(mdData, dfGenBinFY, dfEnvBinFY)
        self.data  = mdData
        self.vBetaGLS = None
        self.mVarGLS = None
    
    def ComputeLogLik(self, bGrad = False, bInfo = False, vNew = None, bSilent = False):
        """
        Author     : Ronald de Vlaming
        Date       : December 3, 2020
        
        Summary    : compute log-likelihood, gradient and information matrix
        
        Input
        bGrad      : boolean: is gradient desired?
        bInfo      : boolean: is average information matrix required?
        
        """
        # for brevity of code, get sample siz and so on
        iN           = self.data.iN
        iT           = self.data.iT
        iFG          = self.model.genmod.iF
        iFE          = self.model.envmod.iF
        bCovs        = self.data.bCovs
        bSameCovs    = self.data.bSameCovs
        iParamsG     = self.model.iParamsG
        iParamsE     = self.model.iParamsE
        iParams      = self.model.iParams
        iParamCombos = self.model.iParamCombos
        if not(bGrad): # if gradient is not desired
            # do not return identity as information matrix
            bInfo = False
        if bCovs: # if we have covariates
            # get number of covariates
            iK = self.data.iK
            if not(bSameCovs): # but different covariates per trait
                # get indices of these covariates in terms of Z matrix etc.
                vIndCovs = np.array(np.where(np.array(self.data.mBinXY).ravel()==1)).ravel()
                # count total number of covariates across traits
                iKtotal = self.data.mBinXY.sum()
        # convert parameters to variance components
        (mVG,mCG,mVE,mCE) = self.model.GetVandC(vNew)
        # stabilize mVE and mVG
        mVG = (mVG + mVG.T)/2
        mVE = (mVE + mVE.T)/2
        # get EVD of matrix of environment covariance matrix
        (vPhi,mQ) = np.linalg.eigh(mVE)
        # abort if eigenvalues of environment variance matrix too low
        if min(vPhi) <= MgremlModel.dMinEigVal:
            print('The environment covariance matrix is rank deficient. Possible reasons: (1) multicollinearity between phenotypes, (2) a poorly specified model, and/or (3) poor starting values.')
            self.model.envmod.DiagnoseProblem()
            raise ValueError('Rank deficient environment covariance matrix')
        # set the outer product of the vector of one over EVs
        vPhiNegSqrt  = np.power(vPhi,-0.5)
        mOuterPhiInv = np.outer(vPhiNegSqrt,vPhiNegSqrt)
        # rescale genetic covariance components matrix based on EVD
        mVGadj = np.multiply(mOuterPhiInv,np.matmul(np.matmul(mQ.T,mVG),mQ))
        # stabilise mVGadj
        mVGadj = (mVGadj + mVGadj.T)/2
        # get EVs of rescaled genetic components and rescale Q
        (vLambda,mL) = np.linalg.eigh(mVGadj)
        mF           = np.matmul(np.multiply(mQ,repmat(vPhiNegSqrt,iT,1)),mL)
        mFT          = mF.T
        # get pairwise product of EVs in vLambda and vD and transform
        mDS = 1/(np.outer(self.data.vD,vLambda)+1)
        if bGrad: # if gradient is desired
            # multiply one over paired EVs by EVs of GRM and store
            mDDS = np.multiply(mDS,repmat(self.data.vD,iT,1).T)
        # compute first terms for log-lik
        dLogDetV = iN*np.log(vPhi).sum() - np.log(mDS).sum()
        dCons    = iN*iT*np.log(2*np.pi)
        # initialise log|X'inv(V)X| at minus iK*log|V_E| (=0 if no covariates)
        dLogDetZ = -iK*np.log(vPhi).sum()
        if bCovs: # if we have covariates
            if bSameCovs: # if same covs across traits
                # set precursor for inv(X'inv(V)X)
                mJ = np.zeros((iK*iT,iK*iT))
                # if gradient desired: set precursors for mUG and mUE
                if bGrad:
                    vDiagUG = np.zeros(iT)
                    vDiagUE = np.zeros(iT)
            else: # if different covs
                # set mZ matrix and grand BG and BE matrices
                mZ       = np.zeros((iK*iT,iK*iT))
                mGrandBG = np.zeros((iK*iT,iK*iT))
                mGrandBE = np.zeros((iK*iT,iK*iT))
            # for each trait
            for i in range(0,iT):
                # compute the corresponding matrix X'WX, with W = diagonal
                mB = np.matmul(np.multiply(self.data.mXT,repmat(mDS[:,i],iK,1)),self.data.mX)
                # stabilise mB
                mB = (mB + mB.T)/2
                if bSameCovs: # if same covs across traits
                    # compute the EVD of mB
                    (vEigValsB,mEigVecsB) = np.linalg.eigh(mB)
                    # update log|X'inv(V)X|
                    dLogDetZ = dLogDetZ + np.log(vEigValsB).sum()
                    # compute and store inverse of mB
                    mThisInvB = np.matmul(np.multiply(mEigVecsB,repmat(1/vEigValsB,iK,1)),mEigVecsB.T)
                    mJ[i*iK:(i+1)*iK,i*iK:(i+1)*iK] = mThisInvB
                    # if gradient desired
                    if bGrad:
                        # compute squared values of relevant column in mDS
                        vDSSquared = np.power(mDS[:,i],2)
                        # compute precursors for mUG and mUE
                        mBG = np.matmul(np.multiply(self.data.mXT,repmat(self.data.vD*vDSSquared,iK,1)),self.data.mX)
                        mBE = np.matmul(np.multiply(self.data.mXT,repmat(vDSSquared,iK,1)),self.data.mX)
                        vDiagUG[i] = np.diag(np.matmul(mThisInvB,mBG)).sum()
                        vDiagUE[i] = np.diag(np.matmul(mThisInvB,mBE)).sum()
                else:
                    # update mZ matrix
                    mZ[i*iK:(i+1)*iK,i*iK:(i+1)*iK] = mB
                    # if gradient desired
                    if bGrad:
                        # for each other trait
                        for k in range(i,iT):
                            # compute product of values in relevant column in mDS
                            vDSDS = mDS[:,i]*mDS[:,k]
                            # compute precursors for second trace term in grad
                            mBG = np.matmul(np.multiply(self.data.mXT,repmat(self.data.vD*vDSDS,iK,1)),self.data.mX)
                            mBE = np.matmul(np.multiply(self.data.mXT,repmat(vDSDS,iK,1)),self.data.mX)
                            # store
                            mGrandBG[i*iK:(i+1)*iK,k*iK:(k+1)*iK] = mBG
                            mGrandBE[i*iK:(i+1)*iK,k*iK:(k+1)*iK] = mBE
                            if k > i:
                                mGrandBG[k*iK:(k+1)*iK,i*iK:(i+1)*iK] = mBG
                                mGrandBE[k*iK:(k+1)*iK,i*iK:(i+1)*iK] = mBE
            if bSameCovs: # if same covariates across traits
                if bGrad: # if gradient desired
                    # compute weighted sums mUG,mUE
                    mUG = np.matmul(np.multiply(mF,repmat(vDiagUG,iT,1)),mFT)
                    mUE = np.matmul(np.multiply(mF,repmat(vDiagUE,iT,1)),mFT)
                # compute inverse of mF
                mInvF = np.linalg.inv(mF)
                # compute Kronecker of mInvF and identity
                mKronInvFI = np.zeros((iT, iK, iT, iK), mInvF.dtype)
                for j in range(0,iK): mKronInvFI[:,j,:,j] = mInvF
                mKronInvFI = mKronInvFI.reshape(iK*iT,iK*iT)
                # compute inv(X'inv(V)X) and store
                mInvZ = np.matmul(np.matmul(mKronInvFI.T,mJ),mKronInvFI)
                self.mVarGLS = mInvZ
            else: # if different covariates
                # compute Kronecker of mF and identity
                mKronFI = np.zeros((iT, iK, iT, iK), mF.dtype)
                for j in range(0,iK): mKronFI[:,j,:,j] = mF
                mKronFI = mKronFI.reshape(iK*iT,iK*iT)
                # select appropriate submatrix of mKronFI
                mKronFI = mKronFI[vIndCovs,:]
                # finalise mZ
                mZ = np.matmul(np.matmul(mKronFI,mZ),mKronFI.T)
                # stabilise mZ and compute EVD
                mZ = (mZ + mZ.T)/2
                (vEigValsZ,mEigVecsZ) = np.linalg.eigh(mZ)
                # compute log|Z|
                dLogDetZ = np.log(vEigValsZ).sum()
                # compute inv(Z) and store
                mInvZ = np.matmul(np.multiply(mEigVecsZ,repmat(1/vEigValsZ,iKtotal,1)),mEigVecsZ.T)
                self.mVarGLS = mInvZ
                # compute mJ matrix
                mJ = np.matmul(mKronFI.T,np.matmul(mInvZ,mKronFI))
                # compute CGTF and CETF
                mCGTF = np.matmul(mCG.T,mF)
                mCETF = np.matmul(mCE.T,mF)
                # compute appropriate kroneckers of M with vectors of ones
                mKronCGTFiota = np.empty((iFG, iT, iK), mCGTF.dtype)
                mKronCETFiota = np.empty((iFE, iT, iK), mCETF.dtype)
                mKronFTiota   = np.empty((iT, iK, iT), mFT.dtype)
                mKronCGTFiota[...] = mCGTF[:, :, None]
                mKronCETFiota[...] = mCETF[:, :, None]
                mKronFTiota[...]   = mFT[:, None, :]
                mKronCGTFiota = mKronCGTFiota.reshape(iFG, iT*iK)            
                mKronCETFiota = mKronCETFiota.reshape(iFE, iT*iK)
                mKronFTiota   = mKronFTiota.reshape(iT*iK, iT)
                # compute contribution of second trace term to grad
                mTrG2 = np.matmul(np.matmul(mKronCGTFiota,np.multiply(mGrandBG,mJ)),mKronFTiota)
                mTrE2 = np.matmul(np.matmul(mKronCETFiota,np.multiply(mGrandBE,mJ)),mKronFTiota)
        # if the gradient is desired: compute weighted sums mTG,mTE
        if bGrad:
            mTG = np.matmul(np.multiply(mF,repmat(mDDS.sum(axis=0),iT,1)),mFT)
            mTE = np.matmul(np.multiply(mF,repmat(mDS.sum(axis=0),iT,1)),mFT)
        # if information matrix is desired: compute precursors
        if bInfo:
            if not(bSilent):
                # print statement
                print("Preparing for efficient calculation of AI matrix.")
                tAIprogress = tqdm(total=iN)
            if bCovs and not(bSameCovs):
                mFTCG = mCGTF.T
                mFTCE = mCETF.T
            else:
                mFTCG = np.matmul(mFT,mCG)
                mFTCE = np.matmul(mFT,mCE)
            # set 3d-arrays for inv(Vj)*C etc.
            m3InvVj      = np.zeros((iT,iT,iN))
            m3InvVjCG    = np.zeros((iT,iFG,iN))
            m3InvVjCE    = np.zeros((iT,iFE,iN))
            m3CGTInvVjCG = np.zeros((iFG,iFG,iN))
            m3CETInvVjCG = np.zeros((iFE,iFG,iN))
            m3CETInvVjCE = np.zeros((iFE,iFE,iN))
            # for each observation
            for j in range(0,iN):
                mFInvEVs      = np.multiply(mF,repmat(mDS[j,:],iT,1))
                mInvVjCG      = np.matmul(mFInvEVs,mFTCG)
                mInvVjCE      = np.matmul(mFInvEVs,mFTCE)
                m3InvVjCG[:,:,j] = mInvVjCG
                m3InvVjCE[:,:,j] = mInvVjCE
                m3InvVj[:,:,j]   = np.matmul(mFInvEVs,mFT)
                m3CGTInvVjCG[:,:,j] = np.matmul(mCG.T,mInvVjCG)
                m3CETInvVjCG[:,:,j] = np.matmul(mCE.T,mInvVjCG)
                m3CETInvVjCE[:,:,j] = np.matmul(mCE.T,mInvVjCE)
                if not(bSilent):
                    tAIprogress.update(1)
            if not(bSilent):
                tAIprogress.close()
        # compute YjTilde = Vj^(-1)Yj in one go
        mYjTilde = np.matmul(mF,np.multiply(mDS,np.matmul(self.data.mY,mF)).T)
        if bCovs: # if we have covariates
            if bSameCovs: # if identical across traits
                # compute the fixed effects and store
                vB = np.array(np.matmul(mInvZ,np.array(np.matmul(mYjTilde,self.data.mX)).ravel())).ravel()
                self.vBetaGLS = vB
            else:
                # set vector with fixed effects for all combinations of traits and params
                vB = np.zeros(iK*iT)
                # compute the fixed effects only for covariates that matter, and store
                vB[vIndCovs] = np.array(np.matmul(mInvZ,np.array(np.matmul(mYjTilde,self.data.mX)).ravel()[vIndCovs])).ravel()
                self.vBetaGLS = vB[vIndCovs]
            # set matrix of GLS residuals
            mR = np.zeros((iT,iN))
            # for each trait
            for i in range(0,iT):
                # get fixed-effects to start working on GLS resid
                mR[i,:] = np.array(np.matmul(self.data.mX,vB[i*iK:(i+1)*iK])).ravel()
            # finalise calculating of residuals
            mR = mYjTilde - (np.matmul(mF,np.multiply(mDS,np.matmul(mR.T,mF)).T))
        else:
            # set YjTilde as residual
            mR = mYjTilde
        # compute SSR in one go
        dSSR = np.multiply(mR,self.data.mY.T).sum()
        # compute log-likelihood per observation
        dLogL = (-0.5*(dCons + dLogDetV + dLogDetZ - self.data.dLogDetXTX + dSSR))/iN
        if bGrad: # if gradient is desired
            if bCovs: # if we have covariates
                if bSameCovs: # if identical across traits 
                    # compute trace matrices
                    mTrG = np.matmul((mTG-mUG),mCG)
                    mTrE = np.matmul((mTE-mUE),mCE)
                else:
                    mTrG = np.matmul(mTG,mCG) - mTrG2.T
                    mTrE = np.matmul(mTE,mCE) - mTrE2.T
            else:
                # compute trace matrices
                mTrG = np.matmul(mTG,mCG)
                mTrE = np.matmul(mTE,mCE)
            # compute SSR
            mCGTR = np.matmul(mCG.T,mR)
            mCETR = np.matmul(mCE.T,mR)
            mSSRG = np.matmul(np.multiply(mR,repmat(self.data.vD,iT,1)),mCGTR.T)
            mSSRE = np.matmul(mR,mCETR.T)
            # construct SSR vector for free parameters
            vSSRG = np.array(2*mSSRG[self.model.genmod.vIndT,self.model.genmod.vIndF]).T
            vTrG  = np.array(2*mTrG[self.model.genmod.vIndT,self.model.genmod.vIndF]).T
            vSSRE = np.array(2*mSSRE[self.model.envmod.vIndT,self.model.envmod.vIndF]).T
            vTrE  = np.array(2*mTrE[self.model.envmod.vIndT,self.model.envmod.vIndF]).T
            # compute gradient
            vGrad = (0.5*np.hstack((vSSRG-vTrG,vSSRE-vTrE)))/iN
            if bInfo: # if info matrix is desired
                # set info matrix
                mInfo = np.zeros((iParams,iParams))
                if bCovs: # if we have covariates
                    # set correction weight matrix
                    mW = np.zeros((iK*iT,iParams))
                if not(bSilent):
                    # print statement
                    print("Computing AI matrix.")
                    tAIprogress = tqdm(total=iParamCombos)
                # for each genetic parameter
                for g1 in range(0,iParamsG):
                    # get row and column index for parameter
                    iX = self.model.genmod.vIndT[g1]
                    iA = self.model.genmod.vIndF[g1]
                    # get fixed stuff
                    vR_X = np.array(mR[iX,:]).ravel()
                    vCGTR_A = np.array(mCGTR[iA,:]).ravel()
                    if bCovs: # if we have covariates
                        # get fixed stuff
                        vCTGQTilde_A = mFTCG[:,iA]
                        vQTilde_X = mF[iX,:]
                        # update correction weight matrix
                        mThisW = np.matmul(self.data.mXT,np.multiply(mDDS,(np.outer(vR_X,vCTGQTilde_A) + np.outer(vCGTR_A,vQTilde_X))))
                        # for each trait
                        for i in range(0,iT):
                            mW[i*iK:(i+1)*iK,g1] = np.array(mThisW[:,i]).ravel()
                    # for each other genetic parameter
                    for g2 in range(g1,iParamsG):
                        # get row and column index for other parameter
                        iY = self.model.genmod.vIndT[g2]
                        iB = self.model.genmod.vIndF[g2]
                        # get fixed stuff
                        vR_Y = np.array(mR[iY,:]).ravel()
                        vCGTR_B = np.array(mCGTR[iB,:]).ravel()
                        # get info from precursors
                        vInvVj_YX = m3InvVj[iX,iY,:]
                        vInvVjCG_YA = m3InvVjCG[iY,iA,:]
                        vInvVjCG_XB = m3InvVjCG[iX,iB,:]
                        vCGTInvVjCG_BA = m3CGTInvVjCG[iB,iA,:]
                        # compute element AI matrix and store
                        dElementAI = (self.data.vDSq*(vCGTInvVjCG_BA*vR_Y*vR_X + vCGTR_B*vInvVjCG_YA*vR_X + vCGTR_A*vInvVjCG_XB*vR_Y + vCGTR_B*vInvVj_YX*vCGTR_A)).sum()
                        mInfo[g1,g2] = dElementAI
                        if (g2 != g1):
                            mInfo[g2,g1] = dElementAI
                        if not(bSilent):
                            tAIprogress.update(1)
                    # for each other environment parameter
                    for g2 in range(0,iParamsE):
                        # get row and column index for other parameter
                        iY = self.model.envmod.vIndT[g2]
                        iB = self.model.envmod.vIndF[g2]
                        # get fixed stuff
                        vR_Y = np.array(mR[iY,:]).ravel()
                        vCETR_B = np.array(mCETR[iB,:]).ravel()
                        # get info from precursors
                        vInvVj_YX = m3InvVj[iX,iY,:]
                        vInvVjCG_YA = m3InvVjCG[iY,iA,:]
                        vInvVjCE_XB = m3InvVjCE[iX,iB,:]
                        vCETInvVjCG_BA = m3CETInvVjCG[iB,iA,:]
                        # compute element AI matrix and store
                        dElementAI = (self.data.vD*(vCETInvVjCG_BA*vR_Y*vR_X + vCETR_B*vInvVjCG_YA*vR_X + vCGTR_A*vInvVjCE_XB*vR_Y + vCETR_B*vInvVj_YX*vCGTR_A)).sum()
                        mInfo[g1,iParamsG+g2] = dElementAI
                        mInfo[iParamsG+g2,g1] = dElementAI
                        if not(bSilent):
                            tAIprogress.update(1)
                # for each environment parameter
                for g1 in range(0,iParamsE):
                    # get row and column index for parameter
                    iX = self.model.envmod.vIndT[g1]
                    iA = self.model.envmod.vIndF[g1]
                    # get fixed stuff
                    vR_X = np.array(mR[iX,:]).ravel()
                    vCETR_A = np.array(mCETR[iA,:]).ravel()
                    if bCovs: # if we have covariates
                        # get fixed stuff
                        vCETQTilde_A = mFTCE[:,iA]
                        vQTilde_X = mF[iX,:]
                        # update correction weight matrix
                        mThisW = np.matmul(self.data.mXT,np.multiply(mDS,(np.outer(vR_X,vCETQTilde_A) + np.outer(vCETR_A,vQTilde_X))))
                        # for each trait
                        for i in range(0,iT):
                            mW[i*iK:(i+1)*iK,iParamsG+g1] = np.array(mThisW[:,i]).ravel()
                    # for each other environment parameter
                    for g2 in range(g1,iParamsE):
                        # get row and column index for other parameter
                        iY = self.model.envmod.vIndT[g2]
                        iB = self.model.envmod.vIndF[g2]
                        # get fixed stuff
                        vR_Y = np.array(mR[iY,:]).ravel()
                        vCETR_B = np.array(mCETR[iB,:]).ravel()
                        # get info from precursors
                        vInvVj_YX = m3InvVj[iX,iY,:]
                        vInvVjCE_YA = m3InvVjCE[iY,iA,:]
                        vInvVjCE_XB = m3InvVjCE[iX,iB,:]
                        vCETInvVjCE_BA = m3CETInvVjCE[iB,iA,:]
                        # compute element AI matrix and store
                        dElementAI = (vCETInvVjCE_BA*vR_Y*vR_X + vCETR_B*vInvVjCE_YA*vR_X + vCETR_A*vInvVjCE_XB*vR_Y + vCETR_B*vInvVj_YX*vCETR_A).sum()
                        mInfo[iParamsG+g1,iParamsG+g2] = dElementAI
                        if (g2 != g1):
                            mInfo[iParamsG+g2,iParamsG+g1] = dElementAI
                        if not(bSilent):
                            tAIprogress.update(1)
                if not(bSilent):
                    tAIprogress.close()
                if bCovs: # if we have covariates
                    # finalise info matrix
                    mInfo = (0.5*(mInfo - np.matmul(mW.T,np.matmul(mJ,mW))))/iN
                else: # if we do not have covariates
                    # finalise info matrix
                    mInfo = (0.5*mInfo)/iN
                # stabilise info matrix
                mInfo = (mInfo + mInfo.T)/2
                # return logl, gradient, info matrix
                return dLogL, vGrad, mInfo
            else: # if info matrix is not desired
                # return logl, gradient
                return dLogL, vGrad
        else: # if grad not desired
            # return logl
            return dLogL
