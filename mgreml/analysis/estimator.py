import numpy as np
import pickle
import math
from numpy.matlib import repmat
from model import model

class MgremlEstimator:
    
    # set parameter-change threshold for ending golden section
    # and switching from Newton to gradient descent for one iteration
    dEps = 1E-9
    # set maximum number of iterations for golden section
    iMaxIterGold = 60
    # set the golden section fraction
    dTheta   = 2/(1+math.sqrt(5)) # golden ratio
    # newton requires gradient and info
    bGradNewton = True
    bInfoNewton = True
    # number of active parameter guesses golden section requires
    iGoldGuesses = 6
    # bfgs requires gradient
    bGradBFGS = True
    # threshold for dot products used for updating approximated inverse Hessian
    dTolDot = 1E-30
    # when doing Newton, let it be quiet about calculating AI matrix in
    # each iter
    bSilentNewton = True
    # set maximum number of iterations for Newton and BFGS
    iMaxIter = 2000
    
    def __init__(self, mdData, bNested = False):
        # if null model, read out appropriate attributes
        if bNested:
            dfGenBinFY = mdData.dfGenBinFY0
            dfEnvBinFY = mdData.dfEnvBinFY0
        else: # if alternative model, read out appropriate attributes
            dfGenBinFY = mdData.dfGenBinFY
            dfEnvBinFY = mdData.dfEnvBinFY
        # store the logger
        self.logger = mdData.logger
        # initialise mgreml model
        self.mgreml_model = model.MgremlModel(mdData, dfGenBinFY, dfEnvBinFY, bNested)
        # set iteration counter
        self.iIter = 0
        # indicate convergence has not occurred yet
        self.bNotConverged = True
        # indicate that estimates are still changing from one iter to the next
        self.bEstimatesChanged = True
        # indicate that estimates and statistics have not been finalised
        self.bDone = False
        # set whether to perform BFGS or Netwon
        self.bBFGS = mdData.bBFGS
        # set whether to store results from each iter
        self.bStoreIter = mdData.bStoreIter
        # set frequency with which to store results from iters
        self.iStoreIterFreq = mdData.iStoreIterFreq
        # set whether to compute standard errors
        self.bSEs = mdData.bSEs
        # set whether we are estimating a nested model
        self.bNested = bNested
        # set whether to return all parameters estimates
        # and sampling variance when done
        self.bAllCoeffs = mdData.bAllCoeffs
        # set whether to return all variance components
        # and their sampling variance when done
        self.bVarComp = mdData.bVarComp
        # set convergence threshold
        self.dGradTol = mdData.dGradTol
        # set prefix for storing iteration results
        self.sPrefix = mdData.sPrefix
        self.logger.info('Model initialised\n')
        
    def PerformEstimation(self):
        if self.bBFGS:
            self.PerformBFGS()
        else:
            self.PerformNewton()
        # compute statistics
        self.ComputeStatistics()
        self.logger.info('Model estimation complete\n')
    
    def ComputeGradLength(self):
        vIndT = self.mgreml_model.model.GetTraitIndices()
        vGradRescaled = self.vGrad.ravel()*((np.diag(self.mgreml_model.data.mCovY)**0.5)[vIndT])
        self.dGradLengthPerParam = np.sqrt(np.power(vGradRescaled,2).mean())
        
    def PerformBFGS(self):
        # use two strikes out system to reset inverse hessian to -I,
        # if two subsequent iterations produce to little change
        bStrike1 = False
        bStrike2 = False
        # print statement that BFGS starts now
        self.logger.info('Performing BFGS algorithm to find coefficients that maximise the MGREML log-likelihood')
        # compute logL and gradient for current estimates
        (self.dLogL, self.vGrad) = self.mgreml_model.ComputeLogLik(MgremlEstimator.bGradBFGS)
        # computer value of convergence criterion
        self.ComputeGradLength()
        # assess convergence
        self.bNotConverged = self.dGradLengthPerParam > self.dGradTol
        # initialise approx. inverse Hessian as -I, so 1st step is gradient descent
        self.mIH = -np.eye(self.mgreml_model.model.iParams)
        # while convergence has not occurred
        while self.bNotConverged:
            # update iteration counter
            self.iIter += 1
            if self.iIter > MgremlEstimator.iMaxIter:
                raise RuntimeError('Aborting BFGS algorithm after ' + str(MgremlEstimator.iMaxIter) + ' iterations. Have you set a too stringent convergence threshold using --grad-tol? We recommend using the default value of 1E-5.')
            self.logger.info('BFGS iteration ' + str(self.iIter) + ': log-likelihood per observation = ' + str(self.dLogL) + ', length of gradient per parameter = ' + str(self.dGradLengthPerParam))
            # if results are stored in each iteration
            if self.bStoreIter:
                if (self.iIter % self.iStoreIterFreq) == 0:
                    # print status and set filename for output
                    if self.bNested:
                        sFileName = self.sPrefix + 'estimates0.iter.' + str(self.iIter) + '.bfgs.pkl'
                    else:
                        sFileName = self.sPrefix + 'estimates.iter.' + str(self.iIter) + '.bfgs.pkl'
                    # construct output dict
                    dictOut = {'currentCombinedModel': self.mgreml_model.model, 'iIter': self.iIter, 'dLogL': self.dLogL, 'vGrad': self.vGrad, 'vBetaGLS': self.mgreml_model.vBetaGLS, 'mVarGLS': self.mgreml_model.mVarGLS}
                    # write output to pickled file
                    with open(sFileName, 'wb') as handle:
                        pickle.dump(dictOut, handle)
            # store the old gradient and parameters
            vGradOld  = self.vGrad
            vParamOld = self.mgreml_model.model.GetParams()
            # compute suggested new set of parameters
            vNew = vParamOld - np.array(np.matmul(self.mIH,self.vGrad)).ravel()
            # perform golden section to obtain new parameter estimates
            self.GoldenSection(vNew, MgremlEstimator.bGradBFGS)
            # computer value of convergence criterion
            self.ComputeGradLength()
            # assess convergence
            self.bNotConverged = self.dGradLengthPerParam > self.dGradTol
            # if not converged: update inverse Hessian
            if self.bNotConverged:
                # compute difference in gradient and parameters
                vY  = self.vGrad - vGradOld
                vS  = self.mgreml_model.model.GetParams() - vParamOld
                # compute dot product
                dSTY = np.dot(vS,vY)
                # if magnitude dot product is too low
                if abs(dSTY) < MgremlEstimator.dTolDot:
                    # if first strike has already occurred
                    if bStrike1:
                        # indicate this is second strike
                        bStrike2 = True
                        # print warning
                        self.logger.warning("Warning: too small change in two subsequent BFGS iterations; reinitialising approximate inverse Hessian")
                    else: # if first strike has not yet occurred
                        # indicate this is first strike
                        bStrike1 = True
                        # replace dot product by tolerance
                        dSTY = -MgremlEstimator.dTolDot
                        # print warning
                        self.logger.warning("Warning: too small change in BFGS iteration; may have disproportionate effect on approximate inverse Hessian")
                else:
                    # if change has occurred in this iteration,
                    # set both strike booleans to False
                    bStrike1 = False
                    bStrike2 = False
                # if no two strikes have been reached
                if not(bStrike2):
                    # compute rho
                    dR = 1/dSTY
                    # compute some intermediate vectors
                    vV = vS*dR
                    vW = np.matmul(self.mIH,vY)
                    # compute updated approximated inverse Hessian
                    self.mIH = self.mIH - np.outer(vV,vW) - np.outer(vW,vV) + np.outer(vV,vV)*np.dot(vW,vY) + np.outer(vV,vS)
                    # stabilise approximated inverse hessian
                    self.mIH = (self.mIH + self.mIH.T)/2
                else: # if this is the second strike
                    # reinitalise approximate inverse Hessian as -I
                    self.mIH = -np.eye(self.mgreml_model.model.iParams)
                    # set both strike booleans to False
                    bStrike1 = False
                    bStrike2 = False
        if self.bSEs:
            # print update
            self.logger.info("BFGS converged. Now computing covariance matrix of estimates.")
            self.logger.warning("Warning! This may take hours for a large sample size and a large number of traits")
            # compute information matrix
            (self.dLogL, self.vGrad, self.mInfo) = self.mgreml_model.ComputeLogLik(MgremlEstimator.bGradBFGS, self.bSEs)
            # scale information matrix back to normal scale
            self.mInfo = self.mInfo*self.mgreml_model.data.iN
        # scale dLogL and vGrad back to normal scale
        self.dLogL = self.dLogL*self.mgreml_model.data.iN
        self.vGrad = self.vGrad*self.mgreml_model.data.iN
        
    def PerformNewton(self):
        # set lowest eigenvalue permitted in info matrix per N to 1E-9
        dMinEigVal = 1E-9
        # print statement that BFGS starts now
        self.logger.info("Performing Newton algorithm to find coefficients that maximise the MGREML log-likelihood")
        # while convergence has not occurred
        while self.bNotConverged:
            # update iteration counter
            self.iIter += 1
            if self.iIter > MgremlEstimator.iMaxIter:
                raise RuntimeError('Aborting Newton algorithm after ' + str(MgremlEstimator.iMaxIter) + ' iterations. Have you set a too stringent convergence threshold using --grad-tol? We recommend using the default value of 1E-5.')
            if self.bEstimatesChanged: # if the estimates have changed
                # compute the log likelihoodper observation etc.
                (self.dLogL, self.vGrad, self.mInfo) = self.mgreml_model.ComputeLogLik(MgremlEstimator.bGradNewton,MgremlEstimator.bInfoNewton,bSilent = MgremlEstimator.bSilentNewton)
                # computer value of convergence criterion
                self.ComputeGradLength()
                # assess convergence
                self.bNotConverged = self.dGradLengthPerParam > self.dGradTol
                # if results are stored in each iteration
                if self.bStoreIter:
                    if (self.iIter % self.iStoreIterFreq) == 0:
                        # print status and set filenames for output
                        if self.bNested:
                            sFileName = self.sPrefix + 'estimates0.iter.' + str(self.iIter) + '.newton.pkl'
                        else:
                            sFileName = self.sPrefix + 'estimates.iter.' + str(self.iIter) + '.newton.pkl'
                self.logger.info('Newton iteration ' + str(self.iIter) + ': log-likelihood per observation = ' + str(self.dLogL) + ', length of gradient per parameter = ' + str(self.dGradLengthPerParam))
            else: # if they have not changed, defer testing of convergence for now
                # compute the log likelihood and gradient per observation
                (self.dLogL, self.vGrad) = self.mgreml_model.ComputeLogLik(MgremlEstimator.bGradNewton)
                # set convergence criterion and info matrix tot None
                self.mInfo = None
                self.dGradLengthPerParam = None
                # if results are stored in each iteration
                if self.bStoreIter:
                    if (self.iIter % self.iStoreIterFreq) == 0:
                        # print status and set filenames for output
                        if self.bNested:
                            sFileName = self.sPrefix + 'estimates0.iter.' + str(self.iIter) + '.grad_descent.pkl'
                        else:
                            sFileName = self.sPrefix + 'estimates.iter.' + str(self.iIter) + '.grad_descent.pkl'
                self.logger.warning("Warning: Estimates did not change in previous iteration, while gradient is too large")
                self.logger.info("Carrying out one gradient-descent step")
                self.logger.info('Gradient-descent iteration  ' + str(self.iIter) + ': log-likelihood per observation = ' + str(self.dLogL) + ', length of gradient per parameter = ' + str(self.dGradLengthPerParam))
            # if results are stored in each iteration
            if self.bStoreIter:
                if (self.iIter % self.iStoreIterFreq) == 0:
                    # construct output dict
                    dictOut = {'currentCombinedModel': self.mgreml_model.model, 'iIter': self.iIter, 'dLogL': self.dLogL, 'vGrad': self.vGrad, 'vBetaGLS': self.mgreml_model.vBetaGLS, 'mVarGLS': self.mgreml_model.mVarGLS}
                    # write output to pickled file
                    with open(sFileName, 'wb') as handle:
                        pickle.dump(dictOut, handle)
            # if not converged in this iteration: update parameters
            if self.bNotConverged:
                if self.bEstimatesChanged: # if the estimates have changed do Newton
                    # compute pseudo inverse of unconstrained part
                    (mInvInfo,_) = self.PseudoInvertSymmetricMat(self.mInfo,dMinEigVal)
                    # compute suggested new parameters
                    vNew = self.mgreml_model.model.GetParams() + np.array(np.matmul(mInvInfo,self.vGrad)).ravel()
                else: # if the estimates have not changed do gradient descent
                    # compute suggested new parameters
                    vNew = self.mgreml_model.model.GetParams() + np.array(self.vGrad).ravel()
                # perform golden section
                self.GoldenSection(vNew)
        # scale dLogL, vGrad and information matrix back to normal scale
        self.mInfo = self.mInfo*self.mgreml_model.data.iN
        self.dLogL = self.dLogL*self.mgreml_model.data.iN
        self.vGrad = self.vGrad*self.mgreml_model.data.iN
                
    def GoldenSection(self, vNew, bGradAtEnd = False):
        # set iteration counter to 0
        iIterGold  = 0
        # get current estimates and suggested estimates
        vParam0  = self.mgreml_model.model.GetParams()
        vParam1  = vNew
        # set x2 as convex combination of x0 and x1, leaning more closely to x0
        vParam2  = vParam1 - MgremlEstimator.dTheta * (vParam1 - vParam0)
        # set x3 as convex combination of x0 and x1, leaning more closely to x1
        vParam3  = vParam0 + MgremlEstimator.dTheta * (vParam1 - vParam0)
        # store the starting points of golden section as x4 and x5
        vParam4  = vParam1
        vParam5  = vParam0
        # define matrix of paramter estimates, each column
        # corresponding to a different guess
        mParam = np.vstack((vParam0,vParam1,vParam2,vParam3,vParam4,vParam5)).copy().T
        # define vector of corresponding log likelihoods
        vLogL  = np.zeros(MgremlEstimator.iGoldGuesses)
        # compute the log likelihood for x2 and x3
        vLogL[2] = self.mgreml_model.ComputeLogLik(vNew = mParam[:,2])
        vLogL[3] = self.mgreml_model.ComputeLogLik(vNew = mParam[:,3])
        # while the largest absolute difference between two parameters in vector x3 and x2 is larger than epsilon AND the maximum no. of iterations hasn't passed yet
        while (np.max(abs(mParam[:,2] - mParam[:,3])) >= MgremlEstimator.dEps) and (iIterGold < MgremlEstimator.iMaxIterGold):
            # if f2 > f3: go towards x0
            if (vLogL[2] > vLogL[3]):
                # set x3 as new x1: the upper bound
                mParam[:,1] = mParam[:,3]
                # set x2 as new x3: the point between the upper and lower bound leaning more closely towards to the upper bound
                mParam[:,3] = mParam[:,2]
                # assign the new x3 the proper log likelihood
                vLogL[3] = vLogL[2]
                # set new x2 as a convex combo of lower and upper bound, leaning more closely towards the lower bound
                mParam[:,2] = mParam[:,1] - MgremlEstimator.dTheta * (mParam[:,1] - mParam[:,0])
                # compute the log likelihood
                vLogL[2] = self.mgreml_model.ComputeLogLik(vNew = mParam[:,2])
            else: # if f2 <= f3: go towards x1
                # set x2 as new x0: the lower bound
                mParam[:,0] = mParam[:,2]
                # set x3 as new x2: the point between the upper and lower bound leaning more closely towards to the lower bound
                mParam[:,2] = mParam[:,3]
                # assign the new x2 the proper log likelihood
                vLogL[2] = vLogL[3]
                # set new x3 as a convex combo of lower and upper bound, leaning more closely towards the upper bound
                mParam[:,3] = mParam[:,0] + MgremlEstimator.dTheta * (mParam[:,1] - mParam[:,0])
                # compute the log likelihood
                vLogL[3] = self.mgreml_model.ComputeLogLik(vNew = mParam[:,3])
            # update iteration counter
            iIterGold += 1
        # compute missing log likelihoods
        vLogL[0] = self.mgreml_model.ComputeLogLik(vNew = mParam[:,0])
        vLogL[1] = self.mgreml_model.ComputeLogLik(vNew = mParam[:,1])
        vLogL[4] = self.mgreml_model.ComputeLogLik(vNew = mParam[:,4])
        vLogL[5] = self.mgreml_model.ComputeLogLik(vNew = mParam[:,5])
        # find highest log likelihood
        dMaxLogL  = vLogL.max()
        # find first index where vLogL equals its max
        iIndexMax = (np.array(np.where(vLogL == dMaxLogL)).ravel())[0]
        # commit results
        vNew = mParam[:,iIndexMax]
        self.dLogL = dMaxLogL
        self.bEstimatesChanged = (np.max(abs(vNew - vParam5)) >= MgremlEstimator.dEps)
        self.mgreml_model.model.UpdateParams(vNew)
        self.mInfo = None
        if bGradAtEnd:
            (self.dLogL, self.vGrad) = self.mgreml_model.ComputeLogLik(bGradAtEnd)
        else:
            self.vGrad = None

    def PseudoInvertSymmetricMat(self, mA, dMinEigVal):
        """
        Author     : Ronald de Vlaming
        Date       : December 3, 2020
        
        Summary    : computes pseudo-inverse of symmetric matrix
                     that is (nearly) positive definite, ignoring
                     eigenvalues too close to or just below zero
                     
        Input
        mA         : input matrix
                     
        Output
        mInvA      : pseudo-inverse (=regular inverse when all
                     eigenvalues of mA exceed dMinEigVal)
        bWarning   : boolean indicating if eigenvalues have been
                     ignored when inverting mA
        """
        # get dimensionality and EVD of mA
        iD = mA.shape[0]
        (vD,mP) = np.linalg.eigh(mA)
        # find eigenvalues less than dMinEigVal
        vIndices = np.where(vD<dMinEigVal)
        # count the number of eigenvalues that will be omitted
        iDropped = ((vD < dMinEigVal).astype(int)).sum()
        # eigenvalues below threshold? warning!
        bWarning = iDropped > 0
        if bWarning:
            # print warning
            self.logger.warning('Warning! ' + str(iDropped) + ' eigenvalues ignored in pseudo inverse of matrix')
            self.logger.info('Probably this is related to the AI matrix')
            self.logger.info('Phenotypes may be multicollinear or the current set of estimates may be poor')
            self.logger.info('Pay attention to whether this message persists up until the last iteration;')
            self.logger.info('If so, be very cautious in interpreting standard errors and covariance matrix of estimates!')
        # set the eigenvalues of the pseudo inverse equal to eigenvalues of inverse
        vDinv = 1/vD
        # except for problematic ones; set those eigenvalues in inverse equal to zero
        vDinv[vIndices] = 0
        # compute pseudo-inverse
        mInvA = np.matmul(np.multiply(mP,repmat(vDinv,iD,1)),mP.T)
        # return pseudo-inverse
        return mInvA, bWarning
    
    def IsConverged(self):
        return not(self.bNotConverged)
    
    def IsDone(self):
        return self.bDone
    
    def ComputeStatistics(self):
        # if not converged: raise error
        if self.bNotConverged:
            raise RuntimeError('Trying to calculate final statistics, while estimatese have not converged.')
        # get variance matrices and coefficient matrices
        (mVG, mCG, mVE, mCE) = self.mgreml_model.model.GetVandC()
        # get genetic variances, environment and total variances
        vVG = np.diag(mVG).copy()
        vVE = np.diag(mVE).copy()
        vVY = vVG + vVE
        # compute heritabilities
        self.vHSq = vVG / vVY
        # set variances to nan where zero for further calculations
        vVG[np.where(vVG==0)] = np.nan
        vVE[np.where(vVE==0)] = np.nan
        # compute one over the square root of variances
        vOneOverSqrtVG = np.power(vVG,-0.5)
        vOneOverSqrtVE = np.power(vVE,-0.5)
        mOneOverSqrtVG = np.outer(vOneOverSqrtVG,vOneOverSqrtVG)
        mOneOverSqrtVE = np.outer(vOneOverSqrtVE,vOneOverSqrtVE)
        # compute correlations
        self.mRhoG = np.multiply(mVG,mOneOverSqrtVG)
        self.mRhoE = np.multiply(mVE,mOneOverSqrtVE)
        if self.bSEs:
            # set lowest eigenvalue permitted in info matrix per N to 1E-18
            dMinEigValPerN = 1E-18
            # scale back to normal scale
            dMinEigVal = dMinEigValPerN*self.mgreml_model.data.iN
            # invert the information matrix
            (self.mSamplingV,_) = self.PseudoInvertSymmetricMat(self.mInfo, dMinEigVal)
            # get indices and parameters
            (vIndTG, vIndFG, vParamG, vIndTE, vIndFE, vParamE) = self.mgreml_model.model.GetSplitParamsAndIndices()
            # compute gradient of heritability w.r.t. each parameter
            vGradHSqG = 2*((1-self.vHSq[vIndTG])/vVY[vIndTG])*vParamG
            vGradHSqE = -2*(self.vHSq[vIndTE]/vVY[vIndTE])*vParamE
            vGradHSq  = np.hstack((vGradHSqG,vGradHSqE))
            # scale sampling variance matrix accordingly
            mGradHSq = np.outer(vGradHSq,vGradHSq)
            mSamplingVGradHsq = np.multiply(self.mSamplingV,mGradHSq)
            # define vector for SEs of heritabilities
            self.vHSqSE = np.zeros(self.mgreml_model.data.iT)
            # define matrix for SEs of correlations
            self.mRhoGSE = np.zeros((self.mgreml_model.data.iT,self.mgreml_model.data.iT))
            self.mRhoESE = np.zeros((self.mgreml_model.data.iT,self.mgreml_model.data.iT))
            # for each trait
            for i in range(0,self.mgreml_model.data.iT):
                # get the variances of X
                dVGX = vVG[i]
                dVEX = vVE[i]
                # find all parameter indices that correspond to current trait
                vIndGX = np.array(np.where(vIndTG==i)).ravel()
                vIndEX = np.array(np.where(vIndTE==i)).ravel()
                vIndX  = np.hstack((vIndGX,self.mgreml_model.model.iParamsG+vIndEX))
                # compute standard error of heritability estimate of given trait
                self.vHSqSE[i] = math.sqrt((mSamplingVGradHsq[vIndX,:][:,vIndX]).sum())
                # find all factors that affect X
                vIndFGX = vIndFG[vIndGX]
                vIndFEX = vIndFE[vIndEX]
                # and find corresponding parameters
                vParamGX = vParamG[vIndGX]
                vParamEX = vParamE[vIndEX]
                # get relevant submatrix of sampling variance matrix for SE correlation
                mSamplingVGXX = self.mSamplingV[vIndGX,:][:,vIndGX]
                mSamplingVEXX = self.mSamplingV[self.mgreml_model.model.iParamsG+vIndEX,:][:,self.mgreml_model.model.iParamsG+vIndEX]
                # for each other trait
                for j in range(i+1,self.mgreml_model.data.iT):
                    # get the variance of the other trait
                    dVGY = vVG[j]
                    dVEY = vVE[j]
                    # get the correlations
                    dRhoGXY = self.mRhoG[i,j]
                    dRhoEXY = self.mRhoE[i,j]
                    # for the factors affecting X find corresponding parameters of Y
                    vLambdaG = mCG[j,:][vIndFGX]
                    vLambdaE = mCE[j,:][vIndFEX]
                    # compute gradients of correlations with respect to parameters affecting X
                    vGradRhoGX = (1/math.sqrt(dVGX*dVGY))*vLambdaG - (dRhoGXY/dVGX)*vParamGX
                    vGradRhoEX = (1/math.sqrt(dVEX*dVEY))*vLambdaE - (dRhoEXY/dVEX)*vParamEX
                    # find all parameter indices that correspond to Y
                    vIndGY = np.array(np.where(vIndTG==j)).ravel()
                    vIndEY = np.array(np.where(vIndTE==j)).ravel()
                    # find all factors that affect Y
                    vIndFGY = vIndFG[vIndGY]
                    vIndFEY = vIndFE[vIndEY]
                    # and find corresponding parameters
                    vParamGY = vParamG[vIndGY]
                    vParamEY = vParamE[vIndEY]
                    # for the factors affecting Y find corresponding parameters of X
                    vLambdaG = mCG[i,:][vIndFGY]
                    vLambdaE = mCE[i,:][vIndFEY]
                    # compute gradients of correlations with respect to parameters affecting X
                    vGradRhoGY = (1/math.sqrt(dVGX*dVGY))*vLambdaG - (dRhoGXY/dVGY)*vParamGY
                    vGradRhoEY = (1/math.sqrt(dVEX*dVEY))*vLambdaE - (dRhoEXY/dVEY)*vParamEY
                    # get relevant submatrix of sampling variance matrix for SE of correlation
                    mSamplingVGXY = self.mSamplingV[vIndGX,:][:,vIndGY]
                    mSamplingVEXY = self.mSamplingV[self.mgreml_model.model.iParamsG+vIndEX,:][:,self.mgreml_model.model.iParamsG+vIndEY]
                    mSamplingVGYY = self.mSamplingV[vIndGY,:][:,vIndGY]
                    mSamplingVEYY = self.mSamplingV[self.mgreml_model.model.iParamsG+vIndEY,:][:,self.mgreml_model.model.iParamsG+vIndEY]
                    # compute standard errors of correlations
                    self.mRhoGSE[i,j] = math.sqrt(np.matmul(vGradRhoGX,np.matmul(mSamplingVGXX,vGradRhoGX.T)) + \
                                   2*np.matmul(vGradRhoGX,np.matmul(mSamplingVGXY,vGradRhoGY.T)) + \
                                   np.matmul(vGradRhoGY,np.matmul(mSamplingVGYY,vGradRhoGY.T)))
                    self.mRhoESE[i,j] = math.sqrt(np.matmul(vGradRhoEX,np.matmul(mSamplingVEXX,vGradRhoEX.T)) + \
                                   2*np.matmul(vGradRhoEX,np.matmul(mSamplingVEXY,vGradRhoEY.T)) + \
                                   np.matmul(vGradRhoEY,np.matmul(mSamplingVEYY,vGradRhoEY.T)))
                    self.mRhoGSE[j,i] = self.mRhoGSE[i,j]
                    self.mRhoESE[j,i] = self.mRhoESE[i,j]
        # if variance components are needed
        if self.bVarComp:
            # get indices and parameters
            (vIndTG, vIndFG, vParamG, vIndTE, vIndFE, vParamE) = self.mgreml_model.model.GetSplitParamsAndIndices()
            # compute how many VCs there are
            iSize = int(self.mgreml_model.data.iT*(self.mgreml_model.data.iT+1))
            # set labels for genetic and environment components
            self.lComponents = ['genetic covariance'] * int(iSize/2) + ['environment covariance'] * int(iSize/2)
            # if SEs are needed:
            if self.bSEs:
                # set matrices for covariance of variance components and VCs + SEs
                self.mSamplingVarVCs = np.zeros((iSize,iSize))
                self.mVCs = np.zeros((iSize,4))
            else:
                # set matrix for VCs
                self.mVCs = np.zeros((iSize,3))
            # initialise row indices
            iRowG = 0
            iRowE = int(iSize/2)
            # for each trait
            for i in range(0,self.mgreml_model.data.iT):
                # for each other trait
                for j in range(i,self.mgreml_model.data.iT):
                    # store trait indices and variance components
                    self.mVCs[iRowG,0] = i
                    self.mVCs[iRowG,1] = j
                    self.mVCs[iRowG,2] = mVG[i,j]
                    self.mVCs[iRowE,0] = i
                    self.mVCs[iRowE,1] = j
                    self.mVCs[iRowE,2] = mVE[i,j]
                    # update row indices
                    iRowG += 1
                    iRowE += 1
            # if SEs are needed:
            if self.bSEs:
                # re-initialise row indices
                iRowG = 0
                iRowE = int(iSize/2)
                # for each trait X1
                for i in range(0,self.mgreml_model.data.iT):
                    # find all parameter indices that correspond to X1
                    vIndGX1 = np.array(np.where(vIndTG==i)).ravel()
                    vIndEX1 = np.array(np.where(vIndTE==i)).ravel()
                    # find all factors that affect X1
                    vIndFGX1 = vIndFG[vIndGX1]
                    vIndFEX1 = vIndFE[vIndEX1]
                    # for each other trait Y1
                    for j in range(i,self.mgreml_model.data.iT):
                        # find all parameter indices that correspond to Y1
                        vIndGY1 = np.array(np.where(vIndTG==j)).ravel()
                        vIndEY1 = np.array(np.where(vIndTE==j)).ravel()
                        # find all factors that affect Y1
                        vIndFGY1 = vIndFG[vIndGY1]
                        vIndFEY1 = vIndFE[vIndEY1]
                        # for the factors affecting Y1 find corresponding parameters of X1
                        vGradGY1 = mCG[i,:][vIndFGY1]
                        vGradEY1 = mCE[i,:][vIndFEY1]
                        # for the factors affecting X1 find corresponding parameters of Y1
                        vGradGX1 = mCG[j,:][vIndFGX1]
                        vGradEX1 = mCE[j,:][vIndFEX1]
                        # initialise column indices
                        iColG = iRowG
                        iColE = iRowE
                        # for each other other trait X2
                        for k in range(i,self.mgreml_model.data.iT):
                            # find all parameter indices that correspond to X2
                            vIndGX2 = np.array(np.where(vIndTG==k)).ravel()
                            vIndEX2 = np.array(np.where(vIndTE==k)).ravel()
                            # find all factors that affect X2
                            vIndFGX2 = vIndFG[vIndGX2]
                            vIndFEX2 = vIndFE[vIndEX2]
                            # for each other other other trait Y2
                            for l in range(max(k,j),self.mgreml_model.data.iT):
                                # find all parameter indices that correspond to Y2
                                vIndGY2 = np.array(np.where(vIndTG==l)).ravel()
                                vIndEY2 = np.array(np.where(vIndTE==l)).ravel()
                                # find all factors that affect Y2
                                vIndFGY2 = vIndFG[vIndGY2]
                                vIndFEY2 = vIndFE[vIndEY2]
                                # for the factors affecting Y2 find corresponding parameters of X2
                                vGradGY2 = mCG[k,:][vIndFGY2]
                                vGradEY2 = mCE[k,:][vIndFEY2]
                                # for the factors affecting X2 find corresponding parameters of Y2
                                vGradGX2 = mCG[l,:][vIndFGX2]
                                vGradEX2 = mCE[l,:][vIndFEX2]
                                # get relevant submatrices of covariance matrix
                                mSamplingVGX1GX2 = self.mSamplingV[vIndGX1,:][:,vIndGX2]
                                mSamplingVGX1GY2 = self.mSamplingV[vIndGX1,:][:,vIndGY2]
                                mSamplingVGY1GX2 = self.mSamplingV[vIndGY1,:][:,vIndGX2]
                                mSamplingVGY1GY2 = self.mSamplingV[vIndGY1,:][:,vIndGY2]
                                mSamplingVEX1GX2 = self.mSamplingV[self.mgreml_model.model.iParamsG+vIndEX1,:][:,vIndGX2]
                                mSamplingVEX1GY2 = self.mSamplingV[self.mgreml_model.model.iParamsG+vIndEX1,:][:,vIndGY2]
                                mSamplingVEY1GX2 = self.mSamplingV[self.mgreml_model.model.iParamsG+vIndEY1,:][:,vIndGX2]
                                mSamplingVEY1GY2 = self.mSamplingV[self.mgreml_model.model.iParamsG+vIndEY1,:][:,vIndGY2]
                                mSamplingVGX1EX2 = self.mSamplingV[vIndGX1,:][:,self.mgreml_model.model.iParamsG+vIndEX2]
                                mSamplingVGX1EY2 = self.mSamplingV[vIndGX1,:][:,self.mgreml_model.model.iParamsG+vIndEY2]
                                mSamplingVGY1EX2 = self.mSamplingV[vIndGY1,:][:,self.mgreml_model.model.iParamsG+vIndEX2]
                                mSamplingVGY1EY2 = self.mSamplingV[vIndGY1,:][:,self.mgreml_model.model.iParamsG+vIndEY2]
                                mSamplingVEX1EX2 = self.mSamplingV[self.mgreml_model.model.iParamsG+vIndEX1,:][:,self.mgreml_model.model.iParamsG+vIndEX2]
                                mSamplingVEX1EY2 = self.mSamplingV[self.mgreml_model.model.iParamsG+vIndEX1,:][:,self.mgreml_model.model.iParamsG+vIndEY2]
                                mSamplingVEY1EX2 = self.mSamplingV[self.mgreml_model.model.iParamsG+vIndEY1,:][:,self.mgreml_model.model.iParamsG+vIndEX2]
                                mSamplingVEY1EY2 = self.mSamplingV[self.mgreml_model.model.iParamsG+vIndEY1,:][:,self.mgreml_model.model.iParamsG+vIndEY2]
                                # compute sampling covariances
                                self.mSamplingVarVCs[iRowG,iColG] = \
                                        np.matmul(vGradGX1,np.matmul(mSamplingVGX1GX2,vGradGX2.T)) + \
                                        np.matmul(vGradGX1,np.matmul(mSamplingVGX1GY2,vGradGY2.T)) + \
                                        np.matmul(vGradGY1,np.matmul(mSamplingVGY1GX2,vGradGX2.T)) + \
                                        np.matmul(vGradGY1,np.matmul(mSamplingVGY1GY2,vGradGY2.T))
                                self.mSamplingVarVCs[iColG,iRowG] = self.mSamplingVarVCs[iRowG,iColG]
                                self.mSamplingVarVCs[iRowE,iColG] = \
                                        np.matmul(vGradEX1,np.matmul(mSamplingVEX1GX2,vGradGX2.T)) + \
                                        np.matmul(vGradEX1,np.matmul(mSamplingVEX1GY2,vGradGY2.T)) + \
                                        np.matmul(vGradEY1,np.matmul(mSamplingVEY1GX2,vGradGX2.T)) + \
                                        np.matmul(vGradEY1,np.matmul(mSamplingVEY1GY2,vGradGY2.T))
                                self.mSamplingVarVCs[iColG,iRowE] = self.mSamplingVarVCs[iRowE,iColG]
                                self.mSamplingVarVCs[iRowG,iColE] = \
                                        np.matmul(vGradGX1,np.matmul(mSamplingVGX1EX2,vGradEX2.T)) + \
                                        np.matmul(vGradGX1,np.matmul(mSamplingVGX1EY2,vGradEY2.T)) + \
                                        np.matmul(vGradGY1,np.matmul(mSamplingVGY1EX2,vGradEX2.T)) + \
                                        np.matmul(vGradGY1,np.matmul(mSamplingVGY1EY2,vGradEY2.T))
                                self.mSamplingVarVCs[iColE,iRowG] = self.mSamplingVarVCs[iRowG,iColE]
                                self.mSamplingVarVCs[iRowE,iColE] = \
                                        np.matmul(vGradEX1,np.matmul(mSamplingVEX1EX2,vGradEX2.T)) + \
                                        np.matmul(vGradEX1,np.matmul(mSamplingVEX1EY2,vGradEY2.T)) + \
                                        np.matmul(vGradEY1,np.matmul(mSamplingVEY1EX2,vGradEX2.T)) + \
                                        np.matmul(vGradEY1,np.matmul(mSamplingVEY1EY2,vGradEY2.T))
                                self.mSamplingVarVCs[iColE,iRowE] = self.mSamplingVarVCs[iRowE,iColE]
                                # if same row as column: get SE from sampling variance
                                if iRowG == iColG:
                                    self.mVCs[iRowG,3] = np.sqrt(self.mSamplingVarVCs[iRowG,iColG])
                                    self.mVCs[iRowE,3] = np.sqrt(self.mSamplingVarVCs[iRowE,iColE])
                                # update column indices
                                iColG += 1
                                iColE += 1
                        # update row indices
                        iRowG += 1
                        iRowE += 1
        # indicate estimates are now done
        self.bDone = True
