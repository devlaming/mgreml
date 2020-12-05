import numpy as np
import pickle
import math
from numpy.matlib import repmat
from mgreml.model import model

class MgremlEstimator:
    
    # set the threshold for length gradient per N per parameter
    dGradLengthTolPerParam = 1E-5
    # set lowest eigenvalue permitted in info matrix per N to 1E-9
    dMinEigVal = 1E-9
    # set parameter-change threshold for ending golden section
    # and switching from Newton to gradient descent for one iteration
    dEps = 1E-9
    # set maximum number of iterations for golden section
    iMaxIter = 60
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
    
    def __init__(self, mdData, dfGenBinFY = None, dfEnvBinFY = None, bStoreIters= False):
        self.mgreml_model = model.MgremlModel(mdData, dfGenBinFY, dfEnvBinFY)
        # set iteration counter
        self.iIter = 0
        # indicate convergence has not occurred yet
        self.bNotConverged = True
        # indicate that estimates are still changing from one iter to the next
        self.bEstimatesChanged = True
        # set whether to store results from each iter
        self.bStoreIters = bStoreIters
        
    def PerformBFGS(self, bInfoAtEnd = False):
        # use two strikes out system to reset inverse hessian to -I,
        # if two subsequent iterations produce to little change
        bStrike1 = False
        bStrike2 = False
        # print statement that BFGS starts now
        print("Performing BFGS algorithm to find coefficients that maximise the MGREML log-likelihood")
        # compute logL and gradient for current estimates
        (self.dLogL, self.vGrad) = self.mgreml_model.ComputeLogLik(MgremlEstimator.bGradBFGS)
        # computer value of convergence criterion
        self.dGradLengthPerParam = np.sqrt(np.power(self.vGrad.ravel(),2).mean())
        # assess convergence
        self.bNotConverged = self.dGradLengthPerParam > MgremlEstimator.dGradLengthTolPerParam
        # initialise approx. inverse Hessian as -I, so 1st step is gradient descent
        self.mIH = -np.eye(self.mgreml_model.model.iParams)
        # while convergence has not occurred
        while self.bNotConverged:
            # update iteration counter
            self.iIter += 1
            print("BFGS iteration",self.iIter,": log-likelihood per observation =",self.dLogL,", length of gradient per parameter = ",self.dGradLengthPerParam)
            # if results are stored in each iteration
            if self.bStoreIters:
                # print status and set filename for output
                sFileName = 'ESTIMATES.ITER.' + str(self.iIter) + '.BFGS.pkl'
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
            self.dGradLengthPerParam = np.sqrt(np.power(self.vGrad.ravel(),2).mean())
            # assess convergence
            self.bNotConverged = self.dGradLengthPerParam > MgremlEstimator.dGradLengthTolPerParam
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
                        print("Warning: too small change in two subsequent BFGS iterations; reinitialising approximate inverse Hessian")
                    else: # if first strike has not yet occurred
                        # indicate this is first strike
                        bStrike1 = True
                        # replace dot product by tolerance
                        dSTY = -MgremlEstimator.dTolDot
                        # print warning
                        print("Warning: too small change in BFGS iteration; may have disproportionate effect on approximate inverse Hessian")
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
        if bInfoAtEnd:
            # print update
            print("BFGS CONVERGED. NOW COMPUTING SAMPLING COVARIANCE MATRIX OF BFGS ESTIMATES.")
            print("WARNING! THIS MAY TAKE HOURS FOR LARGE SAMPLE SIZE AND LARGE NUMBER OF TRAITS!")
            # compute information matrix
            (self.dLogL, self.vGrad, self.mInfo) = self.mgreml_model.ComputeLogLik(MgremlEstimator.bGradBFGS, bInfoAtEnd)
        
    def PerformNewton(self):
        # print statement that BFGS starts now
        print("Performing Newton algorithm to find coefficients that maximise the MGREML log-likelihood")
        # while convergence has not occurred
        while self.bNotConverged:
            # update iteration counter
            self.iIter += 1
            if self.bEstimatesChanged: # if the estimates have changed
                # compute the log likelihoodper observation etc.
                (self.dLogL, self.vGrad, self.mInfo) = self.mgreml_model.ComputeLogLik(MgremlEstimator.bGradNewton,MgremlEstimator.bInfoNewton,bSilent = MgremlEstimator.bSilentNewton)
                # computer value of convergence criterion
                self.dGradLengthPerParam = np.sqrt(np.power(self.vGrad.ravel(),2).mean())
                # assess convergence
                self.bNotConverged = self.dGradLengthPerParam > MgremlEstimator.dGradLengthTolPerParam
                # if results are stored in each iteration
                if self.bStoreIters:
                    # print status and set filenames for output
                    sFileName = 'ESTIMATES.ITER.' + str(self.iIter) + '.NEWTON.pkl'
                print("Newton iteration",self.iIter,": log-likelihood per observation =",self.dLogL,", length of gradient per parameter = ",self.dGradLengthPerParam)
            else: # if they have not changed, defer testing of convergence for now
                # compute the log likelihood and gradient per observation
                (self.dLogL, self.vGrad) = self.mgreml_model.ComputeLogLik(MgremlEstimator.bGradNewton)
                # set convergence criterion and info matrix tot None
                self.mInfo = None
                self.dGradLengthPerParam = None
                # if results are stored in each iteration
                if self.bStoreIters:
                    # print status and set filenames for output
                    sFileName = 'ESTIMATES.ITER.' + str(self.iIter) + '.GRAD_DESCENT.pkl'
                print("WARNING: Estimates did not change in previous iteration, while gradient is too large")
                print("Carrying out one gradient-descent step")
                print("Gradient-descent iteration",self.iIter,": log-likelihood per observation =",self.dLogL,", length of gradient per parameter = ",self.dGradLengthPerParam)
            # if results are stored in each iteration
            if self.bStoreIters:
                # construct output dict
                dictOut = {'currentCombinedModel': self.mgreml_model.model, 'iIter': self.iIter, 'dLogL': self.dLogL, 'vGrad': self.vGrad, 'vBetaGLS': self.mgreml_model.vBetaGLS, 'mVarGLS': self.mgreml_model.mVarGLS}
                # write output to pickled file
                with open(sFileName, 'wb') as handle:
                    pickle.dump(dictOut, handle)
            # if not converged in this iteration: update parameters
            if self.bNotConverged:
                if self.bEstimatesChanged: # if the estimates have changed do Newton
                    # compute pseudo inverse of unconstrained part
                    (mInvInfo,_) = MgremlEstimator.PseudoInvertSymmetricMat(self.mInfo)
                    # compute suggested new parameters
                    vNew = self.mgreml_model.model.GetParams() + np.array(np.matmul(mInvInfo,self.vGrad)).ravel()
                else: # if the estimates have not changed do gradient descent
                    # compute suggested new parameters
                    vNew = self.mgreml_model.model.GetParams() + np.array(self.vGrad).ravel()
                # perform golden section
                self.GoldenSection(vNew)
                
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
        while (np.max(abs(mParam[:,2] - mParam[:,3])) >= MgremlEstimator.dEps) and (iIterGold < MgremlEstimator.iMaxIter):
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
        self.bEstimatesChanged = (np.max(abs(vParam4 - vParam5)) >= MgremlEstimator.dEps)
        self.mgreml_model.model.UpdateParams(vNew)
        self.mInfo = None
        if bGradAtEnd:
            (self.dLogL, self.vGrad) = self.mgreml_model.ComputeLogLik(bGradAtEnd)
        else:
            self.vGrad = None

    @staticmethod
    def PseudoInvertSymmetricMat(mA):
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
        vIndices = np.where(vD<MgremlEstimator.dMinEigVal)
        # count the number of eigenvalues that will be omitted
        iDropped = ((vD < MgremlEstimator.dMinEigVal).astype(int)).sum()
        # eigenvalues below threshold? warning!
        bWarning = iDropped > 0
        if bWarning:
            # print warning
            print("WARNING!",iDropped,"eigenvalues ignored in pseudo inverse of matrix; probably this is related to the AI matrix; phenotypes may be multicollinear or the current set of estimates may be poor; pay attention to whether this message persists up until the last iteration; if so, be very cautious in interpreting sampling variance matrix and standard errors of estimates!")
        # set the eigenvalues of the pseudo inverse equal to eigenvalues of inverse
        vDinv = 1/vD
        # except for problematic ones; set those eigenvalues in inverse equal to zero
        vDinv[vIndices] = 0
        # compute pseudo-inverse
        mInvA = np.matmul(np.multiply(mP,repmat(vDinv,iD,1)),mP.T)
        # return pseudo-inverse
        return mInvA, bWarning
