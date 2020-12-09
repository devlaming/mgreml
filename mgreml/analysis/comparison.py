import numpy as np
import pandas as pd
from scipy.stats import chi2
from analysis import estimator

class NestedEstimators:
    
    def __init__(self, mdData, dfGenBinFY0 = None, dfEnvBinFY0 = None, dfGenBinFY = None, dfEnvBinFY = None, bBFGS = True, bStoreIters= False, bSEs = False, bReturnFullModelSpecs = False):
        self.estimatorA = estimator.MgremlEstimator(mdData, dfGenBinFY, dfEnvBinFY, bBFGS, bStoreIters, bSEs, bReturnFullModelSpecs)
        self.estimator0 = estimator.MgremlEstimator(mdData, dfGenBinFY0, dfEnvBinFY0, bBFGS, bStoreIters, bSEs, bReturnFullModelSpecs)
        # initialisation of StructuralModel instances using
        # same mdData guarantees that the these StructuralModels
        # consider the same set of traits, in the same order
        # checking whether models are nest, is thus only about
        # checking factors and their free coefficients
        self.CheckModelsNested()
        
    def CheckModelsNested(self):
        ''' For each factor in the null model,
        find counterpart in alternative model,
        and assert for those factors, if there
        are free coefficients under the null,
        that are constrained under alternative.
        If so: the null model is not nested.'''
        # which coefficients are free for each part
        # of restricted and unrestricted model?
        (mBGA,mBEA) = self.estimatorA.mgreml_model.model.GetFreeCoeffs()
        (mBG0,mBE0) = self.estimator0.mgreml_model.model.GetFreeCoeffs()
        # convert factor labels to pandas indices
        indFGA = pd.Index(self.estimatorA.mgreml_model.model.genmod.lFactors)
        indFEA = pd.Index(self.estimatorA.mgreml_model.model.envmod.lFactors)
        indFG0 = pd.Index(self.estimator0.mgreml_model.model.genmod.lFactors)        
        indFE0 = pd.Index(self.estimator0.mgreml_model.model.envmod.lFactors)
        # find out if all factors in restricted genetic model
        # appear in unrestricted genetic model
        if not(indFG0.isin(indFGA).all()):
            raise ValueError('There is at least one genetic factor in your nested model that does not appear in your alternative model')        
        # find out if all factors in restricted environment model
        # appear in unrestricted environment model
        if not(indFE0.isin(indFEA).all()):
            raise ValueError('There is at least one environment factor in your nested model that does not appear in your alternative model')
        # select submatrices of free coeffs of unrestricted models
        # that align with free coeffs of restricted model
        mBGA_aligned = np.array(pd.DataFrame(mBGA,columns=indFGA)[indFG0])
        mBEA_aligned = np.array(pd.DataFrame(mBEA,columns=indFEA)[indFE0])
        # if there is at least one coefficient there where the altnerative
        # model is restricted, while the null model is free: not nested!
        if ((mBGA_aligned - mBG0) < 0).any():
            raise ValueError('There is at least one genetic coefficient where the nested model is free and the alternative model is not')
        if ((mBEA_aligned - mBE0) < 0).any():
            raise ValueError('There is at least one environment coefficient where the nested model is free and the alternative model is not')
        self.iDF = int((mBGA.sum() + mBEA.sum()) - (mBG0.sum() + mBE0.sum()))
        
    def PerformEstimation(self):
        print('Estimating the nested model (null hypothesis):')
        self.estimator0.PerformEstimation()
        print('Estimating the parent model (alternative hypothesis):')
        self.estimatorA.PerformEstimation()
        # perform likelihood-ratio test
        self.PerformLRT()
    
    def IsConverged(self):
        return (self.estimator0.IsConverged() & self.estimatorA.IsConverged())
    
    def IsDone(self):
        return (self.estimator0.IsDone() & self.estimatorA.IsDone())
    
    def PerformLRT(self):
        print('Performing likelihood-ratio test with ' + str(self.iDF) + ' degrees of freedom:')
        if not(self.IsConverged()):
            raise RuntimeError('Estimates of the models have not converged.')
        self.dTestStat = -2*(self.estimator0.dLogL - self.estimatorA.dLogL)
        print('Chi-square test statistic is ' + str(self.dTestStat))
        self.dPval = 1-chi2.cdf(self.dTestStat, self.iDF)
        print('P-value = ' + str(self.dPval))
        
        
        
        