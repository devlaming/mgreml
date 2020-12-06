import numpy as np
import pandas as pd
from scipy.stats import chi2
from mgreml.analysis import estimator


class NestedEstimators:
    
    def __init__(self, mdData, dfNestedGenBinFY = None, dfNestedEnvBinFY = None, dfGenBinFY = None, dfEnvBinFY = None, bStoreIters= False):
        self.estimator_unres = estimator.MgremlEstimator(mdData, dfGenBinFY, dfEnvBinFY, bStoreIters)
        self.estimator_res = estimator.MgremlEstimator(mdData, dfNestedGenBinFY, dfNestedEnvBinFY, bStoreIters)        
        self.CheckModelsNested()
        
    def CheckModelsNested(self):
        # find out which coefficients are free for genetic and env part
        # for restricted and unrestricted models
        (mBGunres,mBEunres) = self.estimator_unres.mgreml_model.model.GetFreeCoeffs()
        (mBGres,mBEres) = self.estimator_res.mgreml_model.model.GetFreeCoeffs()
        # convert factor labels to pandas indices
        indFGunres = pd.Index(self.estimator_unres.mgreml_model.model.genmod.lFactors)
        indFEunres = pd.Index(self.estimator_unres.mgreml_model.model.envmod.lFactors)
        indFGres = pd.Index(self.estimator_res.mgreml_model.model.genmod.lFactors)        
        indFEres = pd.Index(self.estimator_res.mgreml_model.model.envmod.lFactors)
        # find out if all factors in restricted genetic model
        # appear in unrestricted genetic model
        if not(indFGres.isin(indFGunres).all()):
            raise ValueError('There is at least one genetic factor in your nested model that does not appear in your alternative model')        
        # find out if all factors in restricted environment model
        # appear in unrestricted environment model
        if not(indFEres.isin(indFEunres).all()):
            raise ValueError('There is at least one environment factor in your nested model that does not appear in your alternative model')
        # select submatrices of free coeffs of unrestricted models
        # that align with free coeffs of restricted model
        mBGunres_aligned = np.array(pd.DataFrame(mBGunres,columns=indFGunres)[indFGres])
        mBEunres_aligned = np.array(pd.DataFrame(mBEunres,columns=indFEunres)[indFEres])
        if ((mBGunres_aligned - mBGres) < 0).any():
            raise ValueError('There is at least one genetic coefficient where the nested model is free and the alternative model is not')
        if ((mBEunres_aligned - mBEres) < 0).any():
            raise ValueError('There is at least one environment coefficient where the nested model is free and the alternative model is not')
        self.iDF = int((mBGunres.sum() + mBEunres.sum()) - (mBGres.sum() + mBEres.sum()))
        
    def PerformBFGS(self, bInfoAtEnd = False, bSEs = False, bSamplingV = False):
        print('Estimating the nested model (null hypothesis):')
        self.estimator_res.PerformBFGS(bInfoAtEnd)
        self.estimator_res.ComputeStatistics(bSEs, bSamplingV)
        print('Estimating the alternative model (alternative hypothesis):')
        self.estimator_unres.PerformBFGS(bInfoAtEnd)
        self.estimator_unres.ComputeStatistics(bSEs, bSamplingV)
        self.PerformLRT()
        
    def PerformNewton(self, bSEs = False, bSamplingV = False):
        print('Estimating the nested model (null hypothesis):')
        self.estimator_res.PerformNewton()
        self.estimator_res.ComputeStatistics(bSEs, bSamplingV)
        print('Estimating the alternative model (alternative hypothesis):')
        self.estimator_unres.PerformNewton()
        self.estimator_unres.ComputeStatistics(bSEs, bSamplingV)
        self.PerformLRT()
        
    def PerformLRT(self):
        print('Performing likelihood-ratio test with ' + str(self.iDF) + ' degrees of freedom:')
        if self.estimator_res.bNotConverged:
            raise RuntimeError('Estimation of the nested model has not converged.')
        if self.estimator_unres.bNotConverged:
            raise RuntimeError('Estimation of the alternative model has not converged.')
        self.dTestStat = -2*(self.estimator_res.dLogL - self.estimator_unres.dLogL)
        print('Chi-square test statistic is ' + str(self.dTestStat))
        self.dPval = 1-chi2.cdf(self.dTestStat, self.iDF)
        print('P-value = ' + str(self.dPval))
        
        
        
        