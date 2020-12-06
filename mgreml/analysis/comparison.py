import numpy as np
import pandas as pd
from scipy.stats import chi2
from mgreml.analysis import estimator


class NestedEstimators:
    
    def __init__(self, mdData, dfNestedGenBinFY = None, dfNestedEnvBinFY = None, dfGenBinFY = None, dfEnvBinFY = None, bBFGS = True, bStoreIters= False, bSEs = False, bReturnFullModelSpecs = False):
        self.estimator_unres = estimator.MgremlEstimator(mdData, dfGenBinFY, dfEnvBinFY, bBFGS, bStoreIters, bSEs, bReturnFullModelSpecs)
        self.estimator_res = estimator.MgremlEstimator(mdData, dfNestedGenBinFY, dfNestedEnvBinFY, bBFGS, bStoreIters, bSEs, bReturnFullModelSpecs)
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
        
    def PerformEstimation(self):
        print('Estimating the nested model (null hypothesis):')
        self.estimator_res.PerformEstimation()
        print('Estimating the parent model (alternative hypothesis):')
        self.estimator_unres.PerformEstimation()
        # perform likelihood-ratio test
        self.PerformLRT()
    
    def IsConverged(self):
        return (self.estimator_res.IsConverged() & self.estimator_unres.IsConverged())
    
    def IsDone(self):
        return (self.estimator_res.IsDone() & self.estimator_unres.IsDone())
    
    def PerformLRT(self):
        print('Performing likelihood-ratio test with ' + str(self.iDF) + ' degrees of freedom:')
        if not(self.IsConverged()):
            raise RuntimeError('Estimates of the models have not converged.')
        self.dTestStat = -2*(self.estimator_res.dLogL - self.estimator_unres.dLogL)
        print('Chi-square test statistic is ' + str(self.dTestStat))
        self.dPval = 1-chi2.cdf(self.dTestStat, self.iDF)
        print('P-value = ' + str(self.dPval))
        
        
        
        