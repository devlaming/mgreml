import numpy as np
import pandas as pd
from scipy.stats import chi2
from mgreml.analysis import estimator

class ThreeEstimators:
    
    def __init__(self, mdData):
        # store the logger
        self.logger = mdData.logger
        self.logger.info('INITIALISING RESTRICTED MODEL WITH NO GENETIC VARIANCE FOR MEDIATOR')
        self.estimator1 = estimator.MgremlEstimator(mdData, bMedVarGM0 = True)
        self.logger.info('INITIALISING RESTRICTED MODEL WITH NO EFFECT MEDIATOR ON OUTCOME')
        self.estimator2 = estimator.MgremlEstimator(mdData, bMedBeta0 = True)
        self.logger.info('INITIALISING UNRESTRICTED MODEL WITH GENETIC MEDIATION')
        self.estimator3 = estimator.MgremlEstimator(mdData)
    
    def PerformEstimation(self):
        self.logger.info('Estimating the restricted model with no genetic variance of mediator:')
        self.estimator1.PerformEstimation()
        self.logger.info('Estimating the restricted model with no effect mediator on outcome:')
        self.estimator2.PerformEstimation()
        self.logger.info('Estimating the unrestricted model with genetic mediation:')
        self.estimator3.PerformEstimation()
        # perform likelihood-ratio test
        self.PerformLRT()
    
    def IsConverged(self):
        return (self.estimator1.IsConverged() & self.estimator2.IsConverged() & self.estimator3.IsConverged())
    
    def IsDone(self):
        return (self.estimator1.IsDone() & self.estimator2.IsDone() & self.estimator3.IsDone())
    
    def PerformLRT(self):
        self.logger.info('Performing likelihood-ratio test:')
        if not(self.IsConverged()):
            raise RuntimeError('Estimates of the models have not converged.')
        self.logger.info('Log-likelihood of restricted model with no genetic variance of mediator is ' + str(self.estimator1.dLogL))
        self.logger.info('Log-likelihood of restricted model with no effect mediator on outcome is ' + str(self.estimator2.dLogL))
        if self.estimator1.dLogL > self.estimator2.dLogL:
            self.logger.info('Supremum restricted models is achieved under model with no genetic variance of mediator: test has 2 degrees of freedom')
            self.iDF=2
            self.dTestStat = -2*(self.estimator1.dLogL - self.estimator3.dLogL)
        else:
            self.logger.info('Supremum restricted models is achieved under model with no effect mediator on outcome: test has 1 degree of freedom')
            self.iDF=1
            self.dTestStat = -2*(self.estimator2.dLogL - self.estimator3.dLogL)
        self.logger.info('Log-likelihood of unrestricted model with genetic mediation is ' + str(self.estimator3.dLogL))
        self.logger.info('Chi-square test statistic is ' + str(self.dTestStat))
        self.dPval = 1-chi2.cdf(self.dTestStat, self.iDF)
        self.logger.info('P-value = ' + str(self.dPval) + '\n')
