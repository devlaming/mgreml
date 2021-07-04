import numpy as np
import pandas as pd
from scipy.stats import chi2
from mgreml.data import writer
from mgreml.data import reader
from mgreml.analysis import estimator
from mgreml.analysis import comparison

class PairwiseEstimators:
    
    def __init__(self, mdData):
        # store the data
        self.mdData = mdData
        self.mdData.logger.info('Model initialised')
        self.mdData.logger.info('Current memory usage is ' + str(int((self.mdData.process.memory_info().rss)/(1024**2))) + 'MB\n')
        
    def PerformEstimation(self):
        for i in range(0,self.mdData.iT):
            # get label trait i
            sLabelI = self.mdData.lPhenos[i]
            for j in range(i+1,self.mdData.iT):
                # get label trait i
                sLabelJ = self.mdData.lPhenos[j]
                # print update
                self.mdData.logger.info('4B. PERFORMING ESIMATION FOR ' + sLabelI + ' AND ' + sLabelJ)
                self.mdData.logger.info('Fetching relevant subset of data\n')
                # copy data on all traits
                mdThisData = reader.PairwiseMgremlReader(self.mdData,i,j)
                # if nested analysis
                if mdThisData.bNested:
                    # initialise nested estimators: no check on nestedness required given setup of this model
                    myEstimator = comparison.NestedEstimators(mdThisData, bCheck=False)
                else: # else
                    # initialise main estimator
                    myEstimator = estimator.MgremlEstimator(mdThisData)
                # perform estimation and write results
                myEstimator.PerformEstimation()
                myMgremlWriter = writer.DataWriter(myEstimator, mdThisData)
                myMgremlWriter.WriteResults()
