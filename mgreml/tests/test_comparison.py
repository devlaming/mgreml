import numpy as np
import pandas as pd
import pickle
from mgreml.analysis import comparison

def TestNestedEstimators():
    # compute info at end of BFGS together with standard errors
    bInfoAtEnd = True
    bSEs = True
    # set filenames
    sFilepath = './temp/'
    sFileName = sFilepath + 'MyMgremlData.pkl'
    # read GremlData instance
    with open(sFileName, 'rb') as handle:
        MyMgremlData = pickle.load(handle)
    # enforce zero genetic correlations: one unique factor per unique traits
    mGenBinFYres = np.eye(MyMgremlData.iT).astype(int)
    lFactors = ['factor ' + str(x) for x in np.arange(MyMgremlData.iT)]
    # construct dataframe
    dfGenBinFYres = pd.DataFrame(mGenBinFYres,index=MyMgremlData.lPhenos,columns=lFactors)
    # construct estimators for nested models
    MyNestedEstimators = comparison.NestedEstimators(MyMgremlData, dfNestedGenBinFY = dfGenBinFYres)
    # apply BFGS to nested models
    MyNestedEstimators.PerformBFGS(bInfoAtEnd, bSEs)
    # done testing
    print('Done testing the comparison of nested models using MGREML.')
