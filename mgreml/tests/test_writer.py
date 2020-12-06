import numpy as np
import pandas as pd
import pickle
from mgreml.analysis import comparison
from mgreml.data import writer

def TestMgremlEstimationAndWriter():
    # compute standard errors and use BFGS
    bBFGS = True
    bSEs = True
    # return complete specification
    bReturnFullModelSpecs = True
    # set filenames
    sFilepath = './temp/'
    sFileName = sFilepath + 'MyMgremlData.pkl'
    # read GremlData instance
    with open(sFileName, 'rb') as handle:
        MyMgremlData = pickle.load(handle)
    print('Estimating a single model.')
    # enforce zero genetic correlations: one unique factor per unique traits
    mGenBinFYres = np.eye(MyMgremlData.iT).astype(int)
    lFactors = ['factor ' + str(x) for x in np.arange(MyMgremlData.iT)]
    # construct dataframe
    dfGenBinFYres = pd.DataFrame(mGenBinFYres,index=MyMgremlData.lPhenos,columns=lFactors)
    # construct estimators for nested models
    MyNestedEstimators = comparison.NestedEstimators(MyMgremlData, dfNestedGenBinFY = dfGenBinFYres, bBFGS = bBFGS, bSEs = bSEs, bReturnFullModelSpecs = bReturnFullModelSpecs)
    print('Estimating a nested model versus its parent.')
    # apply BFGS to nested models
    MyNestedEstimators.PerformEstimation()
    # done testing
    print('Done estimating. Now writing results.')
    # create an Mgreml writer
    MyMgremlWriter = writer.DataWriter(MyNestedEstimators, sPrefix = 'myresults.')
    MyMgremlWriter.WriteRho()
    MyMgremlWriter.WriteHSq()
    MyMgremlWriter.WriteLRT()
    MyMgremlWriter.WriteLogLik()
    MyMgremlWriter.WriteEstimatesGLS()
    MyMgremlWriter.WriteModelCoefficients()
