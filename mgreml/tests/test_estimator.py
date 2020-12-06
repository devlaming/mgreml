import pickle
from mgreml.analysis import estimator

def TestMgremlEstimator():
    # compute info at end of BFGS
    bInfoAtEnd = True
    # compute SEs at end of BFGS
    bSEs = True
    # set filenames
    sFilepath = './temp/'
    sFileName = sFilepath + 'MyMgremlData.pkl'
    # read GremlData instance
    with open(sFileName, 'rb') as handle:
        MyMgremlData = pickle.load(handle)
    # create an Mgreml estimator
    MyMgremlEstimator = estimator.MgremlEstimator(MyMgremlData)
    # perform BFGS
    MyMgremlEstimator.PerformBFGS(bInfoAtEnd)
    # calculate final statistics
    MyMgremlEstimator.ComputeStatistics(bSEs)
    # done testing
    print('Done testing the MGREML estimator.')
