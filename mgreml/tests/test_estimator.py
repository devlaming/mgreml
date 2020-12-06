import pickle
from mgreml.analysis import estimator

def TestMgremlEstimator():
    # use BFGS
    bBFGS = True
    # compute SEs at end
    bSEs = True
    # set filenames
    sFilepath = './temp/'
    sFileName = sFilepath + 'MyMgremlData.pkl'
    # read GremlData instance
    with open(sFileName, 'rb') as handle:
        MyMgremlData = pickle.load(handle)
    # create an Mgreml estimator
    MyMgremlEstimator = estimator.MgremlEstimator(MyMgremlData, bBFGS = bBFGS, bSEs = bSEs)
    # perform BFGS
    MyMgremlEstimator.PerformEstimation()
    # calculate final statistics
    MyMgremlEstimator.ComputeStatistics()
    # done testing
    print('Done testing the MGREML estimator.')
