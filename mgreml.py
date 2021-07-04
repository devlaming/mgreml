import argparse
import logging
import time
import traceback
import os, psutil
from mgreml.data import reader
from mgreml.data import writer
from mgreml.data import tools
from mgreml.analysis import estimator
from mgreml.analysis import comparison
from mgreml.analysis import pairwise

def main():
    # store starting time
    dStartTime = time.time()
    # initialise parser
    parser = argparse.ArgumentParser()
    # initialise logger
    logger = logging.getLogger(__name__)
    # initialise memory tracker
    process = psutil.Process(os.getpid())
    # try MGREML
    try:
        # read in data
        myMgremlData = reader.MgremlReader(parser, logger, process)
        # if analysis is implied by reader
        if myMgremlData.bAnalyse:
            if myMgremlData.bPairwise: # if we have a pairwise model
                # initialise pairwise estimator
                logger.info('3. INTIALISING MGREML ESTIMATORS FOR PAIRWISE BIVARIATE MODELS')
                logger.info('Current memory usage is ' + str(int((process.memory_info().rss)/(1024**2))) + 'MB')
                myEstimator = pairwise.PairwiseEstimators(myMgremlData)
                logger.info('4A. PERFORMING PAIRWISE BIVARIATE MGREML ESIMATION')
                logger.info('Current memory usage is ' + str(int((process.memory_info().rss)/(1024**2))) + 'MB\n')
                myEstimator.PerformEstimation()
            else:
                if myMgremlData.bNested: # if we have a nested model
                    # initialise nested estimators
                    logger.info('3. INTIALISING MGREML ESTIMATORS FOR NESTED MODELS')
                    logger.info('Current memory usage is ' + str(int((process.memory_info().rss)/(1024**2))) + 'MB\n')
                    myEstimator = comparison.NestedEstimators(myMgremlData)
                else: # if neither nested nor pairwise
                    # initialise estimator of the main model
                    if myMgremlData.bReinitialise:
                        logger.info('3. REINITIALISING MGREML ESTIMATOR FOR MAIN MODEL') 
                    else:
                        logger.info('3. INTIALISING MGREML ESTIMATOR FOR MAIN MODEL')
                    logger.info('Current memory usage is ' + str(int((process.memory_info().rss)/(1024**2))) + 'MB')
                    myEstimator = estimator.MgremlEstimator(myMgremlData)
                logger.info('4. PERFORMING MGREML ESIMATION')
                logger.info('Current memory usage is ' + str(int((process.memory_info().rss)/(1024**2))) + 'MB')
                myEstimator.PerformEstimation()
                myMgremlWriter = writer.DataWriter(myEstimator, myMgremlData)
                myMgremlWriter.WriteResults()
            logger.info('MGREML HAS FINISHED.')
        else:
            logger.warning('No MGREML analysis will be carried out.')
            logger.info('python ./mgreml.py -h shows all options')
    except Exception:
        # print the traceback
        logger.error(traceback.format_exc())
        # wrap up with final error message
        logger.error('Error: MGREML did not exit properly. Please inspect the log file.')
    finally:
        # print total time elapsed for MGREML analysis
        logger.info('Total time elapsed: {T}'.format(T=tools.sec_to_str(round(time.time() - dStartTime, 2))))
        logger.info('Current memory usage is ' + str(int((process.memory_info().rss)/(1024**2))) + 'MB')

if __name__ == '__main__':
    main()
