import argparse
import logging
import time
import traceback
from mgreml.data import reader
from mgreml.data import aligner
from mgreml.data import writer
from mgreml.data import tools
from mgreml.analysis import estimator
from mgreml.analysis import comparison

def main():
    # store starting time
    dStartTime = time.time()
    # initialise parser
    parser = argparse.ArgumentParser()
    # initialise logger
    logger = logging.getLogger(__name__)
    # try MGREML
    try:
        # read in data
        myMgremlReader = reader.MgremlReader(parser, logger)
        # if analysis is implied by reader
        if myMgremlReader.bAnalyse:
            # align the data
            myMgremlData = aligner.MgremlData(myMgremlReader)
            # if we have a nested model
            if myMgremlData.bNested:
                # initialise nested estimators
                logger.info('4. INTIALISING MGREML ESTIMATORS FOR NESTED MODELS')
                myEstimator = comparison.NestedEstimators(myMgremlData)
            else:
                # initialise estimator of the main model
                if myMgremlData.bReinitialise:
                    logger.info('4. REINITIALISING MGREML ESTIMATOR FOR MAIN MODEL') 
                else:
                    logger.info('4. INTIALISING MGREML ESTIMATOR FOR MAIN MODEL')
                myEstimator = estimator.MgremlEstimator(myMgremlData)
            logger.info('5. PERFORMING MGREML ESIMATION')
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

if __name__ == '__main__':
    main()
