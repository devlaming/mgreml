import argparse
import logging
import time
import traceback
from data import tools
from data import reader
from data import aligner
from data import writer
from analysis import estimator
from analysis import comparison

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
                logger.info('4. INTIALISING MGREML ESTIMATOR FOR MAIN MODEL')
                myEstimator = estimator.MgremlEstimator(myMgremlData)
            logger.info('5. PERFORMING MGREML ESIMATION')
            myEstimator.PerformEstimation()
            logger.info('6. WRITING MGREML RESULTS')
            myMgremlWriter = writer.DataWriter(myEstimator, myMgremlData)
            myMgremlWriter.WriteResults()
            logger.info('MGREML HAS FINISHED.')
        else:
            logger.warning('No MGREML analysis will be carried out.')
            logger.info('mgreml_prepare.py -h describes options')
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
