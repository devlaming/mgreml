"""
perform mgreml with set of input arguments
"""
import argparse
import logging
import time
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
        # if analysis is implied
        if myMgremlReader.bAnalyse:
            # align the data
            logger.info('Aligning data and applying the canonical transformation')
            logger.info('This may take some time...')
            myMgremlData = aligner.MgremlData(myMgremlReader.dfY, myMgremlReader.dfA, myMgremlReader.dfX)
            logger.info('Initialising estimator')
            myMgremlEstimator = estimator.MgremlEstimator(myMgremlData, bSEs = True, bReturnFullModelSpecs = True)
            logger.info('Performing MGREML estimation')
            myMgremlEstimator.PerformEstimation()
            logger.info('Initialising file writer for results')
            myMgremlWriter = writer.DataWriter(myMgremlEstimator, sPrefix = myMgremlReader.sPrefix)
            logger.info('Writing results')
            myMgremlWriter.WriteHSq()
            myMgremlWriter.WriteRho()
            myMgremlWriter.WriteLogLik()
            myMgremlWriter.WriteEstimatesGLS()
            MyMgremlWriter.WriteModelCoefficients()
        else:
            logger.warning('No MGREML analysis will be carried out.')
            logger.info('mgreml_prepare.py -h describes options')
    except Exception:
        logger.error('MGREML failed.')
    finally:
        logger.info('Total time elapsed: {T}'.format(T=tools.sec_to_str(round(time.time() - dStartTime, 2))))

if __name__ == '__main__':
    main()
