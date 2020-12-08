"""
perform mgreml with set of input arguments
"""
import argparse
import logging
import time
from data import reader

def main():
    # store starting time
    dStart = time.time()
    # initialise parser
    parser = argparse.ArgumentParser()
    # initialise logger
    logger = logging.getLogger(__name__)
    # read in data
    myMgremlData = reader.MgremlReader(parser, logger, dStart)

if __name__ == '__main__':
    main()
