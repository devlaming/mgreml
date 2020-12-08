"""
perform mgreml with set of input arguments
"""
import argparse
import logging
from data import reader

def main():
    # initialise parser
    parser = argparse.ArgumentParser()
    # initialise logger
    logger = logging.getLogger(__name__)
    # read in data
    myMgremlData = reader.MgremlReader(parser, logger)

if __name__ == '__main__':
    main()
