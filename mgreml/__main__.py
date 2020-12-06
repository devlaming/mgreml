from mgreml.tests import test_aligner
from mgreml.tests import test_writer

if __name__ == '__main__':
    test_aligner.TestMgremlAligner()
    test_writer.TestMgremlEstimationAndWriter()
    print('Done testing')