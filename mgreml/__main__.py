from mgreml.tests import test_aligner
from mgreml.tests import test_estimator
from mgreml.tests import test_comparison

if __name__ == '__main__':
    test_aligner.TestMgremlAligner()
    test_estimator.TestMgremlEstimator()
    test_comparison.TestNestedEstimators()