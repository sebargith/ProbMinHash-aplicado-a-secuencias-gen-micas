#######################################
# Copyright (C) 2019-2020 Otmar Ertl. #
# All rights reserved.                #
#######################################

from scipy.stats import chisquare
from scipy.stats import kstest
from scipy.stats import binom_test
from numpy import expm1

def testUniformIntegerDistribution(fileName, n, significanceLevel):
    with open(fileName) as f:
        content = f.readlines()
    values = [int(x) for x in content]
    histogram = [0] * n
    for v in values:
        histogram[v] += 1

    result = chisquare(histogram)
    pValue = result[1]
    assert(pValue > significanceLevel)


def testUniformDoubleDistribution(fileName, significanceLevel):
    with open(fileName) as f:
        content = f.readlines()
    values = [float(x) for x in content]
    _, pValue = kstest(values, "uniform")
    assert(pValue > significanceLevel)

def testUniformDoubleHalfDistribution(fileName, significanceLevel):
    with open(fileName) as f:
        content = f.readlines()
    values = [2*float(x) for x in content]
    _, pValue = kstest(values, "uniform")
    assert(pValue > significanceLevel)

def testExponentialDistribution(fileName, significanceLevel):
    with open(fileName) as f:
        content = f.readlines()
    values = [float(x) for x in content]
    _, pValue = kstest(values, "expon")
    assert(pValue > significanceLevel)


def testBernoulliDistribution(fileName, successProbability, significanceLevel):
    with open(fileName) as f:
        content = f.readlines()
    values = [int(x) for x in content]
    x = sum(values)
    pValue = binom_test(x, len(values), successProbability)
    assert(pValue > significanceLevel)

def cdfTruncatedExp(rate, x):
    if (rate > 0):
        return expm1(-x * rate) / expm1(-rate)
    else:
        return x

def testTruncatedExponentialDistribution(fileName, rate, significanceLevel):
    with open(fileName) as f:
        content = f.readlines()
        values = [float(x) for x in content]

        cdf = lambda x : cdfTruncatedExp(rate, x)

        _, pValue = kstest(values, cdf)
        assert(pValue > significanceLevel)


significanceLevel = 0.01

testUniformIntegerDistribution("data/uniformLemire3.txt", 3, significanceLevel)
testUniformIntegerDistribution("data/uniformLemire11.txt", 11, significanceLevel)
testUniformIntegerDistribution("data/uniformLemire29.txt", 29, significanceLevel)
testUniformIntegerDistribution("data/uniformLemire256.txt", 256, significanceLevel)
testUniformIntegerDistribution("data/uniformLumbroso3.txt", 3, significanceLevel)
testUniformIntegerDistribution("data/uniformLumbroso11.txt", 11, significanceLevel)
testUniformIntegerDistribution("data/uniformLumbroso29.txt", 29, significanceLevel)
testUniformIntegerDistribution("data/uniformLumbroso256.txt", 256, significanceLevel)
testUniformIntegerDistribution("data/intPow3.txt", pow(2, 3), significanceLevel)
testUniformIntegerDistribution("data/intPow8.txt", pow(2, 8), significanceLevel)

testUniformDoubleDistribution("data/uniformDouble.txt", significanceLevel)
testUniformDoubleHalfDistribution("data/uniformDoubleHalf.txt", significanceLevel)

testExponentialDistribution("data/expStandard.txt", significanceLevel)
testExponentialDistribution("data/expZiggurat.txt", significanceLevel)

testBernoulliDistribution("data/boolean.txt", 0.5, significanceLevel)
testBernoulliDistribution("data/bernoulliReal0_2.txt", 0.2, significanceLevel)
testBernoulliDistribution("data/bernoulliRatio1_3.txt", 1./3., significanceLevel)

testTruncatedExponentialDistribution("data/truncatedExp0.txt", 0, significanceLevel)
testTruncatedExponentialDistribution("data/truncatedExp0_1.txt", 0.1, significanceLevel)
testTruncatedExponentialDistribution("data/truncatedExp0_5.txt", 0.5, significanceLevel)
testTruncatedExponentialDistribution("data/truncatedExp1.txt", 1, significanceLevel)
testTruncatedExponentialDistribution("data/truncatedExp2.txt", 2, significanceLevel)
