#!/usr/bin/env python

import argparse
import glob
import os
import sys
import unittest

from tests.test_f5class import TestFast5Class


parser = argparse.ArgumentParser(description = """Given path(s) to fast5
                                 file(s) and/or directories of fast5s, 
                                 run all tests in fast5tools.tests.test_f5class.
                                 """, formatter_class = argparse.RawTextHelpFormatter)


parser.add_argument("fast5", metavar="fast5", nargs="+", type=str, 
                    help="""Paths to as many fast5 files and/or directories 
                    filled with fast5 files as you want.""")


if __name__ == "__main__":

    args = parser.parse_args()
    f5list = args.fast5

    suite = unittest.TestSuite()
    fast5class_testnames = unittest.TestLoader().getTestCaseNames(TestFast5Class)

    for pathname in f5list:
        if os.path.isdir(pathname):
            for fn in glob.glob(os.path.join(pathname, "*.fast5")):
                for test in fast5class_testnames:
                    suite.addTest(TestFast5Class(test, filename=fn))
        else:
            for test in fast5class_testnames:
                suite.addTest(TestFast5Class(test, filename=pathname))
    unittest.TextTestRunner().run(suite)
