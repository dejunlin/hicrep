import os
import numpy as np
import math
import sys
import warnings
from hicrep.utils import (
    readMcool, cool2pixels, getSubCoo,
    trimDiags, meanFilterSparse, varVstran,
    resample
    )
from hicrep.hicrep import (
    sccOfDiag, hicrepSCC
    )

def main(*args):
    import argparse
    import subprocess
    import re

    np.random.seed(10)

    parser = argparse.ArgumentParser()
    parser.add_argument("fmcool1", type=str,
                        help="First cooler multiple-binsize contact files")
    parser.add_argument("fmcool2", type=str,
                        help="Second cooler multiple-binsize contact files")
    parser.add_argument("fout", type=str,
                        help="Output results to this file. Output format would be\
                        one column of scc scores for each chromosome")
    parser.add_argument("--binSize", type=int, default=-1,
                        help="Use this to select the bin size from the input mcool\
                        file. Default to -1, meaning that the inputs are treated as\
                        single-binsize .cool files")
    parser.add_argument("--h", type=int, default=0,
                        help="Smooth the input contact matrices using a 2d mean\
                        filter with window size of 1 + 2 * value")
    parser.add_argument("--dBPMax", type=int, default=-1,
                        help="Only consider contacts at most this number of bp away\
                        from the diagonal. Default to -1, meaning the entire\
                        contact matrix is used")
    parser.add_argument("--bDownSample", action='store_true', default=False,
                        help="Down sample the input with more contact counts to\
                        the the same number of counts as the other input with less\
                        contact counts. If turned off, the input matrices will be\
                        normalized by dividing the counts by their respective total\
                        number of contacts.")

    args = parser.parse_args()

    header = "#"+" ".join(sys.argv)+"\n"

    # Check if current script is under revision control
    gitls = subprocess.Popen('cd '+os.path.dirname(os.path.realpath(__file__)) +
                             ' && git ls-files --error-unmatch ' + __file__,
                             shell=True, stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE, encoding='utf-8')
    if not gitls.stderr.read():
        gitrev = subprocess.Popen('cd '+os.path.dirname(os.path.realpath(__file__))
                                  + ' && git rev-parse HEAD --abbrev-ref HEAD',
                                  shell=True, stdout=subprocess.PIPE,
                                  stderr=subprocess.PIPE, encoding='utf-8')
        if not gitrev.stderr.read():
            p = re.compile("\n")
            header += "# @rev " + p.sub("\n# @branch ",
                                        gitrev.stdout.read().strip()) + "\n"

    fmcool1 = args.fmcool1
    fmcool2 = args.fmcool2
    fout = args.fout
    binSize = args.binSize
    h = args.h
    dBPMax = args.dBPMax
    bDownSample = args.bDownSample

    cool1, binSize1 = readMcool(fmcool1, binSize)
    cool2, binSize2 = readMcool(fmcool2, binSize)

    scc = hicrepSCC(cool1, cool2, h, dBPMax, bDownSample)

    np.savetxt(fout, scc, "%30.15e", header=header)
