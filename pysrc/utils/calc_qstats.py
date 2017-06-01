#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# vim:fileencoding=utf-8
#
# This script is meant as a general template for Py3
#
# @author: James B. Pease
# @version: 1.0

import os
import sys
import argparse
import numpy as np


def basic_stats(nums):
    print("Mean (IQR)")
    print("{} ({},{})".format(
        round(np.mean(nums), 2),
        round(np.percentile(nums, 25), 2),
        round(np.percentile(nums, 75), 2)))
    return ''


def main(arguments=None):
    arguments = arguments if arguments is not None else sys.argv[1:]
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-d', '--data', type=os.path.abspath, nargs=1)
    parser.add_argument("-c", "--clade", nargs=1)
    parser.add_argument("-v", "--verbose", action="store_true")
    parser.add_argument("-s", "--startk", type=int, default=0)
    parser.add_argument("-p", "--stopk", type=int)
    parser.add_argument("-o", "--out", type=os.path.abspath, nargs=1)
    args = parser.parse_args(args=arguments)
    data = {}
    with open(args.data[0]) as datafile:
        firstline = True
        for line in datafile:
            if firstline is True:
                hdr = line.rstrip().split(',')[1:]
                firstline = False
                continue
            row = line.rstrip().split(',')
            data[row[0]] = row[1:]
    print("data read")
    data_qc = []
    data_qd = []
    data_qu = []
    data_qf = []
    nterm = 0
    nint = 0
    qc_index = hdr.index('qc')
    qd_index = hdr.index('qd')
    qu_index = hdr.index('qu')
    qf_index = hdr.index('qf')
    for entry in data:
        if data[entry][qc_index] != 'NA':
            data_qc.append(float(data[entry][qc_index]))
        if data[entry][qd_index] != 'NA':
            data_qd.append(float(data[entry][qd_index]))
        if data[entry][qu_index] != 'NA':
            data_qu.append(float(data[entry][qu_index]))
        if data[entry][qf_index] != 'NA':
            nterm += 1
            data_qf.append(float(data[entry][qf_index]))
        else:
            nint += 1
    print(nterm, nint)
    print("QC")
    basic_stats(data_qc)
    print("QD")
    basic_stats(data_qd)
    print("QU")
    basic_stats(data_qu)
    print("QF")
    basic_stats(data_qf)
    return ''


if __name__ == '__main__':
    main()
