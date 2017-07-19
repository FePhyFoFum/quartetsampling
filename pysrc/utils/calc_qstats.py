#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# vim:fileencoding=utf-8
"""Calculate basic statistics on the
   RESULTS.node.score.csv output file
   from quartet_sampling
   """


import os
import sys
import argparse
import numpy as np


LICENSE = """
This file is part of 'quartetsampling'.

'quartetsampling' is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

'quartetsampling' is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with 'quartetsampling'.  If not, see <http://www.gnu.org/licenses/>.
"""


def basic_stats(nums):
    print("Mean (IQR)")
    print("{} ({},{})".format(
        round(np.mean(nums), 2),
        round(np.percentile(nums, 25), 2),
        round(np.percentile(nums, 75), 2)))
    return ''


def generate_argparser():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog=LICENSE)
    parser.add_argument('-d', '--data', type=os.path.abspath, nargs=1,
                        required=True,
                        help=("RESULT.node.score.csv file output from"
                              "quartet_sampling.py"))
    parser.add_argument("-c", "--clade", nargs=1,
                        help=("specify a clade using a comma-separated"
                              "list of 2+ descendant taxa"))
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="verbose screen output")
    parser.add_argument("-s", "--startk", type=int, default=0,
                        help="starting branch numerical index")
    parser.add_argument("-p", "--stopk", type=int,
                        help="stopping branch numerical index")
    parser.add_argument("-o", "--out", type=os.path.abspath, nargs=1,
                        help="output file path for statistics")
    return parser


def main(arguments=None):
    arguments = arguments if arguments is not None else sys.argv[1:]
    parser = generate_argparser()
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
    data_qi = []
    data_qf = []
    nterm = 0
    nint = 0
    qc_index = hdr.index('qc')
    qd_index = hdr.index('qd')
    qi_index = hdr.index('qi')
    qf_index = hdr.index('qf')
    for entry in data:
        if data[entry][qc_index] != 'NA':
            data_qc.append(float(data[entry][qc_index]))
        if data[entry][qd_index] != 'NA':
            data_qd.append(float(data[entry][qd_index]))
        if data[entry][qi_index] != 'NA':
            data_qi.append(float(data[entry][qi_index]))
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
    print("QI")
    basic_stats(data_qi)
    print("QF")
    basic_stats(data_qf)
    return ''


if __name__ == '__main__':
    main()
