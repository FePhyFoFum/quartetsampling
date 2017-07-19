#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# vim:fileencoding=utf-8
"""
Combines RESULT.node.scores.csv files from separate
runs for the same phylogeny into a single set of csv and tree outputs.

http://www.github.com/FePhyFoFum/quartetsampling
"""


import sys
import os
import argparse
from tree_data import TreeData
from rep_data import DataStore


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


def generate_argparser():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog=LICENSE)
    parser.add_argument('-d', '--nodedata', required=True, nargs=1,
                        help=("file containing paths of one or more"
                              "RESULT.node.score.csv files"))
    parser.add_argument('-t', '--tree', required=True, type=open,
                        nargs=1,
                        help="tree file in Newick format")
    parser.add_argument('-o', '--out', required=True,
                        nargs=1,
                        help="new output files prefix")
    parser.add_argument("-v", "--verbose", action="store_true")
    # These args are hidden to pass through to the treedata object
    parser.add_argument("-c", "--clade", nargs=1, help=argparse.SUPPRESS)
    parser.add_argument("-s", "--startk", type=int, default=0,
                        help=argparse.SUPPRESS)
    parser.add_argument("-p", "--stopk", type=int, help=argparse.SUPPRESS)
    return parser


def main(arguments=None):
    """Main method for merge_output"""
    arguments = arguments if arguments is not None else sys.argv[1:]
    parser = generate_argparser()
    args = parser.parse_args(args=arguments)
    scores = {}
    qfscores = {}
    qfreps = {}
    params = {
        'tree_result_file_path': "{}.tre".format(args.out[0]),
        'figtree_file_path': "{}.tre.figtree".format(args.out[0]),
        'freq_file_path': "{}.tre.freq".format(args.out[0]),
        'qc_tree_file_path': "{}.tre.qc".format(args.out[0]),
        'qd_tree_file_path': "{}.tre.qd".format(args.out[0]),
        'qi_tree_file_path': "{}.tre.qi".format(args.out[0]),
        'merged_file_path': "{}.nodes.scores.csv".format(args.out[0]),
        }
    for fpath in params.values():
        if os.path.exists(fpath):
            raise RuntimeError("File {} already exists! Exiting...".format(
                fpath))
    treedata = TreeData(args)
    params.update({
        'starttime': 0,
        'startk': args.startk,
        'stopk': (args.stopk if args.stopk is not None else
                  treedata.nleaves + 100),
        'verbose': args.verbose,
        'temp_wd': os.path.curdir,
        'results_dir': os.path.curdir,
        'suppress_output': True
        })
    mergedata = DataStore(params)
    filepaths = []
    with open(args.nodedata[0]) as nfile:
        for line in nfile:
            filepaths.append(line.strip())
    firstline = True
    for fname in filepaths:
        with open(fname, 'r') as infile:
            for line in infile:
                if firstline is True:
                    hdrs = line.rstrip().split(',')
                    qfindex = hdrs.index('qf')
                    repindex = hdrs.index('num_replicates')
                    firstline = False
                    continue
                arr = line.rstrip().split(',')
                if len(arr) < 2:
                    continue
                if arr[qfindex] != 'NA':
                    qfscores[arr[0]] = qfscores.get(arr[0], 0.0) + (
                        float(arr[qfindex]) * float(arr[repindex]))
                    qfreps[arr[0]] = qfreps.get(
                        arr[0], 0.0) + float(arr[repindex])
                else:
                    if arr[0] in scores:
                        print("WARNING: {} appeared in two"
                              "separate input files".format(
                                  arr[0]))
                    scores[arr[0]] = dict([
                        (hdrs[i], arr[i]) for i in range(1, len(arr))])
    for xlabel in qfscores:
        print(xlabel, qfscores[xlabel], qfreps[xlabel])
        qfscores[xlabel] = "{:.2g}".format(
            qfscores[xlabel] / qfreps[xlabel])
    k = 1
    mergedata.write_headers(params['merged_file_path'])
    for xnode in treedata.tree.iternodes():
        k, _ = treedata.check_node(xnode, k, params)
        xlabel = xnode.label
        if xlabel in scores:
            xnode.label = xlabel[:]
            xnode.data['freq0'] = scores[xlabel]['freq0']
            xnode.data['qc_score'] = scores[xlabel]['qc']
            xnode.data['qd_score'] = scores[xlabel]['qd']
            xnode.data['qi_score'] = scores[xlabel]['qi']
            xnode.data['replicates'] = scores[xlabel]['num_replicates']
            entry = {'node_label': xlabel,
                     'freq0': scores[xlabel]['freq0'],
                     'qc': scores[xlabel]['qc'],
                     'qd': scores[xlabel]['qd'],
                     'qi': scores[xlabel]['qi'],
                     'num_replicates': scores[xlabel]['num_replicates'],
                     'diff': scores[xlabel]['diff'],
                     'notes': scores[xlabel]['notes']}
            mergedata.write_entry(params['merged_file_path'], entry)
            entry = {}
    for xlabel in sorted(qfscores):
        entry = {'node_label': xlabel,
                 'qf': qfscores[xlabel],
                 'num_replicates': qfreps[xlabel]}
        mergedata.write_entry(params['merged_file_path'], entry)
        entry = {}
    treedata.write_scoretrees(params)
    treedata.write_figtree(params['figtree_file_path'], qfscores)
    return ''


if __name__ == '__main__':
    main()
