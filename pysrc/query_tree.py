#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# vim:fileencoding=utf-8
"""
Tree query script to find specific nodes numbers in large trees
when using the post-run annotated trees.

http://www.github.com/FePhyFoFum/quartetsampling
"""

import os
import sys
import argparse
from tree_data import TreeData
from phylo.tree_utils import get_mrca


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
    """Generates the argparsr ArgumentParser
    """
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog=LICENSE)
    parser.add_argument('-t', '--tree', type=open, nargs=1,
                        help="input tree in newick format")
    parser.add_argument('-d', '--data', type=os.path.abspath, nargs=1,
                        help=("CSV output from quartet_sampling"
                              " (RESULT.node.score.csv)"))
    parser.add_argument("-c", "--clade", nargs=1, help=argparse.SUPPRESS)
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="verbose screen output")
    parser.add_argument("-s", "--startk", type=int, default=0,
                        help=argparse.SUPPRESS)
    parser.add_argument("-p", "--stopk", type=int, help=argparse.SUPPRESS)
    return parser


def main(arguments=None):
    """Main method for query_tree.py
    """
    arguments = arguments if arguments is not None else sys.argv[1:]
    parser = generate_argparser()
    args = parser.parse_args(args=arguments)
    treedata = TreeData(args)
    params = {
        'starttime': 0,
        'startk': args.startk,
        'stopk': (args.stopk if args.stopk is not None else
                  treedata.nleaves + 100),
        'verbose': args.verbose,
        'temp_wd': os.path.curdir,
        'results_dir': os.path.curdir,
        'suppress_feedback': args.verbose is not True,
        }
    data = {}
    with open(args.data[0]) as datafile:
        firstline = True
        for line in datafile:
            if firstline is True:
                hdr = line.rstrip().split(',')[1:]
                firstline = False
            row = line.rstrip().split(',')
            data[row[0]] = row[1:]
    print("data read")
    # args.nodes = args.nodes[0].split(',')
    nodes = {}
    parent = {}
    k = 1
    ntotal = 0
    for xnode in treedata.tree.iternodes():
        k, _ = treedata.check_node(xnode, k, params)
        if xnode.label == '':
            xnode.label = str(k)
        xlabel = xnode.label
        nodes[xlabel] = xnode
        for xchild in xnode.children:
            parent[xchild] = xnode
        if xlabel not in data:
            continue
        xnode.data.update([(hdr[i], data[xlabel][i])for i in range(len(hdr))])
        if xnode.istip:
            ntotal += 1
    ret = ''
    print(ntotal, "total nodes in tree")
    print("Get MRCA data based on two nodes, separated by a comma")
    while ret not in ['quit', 'exit']:
        ret = input("Enter node names:")
        if ret in ['quit', 'exit']:
            sys.exit()
        try:
            ret = [x.strip() for x in ret.split(',')]
            if len(ret) == 1:
                mrca = nodes[ret[0]]
            else:
                mrca = get_mrca([nodes[x] for x in ret], treedata.tree)
            print("MRCA FOUND", mrca.label)
            print(mrca.data)

            def naround(num):
                """Round to two digits unless value is NA
                """
                return 'NA' if num == 'NA' else round(float(num), 2)

            print("SCORE={}/{}/{}".format(
                naround(mrca.data['qc']),
                naround(mrca.data['qd']),
                naround(mrca.data['qi']),
                ))
            print("TREE divided in to four groups:")

            def node_profile(xnode, grp):
                """Get profile of a specific node in the tree
                """
                ntip = 0
                tipname = ''
                qfscores = []
                for ynode in xnode.iternodes():
                    if ynode.istip:
                        if tipname == '':
                            tipname = ynode.label[:]
                        ntip += 1
                        if 'qf' in ynode.data:
                            qfscores.append(float(ynode.data['qf']))
                print("GROUP{} has {} tips "
                      "including {}, with meanQF={}".format(
                          grp, ntip, tipname, sum(qfscores)/len(qfscores)))
                return ntip
            ntip1 = node_profile(mrca.children[0], 0)
            ntip2 = node_profile(mrca.children[1], 1)
            mrca_parent = parent[mrca]
            if mrca_parent == treedata.tree:
                print("is root")
            else:
                for child in mrca_parent.children:
                    if child != mrca:
                        ntip3 = node_profile(child, 2)
            mrca_grandparent = parent[mrca_parent]
            if mrca_grandparent == treedata.tree:
                print("is root")
            else:
                ntip4 = ntotal - (ntip1 + ntip2 + ntip3)
                for child in mrca_grandparent.children:
                    if child != mrca_parent:
                        for xnode in child.iternodes():
                            if xnode.istip:
                                print("GROUP3 has {} tips including {}".format(
                                    ntip4, xnode.label))
                                break
        except Exception as exc:
            print(exc, sys.exc_info())
            print("error, try again or enter 'quit' to quit")
    return ''


if __name__ == '__main__':
    main()
