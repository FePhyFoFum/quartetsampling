#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# vim:fileencoding=utf-8
"""
http://www.github.com/FePhyFoFum/quartetsampling

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

import time
from phylo.tree_reader import read_tree_string
from phylo.tree_utils import get_mrca


class TreeData(object):
    """Main Tree Information"""

    def __init__(self, args):
        self.tree = None
        print("reading tree from {}".format(args.tree[0].name))
        self.tree = read_tree_string(args.tree[0].readline())
        if self.tree is None:
            raise RuntimeError(
                "Could not find a tree in the treefile: {}".format(
                    args.tree[0].name))
        args.tree[0].close()
        self.leaves = self.tree.leaves()
        self.nleaves = len(self.leaves)
        if args.verbose:
            print("tree has {} leaves".format(self.nleaves))
        self.rootlabelset = set()
        for leaf in self.leaves:
            self.rootlabelset.add(leaf.label)
        self.numnodes = (len(list(self.tree.iternodes())) -
                         len(list(self.tree.leaves())))
        self.clade = None
        if args.clade is not None:
            names = args.clade.split(",")
            nodes = []
            print("Just calculating for clade {}".format(names))
            for i in names:
                found = False
                for j in self.leaves:
                    if j.label == i:
                        nodes.append(j)
                        found = True
                        break
                if found is False:
                    raise RuntimeError("{} not found. Exiting...".format(i))
            self.clade = get_mrca(nodes, self.tree)

    def __str__(self):
        print("tree:", self.tree.get_newick_repr())
        print("leaves", [x.label for x in self.leaves])
        print("nleaves:", self.nleaves)
        print("clade:", self.clade)
        print("numnodes:", self.numnodes)
        return ''

    def check_node(self, fnode, k, params):
        """Run a series of checks on the node for suitability"""
        root_bipart_label = None
        if fnode.istip or fnode.parent is None:
            if params['verbose']:
                print("\nskipping {}".format(fnode.label is not None and
                                             "root" or
                                             "tip {}".format(fnode.label)))
            return k, False
        #  record the node label in the tree, these are required for user to
        #  match scores with corresponding branches
        fnode.label = "QS{}".format(k)
        if k < params["startk"]:
            #  skip nodes if we have a specified start node (i.e. not the root)
            #  and we haven't hit it yet
            k += 1
            return k, False
        else:
            # provide user feedback before incrementing
            time_string = (user_feedback_time(k, root_bipart_label, params,
                                              self.numnodes)
                           if k > params["startk"] else "")
            if not params.get("suppress_feedback", False) is True:
                print("\nprocessing node {}/{}{}".format(
                    k, self.numnodes, time_string))
        #  require a bifurcating tree
        #  assert(len(fnode.children) == 2)
        if len(fnode.children) != 2:
            print("Node {} does not have exactly 2 children."
                  "It will be skipped.".format(k))
            return k, False
        # get leaf sets for the four connected subtrees
        (leafsets, testsib, root_bipart_label, is_other_side_of_root,
         skip_tip_child_of_root, tip_child_label) = self.set_leaf_sets(fnode)
        # there is some case that the sib has no other
        if testsib is False:
            print("other sib not found")
            return k, False
        # no more user feedback, now we can increment k
        k += 1
        if skip_tip_child_of_root:
            print("not calculating qc for tip child '" + tip_child_label +
                  "' of the root (qc is 1.0, as for all tips).")
            return k, False
        # if we already processed bipart at root and this is other side of that
        if is_other_side_of_root:
            print("\nskipping second instance of root-adjacent bipartition" +
                  "(it was already processed at node " +
                  root_bipart_label + ").")
            fnode.label = root_bipart_label
            return k, False
        # sanity check
        allset = set()
        for leafset in leafsets.values():
            assert len(leafset) > 0
            allset.update(leafset)
        assert len(allset) == self.nleaves
        del allset
        return k, leafsets

    def write_scoretrees(self, params):
        """Various trees annotated with scores"""
        # LEAVING THESE FOR SIMS, BUT REMOVE AFTER THAT AND JUST USE THE
        # MULTILABELED ONE, FREQ IS GOING TO BE QS_FREQ UNTIL AFTER SIMS
        # write the tree with all processed nodes labeled
        with open(params['tree_result_file_path'], "w") as tree_file_path:
            tree_file_path.write(self.tree.get_newick_repr(True) + ";")
        # write a tree with the frequencies at the nodes
        for datakey, datafile in (
                ('freq0', 'freq_file_path'),
                ('qc_score', 'qc_tree_file_path'),
                ('qd_score', 'qd_tree_file_path'),
                ('qi_score', 'qi_tree_file_path')):
            for i in self.tree.iternodes():
                if (len(i.children) > 1 and
                        i is not self.tree and datakey in i.data):
                    i.label = "{}={}".format(datakey.replace('_score', ''),
                                             i.data.get(datakey, ''))
                elif len(i.children) > 1:
                    i.label = "{}={}".format(datakey.replace('_score', ''),
                                             i.data.get(datakey, ''))
            with open(params[datafile], "w") as freq_file_path:
                freq_file_path.write(self.tree.get_newick_repr(True) + ";")
        return ''

    def write_figtree(self, outfname, qfscores):
        """Write output in FigTree format with scores"""
        # change the names (internal only) to be the full figtree style label
        for xnode in self.tree.iternodes():
            xnode.data["lab"] = xnode.label  # so we can use it later
            if "freq0" not in xnode.data:
                if xnode.label in qfscores:
                    lab = ("{}[&qf={}]".format(
                        xnode.label,
                        qfscores[xnode.label]))
                    xnode.label = lab
            #    continue
            elif len(xnode.children) == 0:
                if xnode.label in qfscores:
                    lab = ("{}[&qf={}]".format(
                        xnode.label,
                        qfscores[xnode.label]))
                    xnode.label = lab
            elif len(xnode.data["freq0"]) > 0:
                lab = ("[&label={},freq={},qc={},qd={},qi={},"
                       "reps={},score={}/{}/{}]").format(
                           xnode.label,
                           xnode.data["freq0"],
                           xnode.data["qc_score"],
                           xnode.data["qd_score"],
                           xnode.data["qi_score"],
                           xnode.data["replicates"],
                           xnode.data["qc_score"],
                           xnode.data["qd_score"],
                           xnode.data["qi_score"]
                           )
                xnode.label = lab
            else:
                lab = ("[&label={},reps={}]").format(
                    xnode.label,
                    xnode.data["replicates"])
                xnode.label = lab
        outf = open(outfname, "w")
        outf.write('#NEXUS\nbegin trees;\n\n')
        outf.write('\ttree tree1 = [&R] {};\n\n'.format(
            self.tree.get_newick_repr(True)))
        outf.write('end;\n\n\n\nbegin figtree;\n\n'
                   '\tset appearance.backgroundColour=#-1;\n'
                   '\tset appearance.branchColorAttribute="freq";\n'
                   '\tset appearance.branchLineWidth=1.0;\n'
                   '\tset appearance.foregroundColour=#-16777216;\n'
                   '\tset appearance.selectionColour=#-2144520576;\n'
                   '\tset branchLabels.displayAttribute="Branch times";\n'
                   '\tset branchLabels.fontName="sansserif";\n'
                   '\tset branchLabels.fontSize=8;\n'
                   '\tset branchLabels.fontStyle=0;\n'
                   '\tset branchLabels.isShown=false;\n'
                   '\tset branchLabels.significantDigits=4;\n'
                   '\tset layout.expansion=0;\n'
                   '\tset layout.layoutType="RECTILINEAR";\n'
                   '\tset layout.zoom=0;\n'
                   '\tset nodeBars.barWidth=4.0;\n'
                   '\tset nodeLabels.displayAttribute="label";\n'
                   '\tset nodeLabels.fontName="sansserif";\n'
                   '\tset nodeLabels.fontSize=8;\n'
                   '\tset nodeLabels.fontStyle=0;\n'
                   '\tset nodeLabels.isShown=false;\n'
                   '\tset nodeLabels.significantDigits=4;\n'
                   '\tset polarLayout.alignTipLabels=false;\n'
                   '\tset polarLayout.angularRange=0;\n'
                   '\tset polarLayout.rootAngle=0;\n'
                   '\tset polarLayout.rootLength=100;\n'
                   '\tset polarLayout.showRoot=true;\n'
                   '\tset radialLayout.spread=0.0;\n'
                   '\tset rectilinearLayout.alignTipLabels=false;\n'
                   '\tset rectilinearLayout.curvature=0;\n'
                   '\tset rectilinearLayout.rootLength=100;\n'
                   '\tset scale.offsetAge=0.0;\n'
                   '\tset scale.rootAge=1.0;\n'
                   '\tset scale.scaleFactor=1.0;\n'
                   '\tset scale.scaleRoot=false;\n'
                   '\tset scaleAxis.automaticScale=true;\n'
                   '\tset scaleAxis.fontSize=8.0;\n'
                   '\tset scaleAxis.isShown=false;\n'
                   '\tset scaleAxis.lineWidth=1.0;\n'
                   '\tset scaleAxis.majorTicks=1.0;\n'
                   '\tset scaleAxis.origin=0.0;\n'
                   '\tset scaleAxis.reverseAxis=false;\n'
                   '\tset scaleAxis.showGrid=true;\n'
                   '\tset scaleAxis.significantDigits=4;\n'
                   '\tset scaleBar.automaticScale=true;\n'
                   '\tset scaleBar.fontSize=10.0;\n'
                   '\tset scaleBar.isShown=true;\n'
                   '\tset scaleBar.lineWidth=1.0;\n'
                   '\tset scaleBar.scaleRange=0.0;\n'
                   '\tset scaleBar.significantDigits=4;\n'
                   '\tset tipLabels.displayAttribute="Names";\n'
                   '\tset tipLabels.fontName="sansserif";\n'
                   '\tset tipLabels.fontSize=8;\n'
                   '\tset tipLabels.fontStyle=0;\n'
                   '\tset tipLabels.isShown=true;\n'
                   '\tset tipLabels.significantDigits=4;\n'
                   '\tset trees.order=true;\n'
                   '\tset trees.orderType="decreasing";\n'
                   '\tset trees.rooting=false;\n'
                   '\tset trees.rootingType="User Selection";\n'
                   '\tset trees.transform=false;\n'
                   '\tset trees.transformType="cladogram";\nend;\n\n')
        outf.close()
        for xnode in self.tree.iternodes():
            if "freq0" not in xnode.data:
                continue
            if len(xnode.children) == 0:
                continue
            xnode.label = xnode.data["lab"]
        return ''

    def set_leaf_sets(self, fnode, root_bipart_label=None):
        """Get sets of leaves"""
        leafsets = {}
        # two daughter subtrees
        for i, lbl in enumerate(["R1", "R2"]):
            leafsets[lbl] = set([fnode.children[i].label, ]
                                if fnode.istip
                                else [l.label for l in
                                      fnode.children[i].leaves()])
        # sibling/parent subtrees
        is_other_side_of_root = False  # used when we hit root a second time
        skip_tip_child_of_root = False  # used when a child of root is a tip
        tip_child_label = None
        testsib = False
        for sib in fnode.parent.children:
            if sib != fnode:
                testsib = True
                # if one of the subtrees is the root, skip over it
                if len(sib.leaves()) + len(fnode.leaves()) == self.nleaves:
                    # if we already processed this bipart (on other side of
                    # the root), don't do it again
                    if root_bipart_label is not None:
                        is_other_side_of_root = True
                        break
                    # get the subtrees opposite the root
                    if len(sib.children) == 2:
                        for i, lbl in enumerate(["L1", "L2"]):
                            leafsets[lbl] = set(
                                [sib.children[i].label, ]
                                if sib.children[i].istip
                                else [l.label for l in
                                      sib.children[i].leaves()])
                    elif len(sib.children) == 0:
                        skip_tip_child_of_root = True
                        tip_child_label = sib.label
                    else:
                        print(("Node {} does not have exactly 2 children."
                               "It will be skipped.").format(fnode.label))
                        continue
                    # remember that we've already done the root,
                    # so we can skip it when we hit the other side
                    root_bipart_label = fnode.label
                # otherwise not at root, all connected subtrees have children
                else:
                    # sibling subtree
                    leafsets["L1"] = set([l.label for l in sib.leaves()])
                    # the rest of the tree
                    leafsets["L2"] = self.rootlabelset - \
                        (leafsets["R1"].union(leafsets["R2"]).union(
                            leafsets["L1"]))  # set()
        return leafsets, testsib, root_bipart_label, is_other_side_of_root, \
            skip_tip_child_of_root, tip_child_label


def user_feedback_time(k, root_bipart_label, params, numnodes):
    """User feedback of time left"""
    mean_time_secs = ((time.time() - params['starttime']) /
                      float(k - params["startk"]))
    if mean_time_secs > 60:
        if mean_time_secs > 3600:  # more than one hour
            mean_time_units = "hours"
            mean_time_scalar = 3600
        else:  # between 1 and 60 minutes
            mean_time_units = "minutes"
            mean_time_scalar = 60
    else:  # less than 60 seconds
        mean_time_units = "seconds"
        mean_time_scalar = 1
    # adjust for the duplicate bipart at the root
    # (until we hit it, then stop adjusting)
    adj = -1 if root_bipart_label is None else 0
    est_remaining_time_secs = (mean_time_secs * (numnodes - k + adj))
    if est_remaining_time_secs > 60:
        if est_remaining_time_secs > 3600:
            est_remaining_time_units = "hours"
            est_remaining_time_scalar = 3600
        else:  # between 1 and 60 minutes
            est_remaining_time_units = "minutes"
            est_remaining_time_scalar = 60
    else:  # less than 60 seconds
        est_remaining_time_units = "seconds"
        est_remaining_time_scalar = 1
    time_string = " | average node time {0:.2f} {1:s}".format(
        mean_time_secs / mean_time_scalar, mean_time_units) + \
        " | est. remaining time {0:.2f} {1:s}".format(
            est_remaining_time_secs / est_remaining_time_scalar,
            est_remaining_time_units)
    return time_string


def write_test_trees(temp_wd="."):
    """Writing testing trees"""
    topos = ("(L1,L2,(R1,R2));\n",
             "(L1,R1,(L2,R2));\n",
             "(L1,R2,(L2,R1));\n")
    for i, topo in enumerate(topos):
        with open("{}/test.trees.{}".format(temp_wd, i), "w") as xtrees:
            xtrees.write(topo)
    return ''


if __name__ == "__main__":
    print("This file is a function library, please run quartet_sampling.py")
