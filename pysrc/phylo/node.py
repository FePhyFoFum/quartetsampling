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


class Node:
    """Node phylogeny"""
    def __init__(self):
        self.label = ""
        self.length = 0.0
        self.parent = None
        self.children = []
        self.data = {}
        self.istip = False

    def add_child(self, child):
        # make sure that the child is not already in there
        assert child not in self.children
        self.children.append(child)
        child.parent = self

    def remove_child(self, child):
        # make sure that the child is in there
        assert child in self.children
        self.children.remove(child)
        child.parent = None

    def leaves(self, v=None):
        if v is None:
            v = []
        if len(self.children) == 0:
            v.append(self)
        else:
            for child in self.children:
                child.leaves(v)
        return v

    def leaves_fancy(self):
        return [n for n in self.iternodes() if n.istip]

    def lvsnms(self):
        return [n.label for n in self.iternodes() if n.istip]

    def iternodes(self, order="preorder"):
        if order.lower() == "preorder":
            yield self
        for child in self.children:
            for d in child.iternodes(order):
                yield d
        if order.lower() == "postorder":
            yield self

    def prune(self):
        p = self.parent
        if p is not None:
            p.remove_child(self)
        return p

    def branch_lengths(self):
        return [n.length for n in self.iternodes()]

    def get_newick_repr(self, showbl=False):
        ret = ""
        for i in range(len(self.children)):
            if i == 0:
                ret += "("
            ret += self.children[i].get_newick_repr(showbl)
            if i == len(self.children)-1:
                ret += ")"
            else:
                ret += ","
        if self.label is not None:
            ret += self.label
        if showbl is True:
            ret += ":" + str(self.length)
        return ret

    def _calc_depth(self):
        '''recursively calculate the depth of this node'''
        if len(self.children) < 1:
            return self.length
        else:
            return self.length + max([c.depth for c in self.children])

    @property
    def depth(self):
        '''return the depth of this node in the tree'''
        # currently just calculating the depth for every call.
        # should really store this and just update it when necessary
        return self._calc_depth()


if __name__ == "__main__":
    print("This file is a function library, please run quartet_sampling.py")
