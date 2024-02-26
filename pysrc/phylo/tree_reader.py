#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# vim:fileencoding=utf-8
"""
this takes a newick string as instr
and reads the string and makes the
nodes and returns the root node

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


from .node import Node


def read_tree_string(instr):
    root = None
    index = 0
    nextchar = instr[index]
    start = True
    keepgoing = True
    curnode = None
    while keepgoing is True:
        if nextchar == "(":
            if start is True:
                root = Node()
                curnode = root
                start = False
            else:
                newnode = Node()
                curnode.add_child(newnode)
                curnode = newnode
        elif nextchar == ',':
            curnode = curnode.parent
        elif nextchar == ")":
            curnode = curnode.parent
            index += 1
            nextchar = instr[index]
            name = ""
            while True:
                if nextchar in ',):;[':
                    break
                name += nextchar
                index += 1
                nextchar = instr[index]
            curnode.label = name
            index -= 1
        elif nextchar == ';':
            keepgoing = False
            break
        elif nextchar == ":":
            index += 1
            nextchar = instr[index]
            brlen = ""
            while True:
                if nextchar in ',):;[':
                    break
                brlen += nextchar
                index += 1
                nextchar = instr[index]
            curnode.length = float(brlen)
            index -= 1
        elif nextchar == ' ':
            index += 1
            nextchar = instr[index]
        else:  # this is an external named node
            newnode = Node()
            curnode.add_child(newnode)
            curnode = newnode
            curnode.istip = True
            name = ""
            while True:
                if nextchar in ',):;[':
                    break
                name += nextchar
                index += 1
                nextchar = instr[index]
            curnode.label = name
            index -= 1
        if index < len(instr) - 1:
            index += 1
        nextchar = instr[index]
    return root


if __name__ == "__main__":
    treestring = "(a:3,(b:1e-05,c:1.3)int_|_and_33.5:5)root;"
    xnode = read_tree_string(treestring)
    print(xnode.get_newick_repr(True))
