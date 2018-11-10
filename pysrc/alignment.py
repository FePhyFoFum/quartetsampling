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

import os
from random import sample as rsample


class Alignment(object):
    """Alignment object containing all sequences.
       Read the alignment into a dict,
       assumes phylip format with seqs unbroken on lines
    """

    def __init__(self, params):
        self.length = 0
        self.seqs = {}
        self.partitions = {}
        self.valid_sites = {}
        self.min_overlap = params['min_overlap']
        self.valid_chars = 'AaCcGgTtUu'
        self.invalid_chars = "NnXx-?"
        if params['data_type'] == 'amino':
            self.valid_chars = 'AaCcDdEeFfGgHhIiKkLlMmNnPpQqRrSsTtVvWwYy'
            self.invalid_chars = "Xx-?"
        elif params['data_type'] == 'cat':
            self.valid_chars = (
                '0123456789'
                'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
                'abcdefghijklmnopqrstuvwxyz')
            self.invalid_chars = "-?"
        self.genes = []

    def read_align(self, alnfile, params):
        """Read the full alingment"""
        print("reading alignment from {}".format(alnfile.name))
        firstline = True
        firstentry = True
        self.min_overlap = params['min_overlap']
        for line in alnfile:
            if firstline:
                firstline = False
                continue
            entry = line.split()
            if len(entry) > 1:
                if entry[0] in self.seqs:
                    raise RuntimeError(
                        "Sequence label {} is duplicate".format(
                            entry[0]))
                if params['low_mem'] is True:
                    sfilepath = os.path.join(
                        params['temp_wd'],
                        "seqtemp.{}.phy".format(entry[0]))
                    self.seqs[entry[0]] = sfilepath[:]
                    with open(sfilepath, 'w') as sfile:
                        sfile.write(entry[1])
                else:
                    self.seqs[entry[0]] = entry[1]
                if firstentry:
                    self.length = len(entry[1])
                    firstentry = False
                else:
                    if len(entry[1]) != self.length:
                        raise RuntimeError(
                            "Sequence '{}' is not the same length({})"
                            "as the first sequence ({})".format(
                                entry[0], len(entry[1]), self.length))
                    validchars, invalidchars = (
                        self.count_valid_chars(entry[1]))
                    if params['verbose']:
                        print(entry[0], "has ",
                              validchars, "valid sites and",
                              invalidchars, "invalid sites")
                    if validchars < self.min_overlap:
                        print("WARNING: Sequence {} has {}"
                              "valid sequence characters,"
                              "which is less than the minimum overlap"
                              "all quartets including this taxon"
                              "will be rejected!".format(entry[0],
                                                         validchars))
        alnfile.close()
        return ''

    def read_genes(self, alnfile, params):
        """Read the gene-partitioned alignment"""
        print("reading gene tree boundaries from {}".format(
            params['genetrees_file_path']))
        firstline = True
        firstentry = True
        self.min_overlap = params['min_overlap']
        with open(params['genetrees_file_path']) as gfile:
            for line in gfile:
                entry = [x.strip() for x in line.rstrip().split(',')]
                gname = entry[1].split("=")[0].strip()
                coords = entry[1].split("=")[1].strip().split('-')
                self.partitions[gname] = (int(coords[0]), int(coords[1]))
                self.seqs[gname] = {}
                self.valid_sites[gname] = {}
                self.genes.append(gname)
        print("reading alignment from {}".format(alnfile.name))
        for line in alnfile:
            if firstline:
                firstline = False
                continue
            entry = line.split()
            if len(entry) > 1:
                if firstentry is True:
                    self.length = len(entry[1])
                    firstentry = False
                elif len(entry[1]) != self.length:
                    raise RuntimeError(
                        "Sequence '{}' is not the same length({})"
                        "as the first sequence ({})".format(
                            entry[0], len(entry[1]), self.length))
                validchars, invalidchars = (
                    self.count_valid_chars(entry[1]))
                if params['verbose']:
                    print(entry[0], "has ",
                          validchars, "valid sites and",
                          invalidchars, "invalid sites")
                if validchars < self.min_overlap:
                    print("WARNING: Sequence {} has {}"
                          "valid sequence characters,"
                          "which is less than the minimum overlap"
                          "all quartets including this taxon"
                          "will be rejected!".format(entry[0],
                                                     validchars))
                for gene in self.partitions:
                    self.seqs[gene][entry[0]] = (
                        entry[1][self.partitions[gene][0] - 1:
                                 self.partitions[gene][1]])
        alnfile.close()
        return ''

    def __str__(self):
        print("seqs:", [x for x in self.seqs])
        print("length", self.length)
        return ''

    def count_valid_chars(self, sequence):
        """for getting the overlap between test seqs
        """
        valid_chars = 0
        invalid_chars = 0
        for xchar in sequence:
            if xchar in self.valid_chars:
                valid_chars += 1
            else:
                invalid_chars += 1
        return valid_chars, invalid_chars

    def check_aln_overlap(self, replicate):
        """Check for alignment overlap
        """
        n_overlap = 0
        seqnames = replicate[0:4]
        min_overlap = 1 if self.min_overlap < 1 else self.min_overlap
        genename = replicate[4]
        if genename is None:
            for i in range(len(self.seqs[seqnames[0]])):
                if self.seqs[seqnames[0]][i] not in self.valid_chars:
                    continue
                is_valid_site = True
                for seqname in seqnames[1:]:
                    if self.seqs[seqname][i] not in self.valid_chars:
                        is_valid_site = False
                        break
                if is_valid_site:
                    n_overlap += 1
                if n_overlap >= min_overlap:
                    break
        else:
            for i in range(len(self.seqs[genename][seqnames[0]])):
                if (self.seqs[genename][seqnames[0]][i] not in
                        self.valid_chars):
                    continue
                is_valid_site = True
                for seqname in seqnames[1:]:
                    if (self.seqs[genename][seqname][i] not in
                            self.valid_chars):
                        is_valid_site = False
                        break
                if is_valid_site:
                    n_overlap += 1
                if n_overlap >= min_overlap:
                    break
        return n_overlap >= min_overlap

    def get_random_gene(self):
        """get a random gene from the set of genes in seq
        """
        return rsample(self.genes, 1)[0]


if __name__ == "__main__":
    print("This file is a function library, please run quartet_sampling.py")
