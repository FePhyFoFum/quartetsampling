#!/usr/bin/env python3

"""
Convert fasta alignemnts to relaxed phylip ones in constant memory.
Written by Lucas Sinclair.
Kopimi.
This version edited by: James B. Pease

You can use this script from the shell like this::
$ fasta_to_phylip seqs.fasta seqs.phylip
"""

###############################################################################

import sys
import os
import random
import re


class Sequence(object):
    """The Sequence object has a string *header* and
    various representations."""

    def __init__(self, header, seq):
        self.header = re.findall('^>(\S+)', header)[0]
        self.seq = seq

    def __len__(self):
        return len(self.seq)

    @property
    def phylip(self):
        return self.header + " " + self.seq.replace('.', '-') + "\n"

    @property
    def fasta(self):
        return ">" + self.header + "\n" + self.seq + "\n"


def fasta_parse(path):
    """Reads the file at *path* and yields
       Sequence objects in a lazy fashion"""
    header = ''
    seq = ''
    with open(path) as f:
        for line in f:
            line = line.strip('\n')
            if line.startswith('>'):
                if header:
                    yield Sequence(header, seq)
                header = line
                seq = ''
                continue
            seq += line
    yield Sequence(header, seq)

###############################################################################


# The libraries we need #
# Get the shell arguments #
fa_path = sys.argv[1]
ph_path = sys.argv[2]
# Check that the path is valid #
if not os.path.exists(fa_path):
    raise RuntimeError("No file at: {}".format(fa_path))
# Use our two functions #
seqs = fasta_parse(fa_path)
# Write the output to temporary file #
tm_path = ph_path + '.' + str(random.randint(100000000, 1000000000))
# Count the sequences #
count = 0
with open(tm_path, 'w') as f:
    for seq in seqs:
        f.write(seq.phylip)
        count += 1
# Add number of entries and length at the top #
with open(tm_path, 'r') as old, open(ph_path, 'w') as new:
    new.write(" " + str(count) + " " + str(len(seq)) + "\n")
    new.writelines(old)
# Clean up #
os.remove(tm_path)
