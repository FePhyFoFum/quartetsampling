#!/usr/bin/env python3

import os
import sys
import subprocess


def run(test_name):
    try:
        os.chdir(test_name)
        print('\ncurrently executing test at: ' + os.getcwd() + '\n')
    except OSError:
        raise OSError(
            "the specified test '{}' does not seem to exist.".format(
                test_name))
    args = ['python', '../../pysrc/quartet_sampling.py',
            'test.config']
    print('running command: \n'+' '.join(args)+'\n')
    subprocess.call(' '.join(args), shell=True)


if __name__ == '__main__':
    run(sys.argv[1])
