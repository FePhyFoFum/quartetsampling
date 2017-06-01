#!/usr/bin/env python3

import os
import subprocess
for name in os.listdir('.'):
    if name[0:4] == 'test':
        subprocess.call('python run_single_test.py ' + name, shell=True)
