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


class ParamSet(dict):
    """Parameter container dictionary"""

    def setup(self, args, nleaves):
        """Setup main parameters"""
        # defaults
        default_mrsp = 0.5
        # Basic params
        self['verbose'] = args.verbose
        self['low_mem'] = args.low_mem
        self['retain_temp'] = args.retain_temp
        self['data_type'] = args.data_type[0]
        self['raxml_model'] = 'GTRGAMMA'
        if self['data_type'] == 'amino':
            self['raxml_model'] = 'PROTGAMMAWAG'
        elif self['data_type'] == 'cat':
            self['raxml_model'] = 'BINGAMMA'
        self['calc_qdstats'] = args.calc_qdstats
        self['max_quartet_enumeration_threshold'] = 0.333333
        self['just_clade'] = args.clade is not None
        self['nprocs'] = args.number_of_threads[0]
        if args.partitions is not None and args.genetrees is not None:
            raise RuntimeError("Cannot use -g (--genetrees)"
                               "and -q (--partitions) simultaneously")
        # Partitions and Gene trees
        self['partitions_file_path'] = (os.path.abspath(args.partitions[0])
                                        if args.partitions is not None
                                        else None)
        self['using_partitions'] = self['partitions_file_path'] is not None
        self['genetrees_file_path'] = (os.path.abspath(args.genetrees[0])
                                       if args.genetrees is not None
                                       else None)
        self['using_genetrees'] = self['genetrees_file_path'] is not None
        # RAxML and PAUP Settings
        self['ignore_error'] = args.ignore_errors
        self['raxml_executable'] = (args.raxml_executable[0] if
                                    args.raxml_executable is not None else
                                    'raxml')
        self['paup'] = bool(args.paup)
        self['paup_executable'] = (args.paup_executable[0] if
                                   args.paup_executable is not None else
                                   'paup')
        if self['using_genetrees'] is True and self['low_mem']:
            raise RuntimeError("Cannot use -g (--genetrees) and"
                               "--low-mem simultaneously")
        # This parameter sets a threshold for how close you need to be to the
        # total number of possible quartets to trigger exhaustive
        # sampling

        print("setting the number of threads to {}".format(self['nprocs']))
        self['nreps'] = args.number_of_reps[0]
        print("setting the number of replicates to {}".format(self['nreps']))
        self['min_overlap'] = (args.min_overlap if
                               args.min_overlap is not None else 1)
        print("setting the minimum overlap to {}".format(self['min_overlap']))
        self['startk'] = (args.start_node_number[0]
                          if args.start_node_number is not None
                          else 1)
        self['stopk'] = (args.stop_node_number[0]
                         if args.stop_node_number is not None
                         else nleaves + 100)
        if self['stopk'] < self['startk']:
            raise RuntimeError("The start node number is higher"
                               "than the stop node number, "
                               "designating no nodes for processing.")
        # File Paths
        # print(args.temp_dir)
        self['temp_wd'] = (os.path.abspath("./QuartetSampling") if
                           args.temp_dir is None else
                           os.path.abspath(args.temp_dir[0]))
        print(self['temp_wd'])
        if 'QuartetSampling' not in self['temp_wd']:
            raise RuntimeError("Temporary directory name set by "
                               "-e/--temp-dir must contain 'QuartetSampling' "
                               "in the name to prevent accidental file "
                               "deletion!")
        print("setting the temp working dir to {}".format(self["temp_wd"]))
        self['result_prefix'] = (args.result_prefix[0]
                                 if args.result_prefix is not None
                                 else "RESULT")
        self['results_dir'] = (os.path.abspath(os.path.curdir) if
                               args.results_dir is None else
                               args.results_dir[0])
        self['score_result_file_path'] = os.path.join(
            self['results_dir'], "{}.node.scores.csv".format(
                self['result_prefix']))
        self['run_stats_file_path'] = os.path.join(
            self['results_dir'], "{}.run.stats".format(
                self['result_prefix']))
        self['tree_result_file_path'] = os.path.join(
            self['results_dir'], "{}.labeled.tre".format(
                self['result_prefix']))
        self['figtree_file_path'] = "{}.figtree".format(
            self['tree_result_file_path'])
        self['freq_file_path'] = "{}.freq".format(
            self['tree_result_file_path'])
        self['qc_tree_file_path'] = "{}.qc".format(
            self['tree_result_file_path'])
        self['qd_tree_file_path'] = "{}.qd".format(
            self['tree_result_file_path'])
        self['qi_tree_file_path'] = "{}.qi".format(
            self['tree_result_file_path'])
        self['verbout'] = args.verbout
        if args.verbout is True:
            self['verbout_file_path'] = "{}/{}.verbout".format(
                self['results_dir'], self['result_prefix'])
            with open(self['verbout_file_path'], "w") as verbout:
                verbout.write("topo1,topo2,topo3,topou,qc,qd,qi\n")
        self['max_random_sample_proportion'] = default_mrsp
        if args.max_random_sample_proportion:
            self['max_random_sample_proportion'] = (
                args.max_random_sample_proportion)
            print("setting the maximum proportion of "
                  "possible quartets to be explored to: {}".format(
                      args.max_random_sample_proportion))
            if (args.max_random_sample_proportion >
                    default_mrsp):
                print("WARNING: for some alignments, the quartet "
                      "randomization procedure may take a long time to finish "
                      "(or fail) when max proportion of quartets to sample "
                      "is greater than {}".format(
                          default_mrsp))
        if args.lnlike_thresh is not None:
            if isinstance(args.lnlike_thresh, list):
                self['lnlikethresh'] = args.lnlike_thresh[0]
            else:
                self['lnlikethresh'] = args.lnlike_thresh
            print("setting the minimum lnL thresh to {}".format(
                self['lnlikethresh']))
        else:
            self['lnlikethresh'] = 0
        return ''

        if self['data_type'] == 'cat' and not self.get('paup', False):
            raise RuntimeError("-d/-datatype 'cat' only currently enabled "
                               "for PAUP mode.")

    def __str__(self):
        """Print all parameters in alphabetical order"""
        return ('\n'.join(sorted("{}={}".format(k, v)
                                 for k, v in self.items())))


def read_config(configfilepath):
    """Reads the config file and returns command line arguments
    """
    args = []
    with open(configfilepath, 'r') as cfile:
        for line in [l.strip() for l in cfile.readlines()]:
            if len(line) < 1 or line[0] == '#':
                continue
            if '=' in line:
                elems = [d.strip() for d in line.split('=')]
            else:
                elems = [d.strip() for d in line.split()]
            args.extend(elems)
    return args


if __name__ == "__main__":
    print("This file is a function library, please run quartet_sampling.py")
