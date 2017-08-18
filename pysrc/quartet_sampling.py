#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# vim:fileencoding=utf-8
"""
quartet_samping.py: Quartet Sampling method for
phylogenetic branch support evaluation

<http://www.github.com/FePhyFoFum/quartetsampling>
"""


import argparse
import os
import sys
import time
from multiprocessing import Manager, Pool
from shutil import copyfile
from tree_data import TreeData, write_test_trees
from rep_data import DataStore, process_replicate_raxml
from rep_data import process_replicate_raxml2lk, process_replicate_paup
from rep_data import get_replicates_exhaustive, get_replicates_random
from rep_data import write_run_stats
from paramset import ParamSet, read_config
from alignment import Alignment


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
        prog="quartet_sampling.py",
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog=LICENSE)
    parser.add_argument("-t", "--tree", type=open, nargs=1, required=True,
                        help=("The input tree in Newick "
                              "(parenthetical) format."))
    parser.add_argument("-a", "--alignment", type=open, nargs=1, required=True,
                        help=("Alignment file in \"relaxed phylip\" format, "
                              "as used by RAxML."))
    parser.add_argument("-N", "--number-of-reps", type=int, nargs=1,
                        required=True, default=100,
                        help=("The number of replicate quartet topology "
                              "searches to be performed at each node."))
    parser.add_argument("-T", "--number-of-threads", type=int, nargs=1,
                        required=True, default=1,
                        help=("The number of parallel threads to be used "
                              "by Python for quartet topology searches."))
    parser.add_argument("-L", "--lnlike-thresh", type=float, nargs=1,
                        default=2.0,
                        help=("The lnlike threshhold that is the minimum "
                              "value by which the log-likelihood value "
                              "of the best-likelihood tree must be "
                              "higher than the second-best-likelihood tree "
                              "for the replicate to register as the "
                              "best-likelihood topology rather than "
                              "'uncertain'. If set to zero, this turns off "
                              "likelihood evaluation mode and invokes tree "
                              "inference mode where a tree is simply inferred "
                              "from the alignment without considering "
                              "likelihood (QI values are N/A in this case)."))
    parser.add_argument("-r", "--result-prefix", type=str, nargs=1,
                        help="A prefix to put on the result files.")
    parser.add_argument("-d", "--data-type", choices=('nuc', 'amino', 'cat'),
                        default=["nuc"], nargs=1,
                        help=("(nuc)leotide, (amino) acid, "
                              "or (cat)egorical data"))
    parser.add_argument("-O", "--min-overlap", type=int,
                        help=("The minimum sites required to be sampled for "
                              "all taxa in a given quartet."))
    parser.add_argument("-o", "--results-dir", type=os.path.abspath, nargs=1,
                        help=("A directory to which output files will "
                              "be saved. If not supplied, the current working "
                              "directory will be used. (default is current "
                              "folder)."))
    parser.add_argument("-V", "--verbout", action="store_true",
                        help=("Provide output of the frequencies of each "
                              "topology and QC."))
    parser.add_argument("-q", "--partitions", type=os.path.abspath, nargs=1,
                        help=("Partitions file in RAxML format. If omitted "
                              "then the entire alignment will be treated "
                              "as one partition for all quartet replicate "
                              "topology searches."))
    parser.add_argument("-g", "--genetrees", type=os.path.abspath, nargs=1,
                        help=("Use partitions file (RAxML format) to divide "
                              "the alignment into separate gene tree regions. "
                              "Gene alignments will be sampled random for the "
                              "quartet topology searches."))
    parser.add_argument("-e", "--temp-dir", type=os.path.abspath, nargs=1,
                        help=("A directory to which temporary files will be "
                              "saved. If not supplied, 'QuartetSampling' "
                              "will be created in the current "
                              "working directory. "
                              "When specifying a custom temporary output "
                              "the characters 'QuartetSampling' must appear "
                              "in the directory name to prevent accidental "
                              "file deletion. (default='./QuartetSampling'"))
    parser.add_argument("--retain-temp", action="store_true",
                        help=("Do not remove temporary files"))
    parser.add_argument("-C", "--clade", type=str,
                        help=("Conduct analysis on specific clade identified "
                              "by CSV taxon list"))
    parser.add_argument("-s", "--start-node-number", type=int, nargs=1,
                        help=("An integer denoting the node to which to start "
                              "from. Nodes will be read from topologically "
                              "identical (and isomorphic!) input trees in "
                              "deterministic order, so this argument may be  "
                              "used to restart at an intermediate position "
                              "(in case the previous run was canceled before "
                              "completion, for example)."))
    parser.add_argument("-p", "--stop-node-number", type=int, nargs=1,
                        help=("An integer denoting the node at which to stop. "
                              "Will include nodes with indices <= the stop "
                              "node number. This argument may be used to "
                              "limit the length of a given run in case only "
                              "a certain part of the tree is of interest. "
                              "Nodes will be read from topologically "
                              "identical (and isomorphic!) input trees "
                              "in deterministic order."))
    parser.add_argument("-X", "--raxml-executable", nargs=1,
                        help=("The name (or absolute path) of the raxml "
                              "executable to be used for calculating "
                              "likelihoods on quartet topologies."
                              "(default='raxml')"))
    parser.add_argument("--raxml-model", nargs=1,
                        help=("Advanced: specify a custom RAxML model name "
                              "for the raxml '-m' parameter"))
    parser.add_argument("-P", "--paup", action="store_true",
                        help="Use PAUP instead of RAxML.")
    parser.add_argument("--paup-executable", nargs=1, default=["paup"],
                        help=("The name or path of the PAUP executable to "
                              "be used for calculated quartets."))
    parser.add_argument("--ignore-errors", action="store_true",
                        help=("Ignore RAxML and PAUP erroneous runs"))
    parser.add_argument("--low-mem", action="store_true",
                        help=("Do not store large alignment in memory "
                              "for whole-alignment (non-genetree) mode"))
    parser.add_argument('--max-random-sample-proportion', type=float,
                        help=("The proportion of possible replicates explored "
                              "unsuccessfully by the random generation "
                              "procedure before it gives up. Because this "
                              "generates random replicates, it takes "
                              "progressively longer as it proceeds. To avoid "
                              "long runtimes, the recommended range is < 0.5 "
                              "(which is the default)."))
    parser.add_argument("--calc-qdstats", action="store_true",
                        help=("EXPERIMENTAL: Calculates Chi-square test "
                              "for QD tree frequencies. Use only "
                              " if Scipy is available. "
                              "Will increase running time."))
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="Provide more verbose output if specified.")
    parser.add_argument('--version', action='version',
                        version='%(prog)s version 1.2')
    return parser


def main(arguments=None):
    """Main method for quartet_sampling"""
    if arguments is None:
        if (len(sys.argv) == 2 and
                sys.argv[1] not in ('-h', '--help', '--version')):
            arguments = read_config(sys.argv[1])
            print("Config file used.")
            print("Executing with arguments: ", " ".join(arguments))
        else:
            arguments = sys.argv[1:]
    parser = generate_argparser()
    args = parser.parse_args(arguments)
    treedata = TreeData(args)
    params = ParamSet()
    params.setup(args, treedata.nleaves)
    if args.verbose:
        print("-----------")
        print("PARAMETERS:")
        print(params)
        print("-----------")
    maindata = DataStore(params)
    #  shared object access for multithreading
    manager = Manager()
    lock = manager.RLock()
    aln = Alignment(params)
    if params['using_genetrees']:
        aln.read_genes(args.alignment[0], params)
    else:
        aln.read_align(args.alignment[0], params)
    params['min_overlap'] = aln.min_overlap
    # k is the node counter
    k = 1
    #  if we are starting at the beginning, initialize the results file
    #  (otherwise assume it's already there and don't overwrite it)
    if not params['startk'] > k:
        maindata.write_headers(params['score_result_file_path'])
    # process the nodes in the tree
    params['starttime'] = time.time()
    for fnode in treedata.tree.iternodes():
        if params['verbose'] is True:
            print("testing node", [x.label for x in fnode.leaves()])
        if treedata.clade is not None:
            if fnode is not treedata.clade:
                continue
        os.chdir(params['temp_wd'])
        if k > params['stopk']:
            print("Processed all nodes up to the stop node. Exiting...")
            break
        write_test_trees()
        # skip tips and root
        k, leafsets = treedata.check_node(fnode, k, params)
        if leafsets is False:
            if params['verbose'] is True:
                print("skipping node...")
            continue
        # Begin multiprocessing queue
        results_queue = manager.Queue()
        n_completed = manager.Value("i", 0, "lock")
        # Establish replicates
        n_possible_replicates = 1
        for leafset in leafsets.values():
            n_possible_replicates *= len(leafset)
        if params['using_genetrees']:
            n_possible_replicates *= len(aln.seqs)
            if params['verbose'] is True:
                print('number of possible gene-quartet combos: {}'.format(
                    n_possible_replicates))
        elif params['verbose'] is True:
            print('number of possible quartets: {}'.format(
                n_possible_replicates))
        if (n_possible_replicates *
                params['max_quartet_enumeration_threshold'] < params['nreps']):
            if params['verbose'] is True:
                print('Number of possible quartets is close enough to the '
                      'total number to be sampled, so will generate all'
                      'and do a random draw')
            replicates, repstats = get_replicates_exhaustive(
                n_completed, results_queue, leafsets,
                params, aln, fnode, lock)
        else:
            if params['verbose']:
                print('Generating random quartets...')
            replicates, repstats = get_replicates_random(
                n_completed, results_queue, leafsets,
                params, aln, fnode, lock)
        nreplicates = len(replicates)
        if nreplicates < 1:  # no suitable replicates
            maindata.process_empty_rep_results(fnode, params, nreplicates)
        else:
            # copy original partitions file, should not change throughout run
            if params['partitions_file_path'] is not None:
                copyfile(params['partitions_file_path'], "temp_parts")
            # run the raxml calls in parallel
            # now designate multiprocessing resource pool.
            # important to do outside node loop. garbage collecting does not
            # apply to threads! set maxtasksperchild to release mem and files
            pool = Pool(params['nprocs'], maxtasksperchild=1)
            if params['paup'] is True:
                pool.map(process_replicate_paup, replicates)
            elif params['lnlikethresh'] > 0:
                pool.map(process_replicate_raxml2lk, replicates)
            else:
                pool.map(process_replicate_raxml, replicates)
            pool.close()
            pool.join()
            del pool
            # print("")
            # now process the results. first open a file to hold topologies
            # sending params['just_clade'] = True will give back detailed
            # name results
            maindata.process_rep_results(fnode, results_queue, params,
                                         nreplicates)
        # clean up
        del results_queue
        del n_completed
        # break # Left in place for troubleshooting
    if params['retain_temp'] is False:
        for the_file in os.listdir(params['temp_wd']):
            file_path = os.path.join(params['temp_wd'], the_file)
            try:
                if os.path.isfile(file_path):
                    if "QuartetSampling" not in file_path:
                        print(file_path,
                              " does not contain 'QuartetSampling' "
                              "and will not be deleted for safety")
                    else:
                        os.remove(file_path)
            except FileNotFoundError as exc:
                print(file_path, " not found")
        if 'QuartetSampling' in params['temp_wd']:
            os.rmdir(params['temp_wd'])
    qf_scores = maindata.write_qf_scores(params["score_result_file_path"])
    treedata.write_figtree(params['figtree_file_path'], qf_scores)
    treedata.write_scoretrees(params)
    write_run_stats(repstats, params)
    print(("\ndone.\nscores written to: {}\nlabeled "
           "tree written to: {}\ntotal time {:.2f} hours").format(
               params['score_result_file_path'],
               params['tree_result_file_path'],
               (time.time() - params['starttime']) / 3600))
    return ''


if __name__ == "__main__":
    main()
