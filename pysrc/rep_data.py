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

import sys
import os
import random
import math
import subprocess
from itertools import product
from shutil import copyfile
from phylo import tree_reader
from phylo import tree_utils


class DataStore():
    """Main data storage and output functions"""

    def __init__(self, params):
        self.tree_counts = {}
        self.leaf_counts = {}
        self.node_replicates = {}
        for dirname in (params['results_dir'], params['temp_wd']):
            if not os.path.exists(dirname):
                print("creating {}".format(dirname))
                os.mkdir(dirname)
        self.headers = {"main": ["node_label", "freq0",
                                 "qc", "qd", "qi", "qf",
                                 "qdsig",
                                 "diff", "num_replicates", "notes"],
                        "clade": ["taxon", "tree1", "tree2", "tree3", "treeu",
                                  "qc", "qd", "qi",
                                  "qdsig", "freq0"]}

    def write_headers(self, file_path, restype="main"):
        """Write the headers into the file"""
        with open(file_path, "w") as resultsfile:
            resultsfile.write("{}\n".format(','.join(self.headers[restype])))
        return ''

    def write_entry(self, file_path, entry, restype="main"):
        """Add an entry to the file"""
        with open(file_path, "a") as resultsfile:
            resultsfile.write("{}\n".format(','.join([
                str(entry.get(x, "NA")) for x in self.headers[restype]])))
        return ''

    def process_rep_results(self, fnode, results_queue, params, nreplicates):
        """Process the output from a replicate"""
        record_detail = False
        detail_name_sets = []
        detail_tree_sets = []
        notes = ''
        tree_counts_detailed = 0
        if fnode not in self.tree_counts:
            self.tree_counts[fnode] = {}
        while not results_queue.empty():
            result = results_queue.get()
            rep_info = {'diff_exceeds_cutoff': None,
                        'best_tree': None,
                        'likelihood_diff': None}
            if params['verbose']:
                print("---")
                print("Seqnames: ", result['seq_names'])
            if params['paup']:
                # using paup
                rep_info = result
            else:
                rep_info = result
            if rep_info['diff_exceeds_cutoff'] is False:
                rep_info['best_tree'] = 3
            # case of not exceeding the likelihood = tree 3 for QI calc
            self.tree_counts[fnode][rep_info['best_tree']] = (
                self.tree_counts[fnode].get(rep_info['best_tree'], 0) + 1)
            for seqname in result['seq_names'].values():
                if seqname not in self.leaf_counts:
                    self.leaf_counts[seqname] = {}
                self.leaf_counts[seqname][rep_info['best_tree']] = (
                    self.leaf_counts[seqname].get(
                        rep_info['best_tree'], 0) + 1)
        # calcluate the q scores
        # qc_score, qd_score, qi_score, freq0
        qscores = calc_qc_qd_qi(
            self.tree_counts[fnode], params)
        fnode.data["qc_score"] = na_fmt(qscores['qc'])
        fnode.data["qd_score"] = na_fmt(qscores['qd'])
        fnode.data["qd_sig"] = na_fmt(qscores['qdsig'])
        fnode.data["qi_score"] = na_fmt(qscores['qi'])
        fnode.data['freq0'] = na_fmt(qscores['freq0'])
        fnode.data['replicates'] = na_fmt(nreplicates)
        # write the scores to the file
        self.write_entry(
            params["score_result_file_path"],
            {"node_label": fnode.label,
             "freq0": qscores['freq0'],
             "qc": qscores['qc'], "qd": qscores['qd'],
             'qdsig': qscores['qdsig'],
             "qi": qscores['qi'], "diff": rep_info['likelihood_diff'],
             "num_replicates": nreplicates, "notes": notes})
        # this will record the taxa included and the tree
        tree_counts_detailed = tree_counts_detailed + (
            rep_info['likelihood_diff']
            if rep_info['likelihood_diff'] is not None else 0)
        if record_detail:
            detail_name_sets.append(set(result["seq_names"].values()))
            detail_tree_sets.append(rep_info['best_tree'])
        if params['just_clade']:
            fnode_dict = {}
            for j in detail_name_sets:
                for k in j:
                    if k not in fnode_dict:
                        fnode_dict[k] = {0: 0, 1: 0, 2: 0, 3: 0}
            count = 0
            for entry in detail_name_sets:
                for subentry in entry:
                    fnode_dict[subentry][detail_tree_sets[count]] += 1
                count += 1
                clade_file_path = os.path.join(
                    params['temp_wd'],
                    "{}.clade".format(
                        params['score_result_file_path']))
                self.write_headers(clade_file_path, restype="clade")
                for xnode in fnode_dict:
                    nqscores = calc_qc_qd_qi(
                        fnode_dict[xnode], params)
                    self.write_entry(clade_file_path, {
                        "taxon": entry,
                        "tree1": fnode_dict[xnode][0],
                        "tree2": fnode_dict[xnode][1],
                        "tree3": fnode_dict[xnode][2],
                        "treeu": fnode_dict[xnode][3],
                        "qc": nqscores['qc'],
                        'qd': nqscores['qd'],
                        "qdsig": nqscores['qdsig'],
                        "qi": nqscores['qi'],
                        "freq0": nqscores['freq0']},
                                     restype="clade")
        return ''

    def process_empty_rep_results(self, fnode, params, nreplicates):
        """Process the results from an empty replicate"""
        fnode.data.update({
            "replicates": "{:.2g}".format(0),
            "freq0": "{:.2g}".format(0),
            "qc_score": "NA", "qd_score": "NA", "qi_score": "NA",
            'qd_sig': 'NA'})
        # write the scores to the file
        self.write_entry(
            params["score_result_file_path"],
            {"node_label": fnode.label,
             "freq0": 0,
             "num_replicates": nreplicates,
             "notes": 'found no suitable replicates'})
        return ''

    def write_qf_scores(self, outfile):
        """Writes the Quartet Taxon (Rogue taxon) score"""
        qf_scores = {}
        for fnode in sorted(self.leaf_counts):
            total = float(sum([self.leaf_counts[fnode].get(x, 0)
                               for x in (0, 1, 2, 3)]))
            qf_score = self.leaf_counts[fnode].get(0, 0.) / float(total)
            self.write_entry(
                outfile, {"node_label": fnode, "qf": qf_score,
                          "num_replicates": total, "notes": ''})
            qf_scores[fnode] = qf_score + 0.0
        return qf_scores


def na_fmt(num):
    """Formats float or returns NA"""
    return "NA" if num == 'NA' else "{:.2g}".format(num)


def chi2_test(val0, val1):
    """Calculate Pearson Chi-Squared for the special case of
       two values that are expected to be equal
       Arguments:
           val0: first value
           val1: second value
    """
    from scipy.stats import chi2
    try:
        chisq = float((val0 - val1)**2) / float(val0 + val1)
        if not chisq:
            return (0, 1)
        pval = 1.0 - chi2.cdf(chisq, 1)
        return (chisq, pval)
    except ZeroDivisionError as _:
        return (0, 1)


def calc_qc_qd_qi(counts, params):
    """Calculate the QC, QD, and QI scores"""
    if params['verbose'] is True:
        print(counts)
    total = float(sum(counts.get(x, 0) for x in (0, 1, 2)))
    utotal = total + counts.get(3, 0)
    if utotal == 0:
        qscores = {'qc': 'NA', 'qd': 'NA', 'qi': 'NA',
                   'freq0': 'NA', 'qdsig': 'NA'}
    elif total == 0:
        qscores = {'qc': 0, 'qd': 0, 'qi': 1, 'freq0': 0,
                   'qdsig': 'NA'}
    else:
        freqs = [counts.get(x, 0) / total for x in (0, 1, 2)]
        ufreqs = [counts.get(x, 0) / utotal for x in (0, 1, 2, 3)]
        ncounts = sum([int(x > 0) for x in freqs])
        qc_score = 1
        if ncounts > 1:
            qc_score += sum([freqs[i] * math.log(freqs[i], ncounts)
                             if freqs[i] > 0.0 else 0.0
                             for i in (0, 1, 2)])
        if counts.get(0, 0) < max(counts.get(1, 0), counts.get(2, 0)):
            qc_score *= -1
        qd_score = (1.0 - (abs(counts.get(1, 0) - counts.get(2, 0)) / (
            counts.get(1, 0) + counts.get(2, 0)))
                    if counts.get(1, 0) + counts.get(2, 0) > 0 else "NA")
        if params['calc_qdstats'] is False:
            qd_sig = 'NA'
        else:
            _, qd_sig = chi2_test(counts.get(1, 0), counts.get(2, 0))
        qi_score = ('NA' if params['lnlikethresh'] == 0
                    else (1.0 - (counts.get(3, 0) / utotal)))
        if params['verbout'] is True:
            with open(params['verbout_file_path'], 'a') as vfile:
                vfile.write('{}\n'.format(','.join([
                    str(x) for x in ufreqs + [
                        qc_score, qd_score, qi_score]])))
        qscores = {'qc': qc_score,
                   'qd': qd_score, 'qdsig': qd_sig,
                   'qi': qi_score, 'freq0': freqs[0]}
    return qscores


def empty_rep(j, fnode, lock, results_queue, n_completed, params):
    """Create empty replicate object"""
    rep = params.copy()
    rep["queue"] = results_queue
    rep["lock"] = lock
    rep["node_id"] = fnode.label
    rep["replicate_id"] = str(j)
    rep["n_completed"] = n_completed
    rep["seqs"] = {}
    rep["seq_names"] = {}
    rep["genename"] = None
    return rep


def get_replicates_exhaustive(n_completed, results_queue, leafsets,
                              params, aln, fnode, lock):
    """Get a single sampling replicate set"""
    replicates = []
    repstats = {}
    # need to make sure we don't repeat, so generate all the quartets
    if params['using_genetrees']:
        possible_quartets = list(product(leafsets['L1'], leafsets['L2'],
                                         leafsets['R1'], leafsets['R2'],
                                         aln.genes))
    else:
        possible_quartets = list(product(leafsets['L1'], leafsets['L2'],
                                         leafsets['R1'], leafsets['R2'],
                                         [None]))
    # look through them in random order for suitable ones
    nonoverlapping_count = 0
    while len(replicates) < params['nreps'] and len(possible_quartets) > 0:
        proposed_quartet = random.sample(possible_quartets, 1)[0]
        possible_quartets.remove(proposed_quartet)
        if aln.check_aln_overlap(proposed_quartet) is False:
            nonoverlapping_count += 1
            if params['verbose']:
                print('non-overlap count: {}'.format(nonoverlapping_count))
            continue
        # if we made it here then the proposed rep is acceptable
        rep = empty_rep(len(replicates), fnode, lock,
                        results_queue, n_completed, params)
        for i, subtree_name in enumerate(['L1', 'L2', 'R1', 'R2']):
            leaf_name = proposed_quartet[i]
            if params['using_genetrees']:
                rgenename = proposed_quartet[4]
                rep['genename'] = rgenename[:]
                rep['seqs'][subtree_name] = aln.seqs[rgenename][leaf_name]
                rep['seq_names'][subtree_name] = leaf_name
            else:
                rep['seqs'][subtree_name] = aln.seqs[leaf_name]
                rep['seq_names'][subtree_name] = leaf_name
        replicates.append(rep)
        rep["unique_label"] = "{}.{}".format(
            rep["node_id"], rep["replicate_id"])
        rep["aln_fname"] = os.path.join(
            params['temp_wd'],
            "temp_inseqs.{}".format(rep["unique_label"]))
        # write file for successful rep
        if rep['paup'] is True:
            write_paup(rep["aln_fname"], rep["seqs"])
        else:
            if params["low_mem"] is True:
                # print(rep['seqs'], rep['seq_names'])
                cat_raxml(rep["aln_fname"], rep["seqs"], rep["seq_names"])
            else:
                write_raxml(rep["aln_fname"], rep["seqs"])
        del rep["seqs"]
        if rep["using_partitions"]:
            # make a copy of the partitions file
            newpartfile = os.path.join(
                params['temp_wd'],
                "temp_parts.{}".format(rep["unique_label"]))
            copyfile(params['partitions_file_path'], newpartfile)
            rep["part_fname"] = newpartfile[:]
    if len(replicates) < 1:
        print('WARNING: generated all possible quartets '
              'and did not find a suitable one! If you have the -O '
              'flag enabled, alignment may not have enough data at '
              'this edge (i.e. low partial decisiveness).')
    elif len(replicates) < params['nreps']:
        print('WARNING: only {} suitable replicates for this node.'.format(
            len(replicates)))
    repstats['nonoverlapping_count'] = nonoverlapping_count + 0
    return replicates, repstats


def get_replicates_random(n_completed, results_queue, leafsets,
                          params, aln, fnode, lock):
    """Get random replicate quartets using full sampling"""
    replicates = []
    repstats = {}
    n_possible_replicates = 1
    for leafset in leafsets.values():
        n_possible_replicates *= len(leafset)
    if params['using_genetrees']:
        n_possible_replicates *= len(aln.seqs)
    observed_quartets = set([])
    for j in range(params['nreps']):
        if params['verbose']:
            print('looking for unique replicate {}'.format(j))
        rep = None
        duplicate_count = 0
        nonoverlapping_count = 0
        attempted_count = 0
        observed_count = 0
        # Set the cutoff for replicate attempts at proportion of total
        # with ceiling of 1000000
        max_attempts = min(n_possible_replicates *
                           params['max_random_sample_proportion'],
                           params['nreps'] * 200)
        # loop until we find an acceptable replicate, or we hit the
        # maximum allowed proportion to attempt
        while attempted_count < max_attempts:
            proposed_quartet = set()
            proposed_rep = empty_rep(len(replicates), fnode, lock,
                                     results_queue, n_completed, params)
            # generate a random replicate
            rgenename = None
            if params['using_genetrees']:
                rgenename = aln.get_random_gene()
            for subtree_name, leaf_names in leafsets.items():
                randleaf = random.sample(leaf_names, 1)[0]
                # should really be doing this check when loading files
                if params['using_genetrees']:
                    if randleaf not in aln.seqs[rgenename]:
                        raise RuntimeError(
                            "FATAL ERROR: name {} not in alignment".format(
                                randleaf))
                    if params['verbose']:
                        print(" {} = {}".format(subtree_name, randleaf))
                    proposed_rep['seq_names'][subtree_name] = randleaf
                    proposed_rep['genename'] = rgenename[:]
                    proposed_rep['seqs'][subtree_name] = (
                        aln.seqs[rgenename][randleaf])
                else:
                    if randleaf not in aln.seqs:
                        raise RuntimeError(
                            "FATAL ERROR: name {} not in alignment".format(
                                randleaf))
                    if params['verbose']:
                        print(" {} = {}".format(subtree_name, randleaf))
                    proposed_rep['seqs'][subtree_name] = aln.seqs[randleaf]
                    proposed_rep['seq_names'][subtree_name] = randleaf
                proposed_quartet.add(randleaf)
            proposed_quartet = tuple(list(sorted(proposed_quartet)) +
                                     [rgenename])
            attempted_count += 1
            if proposed_quartet in observed_quartets:
                duplicate_count += 1
                if params['verbose']:
                    print('duplicate count= ', duplicate_count,
                          ', possible=', n_possible_replicates)
                continue
            # quartet is unique, remember that we tried it
            observed_quartets.add(proposed_quartet)
            observed_count += 1
            if aln.check_aln_overlap(proposed_quartet) is False:
                nonoverlapping_count += 1
                if params['verbose']:
                    print('non-overlap count = {}'
                          ', possible = {}'.format(
                              nonoverlapping_count,
                              n_possible_replicates))
                continue
            # if we made it here then the proposed rep is acceptable
            rep = proposed_rep
            if params['verbose']:
                print("passed taxa", ",".join(list(proposed_quartet[:4])))
                print("passed gene {}".format(rgenename))
            # generate labels for temp files
            rep["unique_label"] = "{}.{}".format(
                rep["node_id"], rep["replicate_id"])
            rep["aln_fname"] = os.path.join(
                params['temp_wd'],
                "temp_inseqs.{}".format(rep["unique_label"]))
            # write file for successful rep
            if rep['paup'] is True:
                write_paup(rep["aln_fname"], rep["seqs"],
                           datatype=rep['data_type'])
            else:
                if params["low_mem"] is True:
                    cat_raxml(rep["aln_fname"], rep["seqs"], rep["seq_names"])
                else:
                    write_raxml(rep["aln_fname"], rep["seqs"])
            del rep["seqs"]
            if rep["using_partitions"]:
                # make a copy of the partitions file
                newpartfile = os.path.join(
                    params['temp_wd'],
                    "temp_parts.{}".format(rep["unique_label"]))
                copyfile(params['partitions_file_path'], newpartfile)
                rep["part_fname"] = newpartfile[:]
                # make a copy of the partitions file
            break
        # if there is no rep, then we hit the max number of attempts
        if rep is None:
            print(("WARNING: attempted ~{:0.0f} of possible quartets "
                   " ~{:0.0f} of which were unique "
                   "and did not find a suitable one! If you have the -O "
                   "flag enabled, alignment may not have enough data at "
                   "this edge (i.e. low partial decisiveness).").format(
                       attempted_count, observed_count))
            break
        elif params['verbose']:
            print('constructed replicate {}\n'.format(j))
        repstats['duplicate_count'] = duplicate_count + 0
        repstats['attempted_count'] = attempted_count + 0
        repstats['observed_count'] = observed_count + 0
        repstats['non_overlapping_count'] = nonoverlapping_count + 0
        replicates.append(rep)
    return replicates, repstats


def process_replicate_raxml2lk(replicate):
    """Process individual replicate sampling"""
    os.chdir(replicate["temp_wd"])
    # just alias dictionary elements for convenience
    result = {}
    result["label"] = replicate['unique_label']
    temp_file_paths = [replicate['aln_fname'],
                       "{}.reduced".format(replicate['aln_fname'])]
    result["seq_names"] = replicate["seq_names"].copy()
    # generate a label that will be unique within this run
    #  (but probably not among runs!)
    temp_ml_search_label = "tts.{}".format(replicate["unique_label"])
    # this will test the three topologies
    # test alignment readability by raxml, also filters missing columns
    base_raxml_args = [replicate['raxml_executable'],
                       "-s", replicate['aln_fname'],
                       "-m", replicate['raxml_model'],
                       "-T", "1",
                       "-p", "11341",
                       "--silent",
                       "-F",
                       "-f", "N"]
    if replicate['using_partitions']:
        temp_file_paths.extend([replicate['part_fname'],
                                "{}.reduced".format(replicate['part_fname'])])
        base_raxml_args.extend(["-q", replicate["part_fname"]])
    treelikelihoods = {0: 0, 1: 0, 2: 0}
    likelihood_diff_exceeds_cutoff = False
    # correct = None
    for i in range(3):
        raxml_args = base_raxml_args + [
            "-n", "{}.{}".format(temp_ml_search_label, i),
            "-z", "test.trees.{}".format(i),
            ]
        result["raxml_args"] = " ".join(raxml_args)
        if replicate['verbose']:
            print('calling: {}'.format(result["raxml_args"]))
        proc = subprocess.Popen(raxml_args, stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
        serr, sout = proc.communicate()
        if replicate['verbose'] is True:
            if "fix your data" in serr.decode('utf-8'):
                print("PARTITIONS WARNING OR MISSING DATA FROM RAXML")
        lpath = os.path.join(replicate['temp_wd'],
                             "RAxML_info.tts.{}.{}".format(
                                 result["label"], i))
        tpath = os.path.join(replicate['temp_wd'],
                             "RAxML_result.tts.{}.{}".format(
                                 result["label"], i))
        temp_file_paths.extend([lpath, tpath])
        if not os.path.exists(lpath):
            if replicate['ignore_error'] is False:
                raise RuntimeError("'{}' does not exist".format(lpath))
        with open(lpath, "r") as info_result:
            for line in info_result:
                if "Tree 0 Likelihood " in line:
                    treelikelihoods[i] = (
                        -1 * float(line.split(" ")[3]))
    srt_likelihoods = [(treelikelihoods[x], x) for x in (0, 1, 2)]
    srt_likelihoods.sort()
    likelihood_diff = abs(srt_likelihoods[0][0] - srt_likelihoods[1][0])
    if likelihood_diff > replicate['lnlikethresh']:
        likelihood_diff_exceeds_cutoff = True
    if replicate['verbose']:
        print("-" * 5)
        print("Filename: ", lpath)
        print("Likelihoods:", treelikelihoods)
        print("Lowest L:", srt_likelihoods[0][1], "; Diff = ", likelihood_diff,
              "; Exceeds cutoff = ", likelihood_diff_exceeds_cutoff)
    result["diff_exceeds_cutoff"] = likelihood_diff_exceeds_cutoff
    result["best_tree"] = srt_likelihoods[0][1]
    result["likelihood_diff"] = likelihood_diff
    if replicate['retain_temp'] is False:
        for fpath in temp_file_paths:
            if os.path.exists(fpath):
                if replicate['verbose'] is True:
                    print("REMOVED", fpath)
                os.remove(fpath)
            elif replicate['verbose'] is True:
                print("NOTFOUND", fpath)
    replicate["queue"].put(result)
    # increment counter and update user feedback
    replicate["n_completed"].value += 1
    replicate["lock"].acquire()
    sys.stdout.flush()
    replicate['lock'].release()
    return ''


def process_replicate_raxml(replicate):
    """Process individual replicate sampling"""
    os.chdir(replicate["temp_wd"])
    temp_file_paths = [replicate['aln_fname'],
                       "{}.reduced".format(replicate['aln_fname'])]
    # just alias dictionary elements for convenience
    result = {"label": replicate['unique_label']}
    result["seq_names"] = replicate["seq_names"].copy()
    # generate a label that will be unique within this run
    #  (but probably not among runs!)
    temp_ml_search_label = "tts.{}".format(replicate["unique_label"])
    raxml_args = [replicate['raxml_executable'],
                  "-s", replicate['aln_fname'],
                  "-m", replicate['raxml_model'],
                  "-T", "1",
                  "-p", "11341",
                  "--silent",
                  "-F",
                  "-n", temp_ml_search_label]
    if replicate['using_partitions']:
        temp_file_paths.extend([replicate['part_fname'],
                                "{}.reduced".format(replicate['part_fname'])])
        raxml_args.extend(["-q", replicate["part_fname"]])
    result["raxml_args"] = " ".join(raxml_args)
    proc = subprocess.Popen(raxml_args, stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE)
    serr, sout = proc.communicate()
    if replicate['verbose'] is True:
        if "fix your data" in serr.decode('utf-8'):
            print("PARTITIONS WARNING OR MISSING DATA FROM RAXML")
    if replicate['verbose']:
        print('calling:{}'.format(result["raxml_args"]))
    result["label"] = replicate['unique_label']
    best_tree = None
    lpath = os.path.join(replicate['temp_wd'],
                         "RAXML_info.tts.{}".format(result["label"]))
    tpath = os.path.join(replicate['temp_wd'],
                         "RAxML_result.tts.{}".format(result["label"]))
    temp_file_paths.extend([lpath, tpath])
    with open(tpath, "r") as tfile:
        tline = tfile.readline()
    restree = tree_reader.read_tree_string(tline)
    likelihood1, likelihood2 = tree_utils.calc_biparts(restree)
    if ("L1" in likelihood1[0] and "L2" in likelihood1[0]) or (
            "R1" in likelihood1[0] and "R2" in likelihood1[0]):
        best_tree = 0
    elif ("L1" in likelihood1[0] and "R1" in likelihood1[0]) or (
            "L2" in likelihood1[0] and "R2" in likelihood1[0]):
        best_tree = 1
    elif ("L1" in likelihood1[0] and "R2" in likelihood1[0]) or (
            "L2" in likelihood1[0] and "R1" in likelihood1[0]):
        best_tree = 2
    if replicate['verbose']:
        print("---")
        print(best_tree, likelihood1, likelihood2)
    result['diff_exceeds_cutoff'] = True
    result['best_tree'] = best_tree
    result['likelihood_diff'] = 0
    if replicate['retain_temp'] is False:
        for fpath in temp_file_paths:
            if os.path.exists(fpath):
                os.remove(fpath)
    replicate["queue"].put(result)
    # increment counter and update user feedback
    replicate["n_completed"].value += 1
    replicate["lock"].acquire()
    sys.stdout.flush()
    replicate['lock'].release()
    return ''


def process_replicate_paup(replicate):
    """Process individual replicate sampling"""
    os.chdir(replicate["temp_wd"])
    temp_file_paths = [replicate['aln_fname']]
    result = {"label": replicate['unique_label']}
    # just alias dictionary elements for convenience
    paup_out_file_path = os.path.join(
        replicate['temp_wd'],
        "temp_inseqs.{}.out".format(result['label']))
    result["seq_names"] = replicate["seq_names"].copy()
    # write the alignment    result["label"] = replicate['unique_label']
    # this will test the three topologies
    paup_args = [replicate["paup_executable"], replicate['aln_fname']]
    result["paup_args"] = " ".join(paup_args)
    if replicate['verbose']:
        print('calling: ', result['paup_args'])
    proc = subprocess.Popen(paup_args, stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE)
    proc.communicate()
    treelikelihoods = {0: 0, 1: 0, 2: 0, 3: 0}
    likelihood_diff_exceeds_cutoff = False
    outfile = open(paup_out_file_path, "r")
    resstr = outfile.readlines()
    for elem in resstr:
        if "Tree" in elem:
            continue
        row = elem.strip().split()
        treelikelihoods[int(row[0]) - 1] = float(row[1])
    srt_likelihoods = [(treelikelihoods[x], x) for x in (0, 1, 2)]
    srt_likelihoods.sort()
    likelihood_diff = abs(srt_likelihoods[0][0] - srt_likelihoods[1][0])
    if likelihood_diff > replicate['lnlikethresh']:
        likelihood_diff_exceeds_cutoff = True
    if replicate['verbose']:
        print("-"*5)
        print("Filename: ", paup_out_file_path)
        print("Likelihoods:", treelikelihoods)
        print("Lowest L:", srt_likelihoods[0][1], "; Diff = ", likelihood_diff,
              "; Exceeds cutoff = ", likelihood_diff_exceeds_cutoff)
    result['diff_exceeds_cutoff'] = likelihood_diff_exceeds_cutoff
    result['best_tree'] = srt_likelihoods[0][1]
    result['likelihood_diff'] = likelihood_diff
    temp_file_paths.append(paup_out_file_path)
    if replicate['retain_temp'] is False:
        for fpath in temp_file_paths:
            if os.path.exists(fpath):
                os.remove(fpath)
    replicate["queue"].put(result)
    # increment counter and update user feedback
    replicate["n_completed"].value += 1
    replicate["lock"].acquire()
    sys.stdout.flush()
    replicate['lock'].release()
    return ''


def write_raxml(fname, seqs):
    """Write FASTA for RAxML"""
    with open(fname, "w") as outfile:
        for hdr, seq in seqs.items():
            outfile.write(">{}\n{}\n".format(hdr, seq))
    return ''


def cat_raxml(fname, seqs, seqnames):
    """Concatenate sequences for fasta"""
    with open(fname, 'w') as outfile:
        for repname, _ in seqnames.items():
            outfile.write(">{}\n".format(repname))
            with open(seqs[repname], 'r') as sfile:
                outfile.write(sfile.readline().rstrip())
            outfile.write("\n")
    return ''


def write_paup(fpath, seqs, datatype="nuc"):
    """Write PAUP output"""
    paup_datatype = 'dna'
    if datatype == 'amino':
        paup_datatype = 'protein'
    elif datatype == 'cat':
        paup_datatype = 'standard'
    test_trees = {0: "(R1,R2,(L1,L2));",
                  1: "(R1,L1,(L2,R2));",
                  2: "(R1,L2,(L1,R2));"}
    slen = len(seqs["L1"])
    with open(fpath, "w") as pfile:
        pfile.write("#nexus\n"
                    "begin data;\n"
                    "  dimensions ntax=4 nchar="+str(slen)+";\n"
                    "  format datatype=" + str(paup_datatype) + " missing=-;\n"
                    "  matrix\n")
        for hdr in seqs:
            pfile.write("    {} {}\n".format(hdr, seqs[hdr]))
        pfile.write("    ;\nend;\n\nbegin trees;\n")
        for hdr, val in test_trees.items():
            pfile.write("  utree t{} = {};\n".format(hdr, val))
        pfile.write(
            "end;\n\n"
            "begin paup;\n"
            # "  lset lcollapse=no precision=double nst=6"
            # " rmatrix=estimate basefreq=empirical;\n"
            "  lset lcollapse=no precision=double nst=1 basefreq=equal;\n"
            "  lscores all /scorefile="+fpath+".out replace=yes;\n"
            "  quit;\n"
            "end;\n")
    return ''


def write_run_stats(repstats, params):
    """Write run stats to file"""
    with open(params['run_stats_file_path'], 'w') as outfile:
        for k, val in repstats.items():
            outfile.write("{}\t{}\n".format(k, val))
    return ''


if __name__ == "__main__":
    print("This file is a function library, please run quartet_sampling.py")
