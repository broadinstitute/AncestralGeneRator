#!/usr/bin/env python
import os
import sys
from sys import stderr
from ete3 import Tree
import argparse
import subprocess

try:
    assert (sys.version_info[0] == 3)
except:
    stderr.write("Please use Python-3 to run this program. Exiting now ...\n");
    sys.exit(1)

script_dir = os.path.abspath(os.path.dirname(os.path.realpath(__file__))) + '/'

def is_number(x):
    try:
        float(x); return True
    except:
        return False

def create_nexus(name_map, tree, homolog_array, cmatrix, num_chars, unrooted, num_strains, outdir):
    nex_file = outdir + 'input.nex'

    tree_line = open(tree).readlines()[0]

    nex_file_handle = open(nex_file, 'w')
    nex_file_handle.write('#NEXUS\n\nBEGIN data;\ndimensions ntax=%d nchar=%d;\nformat symbols="0~%d";\nmatrix\n' % (
        num_strains, len(homolog_array) - 1, num_chars - 1))
    for row_split in zip(*homolog_array):
        new_row_split = [name_map[row_split[0]]]
        for val in row_split[1:]:
            if not type(cmatrix) == list:
                if val != '0':
                    new_row_split.append('1')
                else:
                    new_row_split.append('0')
            else:
                if int(val) >= num_chars:
                    new_row_split.append(str(num_chars - 1))
                else:
                    new_row_split.append(val)

        row = ' '.join(new_row_split) + ' '
        nex_file_handle.write(row + '\n')

    cmchars_string = ''.join([str(x) for x in range(0, num_chars)])
    cost_matrix = '0 10\n2 0\n'
    if type(cmatrix) == list:
        cost_matrix = ""
        for event in cmatrix:
            cost_matrix += ' '.join(event) + '\n'

    nex_file_handle.write(
        ';\nendblock;\n\nBEGIN ASSUMPTIONS;\nusertype costmatrix stepmatrix = %d %s\n' % (num_chars, cmchars_string))
    nex_file_handle.write(cost_matrix)
    nex_file_handle.write(';\ntypeset *a = costmatrix:all;\nENDBLOCK;\n\n')
    rooted_spec = '[&R]'
    if unrooted: rooted_spec = '[&U]'
    nex_file_handle.write('BEGIN trees;\nTREE mytree = %s %s\nEND;\n\n' % (rooted_spec, tree_line))
    nex_file_handle.write('BEGIN paup;\nshowuser;\ndescribe 1/xout=both;\nendblock;\n')
    nex_file_handle.close()
    return nex_file


def run_paup(paup_exe, nexus_input, outdir):
    paup_result_file = outdir + 'gene_flux_results.txt'
    oprf = open(paup_result_file, 'w')
    paup_cmd = [paup_exe, '-n', nexus_input]
    subprocess.call(paup_cmd, stdout=oprf)
    oprf.close()
    return paup_result_file


def create_necessary_files(paup_result_file, homolog_names, outdir):
    parent_child_file = outdir + 'parent_child.txt'
    homolog_name_file = outdir + 'homolog_names.txt'

    pcfo = open(parent_child_file, 'w')
    parse_cmd = ["perl", script_dir + "Related_Scripts/read_paup_tree_tri.pl", paup_result_file]
    subprocess.call(parse_cmd, stdout=pcfo)
    pcfo.close()

    onfo = open(homolog_name_file, 'w')
    for on in homolog_names:
        onfo.write(on + '\n')
    onfo.close()

    return [parent_child_file, homolog_name_file]


def create_results(paup_result_file, parent_child_file, homolog_names_file, name_map, outdir):
    final_result = open(outdir + 'parsed_results.txt', 'w')
    ind_res_dir = outdir + "individual_result_files/"
    script = str(script_dir) + 'Related_Scripts/MakeGeneFluxSpreadsheet.py'
    parse_cmd = ['python', script, parent_child_file, paup_result_file, homolog_names_file,
                 outdir + 'name_mapping.txt', ind_res_dir]
    subprocess.call(parse_cmd, stdout=final_result)
    final_result.close()
    return (outdir + 'parsed_results.txt')


def reformat_tree(phylo_tree, parsed_results, outdir):
    child_to_parent = {}
    with open(parsed_results) as opr:
        for i, line in enumerate(opr):
            line = line.rstrip('\n')
            ls = line.split('\t')
            if i > 0:
                child_to_parent[ls[0]] = ls[1]

    t = Tree(phylo_tree)
    for node in t.traverse("postorder"):
        if (node.name.lstrip('node')) in child_to_parent.keys():
            parent_name = child_to_parent[node.name.lstrip('node')]
            if is_number(parent_name):
                parent_name = 'node' + parent_name
            node.up.name = parent_name
        elif node.name in child_to_parent.keys():
            parent_name = child_to_parent[node.name]
            if is_number(parent_name):
                parent_name = 'node' + parent_name
            node.up.name = parent_name
    for node in t.traverse("postorder"):
        if (node.name.lstrip('node')) in child_to_parent.keys():
            parent_name = child_to_parent[node.name.lstrip('node')]
            if is_number(parent_name):
                parent_name = 'node' + parent_name
            node.up.name = parent_name

    t.write(features=["name"], format=3, outfile=outdir + "phylogeny_with_innernode_names.nwk")


def rename_nodes(parsed_results, rename_phylogeny, curr_phylogeny):
    incorrect_tree_names = Tree(curr_phylogeny, format=1)
    correct_tree_names = Tree(rename_phylogeny, format=1)

    correct_names = {};
    renaming_map = {}
    for n in correct_tree_names.traverse("postorder"):
        if not n.is_leaf():
            children = set([])
            for c in n.traverse("postorder"):
                if c.is_leaf(): children.add(c.name)
            children = tuple(sorted(list(children)))
            correct_names[children] = n.name

    for n in incorrect_tree_names.traverse("postorder"):
        if not n.is_leaf():
            children = set([])
            for c in n.traverse("postorder"):
                if c.is_leaf(): children.add(c.name)
            children = tuple(sorted(list(children)))
            renaming_map[n.name.lstrip('node')] = correct_names[children].lstrip('node')
            n.name = correct_names[children]

    incorrect_tree_names.write(features=["name"], format=1, outfile=curr_phylogeny)

    tmp_file = open(parsed_results + '.tmp', 'w')
    for i, line in enumerate(open(parsed_results)):
        if i == 0:
            tmp_file.write(line)
        else:
            ls = line.rstrip('\n').split('\t')
            if is_number(ls[0]) and ls[1] in renaming_map:
                ls[0] = renaming_map[ls[0]]
            if is_number(ls[1]) and ls[1] in renaming_map:
                ls[1] = renaming_map[ls[1]]
            tmp_file.write('\t'.join(ls) + '\n')
    tmp_file.close()
    os.system('mv %s %s' % (parsed_results + '.tmp', parsed_results))


def rename_strains(tree, outdir):
    name_map = {}
    outf = open(outdir + 'name_mapping.txt', 'w')
    new_tree = outdir + 'phylo.tre'
    node_count = 1
    old_tree = Tree(tree)
    for node in old_tree.traverse("postorder"):
        if node.is_leaf():
            new_name = 'N' + str(node_count)
            outf.write(node.name + '\t' + new_name + '\n')
            name_map[node.name] = new_name
            node.name = new_name
            node_count += 1
    outf.close()
    old_tree.write(format=5, outfile=new_tree)
    return [new_tree, name_map]


def create_itol_piechart(parsed_results, outdir):
    itol_piechart_file = open(outdir + "iTol_piechart.txt", 'w')
    itol_piechart_file.write("DATASET_PIECHART\nSEPARATOR TAB\nDATASET_LABEL\tpiechart\n"
                             "COLOR\t#ff0000\nFIELD_COLORS\t#ff4d4d\t#33bbff\t#70db70\t#ffcc66\n"
                             "FIELD_LABELS\tGain\tLoss\tDuplication\tReduction\nDATA\n")
    with open(parsed_results) as opr:
        for i, line in enumerate(opr):
            if i > 0:
                line = line.rstrip('\n')
                ls = line.split('\t')
                s = ls[0]
                sumg = int(ls[4]) + int(ls[6]) + int(ls[8]) + int(ls[10])
                if is_number(s):
                    s = 'node' + s
                itol_piechart_file.write('\t'.join([str(x) for x in [s, 1, sumg, ls[4], ls[6],
                                                                     ls[8], ls[10]]]) + '\n')
    itol_piechart_file.close()


def cleanUp(outdir):
    os.system('rm %s/phylo.tre %s/homolog_names.txt' % (outdir, outdir))
    os.system('gzip %s/input.nex %s/gene_flux_results.txt' % (outdir, outdir))
    os.system('tar -zcvf %sindividual_result_files.tar.gz %sindividual_result_files' % (outdir, outdir))
    os.system('rm -rf %s/individual_result_files/' % outdir)


def compute_gene_flux(tree, homolog_matrix, paup_exe, cost_matrix, num_chars, unrooted, just_nexus, output,
                      rename_innernodes):
    """ Main function which wraps everything together """
    try:
        assert (num_chars > 0 and num_chars <= 10)
    except:
        sys.stderr.write(
            'Error: Number of type of characters in homolog matrix must be between 1 and 10. Exiting now ...\n');
        sys.exit(
            1)
    try:
        assert (os.path.isfile(tree) and os.path.isfile(homolog_matrix))
    except:
        sys.stderr.write(
            'Error: Unable to locate input newick tree and/or homolog matrix file! Please check input and try again. Exiting now ...\n');
        sys.exit(1)

    # load in cost matrix unless default is used

    cmatrix = 'default'

    if not cost_matrix == 'default':
        num_cols = 0
        num_rows = 0
        cmatrix = []
        for i, line in enumerate(open(cost_matrix)):
            line = line.rstrip('\n')
            ls = line.split()
            cmatrix.append(ls)
            if i == 0:
                num_rows = len(ls)
            num_cols += 1
        # check that cost matrix is square-ish and that its height/width matches num_chars input
        try:
            assert (num_chars == num_cols and num_cols == num_rows)
        except:
            sys.stderr.write(
                'Error: Cost matrix is either not a square as expected or num_chars input doesn\'t quite match up with the cost matrix attributes. Please check input and try again. Exiting now ...\n');
            sys.exit(
                1)

    # load tree into ete object and get set of strains in tree.
    ete_tree = Tree(tree)
    strains_in_tree = set([])
    for leaf in ete_tree:
        strains_in_tree.add(str(leaf).strip('\n').lstrip('-'))

    # read in homolog cluster input and store set of strains in data file.
    homolog_strains = set([])
    homolog_array = []
    homolog_names = []
    max_count_flag = False
    with open(homolog_matrix) as ohm:
        for i, line in enumerate(ohm):
            line = line.rstrip('\n')
            ls = line.split('\t')
            if i == 0:
                homolog_strains = set(ls[1:])
                homolog_array.append(ls[1:])
            else:
                homolog_names.append(ls[0])
                homolog_counts = []
                for g in ls[1:]:
                    if int(g) > 9:
                        max_count_flag = True; homolog_counts.append('9')
                    else:
                        homolog_counts.append(g)
                homolog_array.append(homolog_counts)
    if max_count_flag:
        sys.stderr.write(
            "Warning: Homolog group copy count above 9 observed for at least one homolog group. Counts were capped at 9!\n")
    try:
        assert (len(homolog_strains.symmetric_difference(strains_in_tree)) == 0)
    except:
        sys.stderr.write(
            'Error: Strains in homolog matrix file do not match the strains in the newick file. Please check input and try again. Exiting now ...\n');
        sys.exit(
            1)
    num_strains = len(homolog_strains)

    output = os.path.abspath(output) + '/'
    if not os.path.isdir(output):
        os.system('mkdir ' + output)

    # step 0: get new name dictionary and create new tree
    new_tree, name_map = rename_strains(tree, output)

    # step 1: create input nexus file for paup
    nexus_input_file = create_nexus(name_map, new_tree, homolog_array, cmatrix, num_chars, unrooted, num_strains,
                                    output)

    if not just_nexus:
        # step 2: run paup
        paup_result_file = run_paup(paup_exe, nexus_input_file, output)

        # step 3: process PAUP file for tree inner-node labelling using Abigail Manson's program
        parent_child_file, homolog_names_file = create_necessary_files(paup_result_file, homolog_names, output)

        # step 4: summarize gene flux results
        parsed_results = create_results(paup_result_file, parent_child_file, homolog_names_file, name_map, output)

        # step 5: reformat phylogenetic tree to include inner node names
        reformat_tree(tree, parsed_results, output)

        # step 5.5: rename node if necessary
        if rename_innernodes:
            rename_nodes(parsed_results, rename_innernodes, output + 'phylogeny_with_innernode_names.nwk')

        # step 6: create pie chart input for iTol visualization
        create_itol_piechart(parsed_results, output)

        # step 7: clean up any "junk" files
        cleanUp(output)


if __name__ == '__main__':
    # Pull out the arguments.
    parser = argparse.ArgumentParser(description="""
	This program runs a gene flux analysis using PAUP (tested with v. 4.0b) provided a newick tree (should be rooted),
	a homolog group by sample copy-count matrix, and some recommended but optional parameters.
	""")

    parser.add_argument('-t', '--tree', help='Phylogenetic tree in Newick format.', required=True)
    parser.add_argument('-i', '--homolog_matrix',
                        help='Homolog matrix showing homolog copy-count across samples in phylogeny.',
                        required=True)
    parser.add_argument('-p', '--paup_exe', help='Path to the PAUP executable.', required=False, default='paup')
    parser.add_argument('-c', '--cost_matrix', type=str,
                        help='Input file with cost matrix. Default is a gain has a cost of 10 while a loss has a cost of 2.',
                        default="default", required=False)
    parser.add_argument('-n', '--num_chars', type=int,
                        help='Max copy count of any homolog group in homolog matrix, should correspond to the length/width of cost matrix. Max value allowed 9.',
                        default=2, required=False)
    parser.add_argument('-u', '--unrooted', action='store_true',
                        help='Specify that phylogenetic tree is unrooted. Default is off (rooted tree is assumed).',
                        required=False, default=False)
    parser.add_argument('-j', '--just_nexus', action='store_true',
                        help='Specify that only the input nexus file for PAUP should be created. This is recommended for large datasets with many strains (>= 300). In such a case you would run PAUP and the post processing individually, zooming out in your terminal to make sure the tree with inner node information at the end of the analysis is captured properly.')
    parser.add_argument('-o', '--output_dir', help='Output directory where results from analysis can be found.',
                        required=True)
    parser.add_argument('-r', '--rename_innernodes',
                        help='PAUP doesn\'t always retain the same IDs for innernodes with multiple runs.',
                        required=False, default=None)

    args = parser.parse_args()

    compute_gene_flux(args.tree, args.homolog_matrix, args.paup_exe, args.cost_matrix, args.num_chars, args.unrooted,
                      args.just_nexus,
                      args.output_dir, args.rename_innernodes)