#!/usr/bin/env python

import os
import sys
import argparse
from collections import defaultdict
from ete3 import Tree

def read_orthofile(homolog_matrix_file):
	try:
		col_to_sample = {}
		homolog_info = {}
		for i, line in enumerate(open(homolog_matrix_file)):
			line = line.rstrip('\n')
			ls = line.split('\t')
			if i == 0:
				for j, val in enumerate(ls[1:]):
					col_to_sample[j] = val
			else:
				homolog_info[ls[0]] = ls[1:]
		return [col_to_sample, homolog_info]
	except:
		sys.stderr.write("Problem parsing homolog matrix! Please check for formatting. Exiting now ...")
		raise RuntimeError


def determine_core_size(ortho_info, col_to_sample, strains):
	core_size = 0
	for o in ortho_info.items():
		flag = True
		for j, val in enumerate(o[1]):
			if col_to_sample[j] in strains and int(val) != 1: flag = False
		if flag: core_size += 1
	return core_size


def is_innernode(i):
	if i.startswith("node"):
		return True
	try:
		float(i); return True
	except:
		return False


def recursively_get_children(direct_children_map, curr_node):
	""" feeling like fibonacci """
	children = set([])
	direct_children = direct_children_map[curr_node]
	for child in direct_children:
		if is_innernode(child):
			children = children.union(recursively_get_children(direct_children_map, child))
		else:
			children.add(child)
	return children


def parse_phylogeny(tree):
	try:
		t = Tree(tree, format=1)
		direct_children = defaultdict(set)
		for node in t.traverse("postorder"):
			try:
				parent_name = node.up.name
				direct_children[parent_name].add(node.name)
			except: pass
		return direct_children
	except:
		sys.stderr.write("Problem parsing phylogeny! Please check input newick file. Exiting now ...")
		raise RuntimeError


def CorStrictor(tree, homolog_matrix, output):
	try:
		assert (os.path.isfile(tree) and os.path.isfile(homolog_matrix))
	except:
		sys.stderr.write("Either phylogeny or homolog matrix does not exist. Exiting now ..."); raise RuntimeError

	direct_children = parse_phylogeny(tree)
	col_to_sample, homolog_info = read_orthofile(homolog_matrix)

	output = os.path.abspath(output)

	try:
		assert (not os.path.isfile(output))
	except:
		sys.stderr.write("Output file already exists. Please remove/rename."); raise RuntimeError

	out = open(output, 'w')
	out.write(
		'DATASET_PIECHART\nSEPARATOR TAB\nDATASET_LABEL\tCorStrictor\nCOLOR\t#ff0000\nFIELD_COLORS\t#ff0000\nFIELD_LABELS\tCore_Genome_Size_Amongst_Child_Nodes\nDATA\n')
	for par in direct_children:
		all_children = recursively_get_children(direct_children, par)
		core_size = determine_core_size(homolog_info, col_to_sample, all_children)
		if par.strip():
			out.write('\t'.join([par, '1', str(core_size), str(core_size)]) + '\n')
	out.close()

if __name__ == '__main__':
	# Pull out the arguments.
	parser = argparse.ArgumentParser(
		description=""" This program creates an iTol visualization dataset depicting the narrowing of the core genome as one traverses up a phylogeny and how the strict single copy core for the full strain set ultimately forms. Pun based on second non-snake related definition of constrictor.""")
	parser.add_argument('-t', '--tree', help='Phylogenetic tree in Newick format. Inner nodes must be named!',
						required=True)
	parser.add_argument('-i', '--homolog_matrix',
						help='Homolog matrix showing homolog copy-count across samples in phylogeny.', required=True)
	parser.add_argument('-o', '--output', help="Output iTol dataset file.", required=True)
	args = parser.parse_args()

	CorStrictor(args.tree, args.homolog_matrix, args.output)