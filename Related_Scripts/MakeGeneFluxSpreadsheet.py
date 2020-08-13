import os
import sys
from collections import defaultdict

parent_child_file = sys.argv[1]
paup_result_file = sys.argv[2]
ortholog_name_file = sys.argv[3]
node_naming_file = sys.argv[4]
outdir = os.path.abspath(sys.argv[5]) + '/'

if not os.path.isdir(outdir):
    os.system('mkdir ' + outdir)

parent_mapping = {}
for line in open(parent_child_file):
    line = line.rstrip('\n')
    ls = [x.strip() for x in line.split('\t')]
    parent_mapping[ls[0]] = ls[1]

ortholog_names = []
for line in open(ortholog_name_file):
    ortholog_names.append(line.rstrip('\n'))

node_names = {}
for line in open(node_naming_file):
    line = line.rstrip('\n')
    ls = line.split('\t')
    node_names[ls[1]] = ls[0]

flag = False
sample_data = defaultdict(list)
for line in open(paup_result_file):
    line = line.rstrip('\n')
    chars_in_line = set(list(line))
    if len(chars_in_line) == 1 and list(chars_in_line)[0] == '-':
        flag = True
    elif line.strip() == '':
        flag = False
    elif flag:
        ls = line.split()
        sample_data[ls[0].strip()] += list(ls[1])

## 
## Currently assumes '?' characters, representing missing values are absent from the data.
##
print('\t'.join(
    ['node', 'parent', 'num_total_gene_clusters', 'num_unique_gene_clusters', 'num_gene_clusters_gained',
     'num_genes_gained',
     'num_gene_clusters_lost', 'num_genes_lost', 'num_gene_clusters_duplicated', 'num_genes_duplicated',
     'num_gene_clusters_reduced',
     'num_genes_reduced', 'all_gene_clusters', 'gene_clusters_gained', 'gene_clusters_lost', 'gene_clusters_duplicated',
     'gene_clusters_reduced']))

for samp in parent_mapping:
    parent = parent_mapping[samp].strip()
    samp = samp.strip()
    samp_name = samp
    if samp in node_names.keys():
        samp_name = node_names[samp]
    if not parent == 'none':
        samp_data = sample_data[samp]
        parent_data = sample_data[parent]

        allgcs = set([])
        gained = defaultdict(int)
        lost = defaultdict(int)
        duplicated = defaultdict(int)
        reduced = defaultdict(int)
        genes_gained = 0
        genes_lost = 0
        genes_duplicated = 0
        genes_reduced = 0

        out_gain = open(outdir + samp_name + '_gains.txt', 'w')
        out_lost = open(outdir + samp_name + '_losses.txt', 'w')
        out_reds = open(outdir + samp_name + '_reductions.txt', 'w')
        out_dups = open(outdir + samp_name + '_duplications.txt', 'w')
        out_all = open(outdir + samp_name + '_all_clusters.txt', 'w')

        num_tot_genes = 0
        num_uniq_clusts = 0
        for i in range(0, len(samp_data)):
            par_val = int(parent_data[i])
            sam_val = int(samp_data[i])
            num_tot_genes += sam_val
            if sam_val > 0:
                num_uniq_clusts += 1
                allgcs.add(ortholog_names[i] + ':' + str(sam_val))
                out_all.write(ortholog_names[i] + '\t' + str(sam_val) + '\n')
            if par_val == 0 and sam_val > 0:
                genes_gained += abs(sam_val-par_val)
                gained[ortholog_names[i]] = abs(par_val-sam_val)
                out_gain.write(ortholog_names[i] + '\n')
            elif par_val > 0 and sam_val == 0:
                genes_lost += abs(par_val-sam_val)
                lost[ortholog_names[i]] = abs(par_val-sam_val)
                out_lost.write(ortholog_names[i] + '\n')
            elif par_val > 0 and sam_val > 0 and par_val > sam_val:
                genes_reduced += abs(par_val-sam_val)
                reduced[ortholog_names[i]] = abs(par_val - sam_val)
                out_reds.write(ortholog_names[i] + '\t' + str(abs(par_val - sam_val)) + '\n')
            elif par_val > 0 and sam_val > 0 and sam_val > par_val:
                genes_duplicated += abs(par_val-sam_val)
                duplicated[ortholog_names[i]] = abs(par_val - sam_val)
                out_dups.write(ortholog_names[i] + '\t' + str(abs(par_val - sam_val)) + '\n')
        out_gain.close();
        out_lost.close();
        out_reds.close();
        out_dups.close();
        gains = ['']; losses = ['']; dups = ['']; reds = ['']
        if len(gained) > 0: gains = [s + ':' + str(gained[s]) for s in gained]
        if len(lost) > 0: losses = [s + ':' + str(lost[s]) for s in lost]
        if len(duplicated) > 0: dups = [s + ':' + str(duplicated[s]) for s in duplicated]
        if len(reduced) > 0: reds = [s + ':' + str(reduced[s]) for s in reduced]
        print('\t'.join([samp_name,
                         parent,
                         str(num_tot_genes),
                         str(num_uniq_clusts),
                         str(len(gained)),
                         str(genes_gained),
                         str(len(lost)),
                         str(genes_lost),
                         str(len(dups)),
                         str(genes_duplicated),
                         str(len(reds)),
                         str(genes_reduced),
                         ','.join(sorted(allgcs)),
                         ','.join(sorted(gains)),
                         ','.join(sorted(losses)),
                         ','.join(sorted(dups)),
                         ','.join(sorted(reds))]))
