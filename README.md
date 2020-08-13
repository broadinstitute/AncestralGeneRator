# AncestralGeneRator

This repo contains two programs for studying evolutionary gene flux across a phylogeny:

* **GeneFluxAnalysis.py** : Wrapper program to automate gene flux analysis using PAUP and post-process results.  
* **CorStrictor.py** : Program to determine the size of the number of single-copy-core genes at each inner-node of the phylogeny.

These wrappers/programs were developed by Rauf Salamzade and Abigail Manson in the Bacterial Genomics Group at the Broad Institute.

## Documentation

Check out the [wiki](https://github.com/broadinstitute/AncestralGeneRator/wiki) for more information on running the two programs!

![](https://github.com/broadinstitute/AncestralGeneRator/blob/master/Testing_Data/Gene_Flux_Upon_Phylogeny.png)

## Citations

If you use `GeneFluxAnalysis.py`, which uses PAUP (https://paup.phylosolutions.com/) for parsimony based ancestral state reconstruction, please cite:

> Swofford, D. L. 2003. PAUP*. Phylogenetic Analysis Using Parsimony (*and Other Methods).
> Version X. Sinauer Associates, Sunderland, Massachusetts.
> _Note: Because there are a number of beta and test versions of the program_
> _you should mention the specific version of PAUP* somewhere in the methods._

If you use the interactive tree of life (https://itol.embl.de/) for visualization, please cite:

> Letunic I and Bork P (2019) Nucleic Acids Res doi: 10.1093/nar/gkz239 Interactive Tree Of Life (iTOL) v4: recent updates and new developments
