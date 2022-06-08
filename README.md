# CCGG-Coordinate-Mapping
Mapping Coordinates between pair of linear genomes extracted from CCGG


Anchors are conserved, unique and topologically sorted 45-mers in every linear genome in CCGG. Anchor sequences segment large contigs and identify homologies between linear genomes. We first perform pairwise alignment for parallel sequences between anchor pairs and then use the alignment for coordinate mapping. A intermediate string are constrcuted based on alignment and incorporate every mappable position of the two linear genomes.

Input:
Anchor information file for the linear genome
Genome FASTA file for the linear genome
