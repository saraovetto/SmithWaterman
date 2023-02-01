# SmithWaterman - pairwise local sequence alignment

I developed this Bioconductor-compliant package as a project for the Scientific Programming course delivered by the Politecnico of Milan as a part of the MSc in Bioinformatics for Computational Genomics.

The SmithWaterman package provides a function that performs pairwise local sequence alignment to find the optimal local alignment between two nucleotide sequences according to the Smith Waterman greedy algorithm.

### Description
This R package contains the "smith_waterman" function which takes as input 2 nucleic acid sequences and 3 scoring parameters (match, mismatch and gap penalty). As result, the function returns one optimal alignment between the two input sequences. The implementation works smoothly for sequences of length equals to ~500 nucelotides.

The package provides also two example sequences in FASTA file format to test the alignment. The sequences come from two different species of the Saccharomycetaceae family.

vignette link: https://htmlpreview.github.io/?https://github.com/saraovetto/SmithWaterman/blob/main/SmithWaterman.html
