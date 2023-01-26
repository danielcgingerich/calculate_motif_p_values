# Quickly compute motif p values

Deterministic algorithm to compute exact p values for TFBS motifs in the genome.  
Computes p values using methodology of Tan et al (https://cran.r-project.org/web/packages/TFMPvalue/TFMPvalue.pdf), with code rewritten to be vectorized for efficiency.
Computationally outperforms TFMPvalue in calculation of exact p values for large number of motifs.
