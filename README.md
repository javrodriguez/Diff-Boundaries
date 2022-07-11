# Diff-Boundaries

This pipeline performs a pairwise differential analysis of TAD boundary insulation in two conditions. It was designed to use output files produced by the HiC-Bench pipeline. It uses the domains (TADs) identified by the Crane method (params.crane.ins_0500K.tcsh) and the 'ratio' insulation scores obtained for each replicate (matrix.ratio.k=001.tsv).

This tool generates a table with the insulation metrics for each boundary, a volcano-plot, a boxplot and three bed files with the boundary coordinates of the boundaries classified as 'increased', 'stable' and 'decreased' insulation-wise.

RUN:

Rscript diff.boundaries.r <condition1.name> <condition2.name> <condition1.domains.file <condition2.domains.file> <ratio.scores.file>

Notes:
- the ratio.scores.file (matrix.ratio.k=001.tsv) is generated in the 'boundary-scores-pca' step of HiC-Bench.
- the domains bed files (domains.k=001.bed) are generated in the 'domains' step of HiC-Bench.

RUN TEST:

git clone https://github.com/javrodriguez/Diff-Boundaries.git

cd Diff-Boundaries/

Rscript diff.boundaries.r PDX.D PDX.t1416 ./test_data/domains/PDX-D/domains.k\=001.bed ./test_data/domains/PDX-t1416/domains.k\=001.bed ./test_data/boundary-scores-pca/matrix.ratio.k\=001.tsv 
