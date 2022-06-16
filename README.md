# Diff-Boundaries

This pipeline performs a pairwise differential analysis of TAD boundary insulation in two conditions. It was designed to use output files produced by the HiC-Bench pipeline. It uses the domains (TADs) identified by the Crane method and the 'ratio' insulation scores obtained for each replicate.

This tool generates a master table with the insulation metrics for each boundary, a volcano-plot and a boxplot and three bed files with the boundary coordinates of the boundaries classified as 'increased', 'stable' and 'decreased' insulation-wise.

RUN:

Rscript diff.boundaries.r <condition1.name> <condition2.name> <condition1.domains.file <condition2.domains.file> <ratio.scores.file>

Notes:
- the ratio.scores.file (matrix.ratio.k=001.tsv) is generated in the 'boundary-scores-pca' step of HiC-Bench.
- the domains bed files (domains.k=001.bed) are generated in the 'domains' step of HiC-Bench.

RUN EXAMPLE:

Rscript diff.boundaries.r PDX.D PDX.t1416 __10a-domains/by_group/PDX-D/domains.k=001.bed __10a-domains/by_group/PDX-t1416/domains.k=001.bed files/matrix.ratio.k=001.tsv
