1. fastq.procesing.sh: initial data processing with fastq files as input
2. monocle.fxns.R: we initially used monocle3 for analyses and conducted ESR activity quantification using mononcle3 objects. These functions are related to this phase of analysis.
3. byxrm.sp.monocle.aucell.R: processing of cellranger outputs using the above functions as well as monocle3 functions. We only used the ESR activity values from this phase of analysis in the later mapping. We switched to seurat for analyses but kept these numbers.
4. cbsxyjm.salt.monocle.aucell.R: same as the above but with data for the CBSxYJM cross. 
#note - add the byxrm salt 
