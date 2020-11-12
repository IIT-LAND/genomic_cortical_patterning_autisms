# Analysis code for paper on genomic cortical patterning in autism early language outcome subtypes

This repository has all the code for analyses conducted for:

`Lombardo, M. V., Eyler, L., Pramparo, T., Gazestani, V. H., Hagler, D. J., Chen, C.-H., Dale, A. M., Seidlitz, J., Bethlehem, R. A. I., Bertelsen, N., Carter Barnes, C., Lopez, L., Lewis, N. E., Pierce, K., & Courchesne, E. Atypical genomic patterning of the cerebral cortex in autism with poor early language outcome. bioRxiv. doi:10.1101/2020.08.18.253443.`


## Requirements

All analysis is primarily conducted within **R** (R version 3.6.1 (2019-07-05) -- "Action of the Toes"), with the exception of the PLS analysis, which runs in **MATLAB** (R2019b).

**PLS analysis** uses the PLS toolbox in MATLAB that can be found here: https://www.rotman-baycrest.on.ca/index.php?section=84

Genetic cluster parcellations - **GCLUST** - can be found here:  https://github.com/ENIGMA-git/GCLUST.

Plotting of the Chen GCLUST parcellation is implemented with **ggseg**:  https://github.com/LCBC-UiO/ggsegChen.

**acPCA** can found here: https://github.com/linzx06/AC-PCA.

**Gene lists for enrichment analyses** can be found in the .Rdata files (e.g., enrichment_data2.Rdata, gandal_genelists.Rdata, etc).


## Analysis workflow and explanations

The *.Rmd files explained below, show the code, and the accompanying html files show the code along with all the output generated by the code (e.g., statistical results, plots, etc).

`_01_demographics.Rmd` - Runs analysis on the clinical and demographic variables.

`_01_covAdj_gexData.Rmd` - Does covariate adjustment of the gene expression data. Part of this analysis requires use of the CellCODE R library found here (https://github.com/mchikina/CellCODE).

`_01_covAdj_mriData_gentemp_ct.Rmd` and `_01_covAdj_mriData_gentemp_sa.Rmd` - Covariate adjustment of Freesurfer MRI features of surface area (SA) and cortical thickness (CT).

`_02_WGCNA.Rmd` - Runs the WGCNA analyses on gene expression data.

`_03_batch_PLS_all_gentemp.m` - Runs the PLS analysis in MATLAB.

`_03_plotPLSbrainBSR_gradient.Rmd` - Generates plots of the BSRs from the PLS results.

`_04_enrichmentAnalysis_sa_LV1_sexAdj_gentemp.Rmd` - Enrichment analyses on SA LV1 results.
`_04_enrichmentAnalysis_ct_LV1_sexAdj_gentemp.Rmd` - Enrichment analyses on CT LV1 results.
`_04_enrichmentAnalysis_ct_LV2_sexAdj_gentemp.Rmd` - Enrichment analyses on CT LV2 results.

`_05_FSanalysis_GenTemp.Rmd` - Runs analyses on just the MRI data.

`_06_psychENCODE.Rmd` - Analysis of the psychENCODE data.
