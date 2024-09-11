# Coloc_Epi_Sleep
Colocalization analysis for epilepsy and sleep-related trait, as part of Master of Binformatics project.  

Ones need to be aware of the assumptions of Coloc before running analysis. Coloc assume 1 single shared causal variant in a region between two traits.  

Depend on the sensitivity analysis of coloc, ones may need to run Susie, which relax the assumptions of coloc and allow the existence of more than one causal variant per region.   

Also, the shared SNPs and the Allele alignment of those SNPs might cause a lot of issue when running Susie, when generating the LD matrix.

**Useful source: **  
GitHub Coloc: https://github.com/chr1swallace/coloc?tab=readme-ov-file   
GitHub SuSiE: https://chr1swallace.github.io/coloc/articles/a06_SuSiE.html   
R Tutorial: https://cran.r-project.org/web/packages/coloc/vignettes/a06_SuSiE.html 
