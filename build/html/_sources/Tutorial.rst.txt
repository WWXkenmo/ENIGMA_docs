Tutorial
=======================

**ENIGMA** is a method for deconvoluting bulk RNA-seq matrix into cell type fractions and cell type-specific expression matrices. User could freely explore the cell heterogeneity within the bulk RNA-seq samples, study the differential expression gene, cell type-specific gene co-expression module or differentiation trajectory. In this tutorial, I applied ENIGMA on NSCLC bulk RNA-seq and corresponding FACS RNA-seq collected by Gentles et al to illustrate the main steps of ENIGMA analysis.

1. Getting Started 
*************************************

The first step of ENIGMA require bulk RNA-seq and reference matrix to construct ENIGMA object, the ENIGMA object contains parameters and datasets for deconvolution in the following steps

.. code-block:: R

    suppressPackageStartupMessages(library(magrittr))
    suppressPackageStartupMessages(library(scater))
    suppressPackageStartupMessages(library(Biobase))
    suppressPackageStartupMessages(library(ENIGMA))
    
    dataNSCLC <- readRDS("/load/Data/Path/dataNSCLC.rds")
    ref_sc <- readRDS("/load/Data/Path/ref.rds")
    
    #We used the third patients to generate reference
    ref_sc_sub <- ref_sc[,ref_sc$PatientID %in% "3" == TRUE]
    ref_sc_sub <- ref_sc_sub[,ref_sc_sub$CellFromTumor %in% "1"]
    ref_sc_sub <- ref_sc_sub[,ref_sc_sub$main_celltype %in% c("Alveolar","Epi") == FALSE]
    
    Bulk <- dataNSCLC[[5]]
    Tumor <- dataNSCLC[[1]]
    Immune <- dataNSCLC[[2]]
    Endothelial <- dataNSCLC[[3]]
    Fibroblast <- dataNSCLC[[4]]
    pheno <- dataNSCLC[[6]]
    
    # The pheno variable contain the label of each samples (LUSC vs LUAD)
    names(pheno) <- colnames(Tumor)
    
    # Create ENIGMA object
    egm = create_ENIGMA(bulk = Bulk, ref = exprs(ref_sc_sub), ref_type = "single_cell", meta_ref = pData(ref_sc_sub))

ENIGMA would automatically identify the type of reference matrix (FACS RNA-seq or scRNA-seq), and run the corresponding batch effect correction methods to generate reference matrix (S-mode and B-mode)

Batch effect correction
-----------------------

We used scRNA-seq to illustrate the S-mode batch effect correction

.. code-block:: R

    egm = batch_correct(egm, varname_cell_type = "main_celltype", n_pseudo_bulk=100)
    
    #check the reference profile
    head(egm@ref)

We used FACS RNA-seq to illustrate the B-mode batch effect correction, regarding the expression profile of T3 patients as the the reference, we noted that the FACS RNA-seq is not necessarily to be generated from the same cohort with the mixture(bulk RNA-seq) but could from the independent study.

.. code-block:: R

    B_ref = dataNSCLC[1:4] %>%
    lapply(function(x) x[,"T3"]) %>%
    do.call(cbind, .) %>%
    set_colnames( c("Tumor","Immune","Endothelial","Fibroblast") )
    B_ref = B_ref[rownames(egm@ref),]
    
    object_b_mode = create_ENIGMA(bulk = Bulk[rownames(B_ref),], ref = B_ref, ref_type = "bulk",meta_ref = as.matrix(colnames(B_ref)))
    object_b_mode = batch_correct(object_b_mode)

2. Deconvolution through ENIGMA
************************************

ENIGMA provides two type of regularized matrix completion methods to deconvolute bulk RNA-seq profile, L2-max norm and trace norm regularization. ENIGMA used sqrt() or log2() transformation to stablize each gene variance, make each gene expression distribution closer to the gaussian distribution. 

ENIGMA first used robust linear regression to estimation cell type fractions across each sample.

.. code-block:: R

    egm@bulk <- Bulk[rownames(tmp$main_celltype),]
    egm = get_cell_proportion(egm)
    plot_proportion(egm)

Next, we performed RNA-seq deconvolution based on L2-max norm regularized matrix completion.

.. code-block:: R
    
    egm = ENIGMA_L2_max_norm(egm, epsilon=0.001, alpha=0.8,
                            beta=10000,tao_k=0.01,max.iter=1000,verbose=TRUE)

We also could used trace norm regularized matrix completion to perform deconvolution

.. code-block:: R
    
    egm = ENIGMA_trace_norm(egm, epsilon=0.0001,alpha=0.8,beta=200,gamma = 1,verbose=TRUE,max.iter = 300)


