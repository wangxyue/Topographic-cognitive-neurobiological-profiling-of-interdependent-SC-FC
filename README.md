# Topographic-cognitive-neurobiological-profiling-of-interdependent-SC-FC
This repository provides data, core code and relevant toolboxes for data analysis in the article “Topographic, cognitive, and neurobiological profiling of the interdependent structural and functional connectome in the human brain”.

# Overview:
Contents include standalone software, source code, and demo data. All data and code necessary to reproduce our results have been made publicly available at https://github.com/wangxyue/Topographic-cognitive-neurobiological-profiling-of-interdependent-SC-FC.

# Installation guide:
All the scripts (*.m files in the code folder) and the toolboxes can be executed by adding the appropriate environment paths according to the corresponding code language. Use the "add path" function for Matlab and the "install.packages()" function for R. These procedures are not time-consuming.

# Code:
The following analyses were carried out using open-source packages:
1. Multilayer modular detection: The multilayer modular detection algorithm was obtained from an open Matlab code package GenLouvain (https://github.com/GenLouvain/GenLouvain) [1, 2].
2. Heritability analysis: The heritability analysis of multilayer modular variability was performed using the Accelerated Permutation Inference for the ACE model (APACE) method (https://github.com/nicholst/APACE) [3].
3. PLS analysis: The cognition-association analysis of multilayer modular variability was performed using an open Matlab toolbox myPLS (https://github.com/danizoeller/myPLS).
4. Allen Human Brain Atlas (AHBA) datasets: Regional microarray expression data from Allen Human Brain Atlas were preprocessed following the AHBA processing pipeline (https://github.com/BMHLab/AHBAprocessing) [4]. 
5. Gene Ontology Enrichment Analysis: Gene ontology enrichment analysis was performed using online tools including GOrilla (http://cbl-gorilla.cs.technion.ac.il/) [5] and REViGO (http://revigo.irb.hr).
6. Visualization: Brain map visualization was implemented using BrainNet Viewer Matlab toolbox (http://www.nitrc.org/projects/bnv) [6].

The code folder contains all the code used to run the analyses and generate the figures. This folder contains the following files:
1. s1_multilayer_community：Perform multilayer modular detection and calculate the regional multilayer modular variability (a_multilayer_mod_var.m). Calculate the multilayer modular variability values for four hierarchical systems (b_hiesys_MV_spin.m). Plot the correlation between multilayer modular variability and cortical gradient (scatter_MVcorgrad.R) and cortical expansion (scatter_MVexpans.R).
2. s2_reliability_analysis: Analyzing the test-retest reliability of multilayer modular variability (a_reliability_analysis.m). Calculate the test-retest reliability of four hierarchical systems (b_hiesys_ICC_spin.m). Plot the intra-individual and inter-individual similarity of multilayer modular variability (plot_MV_reliability.R).
3. s3_heritability_analysis: Calculate the heritability of multilayer modular variability and compare the similarity of multilayer modular variability of three groups (a_WKParallelOutline.m). Calculate the heritability of four hierarchical systems (b_hiesys_h2_spin.m). Plot the similarity of multilayer modular variability of three groups (plot_MV_h2simi.R).
4. s4_cognition_association: Neurocognitive flexibility association analysis (a_flexibility_association.m). Plot the multilayer modular variability values for the four groups of nodes partitioned by neurocognitive flexibility (plot_4nodes_flexibility.R). Analyzing the cognition-association of the multilayer modular variability of low-order (b_Cognition_PLS_priuni.m) and high-order cortex (c_Cognition_PLS_hetepara.m). Calculate the loadings of variables with significant contributions (d_select_salience_loadings.m).
5. s5_receptor_association: Calculate the correlation between multilayer modular variability and neurotransmitter receptor and transport distribution (a_corr_MV_receptor.m), and plot the correlation (corr_MV_receptor.R). Perform the Elastic net regression analysis (b_Elasticnet_regression.m) and plot the beta values of significant features (c_receptor_beta.m).
6. s6_gene_association: Perform the transcriptome-association analysis (a_gene_PLS_analysis.m) and calculate the weights of genes (b_my_PLS_bootstrap.m).
7. mod_visualization: Visualize the brain maps on the inflated cortical 32K surface [7] using BrainNet Viewer.
8. Spintest: Spatial Permutation Testing.

Custom codes and toolboxes were tested on a 64-bit Windows 10 PC. Data analysis and visualization were performed using Matlab 2019b and R 4.1.3. We thank the authors and developers for providing these data analysis tools.

# Data:
The following data were obtained from publicly available datasets:
1. The HCP dataset, including structural MRI, functional MRI, and diffusion-weighted MRI, is available in the HCP ConnectomeDB (https://db.humanconnectome.org/). 
2. The neurocognitive flexibility data [8] is publicly available at https://surfer.nmr.mgh.harvard.edu/fswiki/BrainmapOntology_Yeo2015. 
3. The neurotransmitter receptor and transport density distribution data [9] are publicly available at https://github.com/netneurolab/hansen_receptors. 
4. The AHBA dataset is publicly available at https://human.brain-map.org/static/download.

The data folder contains all the data used to run the analyses. This folder contains the following files:
1. s1_multilayer_community_data: it associated with the s1_multilayer_community analysis. It includes the results of multilayer modular detection (sub_multi_module.mat, sub_modularity.mat), the individual-level (all_sub_MV.mat) and group-level multilayer modular variability (cross_sub_MV.mat). The cortical gradient data (gradient_glasser360.mat), cortical expansion data (CorticalExpansion_R_glasser360.mat) and surrogated multilayer modular variability maps (Perm_360_MVResults.mat).
2. s2_reliability_data: it associated with the s2_reliability_analysis. It includes the individual-level multilayer modular variability (sub42_MV_test.mat, sub42_MV_retest.mat), the nodal intra-class correlation values (ICC.mat) and the surrogated ICC maps of multilayer modular variability (Perm_360_ICCResults.mat).
3. s3_heritability_data: it associated with the s3_heritability_analysis. It includes the individual-level multilayer modular variability (all_sub_MV1009.mat) and the surrogated heritability maps of multilayer modular variability (Perm_360_h2Results.mat).
4. s4_cognition_data: it associated with the s4_cognition_association analysis. It includes the node’s neurocognitive flexibility (pro_node.mat), the multilayer modular variability of the low-order (priuni_MV.mat) and high-order (hetepara_MV.mat) cortex, loadings of cognition terms (final_sigbeh_loadings.mat, final_sigbeh_name.mat) and brain regions (final_sigMV_loadings.mat) with significant contributions.
5. s5_receptor_data: it associated with the s5_receptor_association analysis. It includes the neurotransmitter receptors and transports data from nine neurotransmitter systems, the data of elastic net regression analysis (XYdata.csv), the correlation r (MV_cor_receptors.csv), the features with significant contributions (final_beta.mat, final_name.mat).
6. s6_gene_data: it associated with the s6_gene_association analysis. It includes the gene expression data (parcelExpression.mat, probeInformation.mat), PLS analysis results (genePLS_stats.csv, MV_corr_gene_permPLS.mat), gene weights (PLS_geneWeights_ascend.csv, PLS_geneWeights_descend.csv), surrogated multilayer modular variability (permMV_res.mat).
7. s7_validation: it associated with the sensitivity and robustness analysis. It includes the group-level multilayer modular variability (cross_sub917_MV.mat, FC005_cross_sub_MV.mat, FC02_cross_sub_MV.mat); the sub37_testMV.mat, sub37_retestMV.mat and ICC.mat for reliability analysis; the sub914_MV.mat for heritability analysis; the sub917hetepara_MV.mat and sub917priuni_MV.mat for PLS cognitive-association analysis; the PLS_geneWeights_ascend.csv and PLS_geneWeights_descend.csv for enrichment analysis;


# References:
1. Mucha, P.J., Richardson, T., Macon, K., Porter, M.A. & Onnela, J.P. Community structure in time-dependent, multiscale, and multiplex networks. Science (New York, N.Y.) 328, 876-878 (2010).
2. Jutla, I., Jeub, L. & Mucha, P. A generalized Louvain method for community detection implemented in MATLAB (computer program). Available at netwiki. amath. unc. edu. Accessed September 1, 2014 (2011).
3. Chen, X., et al. Accelerated estimation and permutation inference for ACE modeling. Human brain mapping 40, 3488-3507 (2019).
4. Arnatkeviciute, A., Fulcher, B.D. & Fornito, A. A practical guide to linking brain-wide gene expression and neuroimaging data. NeuroImage 189, 353-367 (2019).
5. Eden, E., Navon, R., Steinfeld, I., Lipson, D. & Yakhini, Z. GOrilla: a tool for discovery and visualization of enriched GO terms in ranked gene lists. BMC bioinformatics 10, 48 (2009).
6. Xia, M., Wang, J., and He, Y. (2013). BrainNet Viewer: a network visualization tool for human brain connectomics. PLoS One 8, e68910. 10.1371/journal.pone.0068910.
7. Glasser, M.F., Coalson, T.S., Robinson, E.C., Hacker, C.D., Harwell, J., Yacoub, E., Ugurbil, K., Andersson, J., Beckmann, C.F., Jenkinson, M., et al. (2016). A multi-modal parcellation of human cerebral cortex. Nature 536, 171-178. 10.1038/nature18933.
8. Yeo, B.T., et al. Functional Specialization and Flexibility in Human Association Cortex. Cerebral cortex (New York, N.Y. : 1991) 25, 3654-3672 (2015).
9. Hansen, J.Y., et al. Mapping neurotransmitter systems to the structural and functional organization of the human neocortex. Nature neuroscience 25, 1569-1581 (2022).

