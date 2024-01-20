# Topographic-cognitive-neurobiological-profiling-of-interdependent-SC-FC
This repository provides data, core code and relevant toolboxes for data analysis in the article “Topographic, cognitive, and neurobiological profiling of the interdependent structural and functional connectome in the human brain”.
## Overview:
Contents include standalone software, source code, and demo data. All data and code necessary to reproduce our results have been made publicly available at https://github.com/wangxyue/Topographic-cognitive-neurobiological-profiling-of-interdependent-SC-FC.
## Installation guide:
All the scripts (*.m files in the code folder) and the toolboxes can be executed by adding the appropriate environment paths according to the corresponding code language. Use the "add path" function for Matlab and the "install.packages()" function for R. These procedures are not time-consuming.
## Code:
The following analyses were carried out using open-source packages:
1. Multilayer modular detection: The multilayer modular detection algorithm was obtained from an open Matlab code package GenLouvain (https://github.com/GenLouvain/GenLouvain) [1, 2].
2. Heritability analysis: The heritability analysis of multilayer modular variability was performed using the Accelerated Permutation Inference for the ACE model (APACE) method (https://github.com/nicholst/APACE) [3].
3. PLS analysis: The cognition-association analysis of multilayer modular variability was performed using an open Matlab toolbox myPLS (https://github.com/danizoeller/myPLS).
4. Allen Human Brain Atlas (AHBA) datasets: Regional microarray expression data from Allen Human Brain Atlas were preprocessed following the AHBA processing pipeline (https://github.com/BMHLab/AHBAprocessing) [4]. 
5. Gene Ontology Enrichment Analysis: Gene ontology enrichment analysis was performed using online tools including GOrilla (http://cbl-gorilla.cs.technion.ac.il/) [5] and REViGO (http://revigo.irb.hr).
6. Visualization: Brain map visualization was implemented using BrainNet Viewer Matlab toolbox (http://www.nitrc.org/projects/bnv) [6].



