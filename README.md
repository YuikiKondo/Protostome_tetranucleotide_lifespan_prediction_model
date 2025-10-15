# Protostome_tetranucleotide_lifespan_prediction_model
All code used for lifespan prediction from genomic tetranucleotide frequencies and for creating tables and figures in the paper: Kondo, Y., &amp; Budd, A. Lifespan prediction in protostome invertebrates from tetranucleotide frequency.

🔍 Overview of Analyses

1. Dataset preparation: the code used for these processes are summarized in another repository: YuikiKondo/protostome_lifespan_DNA_mechanics
・Download genome files from NCBI
・Calculate di-, tri- and tetranucleotides counts and frequencies for whole-genomes
・Calculate di-, tri- and tetranucleotides counts and frequencies for BUSCO genes using Miniprot as a gene prediction tool
・Gene prediction and extraction of BUSCO regions using Miniprot (or Augustus), followed by motif frequency calculations
・Calculate DNA mechanics index from di-, tri- and tetranucleotides frequencies

2. Lifespan prediction model construction using tetranucleotide frequencies
・01_lifespan_prediction_final_model

3. Comparison of lifespan prediction models using different predictors and genomic regions
・02_model_comparison

4. Investigaton of how much each predictor contribute to lifespan prediction in the complete TMPG model (T: 256 tetranucleotide frequencies, M: 30 DNA mechanical indices, P: 12 phylogenetic eigenvectors, and G: genome size)
・03_variance_decomposition

5. Assessment of how DNA mechanics information is incorporated in the final tetranucleotide-based model 
・04_mechanical_index_and_lifespan_correlation
・05_mechanical_index_parameter_and_weight_correlation

6. Cross-clade comparison of lifespan–tetranucleotide-bias association between two major protostome clades 
・06_cross_clade_models_by_Ecdysozoa_and_Spiralia


💻 Script Types and Execution
All analyses are conducted using Python (.py) and R (.r, .Rmd) scripts.
Some scripts are executed on high-performance computing (HPC) front-end servers as part of continuous workflows originating from dataset preparation in another repository (protostome_lifespan_DNA_mechanics), while others are run locally on a laptop.
Please note that the directory paths differ accordingly between these environments.

⚠️ Note:
Please adjust directory paths, SLURM resource specifications (cores, memory, time), and environment settings according to your system.

📂 Included Files
Each analysis directory includes　input files we used for our analyses and several representative output files for reference.

Author: Yuiki Kondo
