The input files
di_file = "Sharma_Supplemantary_file_dinucleotide.csv"
tri_file = "Sharma_Supplemantary_file_trinuc_revised_20250919.csv"
tetra_file = "Sharma_Supplemantary_file_tetranucleotide.csv"
include 30 DNA mechanical index parameters, originally retrieved from four papers:
Brukner, I., Sanchez, R., Suck, D., & Pongor, S. (1995). Sequence‐dependent bending propensity of DNA as revealed by DNase I: parameters for trinucleotides. The EMBO journal, 14(8), 1812-1818.
Sharma, D., Sharma, K., Mishra, A., Siwach, P., Mittal, A., & Jayaram, B. (2023). Molecular dynamics simulation-based trinucleotide and tetranucleotide level structural and energy characterization of the functional units of genomic DNA. Physical Chemistry Chemical Physics, 25(10), 7323-7337.
Singhal, P., Jayaram, B., Dixit, S. B., & Beveridge, D. L. (2008). Prokaryotic gene finding based on physicochemical characteristics of codons calculated from molecular dynamics simulations. Biophysical journal, 94(11), 4173-4183.
Tolstorukov, M. Y., Ivanov, V. I., Malenkov, G. G., Jernigan, R. L., & Zhurkin, V. B. (2001). Sequence-dependent B↔ A transition in DNA evaluated with dimeric and trimeric scales. Biophysical Journal, 81(6), 3409-3421.

Rescaling method of these values are detailed in: Kondo, Y., & Budd, A. (under review). Lifespan is encoded in DNA structures and sequence motifs: comparative genomics in protostomes. 


The input file "pearson_correlation_results_DNA_mechanics_and_lifespan_global_FDR.csv" is the output from "04_01_scatter_ln_lifespan_and_DNA_mechanics.r".


The input file "combined_model_weights.csv" (in 05_04_combine_average_weight_of_each_tetra_for_each_regions.py) is the output from "01_00_04_elasticnet_glmnetUtils_Train_Test_representative_split_alpha1_1000folds.R".

