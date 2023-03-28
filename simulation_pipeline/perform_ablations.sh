#!/bin/bash

# TODO: modify this to be the same as config.SIMULATED_DATA_PATH
ABLATION_PATH=/home/ema30/udn_data/data/simulated_patients/ablations


#################################################
# Run all combinations of phenotype ablations

python simulation_pipeline/perform_ablations.py --no_phen_noise
python simulation_pipeline/perform_ablations.py --no_phen_dropout
python simulation_pipeline/perform_ablations.py --no_phen_corruption
python simulation_pipeline/perform_ablations.py --no_gene_module_phen

python simulation_pipeline/perform_ablations.py --no_phen_noise --no_phen_dropout
python simulation_pipeline/perform_ablations.py --no_phen_corruption --no_phen_dropout
python simulation_pipeline/perform_ablations.py --no_phen_corruption --no_phen_noise
python simulation_pipeline/perform_ablations.py --no_gene_module_phen --no_phen_noise
python simulation_pipeline/perform_ablations.py --no_gene_module_phen --no_phen_dropout
python simulation_pipeline/perform_ablations.py --no_gene_module_phen --no_phen_corruption

python simulation_pipeline/perform_ablations.py --no_gene_module_phen --no_phen_corruption --no_phen_noise 
python simulation_pipeline/perform_ablations.py --no_gene_module_phen --no_phen_corruption --no_phen_dropout
python simulation_pipeline/perform_ablations.py --no_phen_corruption --no_phen_dropout --no_phen_noise 
python simulation_pipeline/perform_ablations.py --no_gene_module_phen --no_phen_dropout --no_phen_noise 

python simulation_pipeline/perform_ablations.py --no_gene_module_phen --no_phen_corruption --no_phen_noise --no_phen_dropout

# Format patients
for filename in ${ABLATION_PATH}/simulated_patients_no_*.jsonl; do
    if [[ $filename != *formatted.jsonl ]] 
    then
        fname=$(basename "$filename" .txt)
        echo $fname
        python simulate_patients/label_patients_with_novelty_categories.py --input ablations/$fname
    fi
done

#################################################
# Run all gene ablations
python simulation_pipeline/perform_ablations.py  --remove_module non_syndromic_phenotype
python simulation_pipeline/perform_ablations.py  --remove_module common_fp
python simulation_pipeline/perform_ablations.py  --remove_module tissue_distractor
python simulation_pipeline/perform_ablations.py  --remove_module pathogenic_pheno_irrel
python simulation_pipeline/perform_ablations.py  --remove_module insufficient_explanatory
python simulation_pipeline/perform_ablations.py  --remove_module phenotype_distractor

# Format patients
for filename in ${ABLATION_PATH}/simulated_patients_remove_*.jsonl; do
    if [[ $filename != *formatted.jsonl ]] 
    then
        fname=$(basename "$filename" .txt)
        echo $fname
        python simulate_patients/label_patients_with_novelty_categories.py --input ablations/$fname
    fi
done




#################################################
# Run phen + gene ablations

python perform_ablations.py --random_genes
python simulation_pipeline/label_patients_with_novelty_categories.py --input ablations/simulated_patients_rand_genes_no_phengenemod.jsonl

python perform_ablations.py --random_genes --no_phen_noise --no_phen_dropout --no_phen_corruption --no_gene_module_phen
python simulation_pipeline/label_patients_with_novelty_categories.py --input ablations/simulated_patients_rand_genes_no_phencorrupt_no_phendrop_no_phennoise_no_phengenemod.jsonl


