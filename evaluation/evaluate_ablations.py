import pandas as pd
#import pickle
import numpy as np
import argparse
import pickle5 as pickle
from pathlib import Path
from collections import defaultdict, OrderedDict, Counter
import matplotlib.pyplot as plt
plt.switch_backend('agg')
from scipy.stats import rankdata

import sys
import config
from eval_plots import *
from simulate_patients.utils.util import read_simulated_patients
from evaluate_util import *

RESULTS_DIR = config.PROJECT_ROOT / 'gene_prioritization_results'


'''
Evaluate Gene Prioritization Model on Simulated Patients created via Pipeline Ablation

NOTE: running this script requires that the gene prioritization algorithm has been run on ALL combinations of ablations.
If you wish to only run on some of the ablations, either comment them out in the evaluate_ablations or use the evaluate.py script
to evaluate performance on each patient cohort.  
'''

####################################################
# Plot Ablations


def plot_phrank_ablations(results_stem, results_path, ablations, patients, fname, xlabel, category_type='broad'):

    all_results_dict = OrderedDict()
    for i, (title, ablation) in enumerate(ablations):
        filename = Path(results_path) / (ablation + results_stem)
        gene_rankings = read_rankings(filename)
        all_results_dict[title]  = evaluate(filename, gene_rankings, patients, 'All', category_type)
    
    for k,v in all_results_dict.items():
        v['name'] = str(k)
    all_results_df = pd.DataFrame([v for k,v in all_results_dict.items()])
    all_results_df = all_results_df.round(3).set_index('name')
    all_results_df.to_csv(RESULTS_DIR / f'all_model_ablation_results.csv')
    
    plot_ablations(all_results_dict, fname, xlabel)

def evaluate_ablations(results_filename, results_path, patients):
    print('Evaluating Ablations')

    ablations = ablations = [
    ('Complete Pipeline', 'simulated_patients_equal_probs_revised_formatted'), 
    ('No Phenotype Modules', 'simulated_patients_no_phencorrupt_no_phendrop_no_phennoise_ablations_formatted'),
    ('No Gene Modules', 'simulated_patients_rand_genes_no_phengenemod_ablations_formatted'),
    ('No Phenotype or Gene Modules', 'simulated_patients_rand_genes_no_phencorrupt_no_phendrop_no_phennoise_no_phengenemod_ablations_formatted'),
    ]
    fname = str(RESULTS_DIR /  'overall_ablation_results.png')
    plot_phrank_ablations(results_filename, results_path, ablations, patients, fname, xlabel='Ablation')

    gene_replaced_ablations = [
    ('No Gene Modules', 'simulated_patients_equal_probs_revised_formatted'), 
    ('Similarly Expressed Genes', 'simulated_patients_replace_genes_with_rand_remove_tissue_distractor_ablations_formatted'), 
    ('Insufficiently Explanatory Genes', 'simulated_patients_replace_genes_with_rand_remove_insufficient_explanatory_ablations_formatted'), 
    ('Common False Positive Genes', 'simulated_patients_replace_genes_with_rand_remove_common_fp_ablations_formatted'), 
    ('Phenotypically-distinct Disease Genes', 'simulated_patients_replace_genes_with_rand_remove_pathogenic_pheno_irrel_ablations_formatted'), 
    ('Genes Associated with Incidental Phenotypes', 'simulated_patients_replace_genes_with_rand_remove_non_syndromic_phenotype_ablations_formatted'), 
    ('Phenotypically-similar Disease Genes', 'simulated_patients_replace_genes_with_rand_remove_phenotype_and_universal_distractor_ablations_formatted'), 
    ('All Gene Modules', 'simulated_patients_rand_genes_no_phengenemod_ablations_formatted'), 
     ]
    fname = str(RESULTS_DIR /  'gene_module_ablation_results.png')
    plot_phrank_ablations(results_filename, results_path, gene_replaced_ablations, patients, fname, xlabel='Gene Module Removed')


    ablations = [
    ('None', 'simulated_patients_equal_probs_revised_formatted'), 
    ('Corrupt', 'simulated_patients_no_phencorrupt_ablations_formatted'),
    ('Dropout', 'simulated_patients_no_phendrop_ablations_formatted'),
    ('Noise', 'simulated_patients_no_phennoise_ablations_formatted'),
    ('Gene Module (GM) Phenotypes', 'simulated_patients_no_phengenemod_ablations_formatted'),
    ('Dropout & Corrupt', 'simulated_patients_no_phencorrupt_no_phendrop_ablations_formatted'),
    ('Gene Module Phenotypes & Corrupt', 'simulated_patients_no_phencorrupt_no_phengenemod_ablations_formatted'),
    ('Noise & Corrupt', 'simulated_patients_no_phencorrupt_no_phennoise_ablations_formatted'),
    ('Noise & Dropout', 'simulated_patients_no_phendrop_no_phennoise_ablations_formatted'),
    ('Gene Module Phenotypes & Dropout', 'simulated_patients_no_phendrop_no_phengenemod_ablations_formatted'),
    ('Gene Module Phenotypes & Noise', 'simulated_patients_no_phennoise_no_phengenemod_ablations_formatted'),
    ('Noise, Dropout, & Corrupt', 'simulated_patients_no_phencorrupt_no_phendrop_no_phennoise_ablations_formatted'),
    ('GM Phenotypes, Dropout, & Corrupt', 'simulated_patients_no_phencorrupt_no_phendrop_no_phengenemod_ablations_formatted'),
    ('GM Phenotypes, Noise, & Corrupt', 'simulated_patients_no_phencorrupt_no_phennoise_no_phengenemod_ablations_formatted'),
    ('GM Phenotypes, Noise, & Dropout', 'simulated_patients_no_phendrop_no_phennoise_no_phengenemod_ablations_formatted'),
    ('GM Phenotypes, Noise, Dropout, & Corrupt', 'simulated_patients_no_phencorrupt_no_phendrop_no_phennoise_no_phengenemod_ablations_formatted')
    ]
    fname = str(RESULTS_DIR/  'phenotype_ablation_results.png')
    plot_phrank_ablations(results_filename, results_path, ablations, patients, fname, xlabel='Removed Simulation Step')


####################################################
# Main

def main():
    parser = argparse.ArgumentParser(description='Evaluate performance on simulated & UDN patients.')
    parser.add_argument('--patients', type=str, default=f'{config.SIMULATED_DATA_PATH}/ablations/simulated_patients_equal_probs_revised_formatted.jsonl') 
    parser.add_argument('--results_suffix', type=str, default='_phrank_rankgenes_directly=False.pkl', help='Suffix of pkl file with gene prioritization results to evaluate')
    parser.add_argument('--path_to_results', type=str, default=str(RESULTS_DIR/'phrank_ablation'), help='Path containing results of Gene Prioritization models on ablations')

    args = parser.parse_args()

    # Read in simulated patients
    patients = read_simulated_patients(args.patients)

    # Evaluate ablations
    evaluate_ablations(args.results_suffix, args.path_to_results, patients) 

   
if __name__ == "__main__":
    main()
