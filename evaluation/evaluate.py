import pandas as pd
#import pickle
import numpy as np
import argparse
import pickle5 as pickle
from pathlib import Path
from collections import defaultdict, OrderedDict, Counter
from util import read_simulated_patients, read_udn_patients
import matplotlib.pyplot as plt
plt.switch_backend('agg')
from scipy.stats import rankdata

import sys
sys.path.insert(0, '../') # add config to path
import config
from eval_plots import *

PHENOLYZER_DIR = config.PROJECT_ROOT / 'phenolyzer'
PHENOMIZER_DIR = config.PROJECT_ROOT / 'phenomizer'
PHRANK_DIR = config.PROJECT_ROOT / 'phrank'
OUTPUT_DIR = config.PROJECT_ROOT / 'gene_rank_results'


####################################################
# Plot either combined performance or by category

def evaluate_all_methods(filenames, output_base_fname, patients, category_type='broad', is_udn=False, plot_labels=None):
    print('\n----------- evaluate all methods ----------------')
    category = 'All'
    all_results_dict = OrderedDict()

    # Loop through rest of the results & add them
    for filename in filenames:
        print(f'evaluating.... {filename}')
        run_name= str(filename)
        gene_rankings = read_rankings(filename)
        all_results_dict[run_name] = evaluate(filename, gene_rankings, patients, category, category_type)

    fname = str(config.PROJECT_ROOT / 'gene_rank_results' / (output_base_fname + '_' + category + '.png'))
    print(fname)
    plot(all_results_dict, pretty_print_category(category), fname, is_udn, None)  #plot_labels 
    print(all_results_dict)
    for k,v in all_results_dict.items():
        v['name'] = str(k).replace('simulated/simulated_patients_formatted', '')
    all_results_df = pd.DataFrame([v for k,v in all_results_dict.items()])
    all_results_df = all_results_df.round(3).set_index('name')
    print(all_results_df)
    return all_results_dict

def evaluate_all_methods_all_categories(filenames, output_base_fname, category_dict, category_type, is_udn=False, plot_labels=None):
    print('\n----------- evaluate all methods, all categories ----------------')

    # get results for all categories
    categories = list(category_dict.keys())
    all_categories_results_dict = {}
    for category in categories:
        print('category: ', category)
        patients = category_dict[category]
        print(f'There are {len(patients)} patients with {category} category')
        all_results_dict = OrderedDict()

        for filename in filenames:
            print(f'evaluating.... {filename}')
            run_name= str(filename)


            gene_rankings = read_rankings(filename)
            all_results_dict[run_name] = evaluate(filename, gene_rankings, patients, category, category_type)
        for k,v in all_results_dict.items():
            v['name'] = str(k).replace('simulated/simulated_patients_formatted', '')
        all_results_df = pd.DataFrame([v for k,v in all_results_dict.items()])
        all_results_df = all_results_df.round(3).set_index('name')
        print(all_results_df)
        all_results_df.to_csv(OUTPUT_DIR / f'{category}_all_model_results.csv')

        all_categories_results_dict[category] = all_results_dict


    # Plot results for all categories
    categories_fname = str(config.PROJECT_ROOT / 'gene_rank_results' /  (output_base_fname + '_all_categories.png'))

    grouped_plot_all_categories(categories_fname, all_categories_results_dict)


    return all_categories_results_dict


####################################################
# Main

def main():
    plt.style.use('ggplot')

    parser = argparse.ArgumentParser(description='Evaluate performance on simulated & UDN patients.')
    parser.add_argument('--patients', type=str, default='simulated_patients_formatted.jsonl')
    parser.add_argument('--results', type=str, default=None, help='Path to pkl file with additional gene prioritization results to evaluate, if any')
    args = parser.parse_args()


    # Read in simulated patients
    patients = read_simulated_patients(config.SIMULATED_DATA_PATH / args.patients)
    patient_categories_dict, patient_broad_categories_dict = categorize_patients(patients) 
    sim_base_fname = args.patients.split('.jsonl')[0]

    
    if 'ablations/' in sim_base_fname:
        results_dir = PHRANK_DIR / 'ablation_results' 
        filenames = [
            results_dir / (sim_base_fname.split('ablations/')[1] + '_phrank_rankgenes_directly=True.pkl'), 
            results_dir / (sim_base_fname.split('ablations/')[1] + '_phrank_rankgenes_directly=False.pkl'),
            ] 

    else:
        results_dir = PHRANK_DIR / 'simulated'
        filenames = [
            results_dir / (sim_base_fname + '_phrank_rankgenes_directly=True.pkl'), 
            results_dir / (sim_base_fname + '_phrank_rankgenes_directly=False.pkl'),
            PHENOMIZER_DIR / 'results' / (sim_base_fname + '_1000.pkl'),
            PHENOLYZER_DIR / 'results' / (sim_base_fname + '_rank_dict.pkl'),
            config.PROJECT_ROOT / 'random' / (sim_base_fname + '_random_baseline.pkl')
        ] 

    if args.results is not None:
        filenames = filenames + [args.results]

    # Evaluate all patients together 
    sim_all_results_dict = evaluate_all_methods(filenames, sim_base_fname, patients, is_udn=False ) 
    
    # Evaluate each category separately
    sim_all_results_dict_all_cat = evaluate_all_methods_all_categories(filenames, sim_base_fname, patient_broad_categories_dict, category_type='broad', is_udn=False) 
    
        
if __name__ == "__main__":
    main()
