import pandas as pd
import numpy as np
import argparse
import pickle5 as pickle
from pathlib import Path
from collections import defaultdict, OrderedDict, Counter
import matplotlib.pyplot as plt
plt.switch_backend('agg')
from scipy.stats import rankdata, wilcoxon
import sys

import config
from eval_plots import *
from simulate_patients.utils.util import read_simulated_patients
from evaluate_util import *

RESULTS_DIR = config.PROJECT_ROOT / 'gene_prioritization_results'


####################################################
# Plot either combined performance or by category

def evaluate_all_methods(filenames, output_base_fname, patients, category_type='broad', is_udn=False, plot_labels=None):
    print('\n----------- Evaluate All Methods ----------------')
    category = 'All'
    all_results_dict = OrderedDict()

    # Loop through rest of the results & add them
    for filename in filenames:
        print(f'evaluating.... {filename}')
        run_name= str(filename)
        gene_rankings = read_rankings(filename)
        all_results_dict[run_name] = evaluate(filename, gene_rankings, patients, category, category_type)

    fname = str(RESULTS_DIR / (output_base_fname + '_' + category + '.png'))
    plot_top_k_acc(all_results_dict, fname, is_udn, plot_labels)   
    for k,v in all_results_dict.items():
        v['name'] = str(k).replace('simulated/simulated_patients_formatted', '')
    all_results_df = pd.DataFrame([v for k,v in all_results_dict.items()])
    all_results_df = all_results_df.round(3).set_index('name')
    print(all_results_df)
    return all_results_dict

def evaluate_all_methods_all_categories(filenames, output_base_fname, category_dict, category_type, is_udn=False, plot_labels=None):
    print('\n----------- Evaluate All Methods, All Categories ----------------')

    # get results for all categories
    categories = list(category_dict.keys())
    all_categories_results_dict = {}
    for category in categories:
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
        all_results_df.to_csv(RESULTS_DIR / f'{category}_all_model_results.csv')
        all_categories_results_dict[category] = all_results_dict

        # calculate p-values
        for name_1, results_1 in all_results_dict.items():
            for name_2, results_2 in all_results_dict.items():
                if name_1 != name_2:
                    print(name_1, name_2)
                    assert len(results_1['ranks']) == len(results_2['ranks'])
                   
                    if results_1['ranks'] == results_2['ranks']:
                        print(f'Ranks for {name_1} and {name_2} are identical.')
                    else:
                        wilcoxon_test_greater = wilcoxon(results_1['ranks'], results_2['ranks'], alternative='greater')
                        wilcoxon_test_less = wilcoxon(results_1['ranks'], results_2['ranks'], alternative='less')
                        wilcoxon_test_two_sided = wilcoxon(results_1['ranks'], results_2['ranks'])

                        print('wilcoxon greater', wilcoxon_test_greater)
                        print('wilcoxon less', wilcoxon_test_less)
                        print('wilcoxon two-sided',wilcoxon_test_two_sided)


    # Plot results for all categories
    categories_fname = str(RESULTS_DIR /  (output_base_fname + '_all_categories.png'))
    grouped_plot_all_categories(categories_fname, all_categories_results_dict, plot_labels=plot_labels)

    return all_categories_results_dict


####################################################
# Main

def main():
    plt.style.use('ggplot')

    parser = argparse.ArgumentParser(description='Evaluate performance on simulated & UDN patients.')
    parser.add_argument('--patients', type=str, default='simulated_patients_formatted.jsonl')
    parser.add_argument('--results', type=str, default=None, help='Path to pkl file with additional gene prioritization results to evaluate, if any')
    parser.add_argument('--results_name', type=str, default=None, help='Name of Gene Prioritization model. This will appear in the plot as the x-axis label.')
    parser.add_argument('--reproduce_paper_results', action='store_true', help='Specify this option to plot results of the gene prioritization algorithms already run on simulated patients.')
    args = parser.parse_args()
    # NOTE: if reproduce_paper_results is True & an additional results file is passed in, then the user is responsible for making sure that they are all run on the same patient cohort.


    # Read in simulated patients
    patients = read_simulated_patients(config.SIMULATED_DATA_PATH / args.patients)
    patient_categories_dict, patient_broad_categories_dict = categorize_patients(patients) 
    sim_base_fname = args.patients.split('.jsonl')[0]

    if args.reproduce_paper_results:
        filenames = [
            RESULTS_DIR / 'phrank' / (sim_base_fname + '_phrank_rankgenes_directly=True.pkl'), 
            RESULTS_DIR / 'phrank' / (sim_base_fname + '_phrank_rankgenes_directly=False.pkl'),
            RESULTS_DIR / 'phenomizer' / (sim_base_fname + '.pkl'),
            RESULTS_DIR / 'phenolyzer' / (sim_base_fname + '.pkl'),
            RESULTS_DIR / 'random' / (sim_base_fname + '_random_baseline.pkl')
        ] 
        plot_labels = ['Phrank \nPatient-Gene', 'Phrank \nPatient-Disease', 'Phenomizer', 'Phenolyzer', 'Random']
    else:
        filenames = []
        plot_labels = []

    if args.results is not None and args.results_name is not None:
        filenames =  filenames + [Path(args.results)]
        print('filenames', filenames)
        plot_labels = plot_labels + [args.results_name]

    # Evaluate all patients together 
    sim_all_results_dict = evaluate_all_methods(filenames, sim_base_fname, patients, is_udn=False , plot_labels=plot_labels) 
    
    # Evaluate each category separately
    sim_all_results_dict_all_cat = evaluate_all_methods_all_categories(filenames, sim_base_fname, patient_broad_categories_dict, category_type='broad', is_udn=False, plot_labels=plot_labels) 
    
        
if __name__ == "__main__":
    main()
