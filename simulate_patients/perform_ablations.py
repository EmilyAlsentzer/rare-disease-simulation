import json
import pandas as pd
import numpy as np
import jsonlines
from tqdm import tqdm
import random
from collections import defaultdict, Counter
import argparse 

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

import sys

import config
from create_dataset_for_baselines import read_jsonl

'''
We want to be able to measure the contribution of each component of the simulation pipeline. 
To do this, we will measure performance of Phrank on different versions of the simulated cohort
with a single module removed at a time. To minimize randomness between different simulated cohorts, 
we simulate all patients with extra genes and phenotypes then downsample to get the desired counts
using the allowed modules for the given ablation.
'''


def parse_args():
    parser = argparse.ArgumentParser(description='Simulate patients.')

    # parameters for phenotype ablation analyses 
    parser.add_argument('--no_phen_noise', action='store_true', help='Remove phenotypic noise from the pipeline.')
    parser.add_argument('--no_phen_dropout', action='store_true', help='Remove phenotypic dropout from the pipeline.')
    parser.add_argument('--no_phen_corruption', action='store_true', help='Remove phenotype corruption from the pipeline.')
    parser.add_argument('--no_gene_module_phen', action='store_true', help='Remove phenotypes added through gene modules from the pipeline.')

    # parameters for gene ablation analyses 
    parser.add_argument('--random_genes', action='store_true', help='Randomly sample distractor genes instead of using gene modules.')
    parser.add_argument('--remove_module', type=str, default=None, choices= ['non_syndromic_phenotype', 'common_fp', 'tissue_distractor' , 
        'pathogenic_pheno_irrel', 'insufficient_explanatory', 'universal_distractor', 'phenotype_distractor', None], help='Replace the specified removed gene module with random gene modules.')

    args = parser.parse_args()
    return args

def filter_phens(phenotypes, distractor_ids):
    filtered_phenotypes = {}
    for phen, module_source_list in phenotypes.items():
        if 'init_phenotypes' in module_source_list: filtered_phenotypes[phen] = module_source_list
        elif 'noisy_phenotype' in module_source_list: filtered_phenotypes[phen] = module_source_list
        elif 'phenotype_corruption' in module_source_list: filtered_phenotypes[phen] = module_source_list
        elif len(set([m.split('.')[-1] for m in module_source_list]).intersection(set(distractor_ids))) > 0: filtered_phenotypes[phen] = module_source_list
    return filtered_phenotypes

def replace_gene_module(args, patients, random_gene_patients):
    '''
    replace ablated gene modules with randomly sampled gene modules instead of with the next sampled gene modules
    '''
    ablated_patients = []
    for patient, rand_patient in zip(patients, random_gene_patients):
        n_distractors = patient['n_distractor_genes'] 

        # make sure the current patient and patient with random genes have the same gene/disease
        assert patient['disease_id'] == rand_patient['disease_id'], 'Diseases don\'t match'
        assert patient['true_genes'][0] == rand_patient['true_genes'][0], 'True genes don\'t match'
        assert len(rand_patient['distractor_genes']) >= n_distractors, f'There aren\'t enough sampled genes of the specified type. Requested {n_distractors} genes but only have {len(rand_patient["distractor_genes"])}'

        # first subset distractor gene list to n_distractors 
        patient['distractor_genes'] = {gene:module_source_list for i, (gene,module_source_list) in enumerate(patient['distractor_genes'].items()) if i < n_distractors}

        #remove genes from specified module
        module = [args.remove_module]
        if module[0] == 'phenotype_distractor': module.append('universal_distractor')
        patient['distractor_genes'] = {gene:[m for m in module_source_list if m.split('.')[0] not in module] for gene, module_source_list in patient['distractor_genes'].items() if len([m for m in module_source_list if m.split('.')[0] not in module]) > 0} 
        pos_phenotypes = {phen:[m for m in module_source_list if m.split('.')[0] != module] for phen, module_source_list in patient['positive_phenotypes'].items() if len([m for m in module_source_list if m.split('.')[0] not in module]) > 0} 
        neg_phenotypes = {phen:[m for m in module_source_list if m.split('.')[0] != module] for phen, module_source_list in patient['negative_phenotypes'].items() if len([m for m in module_source_list if m.split('.')[0] not in module]) > 0} 

        #add in random genes
        random_genes = list(rand_patient['distractor_genes'].items())
        num_missing_genes = n_distractors - len(patient['distractor_genes'])
        
        # NOTE: added genes were all randomly sampled
        i = 0
        if num_missing_genes > 0:
            while True:
                key, val = random_genes[i]
                if key not in patient['distractor_genes']:
                    patient['distractor_genes'][key] = val
                i += 1
                if len(patient['distractor_genes']) == n_distractors:
                    break

        assert len(patient['distractor_genes']) == n_distractors, f'There aren\'t enough genes {len(patient["distractor_genes"])} vs {n_distractors}.'

        # subset distractor gene list to n_distractors if haven't already
        distractor_ids = [module_source_list[0].split('.')[-1] for gene, module_source_list in patient['distractor_genes'].items()]
        
        # filter pos & neg phenotypes to those added by the subsetted distractor gene list
        patient['positive_phenotypes'] = filter_phens(pos_phenotypes, distractor_ids)
        patient['negative_phenotypes'] = filter_phens(neg_phenotypes, distractor_ids)

        ablated_patients.append(patient)
    return ablated_patients

def remove_phen_module(args, patients):
    '''
    Remove phenotypes from the specified phenotype modules
    '''

    ablated_patients = []
    for patient in patients:

        # subset distractor gene list to n_distractors 
        n_distractors = patient['n_distractor_genes'] 
        patient['distractor_genes'] = {gene:module_source_list for i, (gene,module_source_list) in enumerate(patient['distractor_genes'].items()) if i < n_distractors}
        distractor_ids = [module_source_list[0].split('.')[-1] for gene, module_source_list in patient['distractor_genes'].items()]

        assert len(patient['distractor_genes']) >= n_distractors, 'There aren\'t enough sampled genes of the specified type.'

        # filter pos & neg phenotypes to those added by the subsetted distractor gene list
        patient['positive_phenotypes'] = filter_phens(patient['positive_phenotypes'], distractor_ids)
        patient['negative_phenotypes'] = filter_phens(patient['negative_phenotypes'], distractor_ids)

        # remove phenotypes added via noise
        if args.no_phen_noise:
            patient['positive_phenotypes'] = {phen:[m for m in module_source_list if m != 'noisy_phenotype'] for phen, module_source_list in patient['positive_phenotypes'].items() if len([m for m in module_source_list if m != 'noisy_phenotype']) > 0}
            patient['negative_phenotypes'] = {phen:[m for m in module_source_list if m != 'noisy_phenotype'] for phen, module_source_list in patient['negative_phenotypes'].items() if len([m for m in module_source_list if m != 'noisy_phenotype']) > 0}

        # remove phenotypes added via gene modules
        if args.no_gene_module_phen:
            phen_list = ['init_phenotypes', 'noisy_phenotype', 'phenotype_corruption']
            patient['positive_phenotypes'] = {phen:[m for m in module_source_list if m in phen_list] for phen, module_source_list in patient['positive_phenotypes'].items() if len([m for m in module_source_list if m in phen_list]) > 0}
            patient['negative_phenotypes'] = {phen:[m for m in module_source_list if m in phen_list] for phen, module_source_list in patient['negative_phenotypes'].items() if len([m for m in module_source_list if m in phen_list]) > 0}
        
        # remove corrupted phenotypes & add back in the phenotypes that were removed
        if args.no_phen_corruption:
            patient['positive_phenotypes'] = {phen:[m for m in module_source_list if m != 'phenotype_corruption'] for phen, module_source_list in patient['positive_phenotypes'].items() if len([m for m in module_source_list if m != 'phenotype_corruption']) > 0}
            patient['negative_phenotypes'] = {phen:[m for m in module_source_list if m != 'phenotype_corruption'] for phen, module_source_list in patient['negative_phenotypes'].items() if len([m for m in module_source_list if m != 'phenotype_corruption']) > 0}
            
            if 'positive_phenotypes' in patient['corruption_phenotypes']:
                for p in patient['corruption_phenotypes']['positive_phenotypes']: 
                    if p in patient['positive_phenotypes']: patient['positive_phenotypes'][p].append(['init_phenotypes'])
                    else: patient['positive_phenotypes'][p] = ['init_phenotypes']
            if 'negative_phenotypes' in patient['corruption_phenotypes']:
                for p in patient['corruption_phenotypes']['negative_phenotypes']: 
                    if p in patient['negative_phenotypes']: patient['negative_phenotypes'][p].append(['init_phenotypes'])
                    else: patient['negative_phenotypes'][p] = ['init_phenotypes']

        # add back in phenotypes that were removed with dropout
        if args.no_phen_dropout:
            if 'positive_phenotypes' in patient['dropout_phenotypes']:
                for p in patient['dropout_phenotypes']['positive_phenotypes']: 
                    if p in patient['positive_phenotypes']: patient['positive_phenotypes'][p].append(['init_phenotypes'])
                    else: patient['positive_phenotypes'][p] = ['init_phenotypes']
            if 'negative_phenotypes' in patient['dropout_phenotypes']:
                for p in patient['dropout_phenotypes']['negative_phenotypes']: 
                    if p in patient['negative_phenotypes']: patient['negative_phenotypes'][p].append(['init_phenotypes'])
                    else: patient['negative_phenotypes'][p] = ['init_phenotypes']

        ablated_patients.append(patient)
    
    return ablated_patients

def add_random_gene_modules(patients, random_gene_patients):
    '''
    Add random genes to each patient from the set of random_gene_patients
    '''
    ablated_patients = []
    max_n_distractors = 0
    for patient, rand_patient in zip(patients, random_gene_patients):
        n_distractors = patient['n_distractor_genes']
        if n_distractors > max_n_distractors: max_n_distractors = n_distractors
        assert len(patient['true_genes']) == 1
        assert patient['disease_id'] == rand_patient['disease_id'], 'Diseases don\'t match'
        assert patient['true_genes'][0] == rand_patient['true_genes'][0], 'True genes don\'t match'
        assert len(rand_patient['distractor_genes']) >= n_distractors, f'There aren\'t enough sampled genes of the specified type. Requested {n_distractors} genes but only have {len(rand_patient["distractor_genes"])}'
        
        # filter to the correct number of genes
        patient['distractor_genes'] = {gene:module_source_list for i, (gene,module_source_list) in enumerate(rand_patient['distractor_genes'].items()) if i < n_distractors}
        distractor_ids = [module_source_list[0].split('.')[-1] for gene, module_source_list in patient['distractor_genes'].items()]
        
        # filter pos & neg phenotypes to those added by the subsetted distractor gene list
        patient['positive_phenotypes'] = filter_phens(patient['positive_phenotypes'], distractor_ids)
        patient['negative_phenotypes'] = filter_phens(patient['negative_phenotypes'], distractor_ids)

        ablated_patients.append(patient)
    print(f'Max N distractors: {max_n_distractors}')
    return ablated_patients

def write_patients(patients, filename):
    '''
    write patients to jsonl file
    '''
    with open(config.SIMULATED_DATA_PATH / filename, "w") as output_file:
        for n, patient in enumerate(patients):
            json.dump(patient, output_file)
            output_file.write('\n')

def get_filename(args):
    '''
    Construct filename for ablation
    '''
    rand_genes = '_rand_genes' if args.random_genes else ''
    phen_corrupt = '_no_phencorrupt' if args.no_phen_corruption else ''
    phen_dropout = '_no_phendrop' if args.no_phen_dropout else ''
    phen_noise = '_no_phennoise' if args.no_phen_noise else ''
    phen_genemod = '_no_phengenemod' if args.no_gene_module_phen else ''
    if args.remove_module == 'phenotype_distractor': module = 'phenotype_and_universal_distractor'
    else: module = args.remove_module
    remove_module_type = ('_remove_' + module ) if args.remove_module != None else ''

    filename = "simulated_patients"  + rand_genes + phen_corrupt + phen_dropout + phen_noise + phen_genemod  + remove_module_type + "_ablations.jsonl"

    return filename

def get_equal_n_gene_modules(patients):
    modules = ['phenotype_universal_distractor', 'insufficient_explanatory', 'common_fp', 'tissue_distractor', 'pathogenic_pheno_irrel',  'non_syndromic_phenotype']
    for patient in tqdm(patients):
        n_distractors = patient['n_distractor_genes'] 

        module_types_to_genes = defaultdict(list)
        for gene, mods in patient['distractor_genes'].items():
             for m in mods: 
                mod_name = m.split('.')[0]
                if mod_name == 'phenotype_distractor' or mod_name == 'universal_distractor': mod_name = 'phenotype_universal_distractor'
                module_types_to_genes[mod_name].append({gene:[m]})
        module_types_to_genes_counts = {k: len(v) for k, v in module_types_to_genes.items()}
        avail_modules = [m for m in modules if m in module_types_to_genes_counts.keys()]
        while True:
            sampled_module_types = np.random.choice(avail_modules, n_distractors, p=[1.0/len(avail_modules)] * len(avail_modules))
            sampled_module_types_counts = Counter(sampled_module_types)
            enough_samples = [v <= module_types_to_genes_counts[k] for k, v in sampled_module_types_counts.items()]
            if  all(enough_samples): break 
        
        sampled_modules = [sampled for mod, count in sampled_module_types_counts.items() for sampled in module_types_to_genes[mod][:count]]
        sampled_modules = {list(d.keys())[0]:list(d.values())[0] for d in sampled_modules}
        assert len(np.unique(list(sampled_modules.keys()))) == len(sampled_modules), 'There are duplicate sampled genes'
        patient['distractor_genes'] = sampled_modules
        distractor_ids = [module_source_list[0].split('.')[-1] for gene, module_source_list in patient['distractor_genes'].items()]
        
        # subset distractor gene list to n_distractors 
        #patient['distractor_genes'] = {gene:module_source_list for i, (gene,module_source_list) in enumerate(patient['distractor_genes'].items()) if i < n_distractors}
        #distractor_ids = [module_source_list[0].split('.')[-1] for gene, module_source_list in patient['distractor_genes'].items()]

        assert len(patient['distractor_genes']) == n_distractors, f'There should be {n_distractors} distractors from {sampled_module_types_counts} categories, but only {len(patient["distractor_genes"])} were found: {patient["distractor_genes"]}\n {module_types_to_genes_counts}'

        # filter pos & neg phenotypes to those added by the subsetted distractor gene list
        patient['positive_phenotypes'] = filter_phens(patient['positive_phenotypes'], distractor_ids)
        patient['negative_phenotypes'] = filter_phens(patient['negative_phenotypes'], distractor_ids)

    return patients

def count_gene_modules_in_patients(patients):
    '''
    check how many gene modules exist per patient
    '''

    gene_modules = [Counter(['phenotype_universal_distractor' if (m.split('.')[0] == 'phenotype_distractor' or m.split('.')[0] == 'universal_distractor') else m.split('.')[0] for gene, modules in patient['distractor_genes'].items() for m in modules]) for patient in patients]
    gene_df = pd.DataFrame.from_records(gene_modules)
    gene_df = pd.melt(gene_df, var_name='module', value_name='count')
    gene_df.loc[pd.isnull(gene_df['count']), 'count'] = 0 #replace NaN with 0
    
    print('mean', gene_df.groupby('module').mean())
    print('median', gene_df.groupby('module').median())
    print('max', gene_df.groupby('module').max())
    print('min', gene_df.groupby('module').min())
    
    ax = sns.violinplot(x="module", y="count", data=gene_df)
    ax.set_xticklabels(ax.get_xticklabels(),rotation = 90)
    plt.tight_layout()
    plt.savefig(config.SIMULATED_DATA_PATH / 'ablations' / f'violin_plot_distactor_gene_modules.jpeg')
    plt.clf()
    plt.figure()
    ax = sns.boxplot(x="module", y="count", data=gene_df)
    ax.set_xticklabels(ax.get_xticklabels(),rotation = 90)
    plt.tight_layout()
    plt.savefig(config.SIMULATED_DATA_PATH / 'ablations' / f'box_plot_distactor_gene_modules.jpeg')
    plt.clf()

    custom_params = {"axes.spines.right": False, "axes.spines.top": False}
    sns.set_theme(style="ticks", rc=custom_params)
    palette = sns.color_palette("Paired", n_colors=6)
    gene_df = gene_df.sort_values(by=['module'])
    ax = sns.histplot(
        gene_df, x="count", hue="module", multiple="dodge", stat="percent", discrete=True, shrink = 0.8, common_norm=False) #, palette=palette 
    ax.set_xlabel( "Number of Candidate Genes")
    ax.set_ylabel( "Percent of Simulated Patients")
    plt.tight_layout()
    plt.savefig(config.SIMULATED_DATA_PATH / 'ablations' / f'hist_distactor_gene_modules.jpeg')

def main(args):
    # set seed
    SEED = 3
    random.seed(SEED)
    np.random.seed(SEED)

    # read in patients with many simulated genes
    equal_n_gene_modules_filename = config.SIMULATED_DATA_PATH / 'ablations' / 'simulated_patients_equal_probs_revised.jsonl'
    if equal_n_gene_modules_filename.exists():
        patients = read_jsonl(equal_n_gene_modules_filename) 
    else:
        patients = read_jsonl(config.SIMULATED_DATA_PATH / 'ablations' / 'simulated_patients_equal_probs_many_genes.jsonl') 
        patients = get_equal_n_gene_modules(patients)
        write_patients(patients, equal_n_gene_modules_filename)
    
    count_gene_modules_in_patients(patients)

    # replace existing gene modules with random genes & make sure that all associated phenotypes are removed
    if args.random_genes:
        random_gene_patients = read_jsonl(config.SIMULATED_DATA_PATH / 'ablations' / 'simulated_patients_rand_genes_many_genes.jsonl')
        patients = add_random_gene_modules(patients, random_gene_patients)
        args.no_gene_module_phen = True

    # remove one gene module & replace with randomly sampled genes
    if args.remove_module:
        random_gene_patients = read_jsonl(config.SIMULATED_DATA_PATH / 'ablations' / 'simulated_patients_rand_genes_many_genes.jsonl')
        patients = replace_gene_module(args, patients, random_gene_patients)

    # remove phenotype modules
    if args.no_phen_noise or args.no_gene_module_phen or args.no_phen_dropout or args.no_phen_corruption:
        patients = remove_phen_module(args, patients)

    # make sure each patient has the specified number of genes
    for patient in patients:
        if len(patient['distractor_genes']) != patient['n_distractor_genes']:
            patient['distractor_genes'] = {gene:module_source_list for i, (gene,module_source_list) in enumerate(patient['distractor_genes'].items()) if i < patient['n_distractor_genes']}
            distractor_ids = [module_source_list[0].split('.')[-1] for gene, module_source_list in patient['distractor_genes'].items()]
            assert len(patient['distractor_genes']) >= patient['n_distractor_genes'], 'There aren\'t enough sampled genes of the specified type.'
            
            # filter pos & neg phenotypes to those added by the subsetted distractor gene list
            patient['positive_phenotypes'] = filter_phens(patient['positive_phenotypes'], distractor_ids)
            patient['negative_phenotypes'] = filter_phens(patient['negative_phenotypes'], distractor_ids)

    # write ablations to file
    fname = get_filename(args)
    write_patients(patients, config.SIMULATED_DATA_PATH / 'ablations' / fname)


if __name__ == "__main__":
    args = parse_args()
    main(args)