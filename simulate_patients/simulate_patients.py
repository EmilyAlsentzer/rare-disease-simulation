
import sys
import os
from tqdm import tqdm
from pathlib import Path
import pandas as pd
import numpy as np
import random
import networkx
import obonet
import pickle
import json
from random import shuffle
import argparse
import logging
from collections import defaultdict

import matplotlib.pyplot as plt
plt.switch_backend('agg')
import seaborn as sns

sys.path.insert(0, '../') # add config to path
import config


from modules.disease import Disease
from modules.patient import Patient
from modules.patient_simulator import PatientSimulator

from utils.util import create_disease_dict
 pd.options.mode.chained_assignment = None


SEED = 42


# dict to keep track of the the # of genes added from each gene module
simulation_outcome_counts = {
    'non_syndromic_phenotype' : defaultdict(int),
    'common_false_positive' : defaultdict(int),
    'tissue_distractor' : defaultdict(int),
    'pathogenic_phenotype_irrelevant' : defaultdict(int),
    'insufficient_explainer' : defaultdict(int),
    'universal_distractor' : defaultdict(int),
    'phenotype_distractor' : defaultdict(int)
}

def simulate_patient(args, simulator, patient, disease):
    '''
    simulates a single patient:
        1. intialize phenotypes
        2. for sampled n_distractor_genes, run distractor gene module 
        3. add noisy phenotypes
    '''

    logging.info('--- Initialize phenotypes --- ')
    
    # whether to perform phenotype corruption/dropout
    perform_phenotype_corruption = ~ args.no_phen_corruption
    perform_phenotype_dropout = ~ args.no_phen_dropout

    # initialize phenotypes for patient with disease
    simulator.initialize_phenotypes(patient, disease, \
        config.PROB_DROPOUT_POS, config.PROB_DROPOUT_NEG, config.PROB_CORRUPTION_POS, config.PROB_CORRUPTION_NEG, \
            perform_phenotype_dropout=perform_phenotype_dropout, perform_phenotype_corruption=perform_phenotype_corruption)

    # sample n_distractor_genes
    logging.info('--- Sampling Distractor Genes ---')
    
    # First determine the number of distractor genes for the patient. 
    n_distractor_genes = 1 + np.random.poisson(config.N_DISTRACTORS_LAMBDA - 1) 
    patient.n_distractor_genes = n_distractor_genes
    if args.sim_many_genes: # for the ablation analysis, we want to simulate patients with a large number of genes & then downsample
        n_distractor_genes = 100

    n_distractor_genes_added = 0
    gene_sampler_names = list(simulator.gene_samplers.keys())
    unfinished = False
    
    if args.random_genes: #randomly add genes
        while n_distractor_genes_added < n_distractor_genes:
            did_add_gene = simulator.get_random_gene(patient)
            if did_add_gene == 'gene_added':
                n_distractor_genes_added += 1
    else: # we use the gene modules to add distractors
        n_tries_to_add_any_distractor = 0
        while n_distractor_genes_added < n_distractor_genes:

            # sample false gene generation module 
            if args.equal_probs: #sample each gene module with the same probability
                sampled_name = str(np.random.choice(gene_sampler_names,1)[0])
            else: # sample gene modules with probabilities specified in the config file
                sampled_name = str(np.random.choice(gene_sampler_names,1, p=[config.NON_SYNDROM_PHEN_PROB, \
                config.COMMON_FP_PROB, config.TISSUE_DIST_PROB, config.PATH_PHEN_PROB, config.INSUFF_EXPLAIN_PROB, \
                config.UNIVERSAL_DIST_PROB, config.PHENO_DIST_PROB])[0])
                

            gene_module, _ = simulator.gene_samplers[sampled_name]
            # In some cases, a specific gene module isn't compatible with a patient/gene. We set MAX_ADD_GENE_ATTEMPTS to limit the # of tries.
            # We also log information about whether a gene module was successfully used to add a distractor gene to the patient.
            for i in range(config.MAX_ADD_GENE_ATTEMPTS):  

                did_add_gene = gene_module(patient, n_distractor_genes_added)

                # Tabulate Simulation Outcomes
                simulation_outcome_counts[sampled_name][did_add_gene] += 1
                
                if did_add_gene == 'gene_added': 
                    if args.verbose:
                        logging.warning(f'Took {i + 1} tries to add gene using gene module {sampled_name}')
                    n_distractor_genes_added += 1
                    break
                elif did_add_gene == 'gene_impossible_to_add': 
                    if args.verbose:
                        logging.warning(f'It is impossible for gene module {sampled_name} to add any gene ') 
                    break 
                elif did_add_gene != 'gene_not_added':
                    raise Exception('One of the gene modules is returning an invalid category')
                
                # the gene module was not able to add any genes in MAX_ADD_GENE_ATTEMPTS attempts
                if i == config.MAX_ADD_GENE_ATTEMPTS - 1:
                    if args.verbose:
                        logging.warning(f'Reached maximum number of attempts using module {sampled_name} for patient {patient}')
                    simulation_outcome_counts[sampled_name]["max_attempts"] += 1
            
            if n_tries_to_add_any_distractor > config.MAX_ADD_ANY_DISTRACTOR_ATTEMPTS:
                if args.verbose:
                    logging.warning(f'Failed to add sufficient distractors to patient: {patient} ')
                unfinished = True
                break   
            n_tries_to_add_any_distractor += 1

    # Sample noisy phenotypes & add to the patient
    logging.info('--- Sampling Noisy Phenotypes ---')
    if not args.no_phen_noise:
        simulator.sample_noisy_phenotypes(patient)

    return patient, unfinished

def get_dataset_statistics(patients):
    '''
    Calculate average # of distractor genes, pos & neg phenotypes, and other basic statistics about simulated data
    '''
    n_distractors = []
    n_positive_phenotypes = []
    n_negative_phenotypes = []
    for patient in patients:
        n_distractors.append(len(patient.get_distractor_genes()))
        n_positive_phenotypes.append(len(patient.get_hpo_set(is_positive=True)))
        n_negative_phenotypes.append(len(patient.get_hpo_set(is_positive=False)))
    logging.info('Total number of patients simulated: {}'.format(len(patients)))
    logging.info('Average number of distractor genes: {}'.format(sum(n_distractors)/len(n_distractors)))
    logging.info('Average number of positive phenotypes: {}'.format(sum(n_positive_phenotypes)/len(n_positive_phenotypes)))
    logging.info('Average number of negative phenotypes: {}'.format(sum(n_negative_phenotypes)/len(n_negative_phenotypes)))


def read_data(args):
    '''
    read in orphanet data
    '''
    orphanet_phenotypes = pd.read_csv(config.ORPHANET_PATH / 'orphanet_final_disease_hpo_normalized_2015.tsv', sep='\t', dtype=str)
    orphanet_genes = pd.read_csv(config.ORPHANET_PATH / 'orphanet_final_disease_genes_normalized_2015.tsv', sep='\t', dtype=str)
    orphanet_metadata = pd.read_csv(config.ORPHANET_PATH / 'orphanet_final_disease_metadata_normalized_2015.tsv', sep='\t', dtype=str)

    print(f'There are {len(orphanet_metadata.index)} diseases to simulate')
    return orphanet_phenotypes, orphanet_genes, orphanet_metadata

def parse_args():
    parser = argparse.ArgumentParser(description='Simulate rare disease patients')

    parser.add_argument('--random_n_patients', action='store_true', help='Whether to specify a random # of patients per orphanet disease')
    parser.add_argument('--verbose', action='store_true', help='Additional logging')

    # parameters for phenotype ablation analyses. These remove different sources of phenotypes from the pipeline
    parser.add_argument('--no_phen_noise', action='store_true', help='Remove phenotypic noise from the pipeline.')
    parser.add_argument('--no_phen_dropout', action='store_true', help='Remove phenotypic dropout from the pipeline.')
    parser.add_argument('--no_phen_corruption', action='store_true', help='Remove phenotype corruption from the pipeline.')
    parser.add_argument('--no_gene_module_phen', action='store_true', help='Remove phenotypes added through gene modules from the pipeline.')

    # parameters for ablation analyses 
    parser.add_argument('--sim_many_genes', action='store_true', help='Simulate patients with many distractor genes each. These can be downsampled for ablation analyses.')
    parser.add_argument('--random_genes', action='store_true', help='Randomly sample distractor genes instead of using gene modules.')
    parser.add_argument('--equal_probs', action='store_true', help='Sample gene modules with equal probabilities. This is used for the ablation analysis.')

    args = parser.parse_args()
    return args

def get_filename(args):
    '''
    Returns the filename associated with the simulation run
    '''

    rand_genes = '_rand_genes' if args.random_genes else ''
    phen_corrupt = '_no_phencorrupt' if args.no_phen_corruption else ''
    phen_dropout = '_no_phendrop' if args.no_phen_dropout else ''
    phen_noise = '_no_phennoise' if args.no_phen_noise else ''
    phen_genemod = '_no_phengenemod' if args.no_gene_module_phen else ''
    equal_probs = '_equal_probs' if args.equal_probs else ''
    sim_many_genes = '_many_genes' if args.sim_many_genes else ''

    filename = "simulated_patients" + rand_genes + phen_corrupt + phen_dropout + phen_noise + phen_genemod + equal_probs + sim_many_genes + ".jsonl"
    
    # if we're performing an ablation analysis, add it to the ablation folder
    if args.random_genes or args.no_phen_corruption or args.no_phen_dropout or args.no_phen_noise or args.no_gene_module_phen or args.equal_probs or args.sim_many_genes:
        filename = 'ablations/' + filename

    return filename


def run_simulation(args, filename):
    '''
    Run Rare Diease Patient Simulation
    '''

    logging.basicConfig(format='%(message)s', level=logging.WARNING)

    # set random seed to ensure replicability
    random.seed(SEED)
    np.random.seed(SEED)
    os.environ['PYTHONHASHSEED']=str(SEED)

    # read in orphanet data & filter to diseases from timstamp = args.timstamp
    orphanet_phenotypes, orphanet_genes, orphanet_metadata = read_data(args)

    # create dict mapping orphanet_id -> Disease
    disease_dict = create_disease_dict(orphanet_phenotypes, orphanet_genes, orphanet_metadata)
    
    # initialize simulator
    logging.info('------Initializing Patient Simulator------')
    add_genemodule_distractor_phenotypes = ~ args.no_gene_module_phen
    simulator = PatientSimulator(disease_dict, config.STRONG_PHENOTYPE_THRESH, config.WEAK_PHENOTYPE_THRESH, \
        add_distractor_phenotypes=add_genemodule_distractor_phenotypes, seed=SEED)

    
    # for each disease, simulate X patients 
    logging.info('\n\n------Simulating Patients------')
    patients = []
    unfinished_patients = []
    for orphanet_id, disease in tqdm(disease_dict.items()):
        if args.random_n_patients:
            # NOTE: random_n_patients will simulate PATIENTS_PER_DISEASE on average
            n_sampled_patients = 1 + np.random.poisson(config.PATIENTS_PER_DISEASE - 1) 
        else:
            # NOTE: if we take the constant approach, exactly PATIENTS_PER_DISEASE will be simulated per disease
            n_sampled_patients = config.PATIENTS_PER_DISEASE
        
        for p in range(n_sampled_patients): 
            patient = Patient(disease)
            simulated_patient, unfinished = simulate_patient(args, simulator, patient, disease)
            if unfinished:
                unfinished_patients.append(simulated_patient)
            patients.append(simulated_patient)

    print("Simulation Outcome Rates by Gene Module:")
    print(simulation_outcome_counts)

    print(f'There were {len(unfinished_patients)} patients with incomplete distractor gene sets.')

    # calculate statistics on the newly simulated data
    get_dataset_statistics(patients)

    # randomly order patients
    shuffle(patients)

    # generate save filename based on different configurations
    logging.info('\n\n------Write Patients to File------')

    #save to jsonl file
    with open(config.SIMULATED_DATA_PATH / filename, "w") as output_file:
        for n, patient in enumerate(patients):
            patient_dict = patient.to_dict()
            patient_dict['id'] = n
            
            json.dump(patient_dict, output_file)
            output_file.write('\n')

if __name__ == "__main__":
    args = parse_args()

    # get filename to save to file
    filename = get_filename(args)

    # run simulation
    run_simulation(args, filename)