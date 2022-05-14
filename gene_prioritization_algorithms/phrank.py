import sys
import argparse
from pathlib import Path
import json
from tqdm import tqdm
import pickle
import pandas as pd
import numpy as np
import jsonlines
from collections import defaultdict

import config
from simulation_pipeline.utils.util import read_simulated_patients

sys.path.insert(0, f'{config.PROJECT_ROOT}/phrank/') # add config to path
from phrank import Phrank
from phrank import utils as phrank_utils


DAG = config.PROJECT_ROOT / 'hpo' / '2019' / "hpodag.txt"
PRE_2015_GENE_TO_PHENO = config.KNOWLEDGE_GRAPH_PATH  / 'formatted_phenolyzer_2015/hpo_gene.tsv'
PRE_2015_DX_TO_PHENO = config.KNOWLEDGE_GRAPH_PATH  / 'formatted_phenolyzer_2015/hpo_disease.tsv'
PRE_2015_DX_TO_GENE = config.KNOWLEDGE_GRAPH_PATH  / 'formatted_phenolyzer_2015/gene_disease.tsv'


def run_phrank(args, patients, p_hpo):
    rank_dict = {}

    number_genes_with_0_score = []
    percent_genes_with_0_score = []
    if args.rank_type == 'genes':
        for patient in tqdm(patients):
            patient_id = patient['id']

            all_genes = set(patient['all_candidate_genes'])
            rank_dict[patient_id] = rank_genes(args, p_hpo, all_genes, patient['positive_phenotypes'])
            
            number_genes_with_0_score.append(len([score for score, gene in rank_dict[patient_id] if score == 0]))
            percent_genes_with_0_score.append(len([score for score, gene in rank_dict[patient_id] if score == 0])/len(rank_dict[patient_id]))
            assert len(all_genes) == len(rank_dict[patient_id]), 'Not all genes were ranked'
            

        assert len(patients) == len(rank_dict.keys()), 'Not all patients\' genes were ranked'
    else:
        raise NotImplementedError

    print(f'There are {np.mean(number_genes_with_0_score)} genes or {np.mean(percent_genes_with_0_score)*100:0.2f}% with scores of 0 on average per patient')
    return rank_dict

###################################################
# Phrank Methods

def rank_diseases(p_hpo, genes, phenotypes, normalized=False, baseline=False):
    # sorting the disease by best match
    return p_hpo.rank_diseases(genes, phenotypes, normalized=normalized, baseline=baseline)

def rank_genes(args, p_hpo, genes, phenotypes, normalized=False, baseline=False):
    # sorting the genes by best match
    if args.rank_genes_directly:
        # automatically gives score of 0 to genes not in KG
        return p_hpo.rank_genes_directly(genes, phenotypes, normalized=normalized, baseline=baseline)
    else:
        # only returns scores for genes that appear in KG so we need to add the missing genes
        ranked_genes = p_hpo.rank_genes_using_disease(genes, phenotypes, normalized=normalized, baseline=baseline)
        only_genes = [g for score, g in ranked_genes]
        missing_ranks = [(0,g) for g in genes if g not in only_genes ]
        return ranked_genes + missing_ranks

def main():
    parser = argparse.ArgumentParser(description='Run Phrank Algorithm')
    parser.add_argument('--sim_input', type=str, default=f'{config.SIMULATED_DATA_PATH}/simulated_patients_formatted.jsonl')
    parser.add_argument('--kg', type=str, default='pre_2015', help='Specify which knowledge graph to use')
    parser.add_argument('--rank_type', type=str, default='genes', help='Specify whether to rank "genes" or "diseases"')
    parser.add_argument('--rank_genes_directly', action='store_true', help='Specify whether rank genes directly or rank genes via diseases')
    parser.add_argument('--use_disease_annotations', action='store_true', help='Specify whether to initialize Phrank via G-D and D-P links rather than just G-P links')
    parser.add_argument('--output_dir', type=str, default=f'{config.PROJECT_ROOT}/gene_prioritization_results/phrank', help='Specify path where the results will be outputted')

    args = parser.parse_args()

    if args.kg == 'pre_2015':
        geneannotationsfile  = PRE_2015_GENE_TO_PHENO
        diseaseannotationsfile = PRE_2015_DX_TO_PHENO
        diseasegenefile = PRE_2015_DX_TO_GENE
    else:
        raise NotImplementedError

    # initalize phrank
    if args.rank_genes_directly and not args.use_disease_annotations:
        p_hpo = Phrank(DAG, geneannotationsfile=geneannotationsfile)
    else:
        assert diseaseannotationsfile != None, 'Missing disease annotation file'
        p_hpo = Phrank(DAG, diseaseannotationsfile=diseaseannotationsfile, diseasegenefile=diseasegenefile)

    # read in patients
    patients = read_simulated_patients(args.sim_input)
    
    # get filename
    sim_fname = Path(args.sim_input).stem #args.sim_input.split('.jsonl')[0]
    # if 'ablations/'  in sim_fname: 
    #     sim_fname = sim_fname.split('ablations/')[1]
    use_dx_annot_text = f'_use_dx_annot=True_' if args.use_disease_annotations else ''
    fname =  Path(args.output_dir) / (sim_fname + '_phrank_rankgenes_directly=' + str(args.rank_genes_directly) + use_dx_annot_text + '.pkl')
    
    
    # Run Phrank
    rank_dict = run_phrank(args, patients, p_hpo)

    # Save Results to file
    print(f'Saving to file...{fname}')
    with open(fname, 'wb') as handle:
        pickle.dump(rank_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)

if __name__ == "__main__":
    main()