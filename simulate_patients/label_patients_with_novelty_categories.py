import json
import pandas as pd
import numpy as np
import jsonlines
import random
from collections import defaultdict, Counter
import argparse 
from tqdm import tqdm

import sys
sys.path.insert(0, '../') # add config to path
import config

#orphanet files
orphanet_phenotypes = pd.read_csv(config.ORPHANET_PATH / 'orphanet_final_disease_hpo_normalized_2015.tsv', sep='\t')
orphanet_genes = pd.read_csv(config.ORPHANET_PATH / 'orphanet_final_disease_genes_normalized_2015.tsv', sep='\t')
orphanet_genes['OrphaNumber'] = orphanet_genes['OrphaNumber'].astype(str)
orphanet_disease_metadata = pd.read_csv(config.ORPHANET_PATH / 'orphanet_final_disease_metadata_normalized_2015.tsv', sep='\t')
orphanet_disease_metadata['OrphaNumber'] = orphanet_disease_metadata['OrphaNumber'].astype(str)
orphanet_gene_metadata = pd.read_csv(config.ORPHANET_PATH / 'orphanet_final_gene_metadata_normalized_2015.tsv', sep = '\t')

# KG files (Phrank format)
gene_disease_kg = pd.read_csv(config.KNOWLEDGE_GRAPH_PATH / 'formatted_phenolyzer_2015' / 'gene_disease.tsv', sep='\t', header=None)
gene_disease_kg.columns = ['ensembl_id', 'disease_id']
hpo_gene_kg = pd.read_csv(config.KNOWLEDGE_GRAPH_PATH / 'formatted_phenolyzer_2015' / 'hpo_gene.tsv', sep='\t', header=None)
hpo_gene_kg.columns = ['hpo_id', 'ensembl_id']
all_kg_genes = set(gene_disease_kg['ensembl_id'].tolist() + hpo_gene_kg['ensembl_id'].tolist())

# G-G KG files (Phenolyzer format)
human_geneid = pd.read_csv(config.KNOWLEDGE_GRAPH_PATH / 'raw_normalized_phenolyzer_2015' / 'DB_HUMAN_GENE_ID', sep='\t', dtype=str, skiprows=1, header=None)
human_geneid.columns= ['', 'GENE_ID', 'GENE_LIST']
human_geneid['GENE_LIST'] = human_geneid['GENE_LIST'].str.replace('|$', '')
human_geneid['GENE_LIST'] = human_geneid['GENE_LIST'].str.replace('^,|', '')
human_geneid['GENE_LIST'] = human_geneid['GENE_LIST'].str.split('|')
human_geneid = human_geneid.explode('GENE_LIST')

mentha = pd.read_csv(config.KNOWLEDGE_GRAPH_PATH / 'raw_normalized_phenolyzer_2015' / 'DB_MENTHA_GENE_GENE_INTERACTION', dtype=str,sep='\t')
protein_interaction = pd.read_csv(config.KNOWLEDGE_GRAPH_PATH / 'raw_normalized_phenolyzer_2015' / 'DB_COMPILED_BINARY_PROTEIN_INTERACTION_SCORE', dtype=str,sep='\t')
biosystem = pd.read_csv(config.KNOWLEDGE_GRAPH_PATH / 'raw_normalized_phenolyzer_2015' / 'DB_COMPILED_BIOSYSTEM_SCORE', dtype=str,sep='\t')
hgnc_family = pd.read_csv(config.KNOWLEDGE_GRAPH_PATH / 'raw_normalized_phenolyzer_2015' / 'DB_HGNC_GENE_FAMILY', dtype=str,sep='\t')
transcription_interaction = pd.read_csv(config.KNOWLEDGE_GRAPH_PATH / 'raw_normalized_phenolyzer_2015' / 'DB_HTRI_TRANSCRIPTION_INTERACTION', dtype=str,sep='\t')
all_gg_kg_genes = set(mentha['Gene A'].tolist() + mentha['Gene B'].tolist() \
        + protein_interaction['PROTEIN_1'].tolist() + protein_interaction['PROTEIN_2'].tolist() \
        + biosystem['GENE'].tolist() + hgnc_family['GENE'].tolist() \
        + transcription_interaction['TF'].tolist() + transcription_interaction['TG'].tolist() \
        + human_geneid['GENE_ID'].tolist() + human_geneid['GENE_LIST'].tolist())


def read_jsonl(filename):
    '''
    Read in patients
    '''
    print('filename: ', filename)
    patients = []
    with jsonlines.open(filename) as reader:
        for patient in reader:
            patients.append(patient)    
    return patients

def disease_known_at_KG_time(patient):
    '''
    when disease metadata file col First_Published_post_2015_03_01 == False
    '''
    disease = str(patient['disease_id'])
    disease_not_known_kg_bool = orphanet_disease_metadata.loc[orphanet_disease_metadata['OrphaNumber'] == disease, 'First_Published_post_2015_01']
    assert len(disease_not_known_kg_bool) > 0, f'disease {disease} is missing from disease metadata file'
    assert type(disease_not_known_kg_bool.tolist()[0]) == bool
    return not disease_not_known_kg_bool.tolist()[0]

def disease_gene_assoc_known_at_KG_time(patient):
    '''
    in disease-genes file, DG_Assoc_First_Published_Post_2015_03_01 == False for the specific gene/disease pair
    '''
    disease = patient['disease_id']
    gene = patient['true_genes'][0]
    dg_assoc_bool = orphanet_genes.loc[(orphanet_genes['OrphaNumber'] == disease) & (orphanet_genes['Ensembl_ID'] == gene), 'DG_Assoc_First_Published_Post_2015_01']
    assert len(dg_assoc_bool) > 0, f'disease-gene pair ({disease}, {gene}) is missing from gene file'
    assert type(dg_assoc_bool.tolist()[0]) == bool

    return not dg_assoc_bool.tolist()[0]

def gene_cause_disease_in_KG(patient):
    '''
    whether the gene is known to cause disease in the KG
    '''
    gene = patient['true_genes'][0]
    return gene in gene_disease_kg['ensembl_id'].tolist()

def gene_known_to_cause_dx_at_KG_time(patient):
    '''
    True if D-G link in KG or in Orphanet metadata file

    NOTE: here we consider whether G-D relationship exists in KG because we don't assume that the orphanet annotations
    contain all D-G associations
    '''

    if disease_gene_assoc_known_at_KG_time(patient):
        return True

    if gene_cause_disease_in_KG(patient):
        return True
        
    # in gene metadata file, Gene_In_Any_Orpha_Disease_2015_03 == True for the gene
    gene = patient['true_genes'][0]
    gene_cause_dx_at_kg_bool = orphanet_gene_metadata.loc[orphanet_gene_metadata['Ensembl_ID'] == gene, 'Gene_In_Any_Orpha_Disease_2015_03']
    assert len(gene_cause_dx_at_kg_bool) > 0, 'gene is missing from gene metadata file'
    return gene_cause_dx_at_kg_bool.tolist()[0]

def gene_in_kg(patient):
    '''Returns true if the gene is found in the phenolyzer KG (excluding gene-gene edges)'''
    gene = patient['true_genes'][0]
    return gene in all_kg_genes

def gene_in_gg_kg(patient):
    '''Returns true if the gene is found in the Gene-Gene KG'''
    gene = patient['true_genes'][0]
    return gene in all_gg_kg_genes

def label_with_novelty_categories(patients):
    #new gene, new disease
    new_gene_new_disease_gene_notinany_kg = 0
    new_gene_new_disease_gene_notin_kg = 0
    new_gene_new_disease_gene_in_kg = 0

    #known gene, new disease
    known_gene_new_disease = 0

    # known gene, known disease
    known_gene_disease = 0
  
    # new gene, known disease
    new_gene_known_disease_gnikg = 0
    new_gene_known_disease_gniany_kg = 0
    new_gene_known_disease_gikg = 0

    # known gene, known disease i.e. Gene previously associated w/ another disease
    known_gene_known_disease = 0

    # keep track of patients where their gene is in G-G KG
    gene_in_g_g_kg = 0
    gene_in_dgp_kg = 0
    gene_in_g_g_kg_not_dgp_kg = []
    genes_not_in_any_kg = []
    
    for patient in tqdm(patients):
        
        # hacky fix to map SLC7A2-IT1 gene to SLC7A2. Ideally this should be in the preprocessor code 
        if patient['true_genes'][0] == 'SLC7A2-IT1':
            patient['true_genes'] = ['ENSG00000003989']
        assert len(patient['true_genes']) == 1, 'Patient has more than one correct gene'


        # KNOWN GENE-DISEASE 
        # -----------------------------------------
        # KGD
        if disease_gene_assoc_known_at_KG_time(patient):
            patient['broad_category'] = 'known_gene_disease'
            if gene_in_kg(patient):
                patient['category'] = 'known_gene_disease_gene_in_kg'
            elif gene_in_gg_kg(patient):
                patient['category'] = 'known_gene_disease_gene_only_in_gg_kg'
            else:
                patient['category'] = 'known_gene_disease_gene_notin_kg'
        

        # NEW GENE, NEW DISEASE 
        # ---------------------------------------
        # NG-ND
        elif not disease_known_at_KG_time(patient) and not gene_known_to_cause_dx_at_KG_time(patient):
            patient['broad_category'] = 'new_gene_new_disease'
            if gene_in_kg(patient):
                patient['category'] = 'new_gene_new_disease_gene_in_kg'
            elif gene_in_gg_kg(patient):
                patient['category'] = 'new_gene_new_disease_gene_only_in_gg_kg'
            else:
                patient['category'] = 'new_gene_new_disease_gene_notin_kg'


        # KNOWN GENE, NEW DISEASE 
        # -----------------------------------------
        # KG-ND
        elif not disease_known_at_KG_time(patient) and gene_known_to_cause_dx_at_KG_time(patient):
            patient['broad_category'] = 'known_gene_new_disease'
            if gene_in_kg(patient):
                patient['category'] = 'known_gene_new_disease_gene_in_kg'
            elif gene_in_gg_kg(patient):
                patient['category'] = 'known_gene_new_disease_gene_only_in_gg_kg'
            else:
                patient['category'] = 'known_gene_new_disease_gene_notin_kg'


        # NEW GENE, KNOWN DISEASE 
        # -----------------------------------------

        # Gene never before associated with any disease
        # NG-KD
        elif disease_known_at_KG_time(patient) and not gene_known_to_cause_dx_at_KG_time(patient):
            patient['broad_category'] = 'new_gene_known_disease'
            if gene_in_kg(patient):
                patient['category'] = 'new_gene_known_disease_gene_in_kg'
            elif gene_in_gg_kg(patient):
                patient['category'] = 'new_gene_known_disease_gene_only_in_gg_kg'
            else:
                patient['category'] = 'new_gene_known_disease_gene_not_in_kg'



        # Gene previously associated w/ another disease
        # KG-KD
        elif disease_known_at_KG_time(patient) and gene_known_to_cause_dx_at_KG_time(patient) and not disease_gene_assoc_known_at_KG_time(patient):
            patient['broad_category'] = 'known_gene_known_disease'
            if gene_in_kg(patient):
                patient['category'] = 'known_gene_known_disease_gene_in_kg'
            elif gene_in_gg_kg(patient):
                patient['category'] = 'known_gene_known_disease_gene_only_in_gg_kg'
            else:
                patient['category'] = 'known_gene_known_disease_gene_notin_kg'


        else:
            print('disease known at KG time: ', disease_known_at_KG_time(patient))
            print('gene in kg', gene_in_kg(patient))
            print('gene_known_to_cause_dx_at_KG_time', gene_known_to_cause_dx_at_KG_time(patient))
            print('d-g association known', disease_gene_assoc_known_at_KG_time(patient) )
            print('d-g in train set', disease_gene_assoc_in_dataset(patient, datasets['train']))
            raise Exception(f"patient with gene: {patient['true_gene']} and disease: {patient['disease_id']} doesn't fit into any of the categories")


        if gene_in_gg_kg(patient): 
            patient['in_gene_gene_kg'] = True
            gene_in_g_g_kg += 1
        else: 
            patient['in_gene_gene_kg'] = False
            genes_not_in_any_kg.append(patient['true_genes'][0])

        if gene_in_kg(patient):
            gene_in_dgp_kg += 1
            patient['in_disease_gene_kg'] = True
        else:
            patient['in_disease_gene_kg'] = False
            if gene_in_gg_kg(patient): gene_in_g_g_kg_not_dgp_kg.append(patient['true_genes'][0])


    print("genes not in any KG: ", genes_not_in_any_kg, set(genes_not_in_any_kg))
    print("Number of patients & the unique genes where the causal gene is in G-G KG, but not the D-G-P KG: ", len(gene_in_g_g_kg_not_dgp_kg), len(set(gene_in_g_g_kg_not_dgp_kg)))
    print('Count of each broad category: ', Counter([p['broad_category'] for p in patients]))
    print('Count of each category: ', Counter([p['category'] for p in patients]))
    print('----------------------------------------\n')

    return patients

def write_to_file(dataset, input_filename):
    '''
    Write patients to file
    '''

    output_filename = input_filename.split('.jsonl')[0] + '_formatted.jsonl'
    with open(config.SIMULATED_DATA_PATH / output_filename, "w") as output_file:
        for patient in dataset:

            # create final list of genes
            patient['all_candidate_genes'] = list(patient['distractor_genes'].keys()) + patient['true_genes']
            random.shuffle(patient['all_candidate_genes']) 

            json.dump(patient, output_file)
            output_file.write('\n')


def main():
    parser = argparse.ArgumentParser(description='Annotate patients with their gene-disease association category.')
    parser.add_argument('--input', type=str)
    args = parser.parse_args()

    #read in patients from file
    print('Reading in patients...')
    patients = read_jsonl(config.SIMULATED_DATA_PATH / args.input)

    # label patients with their gene-disease novelty categories
    print('Labeling patients...')
    patients = label_with_novelty_categories(patients)

    # write formatted patients back to jsonl file
    print('Writing patients to file...')
    write_to_file(patients, args.input)

if __name__ == "__main__":
    main()