import numpy as np
import pickle5 as pickle
from pathlib import Path
from collections import defaultdict
from scipy.stats import rankdata


####################################################
# Helper functions for categorizing patients

def categorize_patients(patients):
    patient_categories_dict = defaultdict(list)
    patient_broad_categories_dict = defaultdict(list)
    in_gene_gene_kg = []
    for patient in patients:
        if 'category' in patient and 'broad_category' in patient:
            category = patient['category']
            broad_category = patient['broad_category']
        elif 'udn_categories' in patient:
            broad_category = [c[3] for c in patient['udn_categories']]
            category = [c[2] for c in patient['udn_categories']]
        else:
            category = 'None'            
            broad_category = 'None'
        
        if type(broad_category) == list:
            for cat in broad_category: patient_broad_categories_dict[cat].append(patient)
            for cat in category: patient_categories_dict[cat].append(patient)
        else:
            patient_broad_categories_dict[broad_category].append(patient)
            patient_categories_dict[category].append(patient)


    
    return patient_categories_dict, patient_broad_categories_dict

####################################################
# Helper functions for getting performance metrics

def read_rankings(filename):
    with open(str(filename), 'rb') as handle:
        gene_rankings = pickle.load(handle)
    return gene_rankings

def find_correct_genes(filename, patients, gene_rankings, category, category_type):
    all_ranks = []
    correct_gene_bools = []
    for patient in patients:
        patient_id = patient['id'] 
        
        correct_genes = patient['true_genes'] 

        
        
        #filter correct genes to only consider those from a specific category
        if 'udn_categories' in patient and category != 'All':
            category_ind = 3 if category_type == 'broad' else 2
            correct_genes = [c[0] for c in patient['udn_categories'] if c[category_ind] == category]
        

        assert patient_id in gene_rankings, f'patient {patient_id} \'s genes, (correct={correct_genes}) were not ranked'
        gene_scores_orig = gene_rankings[patient_id]
        
        if 'phenomizer' in str(filename):
            if len(gene_scores_orig[0]) == 3:
                gene_scores = [(pval, score, 'ENSG00000003989') if gene == 'SLC7A2-IT1' else (pval, score, gene) for pval, score, gene in gene_scores_orig]

                ranked_genes = rankdata([score * -1 for pval, score, gene in gene_scores], method='average')
                pvals = [pval for pval, score, gene in gene_scores]
                rank_pval_added = ranked_genes + pvals
                reranked_genes = rankdata([add_score for add_score in rank_pval_added], method='average')
                ranks = [rank for rank, (pval, score, gene) in zip(reranked_genes, gene_scores) if gene in correct_genes]
                correct_gene_bool = [1 if gene in correct_genes else 0 for pval, score, gene in gene_scores]
            else:
                gene_scores = [(pval, score, 'ENSG00000003989', disease) if gene == 'SLC7A2-IT1' else (pval, score, gene, disease) for pval, score, gene, disease in gene_scores_orig]
                ranked_genes = rankdata([score * -1 for pval, score, gene, disease in gene_scores], method='average')
                pvals = [pval for pval, score, gene, disease in gene_scores]
                rank_pval_added = ranked_genes + pvals
                reranked_genes = rankdata([add_score for add_score in rank_pval_added], method='average')
                ranks = [rank for rank, (pval, score, gene, disease) in zip(reranked_genes, gene_scores) if gene in correct_genes]
                correct_gene_bool = [1 if gene in correct_genes else 0 for pval, score, gene, disease in gene_scores]

        else:
            # Remove genes that might have received two scores
            # This occurs in phenolyzer because two gene names mapped to the same ensembl ID. We take the best score.
            gene_scores_orig = sorted(gene_scores_orig, key=lambda x: x[0])
            seen = set()
            gene_scores = []
            for scores in gene_scores_orig:
                gene_score, gene = scores
                if gene not in seen:
                    gene_scores.append((gene_score, gene))
                    seen.add(gene)

            ranked_genes = rankdata([score * -1 for score, gene in gene_scores], method='average')
            ranks = [rank for rank, (score, gene) in zip(ranked_genes, gene_scores) if gene in correct_genes]
        
            correct_gene_bool = [1 if gene in correct_genes else 0 for score, gene in gene_scores]

        
        all_ranks.append(ranks)
        correct_gene_bools.append(correct_gene_bool)

    return all_ranks, correct_gene_bools

def mean_reciprocal_rank(ranks):
    return np.mean([1.0/r for r in ranks])

def average_rank(ranks):
    return np.mean(ranks)

def top_k_accuracy(ranks, k):
    if len(ranks) > 0:
        return len([r for r in ranks if r <= k])/len(ranks)
    else:
        return None

def mean_average_precision(correct_gene_bools):
    all_map = []
    for bools in correct_gene_bools:
        all_map.append(np.mean([sum(bools[:i+1])/len(bools[:i+1]) for i, b in enumerate(bools) if b==1]))
    return np.mean(all_map)


def evaluate(filename, gene_rankings, patients, category, category_type, eval_type = 'duplicate_patients' ): 
    results_dict = defaultdict(list)
    intermed_results_dict = defaultdict(list)
    ranks, correct_gene_bools = find_correct_genes(filename, patients, gene_rankings, category, category_type)

    # change ranks of patients with multiple correct genes so that all subsequent genes 
    # don't consider previously ranked correct genes

    if eval_type == 'top_gene': #only consider gene with best rank
        ranks = [min(rank) for rank in ranks]
    elif eval_type == 'duplicate_patients': #include ranks of all correct genes
        ranks = [r for rank in ranks for r in rank]
    else:
        raise NotImplementedError
    
    results_dict['avg_rank'] = average_rank(ranks)
    results_dict['mrr'] = mean_reciprocal_rank(ranks)
    results_dict['top_1_acc'] = top_k_accuracy(ranks, 1)
    results_dict['top_3_acc'] = top_k_accuracy(ranks, 3)
    results_dict['top_5_acc'] = top_k_accuracy(ranks, 5)
    results_dict['top_10_acc'] = top_k_accuracy(ranks, 10)
    results_dict['ranks'] = ranks

    return results_dict

