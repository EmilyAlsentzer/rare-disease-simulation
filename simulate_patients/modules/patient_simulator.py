import typing
from typing import Set, Dict, Tuple, List
import numpy as np
import pandas as pd
import networkx
import obonet
import random
import sys
import pickle as pkl
import copy
from collections import defaultdict
import logging
import statistics
from scipy.special import softmax

from sklearn.metrics.pairwise import cosine_similarity
from scipy.special import softmax

sys.path.insert(0, '../') # add config to path
import config

from modules.disease import Disease, phenotype_frequency_dict
from modules.patient import Patient

from utils.util import get_nonspecific_phenotype_list



class PatientSimulator():

###################################
#Gene Sampler

    def _is_intersect(self, a,b) -> bool :
        '''Returns true if there is an intersection between the two sets'''
        return len(a.intersection(b)) > 0

    def _is_phenotype_distractor(self, disease, distractor_disease) -> bool:
        '''
        Return true if distractor_disease can be a phenotype distractor for disease.
        Phenotype distractors have overlapping weakly associated phenotypes & non-overlapping
        strongly associated phenotypes
        '''
        disease_strong_phenotypes = disease.get_phenotype_set(min_freq=config.STRONG_PHENOTYPE_THRESH)
        disease_weak_phenotypes = disease.get_phenotype_set(min_freq="Very rare (<4-1%)", max_freq=config.WEAK_PHENOTYPE_THRESH) 
        distractor_strong_phenotypes = distractor_disease.get_phenotype_set(min_freq=config.STRONG_PHENOTYPE_THRESH)
        distractor_weak_phenotypes = distractor_disease.get_phenotype_set(min_freq="Very rare (<4-1%)", max_freq=config.WEAK_PHENOTYPE_THRESH)
        return (not self._is_intersect(disease_strong_phenotypes,distractor_strong_phenotypes) and 
            self._is_intersect(disease_weak_phenotypes, distractor_weak_phenotypes))                    

    def read_pkl(self, fname):
        with open(str(fname), "rb") as f:
            pkl_file =  pkl.load(f)
        return pkl_file

    def write_pkl(self, obj, fname):
        with open(str(fname), "wb") as f:
            pkl.dump(obj, f) 

    def _initialize_phenotype_distractor_dict(self) -> Dict:
        '''
        Initializes a dict from orphanet id -> set(orphanet id) that contains all 
        of the orphanet diseases that can be a phenotype distractor for each disease

        Returns: 
            - phenotype_distractor_dict: ``Dict[str, Set(str)]``
        '''
        phenotype_distractor_loc = config.SIMULATOR_DIR / 'phenotype_distractor_dict.pkl'
        if phenotype_distractor_loc.exists() and not self.override:
            # load dict
            phenotype_distractor_dict = self.read_pkl(phenotype_distractor_loc)
        else:
            phenotype_distractor_dict = defaultdict(set)
            for orphanet_id, disease in self.diseases_dict.items():
                for distractor_id, distractor_disease in self.diseases_dict.items():
                    if orphanet_id == distractor_id: continue
                    if self._is_phenotype_distractor(disease, distractor_disease):
                            phenotype_distractor_dict[orphanet_id].add(distractor_id)
                
            # calculate statistics on coverage
            n_distractors = [len(distractor_ids) for orphanet_id, distractor_ids in phenotype_distractor_dict.items()]
            logging.info('{} percent of orphanet diseases have at least one candidate phenotype distractor disease.' \
                .format(100*len(phenotype_distractor_dict.keys())/len(self.diseases_dict.keys())))
            logging.info('There are {:0.1f} phenotype distractor diseases on average.'.format(sum(n_distractors)/len(n_distractors)))

            phenotype_distractor_dict = {k:list(v) for k, v in phenotype_distractor_dict.items()}

            #write dict to file
            self.write_pkl(phenotype_distractor_dict, phenotype_distractor_loc)
        
        return phenotype_distractor_dict

    def _initialize_universal_distractor_dict(self) -> Dict:
        '''
        Initializes a dict from orphanet id -> set(orphanet id) that contains all of the 
        orphanet diseases that can be a universal distractor for each disease. 
        Universal distractors must have at least 1 excluded or obligate phenotype.

        Furthermore, there are two categories of universal distractors:
        (1) diseases with obligate phenotypes whose obligate phenotypes don't overlap 
            with the disease's obligate phenotypes
        (2) diseases with excluded phenotypes whose excluded phenotypes don't overlap 
            with the disease's excluded phenotypes

        Returns:
             - universal_distractor_dict: Dict[str,set(str)]
        '''
        universal_distractor_loc = config.SIMULATOR_DIR / 'universal_distractor_dict.pkl'
        if universal_distractor_loc.exists() and not self.override:
            universal_distractor_dict = self.read_pkl(universal_distractor_loc)
        else:
            universal_distractor_dict = defaultdict(set)
            for orphanet_id, disease in self.diseases_dict.items():
                for distractor_id, distractor_disease in self.diseases_dict.items():
                    if orphanet_id == distractor_id: continue
                    distractor_obligates = distractor_disease.get_obligate_phenotypes()
                    distractor_excluded = distractor_disease.get_excluded_phenotypes()
                    if len(distractor_obligates) > 0:
                        #make sure obligate phenotypes don't overlap
                        if not self._is_intersect(disease.get_obligate_phenotypes(), distractor_obligates):
                            universal_distractor_dict[orphanet_id].add(distractor_id)

                    if len(distractor_excluded) > 0:
                        #make sure excluded phenotypes don't overlap
                        if not self._is_intersect(disease.get_excluded_phenotypes(), distractor_excluded):
                            universal_distractor_dict[orphanet_id].add(distractor_id)
            
            universal_distractor_dict = {k:list(v) for k, v in universal_distractor_dict.items()}

            self.write_pkl(universal_distractor_dict, universal_distractor_loc)

        # calculate statistics on coverage
        n_distractors = [len(distractor_ids) for orphanet_id, distractor_ids in universal_distractor_dict.items()]
        logging.info('{} percent of orphanet diseases have at least one candidate universal distractor disease.' \
            .format(100*len(universal_distractor_dict.keys())/len(self.diseases_dict.keys())))
        logging.info('There are {} universal distractor diseases on average.'.format(sum(n_distractors)/len(n_distractors)))
        
        return universal_distractor_dict # orphanet_id -> set(orphanet_ids)
    
    def get_disgenet_exclusive_gene_hpo(self):
        '''
        The goal is to find genes associated with phenotypes, but that are not associated with monogenic diseases. 

        To do this, we create a pandas dataframe containing gene-HPO edges from disegenet,
        excluding any genes listed as monogenic disease-causing in Orphanet or HPOA.

        These non-disease genes are used in the insufficiently explanatory or gene causing causing non-syndromic
        phenotytes gene modules.
        '''

        # Read in HPOA, Orpha, and Disgenet Data
        hpoa_genes_phenotypes = pd.read_csv(config.HPOA_PATH / "ALL_SOURCES_ALL_FREQUENCIES_diseases_to_genes_to_phenotypes_normalized.txt",
                                        delimiter= "\t",
                                        names = ['DiseaseID', "GeneID", "Gene_Symbol", "HPO_Name", "HPO_ID", 'ensembl_ids', "HPO_ID_2019"],
                                       skiprows=1)
        orpha_disease_genes = pd.read_csv(config.ORPHANET_PATH / "orphanet_final_disease_genes_normalized.tsv", delimiter="\t")

        disegenet_pheno_full = pd.read_csv(config.DISGENET_PATH / "all_gene_disease_associations_normalized.tsv", delimiter="\t").\
                            query("diseaseType == 'phenotype'")
        # Find Genes Listed in Disgenet, but not Orpha or HPOA
        genes_hpo_orpha = set(orpha_disease_genes.Ensembl_ID.values).\
                            union(set(hpoa_genes_phenotypes.ensembl_ids.values))
        genes_disegenet_only = set(disegenet_pheno_full.ensembl_ids.values).\
                            difference(genes_hpo_orpha)

        # Subset disegenet to non-orpha/hpo genes
        disegenet_exclusive_pheno_full = disegenet_pheno_full.query("ensembl_ids in @genes_disegenet_only")
        disegenet_exclusive_pheno_full = disegenet_exclusive_pheno_full[['ensembl_ids', 'diseaseId', 'diseaseName']]
        
        disegenet_mappings = pd.read_csv(config.DISGENET_PATH / "hpo_disease_mappings_normalized.tsv", delimiter="\t")[['diseaseId', 'HPO_ID_2019', 'vocabularyName']]
        disegenet_mappings.set_index('diseaseId', inplace=True)

        # Map Disegenet diseaseIDs onto HPO Codes and format column names
        disegenet_exclusive_pheno_full.set_index('diseaseId', inplace=True)
        disgenet_exclusive_gene_hpo = disegenet_exclusive_pheno_full.join(disegenet_mappings, how = "inner")
        disgenet_exclusive_gene_hpo.rename(columns={"ensembl_ids":"ensembl_ids",
                                        "diseaseName":"Pheno_String",
                                        "HPO_ID_2019":"HPO_Code",
                                        "vocabularyName":"HPO_String"},
                                inplace=True)
        disgenet_exclusive_gene_hpo = disgenet_exclusive_gene_hpo[["ensembl_ids", "HPO_Code", "HPO_String"]]
        return disgenet_exclusive_gene_hpo

    def _initialize_nondisease_gene_phenotype_dicts(self):
        '''
        Load the dataframe of disgenet-exclusive gene-HPO edges,
        and turn into dictionaries.
        '''
        nondisease_gene_to_hpo_loc = config.SIMULATOR_DIR / 'nondisease_gene_to_hpo_dict.pkl'
        nondisease_hpo_to_genes_loc = config.SIMULATOR_DIR / 'nondisease_hpo_to_genes_dict.pkl'

        if nondisease_gene_to_hpo_loc.exists() and not self.override:
            # load dict
            nondisease_gene_to_hpo = self.read_pkl(nondisease_gene_to_hpo_loc)
            nondisease_hpo_to_genes = self.read_pkl(nondisease_hpo_to_genes_loc)

        else:
            disgenet_exclusive_gene_hpo = self.get_disgenet_exclusive_gene_hpo()
            logging.info(f"There are {len(disgenet_exclusive_gene_hpo.index)} nondisease gene-HPO edges in disgenet")
            nondisease_gene_to_hpo = defaultdict(set)
            nondisease_hpo_to_genes = defaultdict(set)
            for index, row  in disgenet_exclusive_gene_hpo.iterrows():
                gene, hpo = row.ensembl_ids, row.HPO_Code
                nondisease_gene_to_hpo[gene] = nondisease_gene_to_hpo[gene].union({hpo})
                nondisease_hpo_to_genes[hpo] = nondisease_hpo_to_genes[hpo].union({gene})

            nondisease_gene_to_hpo = {k:list(v) for k, v in nondisease_gene_to_hpo.items()}
            nondisease_hpo_to_genes = {k:list(v) for k, v in nondisease_hpo_to_genes.items()}


            self.write_pkl(nondisease_gene_to_hpo, nondisease_gene_to_hpo_loc)
            self.write_pkl(nondisease_hpo_to_genes, nondisease_hpo_to_genes_loc)

        return nondisease_gene_to_hpo, nondisease_hpo_to_genes

    def _initialize_non_syndromic_phenotype_distractor_dict(self):
        '''
        Initializes a dict from orphanet_id to all non-disease genes that could serve as 
        non syndromic phenotype distractors for the disease. Non syndromic phenotype distractors must 
        not have any overlapping phenotypes with the disease.

        Returns:
            - non_syndromic_phenotype_distractor_dict: Dict[str, set(str)]
        '''
        nonsyndromic_distractor_loc = config.SIMULATOR_DIR / 'non_syndromic_distractor_dict.pkl'
        if nonsyndromic_distractor_loc.exists() and not self.override:
            non_syndromic_phenotype_distractor_dict = self.read_pkl(nonsyndromic_distractor_loc)
        else:

            non_syndromic_phenotype_distractor_dict = defaultdict(set)

            for orphanet_id, disease in self.diseases_dict.items():
                for gene_id, nondisease_hpos in self.nondisease_gene_to_hpo.items():
                    disease_hpos = disease.get_phenotype_set(min_freq="Very rare (<4-1%)")
                    # if the disease & candidate distractor don't have overlapping hpos
                    if not self._is_intersect(disease_hpos, set(nondisease_hpos)):
                        non_syndromic_phenotype_distractor_dict[orphanet_id].add(gene_id)
            
            non_syndromic_phenotype_distractor_dict = {k:list(v) for k, v in non_syndromic_phenotype_distractor_dict.items()}

            self.write_pkl(non_syndromic_phenotype_distractor_dict, nonsyndromic_distractor_loc)

        n_distractors = [len(genes) for orphanet_id, genes in non_syndromic_phenotype_distractor_dict.items()]
        logging.info('{} percent of orphanet diseases have at least one candidate non-syndromic phenotype distractor gene.' \
            .format(100*len(non_syndromic_phenotype_distractor_dict.keys())/len(self.diseases_dict.keys())))
        logging.info('There are {:0.1f} non-syndromic phenotype distractor genes on average.' \
            .format(sum(n_distractors)/len(n_distractors)))
        
        return non_syndromic_phenotype_distractor_dict

    def _initialize_insufficient_explanatory_dict(self):
        '''
        Initializes a dict mapping from disease id to a set of all genes that can 
        serve as insufficient explanatory distractors. These genes are not associated 
        with disease (i.e. they are genes in nondisease_gene_phenotype_dict), 
        but whose associated phenotypes are also associated with the disease's
        positive (non-excluded) phenotypes. These genes only partially explain the disease's
        phenotypes.  

        Returns:
            - insufficient_explanatory_distractor_dict: Dict[str, Set(str)]
        '''
        insufficient_explanatory_loc = config.SIMULATOR_DIR / 'insufficient_explanatory_dict.pkl'
        if insufficient_explanatory_loc.exists() and not self.override:
            insufficient_explanatory_distractor_dict = self.read_pkl(insufficient_explanatory_loc)
        else:
            insufficient_explanatory_distractor_dict = defaultdict(set)

            for orphanet_id, disease in self.diseases_dict.items():
                for gene_id, nondisease_hpos in self.nondisease_gene_to_hpo.items():
                    # don't consider excluded or obligate phenotypes when looking at intersection
                    disease_hpos = disease.get_phenotype_set(min_freq="Very rare (<4-1%)", max_freq="Very frequent (99-80%)")
                    n_overlap_hpos = len(disease_hpos.intersection(set(nondisease_hpos)))
                    
                    # to qualify as a insuff. explanaatory distractor, genes must have phenotypes overlap of (2,3) with the disease
                    if  n_overlap_hpos > 0:
                        if  n_overlap_hpos/ len(disease_hpos) <= config.INSUFF_EXPLAINER_HPO_OVERLAP: #check to make sure overlap isn't > 25%
                            insufficient_explanatory_distractor_dict[orphanet_id].add(gene_id)
            
            insufficient_explanatory_distractor_dict = {k:list(v) for k, v in insufficient_explanatory_distractor_dict.items()}

            self.write_pkl(insufficient_explanatory_distractor_dict, insufficient_explanatory_loc)

        n_distractors = [len(genes) for orphanet_id, genes in insufficient_explanatory_distractor_dict.items()]
        logging.info('{} percent of orphanet diseases have at least one candidate insufficient explanatory distractor gene.' \
            .format(100*len(insufficient_explanatory_distractor_dict.keys())/len(self.diseases_dict.keys())))
        logging.info('There are {:0.1f} insufficient explanatory distractor genes per disease on average.' \
            .format(sum(n_distractors)/len(n_distractors)))
        logging.info('There is a max of {:0.1f} and median of {:0.1f} insufficient explanatory distractor genes per disease.' \
            .format(max(n_distractors), statistics.median(n_distractors)))

        return insufficient_explanatory_distractor_dict

    def _load_gtex_data(self):

        # Load Gtex
        gtex_tpm_by_tissue = pd.read_csv(config.GTEX_PATH / "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.tsv", sep='\t')

        # Average rows with same Ensembl_ID/Description (i.e. different versions of an Ensembl ID)
        gtex_tpm_by_tissue["Ensembl_ID"] = [x.split(".")[0] for x in gtex_tpm_by_tissue.Name.values]
        gtex_tpm_by_tissue = gtex_tpm_by_tissue.groupby("Ensembl_ID").mean()

        # Normalize Columns (i.e. within tissue type)
        gtex_tpm_by_tissue = (gtex_tpm_by_tissue-gtex_tpm_by_tissue.min())/(gtex_tpm_by_tissue.max() - gtex_tpm_by_tissue.min())

        # return df where index of df are Ensembl IDs
        return gtex_tpm_by_tissue

    def _initialize_tissue_driven_distractor_dict(self):
        
        # For each gene in the nondisease_gene_phenotype_dict:
            # calculate tissue expression profile simiarlity score with all other genes in the dict
            # output: dict from ensembl_id -> ([ensembl_id], [probabilities])
        tissue_distractor_loc = config.SIMULATOR_DIR / 'tissue_distractor_dict.pkl'
        tissue_distractor_dx_loc = config.SIMULATOR_DIR / 'tissue_distractor_disease_genes_dict.pkl'
        if tissue_distractor_loc.exists() and not self.override:
            tissue_distractor_nondisease_genes_dict = self.read_pkl(tissue_distractor_loc)
            tissue_distractor_disease_genes_dict = self.read_pkl(tissue_distractor_dx_loc)

        else:
            # Get Gtex
            gtex_tpm_by_tissue = self._load_gtex_data()
            
            # Calculate Similarity Matrix

            # get genes associated with disease
            disease_genes = [set(disease.get_gene_list()) for orphanet_id, disease in self.diseases_dict.items()]
            disease_genes = set.union(*disease_genes)


            nondisease_genes = set(self.nondisease_gene_to_hpo.keys())
            
            gtex_nondisease = gtex_tpm_by_tissue.query("index in @nondisease_genes")
            gtex_disease = gtex_tpm_by_tissue.query("index in @disease_genes")

            print(f'There are {len(gtex_nondisease.index)} nondisease genes and {len(gtex_disease.index)} disease genes in gtex')

            cosine_sim_nondisease = pd.DataFrame(data = cosine_similarity(gtex_disease, gtex_nondisease),
                        index=gtex_disease.index,
                        columns=gtex_nondisease.index)

            cosine_sim_disease = pd.DataFrame(data = cosine_similarity(gtex_disease, gtex_disease),
                        index=gtex_disease.index,
                        columns=gtex_disease.index)


            tissue_distractor_nondisease_genes_dict = {}
            tissue_distractor_disease_genes_dict = {}
            for i in range(cosine_sim_nondisease.shape[0]):
                gene_key = cosine_sim_nondisease.index.get_level_values("Ensembl_ID")[i]
                disease_gene_values = cosine_sim_disease.iloc[i,:].sort_values(ascending=False)
                nondisease_gene_values = cosine_sim_nondisease.iloc[i,:].sort_values(ascending=False)
                nondisease_gene_values_names = nondisease_gene_values.index.values
                disease_gene_values_names = disease_gene_values.index.values

                tissue_distractor_nondisease_genes_dict[gene_key] = (nondisease_gene_values_names, nondisease_gene_values.values)
                tissue_distractor_disease_genes_dict[gene_key] =   (disease_gene_values_names, disease_gene_values.values)
            self.write_pkl(tissue_distractor_nondisease_genes_dict, tissue_distractor_loc)
            self.write_pkl(tissue_distractor_disease_genes_dict, tissue_distractor_dx_loc)

                

        return tissue_distractor_nondisease_genes_dict, tissue_distractor_disease_genes_dict #ensembl_id -> ([ensembl_id], [probabilities])

    def _initialize_pathogenic_but_phenotypically_irrelevant_dict(self):
        '''Loop through pairs of diseases in orphanet, if d1 and d2 have no phenotypes in common,
        add d2's genes to d1's entry in the dict
        Returns: dict from orphanet_id -> set(genes)'''
        pathogenic_pheno_irrel_loc = config.SIMULATOR_DIR / 'pathogenic_pheno_irrel_dict.pkl'
        if pathogenic_pheno_irrel_loc.exists() and not self.override:
            pathogenic_pheno_irrel_dict = self.read_pkl(pathogenic_pheno_irrel_loc)
        else:

            pathogenic_pheno_irrel_dict = defaultdict(set)
            copy_diseases_dict = self.diseases_dict.copy()
            for orphanet_id_1, disease_1 in self.diseases_dict.items():
                for orphanet_id_2, disease_2 in copy_diseases_dict.items():
                    if orphanet_id_1 == orphanet_id_2: continue
                    phenotype_overlap = disease_1.get_phenotype_set(min_freq="Very rare (<4-1%)") \
                        .intersection(disease_2.get_phenotype_set(min_freq="Very rare (<4-1%)"))
                    if len(phenotype_overlap) == 0: #if no phenotype overlap
                        # add genes from disease_2 to the list of candidate genes for disease 1
                        candidate_genes = set(disease_2.get_gene_list())
                        # but first make sure to remove causal genes of disease 1
                        disease_1_causal_genes = set(disease_1.get_gene_list())
                        candidate_genes = candidate_genes.difference(disease_1_causal_genes)
                        pathogenic_pheno_irrel_dict[orphanet_id_1] = pathogenic_pheno_irrel_dict[orphanet_id_1].union(candidate_genes)
            
            # convert sets to list
            pathogenic_pheno_irrel_dict = {key: list(val) for key, val in pathogenic_pheno_irrel_dict.items()}
            
            #write to file
            self.write_pkl(pathogenic_pheno_irrel_dict, pathogenic_pheno_irrel_loc)

        n_distractors = [len(genes) for orphanet_id, genes in pathogenic_pheno_irrel_dict.items()]
        logging.info('{} percent of orphanet diseases have at least one candidate pathogenic but phenotypically irrelevant distractor gene.' \
            .format(100*len(pathogenic_pheno_irrel_dict.keys())/len(self.diseases_dict.keys())))
        logging.info('There are {:0.1f} pathogenic but phenotypically irrelevant distractor genes on average.' \
            .format(sum(n_distractors)/len(n_distractors)))
        
        return pathogenic_pheno_irrel_dict

 
    def _initialize_false_positive_gene_list(self):
        # The flags paper provides frequency of genes from a bunch of exomes
        # This function loads the frequencies and converts to probability distribution
        fp_genes_prob = pd.read_csv(config.FLAGS_PATH / "flags_s3_fp_genes_freqs_normalized.txt", delimiter="\t", index_col=0)
        
        # Limit to top common FP genes
        fp_genes_prob = fp_genes_prob.iloc[0:config.N_COMMON_FP_GENES]
        
        # Count is the number of times exome was mutated in two databases
        fp_genes_prob['Count'] = fp_genes_prob['Count'].astype(float)
        fp_genes_prob.Count = fp_genes_prob.Count/sum(fp_genes_prob.Count) 
        fp_genes_prob.rename(columns={"Count":"Prop"}, inplace=True)

        return (fp_genes_prob.Ensembl_ID.tolist(), fp_genes_prob.Prop.tolist())

    def temperature_softmax(self, x, T=1.0):
        return np.exp(x/T)/sum(np.exp(x/T))

    def _initialize_age_to_hpo_probs_dict(self):
        '''
        Create a dictionary of age-stratified HPO_ID incidences, computed from the Finlayson et al 2014 dataset.
        Keys:  Patient ages (e.g Onset_Infant, Onset_Child, Onset_Adolescent, Onset_Adult, Onset_Elderly)
        Values: Tuple of (hpo_ids, probs) where probs[i] is the sampling probability of hpo_ids[i] in the given agegroup.
        Note that right now the frequencies of sampling each HPO code is identical across all age groups.
        In the V2 version of this pipeline, we plan to replace this with truly age-stratified sampling of noisy phenotypes.
        '''
        hpo_infant = pd.read_csv(config.CLAIMS_PATH / "HPO_Pheno_Counts_Infants_normalized.tsv", delimiter = "\t")
        hpo_child = pd.read_csv(config.CLAIMS_PATH / "HPO_Pheno_Counts_Children_normalized.tsv", delimiter = "\t")
        hpo_adolescent = pd.read_csv(config.CLAIMS_PATH / "HPO_Pheno_Counts_Adolescent_normalized.tsv", delimiter = "\t")
        hpo_adult = pd.read_csv(config.CLAIMS_PATH / "HPO_Pheno_Counts_Adult_normalized.tsv", delimiter = "\t")
        hpo_elderly = pd.read_csv(config.CLAIMS_PATH / "HPO_Pheno_Counts_Elderly_normalized.tsv", delimiter = "\t")

        forbidden_hpos = list(self.hpo_ontology.predecessors('HP:0000118')) + ["HP:0011024", "HP:0011025", "HP:0030680", "HP:0010978"]

        age_to_hpo_probs_dict = {}
        hpo_list_str = ["Infant", "Child", "Adolescent", "Adult", "Elderly"]
        hpo_list = [hpo_infant, hpo_child, hpo_adolescent, hpo_adult, hpo_elderly]
        for i, hp in enumerate(hpo_list):
            hp = hp.query("HPO not in @forbidden_hpos")
            hp.MemberCounts = hp.MemberCounts/sum(hp.MemberCounts)
            hp.rename(columns={"MemberCounts":"Prop"}, inplace=True)
            
            hp['Prop'] = self.temperature_softmax(hp['Prop'], 0.5) #TESTING
            age_to_hpo_probs_dict["Onset_" + hpo_list_str[i]] = (hp.HPO.values.tolist(), hp.Prop.values.tolist())

        return age_to_hpo_probs_dict

    def _initialize_all_genes_list(self):
        all_genes_list_loc = config.SIMULATOR_DIR / 'all_genes_list.pkl'

        if all_genes_list_loc.exists():
            genes_hpo_orpha_disgenet_FPgenes = self.read_pkl(all_genes_list_loc)
        else:
            # Read in HPOA, Orpha, and Disgenet Data & FP gene sources
            fp_genes_prob = pd.read_csv(config.FLAGS_PATH / "flags_s3_fp_genes_freqs_normalized.txt", delimiter="\t")

            hpoa_genes_phenotypes = pd.read_csv(config.HPO_PATH / "ALL_SOURCES_ALL_FREQUENCIES_diseases_to_genes_to_phenotypes_normalized.txt",
                                            delimiter= "\t",
                                            names = ['DiseaseID', "GeneID", "Gene_Symbol", "HPO_Name", "HPO_ID", 'ensembl_ids', "HPO_ID_2019"],
                                        skiprows=1)
            orpha_disease_genes = pd.read_csv(config.ORPHANET_PATH / "orphanet_final_disease_genes_normalized_2015.tsv", delimiter="\t")

            disegenet_pheno_full = pd.read_csv(config.DISGENET_PATH / "all_gene_disease_associations_normalized.tsv", delimiter="\t").\
                                query("diseaseType == 'phenotype'")
            
            # take the union of all gene sources
            genes_hpo_orpha = set(orpha_disease_genes.Ensembl_ID.values).\
                                union(set(hpoa_genes_phenotypes.ensembl_ids.values))
            genes_hpo_orpha_disgenet_FPgenes = set(disegenet_pheno_full.ensembl_ids.values).\
                                union(genes_hpo_orpha).union(fp_genes_prob.Ensembl_ID.values)
            
            self.write_pkl(list(genes_hpo_orpha_disgenet_FPgenes), all_genes_list_loc)

        return genes_hpo_orpha_disgenet_FPgenes

##########################################################################################################

    def __init__(self, diseases_dict, strong_phenotype_thresh, weak_phenotype_thresh, add_distractor_phenotypes=True, override=False, seed=42):
        self.seed = seed
        
        logging.info('Initializing....')
        self.diseases_dict = diseases_dict
        self.override = override
        self.add_distractor_phenotypes = add_distractor_phenotypes
        self.hpo_ontology = obonet.read_obo(config.HPO_FILE) 
        self.broad_phen_list, self.nonspecific_phen_list = get_nonspecific_phenotype_list(self.hpo_ontology)

        self.all_gene_options = self._initialize_all_genes_list() #use this for randomly sampling genes

        logging.info('...Initializing age_to_hpo_probs_dict')
        self.age_to_hpo_probs_dict = self._initialize_age_to_hpo_probs_dict()

        logging.info('...Initializing nondisease_gene_phenotype_dicts')
        self.nondisease_gene_to_hpo, self.nondisease_hpo_to_genes = self._initialize_nondisease_gene_phenotype_dicts()
        
        logging.info('...Initializing non_syndromic_phenotype_distractor_dict')
        self.non_syndromic_phenotype_distractor_dict = self._initialize_non_syndromic_phenotype_distractor_dict()
        
        logging.info('...Initializing insufficient_explanatory_distractor_dict')
        self.insufficient_explanatory_distractor_dict = self._initialize_insufficient_explanatory_dict()

        logging.info('...Initializing pathogenic_but_pheno_irrel_dict')
        self.pathogenic_but_pheno_irrel_dict = self._initialize_pathogenic_but_phenotypically_irrelevant_dict()
        
        logging.info('...Initializing tissue_driven_distractor_dict')
        self.tissue_distractor_nondisease_genes_dict, self.tissue_distractor_disease_genes_dict = self._initialize_tissue_driven_distractor_dict()

        logging.info('...Initializing phenotype_distractor_dict')
        self.phenotype_distractor_dict = self._initialize_phenotype_distractor_dict()

        logging.info('...Initializing universal_distractor_dict')
        self.universal_distractor_dict = self._initialize_universal_distractor_dict()

        logging.info('...Initializing false_positive_gene_list')
        self.common_fp_genes, self.common_fp_gene_probabilities = self._initialize_false_positive_gene_list()

        self.phenotype_to_ancestors_dict = {}
        self.phenotype_to_descendents_dict = {}

        self.gene_samplers = {
            'non_syndromic_phenotype' : (self.get_non_syndromic_phenotype_gene, config.NON_SYNDROM_PHEN_PROB),
            'common_false_positive' : (self.get_common_false_positive_gene, config.COMMON_FP_PROB),
            'tissue_distractor' : (self.get_tissue_distractor_gene, config.TISSUE_DIST_PROB),
            'pathogenic_phenotype_irrelevant' : (self.get_pathogenic_phenotypically_irrelevant_gene, config.PATH_PHEN_PROB),
            'insufficient_explainer': (self.get_insufficient_explanatory_gene, config.INSUFF_EXPLAIN_PROB),
            'universal_distractor': (self.get_universal_distractor_gene, config.UNIVERSAL_DIST_PROB),
            'phenotype_distractor': (self.get_phenotype_distractor_gene, config.PHENO_DIST_PROB)
        }

    def get_children(self, p, descendents):
        '''
        Recursive function that gets the children of a given HPO code and adds to decendents set
        '''
        children = list(self.hpo_ontology.predecessors(p))
        if len(children) == 0: return set()       
        for c in children:
            descendents = descendents.union(self.get_children(c, descendents))
        return descendents.union(set(children))

    def get_phenotype_descendents(self, phenotypes: Set) -> Set:
        '''
        Returns a set containing all of the descendents in the HPO hierarchy of all of the phenotypes

        First tries to look up whether each phenotype is in phenotype_to_descendents_dict. If not, calculates descendents
        '''
        all_descendents = set()
        for hpo_id in phenotypes:
            if hpo_id in self.phenotype_to_descendents_dict:
                all_descendents = all_descendents.union(self.phenotype_to_descendents_dict[hpo_id])
            else:
                ancestors = self.get_children(hpo_id, all_descendents)
                self.phenotype_to_descendents_dict[hpo_id] = ancestors
                all_descendents = all_descendents.union(ancestors)
        return all_descendents

    def get_phenotype_ancestors(self, phenotypes: Set) -> Set:
        '''
        Returns a set containing all of the ancestors in the HPO hierarchy of all of the phenotypes

        First tries to look up whether each phenotype is in phenotype_to_ancestors_dict. If not, calculates ancestors

        '''
        all_ancestors = set()
        for hpo_id in phenotypes:
            if hpo_id in self.phenotype_to_ancestors_dict:
                all_ancestors = all_ancestors.union(self.phenotype_to_ancestors_dict[hpo_id])
            else:
                paths_to_root = list(networkx.all_simple_paths(self.hpo_ontology, source=hpo_id, target="HP:0000001"))
                hpo_in_path_to_root = [p for path in paths_to_root for p in path]
                self.phenotype_to_ancestors_dict[hpo_id] = set(hpo_in_path_to_root)
                all_ancestors = all_ancestors.union(set(hpo_in_path_to_root))
        return all_ancestors

    def can_add_phenotype(self, patient: Patient, hpo_id: str, is_positive: bool, allow_switch: bool = False) -> bool:
        '''
        Returns true if the phenotype can successfully be added to the patient
        i.e. no excluded phenotypes are added to the patient's positive phenotypes &
        no obligate phenotypes are added to the patient's negative phenotypes

        Phenotype can't be added to positive list if it's already in negative list & vice versa

        NOTE: this function does NOT check to see if a positive phenotype is already in the list of positive phenotypes 
        or vice versa. Because of this, a phenotype MAY end up in the true positives (negatives) AND distractor positives (negatives)
        '''

        disease = self.diseases_dict[patient.orphanet_id]
        if is_positive: #check to make sure the pos phenotype is not in excluded list or already in neg list
        # can't be negative, exluded, or descendent of negative or descendent of excluded
            excluded_phenotypes = disease.get_excluded_phenotypes() 
            patient_neg_phenotypes = patient.get_hpo_set(is_positive=False) # includes both true & distractor phenotypes
            decendents_patient_neg_phenotypes = self.get_phenotype_descendents(patient_neg_phenotypes)
            descendents_excluded_phenotypes = self.get_phenotype_descendents(excluded_phenotypes) 
            if allow_switch:
                return hpo_id not in excluded_phenotypes and hpo_id not in descendents_excluded_phenotypes
            else:
                return hpo_id not in excluded_phenotypes and hpo_id not in patient_neg_phenotypes \
                and hpo_id not in decendents_patient_neg_phenotypes and hpo_id not in descendents_excluded_phenotypes

        else: #check to make sure the neg phenotype is not in obligate list or already in pos list
            # can't be positive, obligate or ancestor of either positive or obligate

            obligate_phenotypes = disease.get_obligate_phenotypes()
            patient_pos_phenotypes = patient.get_hpo_set(is_positive=True) 
            ancestors_obligate_phenotypes = self.get_phenotype_ancestors(obligate_phenotypes)
            ancestors_pos_phenotypes = self.get_phenotype_ancestors(patient_pos_phenotypes)
            if allow_switch:
                return hpo_id not in obligate_phenotypes \
                    and hpo_id not in ancestors_obligate_phenotypes 
            else:
                return hpo_id not in obligate_phenotypes and hpo_id not in patient_pos_phenotypes \
                    and hpo_id not in ancestors_obligate_phenotypes and hpo_id not in ancestors_pos_phenotypes

    def can_add_gene(self, patient:Patient, gene: str) -> bool:
        '''
        check to see if the gene to add is a causal gene associated with the patient. If so, return False
        check to see if gene is already in the list of distractor genes. If so, return False.
        '''
        disease = self.diseases_dict[patient.orphanet_id]
        genes_assoc_with_disease = set(disease.get_gene_list())
        patient_distractor_genes = patient.get_distractor_genes()
        if gene in genes_assoc_with_disease:
            return False
        if gene in patient_distractor_genes:
            return False
        return True

    def get_random_gene(self, patient: Patient):
        '''
        Randomly adds a gene to the patient
        '''
        #randomly sample from gene sources
        
        sampled_gene = np.random.choice(self.all_gene_options, 1)[0]
        if self.can_add_gene(patient, sampled_gene):
            patient.add_distractor_gene(sampled_gene, 'random')
            return 'gene_added'
        else:
            return 'gene_not_added'

    def get_phenotype_distractor_gene(self, patient: Patient, current_gene_count: int):
        logging.info('Sampling phenotype distractor gene')

        ######  sample gene  ######
        
        distractor_disease_ids = self.phenotype_distractor_dict[patient.orphanet_id] if patient.orphanet_id in self.phenotype_distractor_dict else []
        if len(distractor_disease_ids) > 0:
            sampled_disease_id = np.random.choice(distractor_disease_ids)
            distractor_disease = self.diseases_dict[sampled_disease_id]
            sampled_gene = np.random.choice(distractor_disease.get_gene_list())
            if self.can_add_gene(patient, sampled_gene):

                #############################################################
                ######  add phenotypes  ######
                patient_disease = self.diseases_dict[patient.orphanet_id]
                disease_weak_phenotypes = patient_disease.get_phenotype_set(min_freq="Very rare (<4-1%)", max_freq=config.WEAK_PHENOTYPE_THRESH) 
                distractor_weak_phenotypes = distractor_disease.get_phenotype_set(min_freq="Very rare (<4-1%)", max_freq=config.WEAK_PHENOTYPE_THRESH)
                distractor_strong_phenotypes = list(distractor_disease.get_phenotype_set(min_freq=config.STRONG_PHENOTYPE_THRESH))
                
                # determine how many phenotypes to add (1 to MAX_PHENOTYPES)
                n_strong = 1 + np.random.poisson(config.STRONG_PHENOTYPES_LAMBDA - 1)
                n_weak = 1 + np.random.poisson(config.WEAK_PHENOTYPES_LAMBDA - 1)

                weak_intersection = list(disease_weak_phenotypes.intersection(distractor_weak_phenotypes))
                if len(weak_intersection) > 0:
                    # check to see if weak phenotypes are already in the patient
                    n_weak = min(n_weak, len(weak_intersection))
                    # sample phenotypes from intersection and add to patient
                    list_weak_distractors_to_add = []
                    np.random.shuffle(weak_intersection) # shuffle in place
                    for weak in weak_intersection:
                        if self.can_add_phenotype(patient, weak, is_positive=True):
                            list_weak_distractors_to_add.append(weak)
                        if len(list_weak_distractors_to_add) == n_weak: break
                    if len(list_weak_distractors_to_add) == 0:
                        return('gene_not_added')
                else:
                    raise Exception('There is no intersection between the weak phenotypes of the disease & distractor disease')

                np.random.shuffle(distractor_strong_phenotypes)
                list_strong_distractors_to_add = []
                for strong_distractor in distractor_strong_phenotypes:
                    if self.can_add_phenotype(patient, strong_distractor, is_positive=False):
                        list_strong_distractors_to_add.append(strong_distractor)
                    if len(list_strong_distractors_to_add) == n_strong: break
                if len(list_strong_distractors_to_add) == 0:
                    return('gene_not_added')

                # modify patient
                if self.add_distractor_phenotypes:
                    for strong_distractor in list_strong_distractors_to_add:                 
                        patient.add_phenotype(strong_distractor, f'phenotype_distractor.{current_gene_count}', is_positive=False, is_distractor=True)
                    for weak_distractor in list_weak_distractors_to_add:
                        patient.add_phenotype(weak_distractor, f'phenotype_distractor.{current_gene_count}', is_positive=True, is_distractor=True)
                patient.add_distractor_gene(sampled_gene, f'phenotype_distractor.{current_gene_count}')
                return('gene_added')
            else:
                return('gene_not_added')
        
        else:
            return("gene_impossible_to_add")

    def get_universal_distractor_gene(self, patient: Patient, current_gene_count: int):
        logging.info('Sampling universal distractor gene')
        ###### sample gene ######
        distractor_disease_ids = self.universal_distractor_dict[patient.orphanet_id]
        if len(list(distractor_disease_ids)) > 0:
            sampled_disease_id = np.random.choice(distractor_disease_ids)
            distractor_disease = self.diseases_dict[sampled_disease_id]
            sampled_gene = np.random.choice(distractor_disease.get_gene_list())
            if self.can_add_gene(patient, sampled_gene):

                ###### add phenotypes ######
                distractor_obligates = distractor_disease.get_obligate_phenotypes()
                distractor_excluded = distractor_disease.get_excluded_phenotypes()

                if len(distractor_obligates) == 0 and len(distractor_excluded) == 0:
                    return('gene_not_added')

                #determine how many phenotypes to add
                while True:
                    n_obligates = np.random.poisson(config.OBLIGATE_PHENOTYPES_LAMBDA) if len(distractor_obligates) > 0 else 0
                    n_excluded = np.random.poisson(config.EXCLUDED_PHENOTYPES_LAMBDA) if len(distractor_excluded) > 0 else 0
                    if n_excluded + n_obligates > 0: break
                
                #add obligates
                obligate_phenotypes_to_add = []
                if len(distractor_obligates) > 0:
                    n_obligates = min(n_obligates, len(list(distractor_obligates)))
                    obligate_samples = list(distractor_obligates)
                    np.random.shuffle(obligate_samples)
                    for obligate in obligate_samples:
                        if self.can_add_phenotype(patient, obligate, is_positive=False):
                            obligate_phenotypes_to_add.append(obligate)
                        if len(obligate_phenotypes_to_add) == n_obligates: break
                    
                #add excluded
                excluded_phenotypes_to_add = []
                if len(distractor_excluded) > 0:
                    n_excluded = min(n_excluded, len(list(distractor_excluded)))
                    excluded_samples = list(distractor_excluded)
                    np.random.shuffle(excluded_samples)
                    for excluded in excluded_samples:
                        if self.can_add_phenotype(patient, excluded, is_positive=True):
                            excluded_phenotypes_to_add.append(excluded)
                        if len(excluded_phenotypes_to_add) == n_excluded: break                      

                if len(obligate_phenotypes_to_add) == 0 and len(excluded_phenotypes_to_add) == 0:
                    return('gene_not_added')
            
                #add phenotypes in the intersection of true and distractor disorder phenotypes (excluding obligate & excluded)
                if self.add_distractor_phenotypes:
                    n_intersect_phenotypes = 1 + np.random.poisson(config.INTERSECT_PHENOTYPES_LAMBDA - 1)
                    distractor_phenotypes = distractor_disease.get_phenotype_set(min_freq="Very rare (<4-1%)", max_freq="Very frequent (99-80%)")
                    true_phenotypes = self.diseases_dict[patient.orphanet_id].get_phenotype_set(min_freq="Very rare (<4-1%)", max_freq="Very frequent (99-80%)")
                    intersect_phenotypes = list(distractor_phenotypes.intersection(true_phenotypes))
                    if len(intersect_phenotypes) > 0: # it's okay if we're not able to add any intersecting phenotypes
                        n_intersect_phenotypes = min(n_intersect_phenotypes, len(intersect_phenotypes))
                        np.random.shuffle(intersect_phenotypes)
                        n_phenotypes_added = 0
                        for intersect in intersect_phenotypes:
                            if self.can_add_phenotype(patient, intersect, is_positive=True):
                                patient.add_phenotype(intersect, sampler_source=f'universal_distractor.{current_gene_count}', is_positive=True, is_distractor=True) 
                                n_phenotypes_added += 1
                            if n_phenotypes_added == n_intersect_phenotypes: break
                    
                    
                # Modify patient
                # add phenotypes
                if self.add_distractor_phenotypes:
                    for excluded in excluded_phenotypes_to_add:
                        patient.add_phenotype(excluded, sampler_source=f'universal_distractor.{current_gene_count}', is_positive=True, is_distractor=True) 
                    for obligate in obligate_phenotypes_to_add:                          
                        patient.add_phenotype(obligate, sampler_source=f'universal_distractor.{current_gene_count}', is_positive=False, is_distractor=True) 
                
                # add gene
                patient.add_distractor_gene(sampled_gene, f'universal_distractor.{current_gene_count}')
                return 'gene_added'
            else:
                return 'gene_not_added'
        else:
            return 'gene_impossible_to_add'

    def get_non_syndromic_phenotype_gene(self, patient: Patient, current_gene_count: int):
        logging.info('Sampling non-syndromic phenotype gene')
        candidate_gene_ids = self.non_syndromic_phenotype_distractor_dict[patient.orphanet_id]
        if len(candidate_gene_ids) > 0:
            sampled_gene = np.random.choice(candidate_gene_ids, 1)[0]
            if self.can_add_gene(patient, sampled_gene):

                #sample phenotypes
                if self.add_distractor_phenotypes:
                    hpos = self.nondisease_gene_to_hpo[sampled_gene] 

                    # limit list to phenotypes not caused by patient's true disease
                    disease_hpos = self.diseases_dict[patient.orphanet_id].get_phenotype_set(min_freq="Very rare (<4-1%)", max_freq="Obligate (100%)")
                    hpos = [ h for h in hpos if h not in disease_hpos]

                    n_phenotypes = 1 + np.random.poisson( config.NON_SYNDROMIC_PHENOTYPES_LAMBDA - 1)
                    n_phenotypes = min(n_phenotypes, len(hpos)) #can't sample more hpos than exist in the list
                    np.random.shuffle(hpos)
                    n_phenotypes_added = 0
                    for p in hpos:
                        if self.can_add_phenotype(patient, p, is_positive=True):
                            patient.add_phenotype(p, f'non_syndromic_phenotype.{current_gene_count}', is_positive=True, is_distractor=True)
                            n_phenotypes_added += 1
                        if n_phenotypes_added == n_phenotypes: break

                if n_phenotypes_added > 0 or not self.add_distractor_phenotypes:
                    # add gene
                    patient.add_distractor_gene(sampled_gene, f'non_syndromic_phenotype.{current_gene_count}')
                    return 'gene_added'
                else:
                    return 'gene_not_added'
            else:
                return 'gene_not_added'
        else:
            logging.info('Unable to add non syndromic distractor gene')
            return 'gene_impossible_to_add'

    def disease_is_better_explainer(self, patient: Patient, phenotypes_to_add: List, candidate_gene_hpos: Set, disease_hpos: Set):
        '''
        Checks to see whether the disease explains more of the patient's phenotypes than the candidate distractor gene
        '''
    
        proposed_patient_phenotypes = patient.get_hpo_set(is_positive=True).union(set(phenotypes_to_add))

        patient_disease_intersect = proposed_patient_phenotypes.intersection(disease_hpos)
        patient_distractor_intersect = proposed_patient_phenotypes.intersection(candidate_gene_hpos)
        if len(patient_disease_intersect) > len(patient_distractor_intersect):
            return True
        return False

    def get_insufficient_explanatory_gene(self, patient: Patient, current_gene_count: int):
        logging.info('Sampling insufficent explainatory gene')
        if patient.orphanet_id not in self.insufficient_explanatory_distractor_dict:
            return 'gene_impossible_to_add'
        candidate_gene_ids = self.insufficient_explanatory_distractor_dict[patient.orphanet_id]
        if len(candidate_gene_ids) > 0:
            sampled_gene = np.random.choice(candidate_gene_ids,1)[0]
            
            if self.can_add_gene(patient, sampled_gene):

                true_pos_phenotypes = patient.get_true_phenotypes(is_positive=True)
                candidate_gene_hpos = self.nondisease_gene_to_hpo[sampled_gene]
                intersect = set(true_pos_phenotypes).intersection(set(candidate_gene_hpos))

                # if there's already an overlap between the patient's HPO & HPO associated with distractor gene, just add the gene
                if len(intersect) > 0:
                    patient.add_distractor_gene(sampled_gene, f'insufficient_explanatory.{current_gene_count}')
                    return 'gene_added' 
                
                # otherwise, add HPOs to the patient
                else:
                    # check to see if overlap >= 1 between patient's true pos phenotypes & phenotypes caused by insuff explanatory gene
                    # if not, add n_phenotypes
                    # check to make sure patient has more phenotypes explained by disease than by insuff explanatory gene. If so, add gene. 
                    
                    disease_hpos = self.diseases_dict[patient.orphanet_id].get_phenotype_set(min_freq="Very rare (<4-1%)", max_freq="Very frequent (99-80%)")
                    intersect_hpos = [ h for h in disease_hpos if h in candidate_gene_hpos] 
                    n_phenotypes = 1 + np.random.poisson(config.INSUFFICIENT_EXPLAINER_PHENOTYPE_LAMBDA - 1) 
                    n_phenotypes = min(n_phenotypes, len(intersect_hpos)) # can't sample more phenotypes that intersect than exist
                    np.random.shuffle(intersect_hpos)
                    phenotypes_to_add = []
                    for p in intersect_hpos:
                        if self.can_add_phenotype(patient, p, is_positive=True):
                            phenotypes_to_add.append(p)
                        if len(phenotypes_to_add) == n_phenotypes: break

                    if len(phenotypes_to_add) > 0 and self.disease_is_better_explainer(patient, phenotypes_to_add, candidate_gene_hpos, disease_hpos):
                        if self.add_distractor_phenotypes:
                            for p in phenotypes_to_add: patient.add_phenotype(p, f'insufficient_explanatory.{current_gene_count}', is_positive=True, is_distractor=True) 
                        patient.add_distractor_gene(sampled_gene, f'insufficient_explanatory.{current_gene_count}')
                        return 'gene_added' 
                    else:
                        return 'gene_not_added' 
            else:
                return 'gene_not_added'
        else:
            logging.info('Unable to add insufficient explanatory distractor gene')
            return 'gene_impossible_to_add'


        
        #use dict from _initialize_insufficient_explanatory_dict
        # loop through all genes in dict & see if phenotypes of gene
        # overlap with phenotypes of the patient's realization of the disease
        # (only the true positive phenotypes)

        # if that set is empty, then uniformly sample from the 
        # insufficient_explanatory_dict

    def get_tissue_distractor_gene(self, patient: Patient, current_gene_count: int):
        logging.info('Sampling tissue distractor gene')
        patient_gene = patient.get_true_gene()

        # consider disease genes with no phenotypic overlap
        if patient_gene in self.tissue_distractor_disease_genes_dict:
            candidate_dx_gene_distractors, dx_gene_probabilities = self.tissue_distractor_disease_genes_dict[patient_gene]

            #filter to only include disease genes that don't have any phenotypic overlap with patient's disease
            candidate_path_phen_irrel_genes = self.pathogenic_but_pheno_irrel_dict[patient.orphanet_id]
            filtered_candidates = [(distractor, prob) for distractor, prob in zip (candidate_dx_gene_distractors, dx_gene_probabilities) if distractor in candidate_path_phen_irrel_genes]
            candidate_distractors, probabilities = zip(*filtered_candidates)

            # merge with nondisease genes
            if patient_gene in self.tissue_distractor_nondisease_genes_dict:
                candidate_nondisease_distractors, nondisease_gene_probabilities = self.tissue_distractor_nondisease_genes_dict[patient_gene]
                
                candidate_distractors = list(candidate_distractors) + list(candidate_nondisease_distractors)
                probabilities = list(probabilities) + list(nondisease_gene_probabilities)

                # sort by probability
                probabilities, candidate_distractors =  zip(*sorted(zip(probabilities, candidate_distractors), reverse=True))

        # only nondisease genes to consider
        elif patient_gene in self.tissue_distractor_nondisease_genes_dict:
                candidate_distractors, probabilities = self.tissue_distractor_nondisease_genes_dict[patient_gene]     
        else:
            logging.info('The patient\'s gene was not in GTEX')
            return 'gene_impossible_to_add'
        
        # filter to top X genes & normalize via softmax
        candidate_distractors = list(candidate_distractors[0:config.N_TISSUE_GENES])

        if np.sum(list(probabilities[0:config.N_TISSUE_GENES])) == 0:
            probabilities = [1/config.N_TISSUE_GENES] * config.N_TISSUE_GENES
        probabilities = list(probabilities[0:config.N_TISSUE_GENES])/np.sum(list(probabilities[0:config.N_TISSUE_GENES]))

        sampled_gene = np.random.choice(candidate_distractors, replace=False,p = probabilities)
        if self.can_add_gene(patient, sampled_gene):
            patient.add_distractor_gene(sampled_gene, f'tissue_distractor.{current_gene_count}')
            logging.info(f'Adding tissue_distractor gene: {sampled_gene}')
            return 'gene_added'
        else:
            return 'gene_not_added' 

    def get_pathogenic_phenotypically_irrelevant_gene(self, patient: Patient, current_gene_count: int):
        logging.info('Sampling pathogenic but phenotypically irrelevant gene')
        # probability of sampling is proportionate to similarity of tissue expression
        candidate_genes = self.pathogenic_but_pheno_irrel_dict[patient.orphanet_id]
        if len(candidate_genes) > 0:
            sampled_gene = np.random.choice(candidate_genes, replace=False)
            if self.can_add_gene(patient, sampled_gene):
                patient.add_distractor_gene(sampled_gene, f'pathogenic_pheno_irrel.{current_gene_count}')
                logging.info(f'Adding pathogenic but phenotypically irrelevant gene: {sampled_gene}')
                return 'gene_added'
            else:
                return 'gene_not_added'
        else:
            return 'gene_impossible_to_add'

    def get_common_false_positive_gene(self, patient: Patient, current_gene_count: int) -> str:
        logging.info('Sampling common false positive gene')
        sampled_gene = np.random.choice(self.common_fp_genes, replace=False, p = self.common_fp_gene_probabilities)
        if self.can_add_gene(patient, sampled_gene):
            patient.add_distractor_gene(sampled_gene, f'common_fp.{current_gene_count}')
            logging.info(f'Adding FP gene: {sampled_gene}')
            return 'gene_added'
        return 'gene_not_added'
    

##############################################################
  
# Phenotype Sampler
    def _sample_initial_phenotypes(self, patient: Patient, disease: Disease) -> None:
        '''
        Loop through all phenotypes associated with the patient's disease
        with prob = freq of phenotype in the disease, add to patient's positive phenotypes
        Otherwise, add to patient's negative phenotypes

        At this point, all phenotypes will be added to either the pos or neg phenotype list
        '''
        logging.info('---- Sampling Initial Phenotypes ----')


        for hpo_id, freq_str in disease.get_phenotypes().items():

            freq = phenotype_frequency_dict[freq_str] #frequency of phenotype in the disease
            is_sampled = np.random.binomial(1, freq) == 1 #true if phenotype was sampled
            
            # NOTE: is_distractor is False because the phenotypes were not added during a distractor module
            if is_sampled and self.can_add_phenotype(patient, hpo_id, is_positive=True):
                patient.add_phenotype(hpo_id, sampler_source='init_phenotypes', is_positive=True, is_distractor=False)
                logging.info(f'Adding true positive phenotype {hpo_id}')
            elif not is_sampled and self.can_add_phenotype(patient, hpo_id, is_positive=False):
                patient.add_phenotype(hpo_id, sampler_source='init_phenotypes', is_positive=False, is_distractor=False) 
                logging.info(f'Adding true negative phenotype {hpo_id}')
            else:
                logging.info('ERROR: trying to add an obligate to list of negative phenotypes OR add an excluded to list of positive phenotypes')
        #print('true positive phenotypes: ', patient.get_true_phenotypes(is_positive=True))
        #print('true negative phenotypes: ', patient.get_true_phenotypes(is_positive=False))

    def _phenotype_dropout(self, patient: Patient, p_dropout_pos: float = 0, p_dropout_neg: float = 0) -> None:
        '''
        Loops through all phenotypes that the patient has (which right now is equivalent to all 
        phenotypes associated with disease) & with some probability p, dropout each phenotype.

        This simulates the imperfect diagnosis pipeline, as clinicians may not record all symptoms a patient
        does or does not have. 

        '''
        logging.info('---- Performing Phenotype Dropout ----')

        # dropout for true positive phenotypes
        true_phenotypes = list(patient.get_true_phenotypes(is_positive=True)).copy() #make a copy because we're removing in a loop
        for hpo_id in true_phenotypes:
            should_dropout = np.random.binomial(1, p_dropout_pos) == 1 #true if should drop phenotype
            if should_dropout:
                patient.remove_phenotype(hpo_id, is_positive=True)
                patient.track_dropout_phenotype(hpo_id, is_positive=True)
                logging.info(f'Dropping true positive phenotype {hpo_id}')

                
        # dropout for true negative phenotypes
        true_phenotypes = list(patient.get_true_phenotypes(is_positive=False)).copy()
        for hpo_id in true_phenotypes:
            should_dropout = np.random.binomial(1, p_dropout_neg) == 1 #true if should drop phenotype
            if should_dropout:
                patient.remove_phenotype(hpo_id, is_positive=False)
                patient.track_dropout_phenotype(hpo_id, is_positive=False)
                logging.info(f'Dropping true negative phenotype {hpo_id}')
        
        #print('true positive phenotypes: ', patient.get_true_phenotypes(is_positive=True))
        #print('true negative phenotypes: ', patient.get_true_phenotypes(is_positive=False))

    def _get_ontology_ancestor(self, hpo_id: str, n_generations: int = 1) -> str:
        '''
        Gets the parent of the hpo_id according to the HPO hierarchy
        '''
        try:
            ancestors  = list(self.hpo_ontology.successors(hpo_id))
            assert len(ancestors) > 0, 'The HPO code has no ancestors'
            if len(ancestors) > 1:
                logging.info(f'The HPO code {hpo_id} has {len(ancestors)} ancestors.')
            ancestor = np.random.choice(ancestors,1)[0] # select random ancestor from list in case there are multiple
            return ancestor
        except:
            logging.info(f'ERROR: Could not find ancestors of HPO ID: {hpo_id}')
            return None

    def _phenotype_corruption(self, patient: Patient, p_corrupt_pos: float = 0.0, p_corrupt_neg: float = 0.0) -> None:
        '''
        With some probability, corrupt each of the patient's phenotypes. Corrupted phenotypes
        are replaced with phenotypes that are parents of the current phenotypes in the HPO hierarchy 
        (i.e. they are less specific that the original).

        Note that when corrupting negative phenotypes, we must make sure that the replaced phenotype doesn't 
        invalidate one of the patient's positive or obligate phenotypes
        '''

        logging.info('---- Performing Phenotype Corruption ----')
        #positive phenotypes
        patient_pos_phenotypes = list(patient.get_true_phenotypes(is_positive=True)).copy()
        for hpo_id in patient_pos_phenotypes:
            should_corrupt = np.random.binomial(1, p_corrupt_pos) == 1 #true if should corrupt phenotype
            if should_corrupt:
                ancestor = self._get_ontology_ancestor(hpo_id, n_generations=1)

                #if the patient doesn't have an ancestor, we leave it as is
                # if ancestor is too nonspecific, we also leave as is
                if ancestor and ancestor not in self.nonspecific_phen_list and self.can_add_phenotype(patient, ancestor, is_positive=True): 
                    patient.remove_phenotype(hpo_id, is_positive=True)
                    patient.track_corruption_phenotype(hpo_id, is_positive=True)
                    patient.add_phenotype(ancestor, 'phenotype_corruption', is_positive=True, is_distractor=False)
                    logging.info(f'Replacing positive phenotype {hpo_id} with {ancestor}')

        #negative phenotypes
        patient_neg_phenotypes = list(patient.get_true_phenotypes(is_positive=False)).copy()
        for hpo_id in patient_neg_phenotypes:
            should_corrupt = np.random.binomial(1, p_corrupt_neg) == 1 #true if should corrupt phenotype
            if should_corrupt:
                ancestor = self._get_ontology_ancestor(hpo_id, n_generations=1)

                #if the patient doesn't have an ancestor, we leave it as is
                if ancestor and self.can_add_phenotype(patient, ancestor, is_positive=False): 
                    patient.remove_phenotype(hpo_id, is_positive=False)
                    patient.track_corruption_phenotype(hpo_id, is_positive=False)
                    patient.add_phenotype(ancestor, 'phenotype_corruption', is_positive=False, is_distractor=False)
                    logging.info(f'Replacing negative phenotype {hpo_id} with {ancestor}')

    def initialize_phenotypes(self, patient: Patient, disease: Disease, 
        p_dropout_pos: float = 0, p_dropout_neg: float = 0, p_corrupt_pos: float = 0, p_corrupt_neg: float = 0,
        perform_phenotype_dropout:bool = True, perform_phenotype_corruption:bool = True) -> None:
        '''
        The first step to simulating a patient is to initialize their positive and negative phenotypes.
        We do this by first independently sampling whether the patient expresses each phenotype associated
        with the disease. If so, the phenotype is added to the list of positive phenotypes. Otherwise, it is 
        added to the list of negative phenotypes.

        Then we perform phenotype dropout, which with some probability drops out the pos & neg phenotypes from the
        patient.

        Finally, we perform phenotype corruption in which the phenotypes are replaced with less specific phenotypes 
        that are higher up in the HPO hierarchy.
        '''

        self._sample_initial_phenotypes(patient, disease)
        if perform_phenotype_dropout:
            self._phenotype_dropout(patient, p_dropout_pos, p_dropout_neg)
        if perform_phenotype_corruption:
            self._phenotype_corruption(patient, p_corrupt_pos, p_corrupt_neg)

    def sample_noisy_phenotypes(self, patient: Patient) -> None:
        '''
        Sample random number of noisy phenotypes for the patient
        For each noisy phenotype, sample with bernoulli whether to add to the pos or neg phenotype list.

        Phenotypes added during the noisy phenotype sampling are commonly observed in the general population. These phenotypes
        may appear in the patient's positive or negative phenotypes. 
        '''
        
        logging.info('---- Sampling Noisy Phenotypes ----')
        age = patient.get_age()
        hpo_ids, probabilities = self.age_to_hpo_probs_dict[age]
        n_to_sample = np.random.poisson(config.NOISY_POS_PHEN_SAMPLES_LAMBDA) #it's okay if no noisy phenotypes are added
        
        # add noisy positives by sampling according to prevalence
        sampled_pos_phenotypes = np.random.choice(hpo_ids, n_to_sample, replace=False, p=probabilities)
        for p in sampled_pos_phenotypes:
            if self.can_add_phenotype(patient, p, is_positive=True):
                patient.add_phenotype(p, 'noisy_phenotype', is_positive=True, is_distractor=True)
                logging.info(f'Adding noisy positive phenotype {p}')
    
        #add noisy negatives by sampling randomly
        n_neg_to_sample = np.random.poisson(config.NOISY_NEG_PHEN_SAMPLES_LAMBDA) #it's okay if no noisy phenotypes are added
        sampled_neg_phenotypes = np.random.choice(hpo_ids, n_neg_to_sample, replace=False)
        for p in sampled_neg_phenotypes:
            if self.can_add_phenotype(patient, p, is_positive=False):
                patient.add_phenotype(p, 'noisy_phenotype', is_positive=False, is_distractor=True)
                logging.info(f'Adding noisy negative phenotype {p}')
    

##############################################################
# Check if genes quality for different categories

    def is_tissue_distractor_gene(self, orphanet_id: str, true_gene: str, gene: str) -> bool:
        """
        Gene is a tissue distractor gene if the patient's gene is in the top k nondisease & phenotypically distinct genes 
        with similar tissue expression to the patient's true gene
        NOTE: this module doesn't add phenotypes
        """

         # consider disease genes with no phenotypic overlap
        if true_gene in self.tissue_distractor_disease_genes_dict:
            candidate_dx_gene_distractors, dx_gene_probabilities = self.tissue_distractor_disease_genes_dict[true_gene]

            #filter to only include disease genes that don't have any phenotypic overlap with patient's disease
            candidate_path_phen_irrel_genes = self.pathogenic_but_pheno_irrel_dict[orphanet_id]
            filtered_candidates = [(distractor, prob) for distractor, prob in zip (candidate_dx_gene_distractors, dx_gene_probabilities) if distractor in candidate_path_phen_irrel_genes]
            candidate_distractors, probabilities = zip(*filtered_candidates)

            # merge with nondisease genes
            if patient_gene in self.tissue_distractor_nondisease_genes_dict:
                candidate_nondisease_distractors, nondisease_gene_probabilities = self.tissue_distractor_nondisease_genes_dict[patient_gene]
                
                candidate_distractors = list(candidate_distractors) + list(candidate_nondisease_distractors)
                probabilities = list(probabilities) + list(nondisease_gene_probabilities)

                # sort by probability
                probabilities, candidate_distractors =  zip(*sorted(zip(probabilities, candidate_distractors), reverse=True))

        # only nondisease genes to consider
        elif true_gene in self.tissue_distractor_nondisease_genes_dict:
            candidate_distractors, probabilities = self.tissue_distractor_nondisease_genes_dict[true_gene]     
        else: 
            return False
        
        if gene in candidate_distractors[0:config.N_TISSUE_GENES]: return True
        return False

    def is_pathogenic_phenotypically_irrelevant_gene(self, orphanet_id: str, gene: str) -> bool:
        """
        Gene is a pathogenic but phenotypically irrelevant gene if the patient's orphanet id is in the pathogenic_but_pheno_irrel_dict
        NOTE: this module doesn't add phenotypes
        """
        if orphanet_id in self.pathogenic_but_pheno_irrel_dict:
            candidate_genes, probabilities = self.pathogenic_but_pheno_irrel_dict[orphanet_id]
            if gene in candidate_genes: return True
        return False

    def is_false_positive_gene(self, gene: str) -> bool:
        """
        The Gene is a FP gene if it's in the list of top k common FP genes from FLAGS
        NOTE: this module doesn't add phenotypes
        """
        return gene in self.common_fp_genes

    def is_insufficient_explanatory_gene(self, orphanet_id: str, positive_phenotypes: Dict, gene: str) -> bool:
        """
        Gene is an insufficient explainor if the patient's orphanet id is in the insufficient_explanatory_distractor_dict
        and if the patient's true positive phenotypes contains one of the phenotypes associated with the gene

        Note that because we're limited ourselves to genes that don't cause disease, we're not worried about all phenotype distractors
        falling into this category
        """
        if orphanet_id in self.insufficient_explanatory_distractor_dict:
            candidate_gene_ids = self.insufficient_explanatory_distractor_dict[orphanet_id]

            if gene not in candidate_gene_ids:
                return False

            disease_hpos = self.diseases_dict[orphanet_id].get_phenotype_set(min_freq="Very rare (<4-1%)", max_freq="Obligate (100%)")
            gene_hpos = self.nondisease_gene_to_hpo[gene]

            # get patient phenotypes associated with disease that are also associated with the gene
            intersect = set(disease_hpos).intersection(set(positive_phenotypes.keys())).intersection(set(gene_hpos))

            if len(intersect) > 0:
                return True
        
        return False

    def is_non_syndromic_phenotype_gene(self, orphanet_id: str, positive_phenotypes: Dict, gene:str) -> bool:
        """
        Gene is a non syndromic phenotype gene if the patient's orphanet id is in non_syndromic_phenotype_distractor_dict
        and if at least one phenotype associated with the gene is in the patient's distractor pos phenotypes

        Note that because we're limited ourselves to genes that don't cause disease, we're not worried about all phenotype distractors
        falling into this category
        """
        if orphanet_id in self.non_syndromic_phenotype_distractor_dict:
            candidate_gene_ids = self.non_syndromic_phenotype_distractor_dict[orphanet_id]

            if gene not in candidate_gene_ids:
                return False

            gene_hpos = self.nondisease_gene_to_hpo[gene]
            disease_hpos = self.diseases_dict[orphanet_id].get_phenotype_set(min_freq="Very rare (<4-1%)", max_freq="Obligate (100%)")

            # get patient phenotypes that are not associated with the patient's disease
            distractor_pos_phenotypes = set(positive_phenotypes.keys()).difference(set(disease_hpos))

            # get distractor positive phenotypes associated with gene
            intersect = set(distractor_pos_phenotypes).intersection(set(gene_hpos))

            if len(intersect) > 0:
                return True

        return False

    def is_phenotype_distractor_gene(self, orphanet_id: str, positive_phenotypes:Dict, negative_phenotypes:Dict, gene: str) -> bool:
        if orphanet_id in self.phenotype_distractor_dict:

            for d in self.phenotype_distractor_dict[orphanet_id]:

                # check if gene is in the list of genes assoc. with the distractor disease
                distractor_disease = self.diseases_dict[d]
                if gene in distractor_disease.get_gene_list():

                    patient_disease = self.diseases_dict[orphanet_id]
                    disease_weak_phenotypes = patient_disease.get_phenotype_set(min_freq="Very rare (<4-1%)", max_freq=config.WEAK_PHENOTYPE_THRESH) 
                    distractor_weak_phenotypes = distractor_disease.get_phenotype_set(min_freq="Very rare (<4-1%)", max_freq=config.WEAK_PHENOTYPE_THRESH)
                    distractor_strong_phenotypes = distractor_disease.get_phenotype_set(min_freq=config.STRONG_PHENOTYPE_THRESH)

                    weak_intersection = disease_weak_phenotypes.intersection(distractor_weak_phenotypes)
                    # check if patient's pos phenotypes contain "weak phenotypes" found in distractor & true disease
                    intersect_patient_weak = weak_intersection.intersection(set(positive_phenotypes.keys()))
                   
                    # check if patient's neg phenotypes contain "strong phenotypes" found in distractor disease
                    intersect_patient_strong = set(distractor_strong_phenotypes).intersection(set(negative_phenotypes.keys()))
                    if len(intersect_patient_weak) > 0 and len(intersect_patient_strong) > 0:
                        return True
                        
        return False

    def is_universal_distractor_gene(self, orphanet_id: str, positive_phenotypes: Dict, negative_phenotypes: Dict,  gene: str) -> bool:
        if orphanet_id in self.universal_distractor_dict:

            for d in self.universal_distractor_dict[orphanet_id]:
                
                # check if gene is in the list of genes assoc. with the distractor disease
                distractor_disease = self.diseases_dict[d]
                if gene in distractor_disease.get_gene_list():

                    distractor_obligates = distractor_disease.get_obligate_phenotypes()
                    distractor_excluded = distractor_disease.get_excluded_phenotypes()

                    # make sure that an obligate is in the patient's neg phen or that an excluded is in the patient's pos phen
                    if len(set(distractor_excluded).intersection(set(positive_phenotypes.keys()))) > 0 or \
                        len(set(distractor_obligates).intersection(set(negative_phenotypes.keys()))) > 0:

                        # NOTE: the universal distractor optionally adds phenotypes in the intersect, but this is
                        # optional so we don't check for it here

                        return True
        return False





    

    
