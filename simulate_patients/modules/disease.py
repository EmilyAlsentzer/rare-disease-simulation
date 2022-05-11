import typing
from typing import Dict, List
import logging

phenotype_frequency_dict = {"Obligate (100%)":1.0, 
"Very frequent (99-80%)":0.895,
"Frequent (79-30%)":0.545,
"Occasional (29-5%)":0.17,
"Very rare (<4-1%)":0.025,
"Excluded (0%)":0.0}

class Disease():

    def __init__(self, orphanet_id: str, disease_name: str, gene_list: List, \
        phenotype_dict: Dict[str,str], age_ranges: List = []):

        self.gene_list = gene_list
        self.phenotype_dict = phenotype_dict #{hpo_id -> freq}
        self.age_range_list = age_ranges
        self.orphanet_id = orphanet_id
        self.disease_name = disease_name


    def get_gene_list(self) -> set:
        '''
        Return a list containing the list of genes known to cause the disease
        '''
        return self.gene_list

    def get_age_range_set(self) -> set:
        '''
        Returns age ranges associated with the disease as a set
        '''
        return set(self.age_range_list)

    def get_phenotypes(self, min_freq='Excluded (0%)', \
        max_freq='Obligate (100%)'): 
        '''
        Return all phenotypes between min-max frequency as a dict {hpo_id -> freq}
        '''

        min_freq = phenotype_frequency_dict[min_freq]
        max_freq = phenotype_frequency_dict[max_freq]
        if min_freq == 0 and max_freq == 1:
            return self.phenotype_dict

        # if min is supplied, get min & above
        # if max is supplied, get max & below
        phenotypes_in_range_dict = {}
        for hpo_id, freq_str in self.phenotype_dict.items():
            freq = phenotype_frequency_dict[freq_str]
            if freq >= min_freq and freq <= max_freq:
                phenotypes_in_range_dict[hpo_id] = freq_str
        return phenotypes_in_range_dict
        


    def get_phenotype_set(self, min_freq='Excluded (0%)', \
        max_freq='Obligate (100%)'):
        '''
        Return all phenotypes between min-max frequency as a set
        '''
        phenotypes_in_range_dict = self.get_phenotypes(min_freq, max_freq)
        return set(phenotypes_in_range_dict.keys())

    def get_obligate_phenotypes(self):
        '''
        Returns a set of all obligate phenotypes (i.e. those phenotypes always associated with the disease)
        '''
        return set(hpo_id for hpo_id, freq in self.phenotype_dict.items() \
            if phenotype_frequency_dict[freq] == 1.0)

    def get_excluded_phenotypes(self):
        '''
        Returns a set of all excluded phenotypes (i.e. those phenotypes never associated with the disease)
        '''
        return set(hpo_id for hpo_id, freq in self.phenotype_dict.items() \
            if phenotype_frequency_dict[freq] == 0.0)

    def __str__(self):
        '''
        string representation of the disease
        '''
        s = '''Orphanet id: {}, Disease name: {}, Phenotypes: {}, Genes: {}''' \
            .format(self.orphanet_id, self.disease_name, self.phenotype_dict, self.gene_list)
        return s
