#!/bin/bash

# change this to the path to the simulated patient jsonl file that you wish to run phrank on
SIM_INPUT='simulated_patients_formatted.jsonl'

python phrank.py \
--simulated \
--sim_input $SIM_INPUT \
--rank_genes_directly

python phrank.py \
--simulated \
--sim_input $SIM_INPUT 

python phrank.py \
--simulated \
--sim_input $SIM_INPUT \
--rank_genes_directly \
--use_disease_annotations


