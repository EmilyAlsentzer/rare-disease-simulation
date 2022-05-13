#!/bin/bash


# RUN SIMULATED PATIENTS - ABLATIONS
PATIENT_PATH=/PATH/TO/SIMULATED/PATIENT/ABLATIONS # TODO: Change this to point to the folder where the ablated simulated patients are found (e.g. config.PROJECT_ROOT/simulated_patients/ablations)
for filename in ${PATIENT_PATH}/simulated_patients_*_formatted.jsonl; do

    # Phrank Patient-Gene
    python phrank.py \
    --sim_input $filename \
    --rank_genes_directly

    # Phrank Patient-Disease
    python phrank.py \
    --sim_input $filename 

done
