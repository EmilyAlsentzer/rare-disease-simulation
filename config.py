from pathlib import Path

PROJECT_ROOT = Path('/home/ema30/zaklab/udn_data/data') # Change this to be the path to the data directory downloaded from Harvard Dataverse

#PROJECT_ROOT = Path('/PATH/TO/DATA/DOWNLOAD') # Change this to be the path to the data directory downloaded from Harvard Dataverse

SIMULATED_DATA_PATH = PROJECT_ROOT / 'simulated_patients' # path to simulated patients
KNOWLEDGE_GRAPH_PATH = PROJECT_ROOT / 'knowledge_graph' # path to knowledge graph
UDN_PATH = PROJECT_ROOT / 'processed_udn'
SIMULATOR_DIR = PROJECT_ROOT / 'simulation_resources' / 'simulator_preprocess'

############
# SIMULATION SOURCES
ORPHANET_PATH = Path(PROJECT_ROOT / 'orphanet' )
HPO_FILE = Path(PROJECT_ROOT / 'hpo' / '2019'/ 'hp.obo')
HPOA_PATH = Path(PROJECT_ROOT / 'simulation_resources' / 'hpoa')
GTEX_PATH = Path(PROJECT_ROOT / 'simulation_resources' / 'gtex' )

DISGENET_PATH = Path(PROJECT_ROOT / 'simulation_resources' /'disgenenet')
FLAGS_PATH = Path(PROJECT_ROOT / 'simulation_resources'/ 'flags') # FLAGS, frequently mutated genes in public exomes - https://bmcmedgenomics.biomedcentral.com/articles/10.1186/s12920-014-0064-y
CLAIMS_PATH = Path(PROJECT_ROOT / 'simulation_resources' / 'claims_database') 

############
# PREPROCESS SOURCES
PREPROCESS_PATH = Path(PROJECT_ROOT / 'preprocess') # path to intermediate files used to mapping to ensembl ids
HGNC_PATH = Path(PROJECT_ROOT / 'preprocess'/ 'hgnc' )
BIOMART_PATH = Path(PROJECT_ROOT / 'preprocess'/ 'biomart')
NCBI_PATH = Path(PROJECT_ROOT / 'preprocess'/  'ncbi')

############
# PARAMETERS FOR SIMULATING PATIENTS
PATIENTS_PER_DISEASE = 20 # Number of patients to simulate per Orphanet Disease
N_DISTRACTORS_LAMBDA = 13 
PROB_DROPOUT_POS = 0.2  
PROB_DROPOUT_NEG = 0.9
PROB_CORRUPTION_POS = 0.25 
PROB_CORRUPTION_NEG = 0.25 

# PHENOTYPE DRIVEN DISTRACTORS
STRONG_PHENOTYPE_THRESH = "Very frequent (99-80%)"
WEAK_PHENOTYPE_THRESH = "Occasional (29-5%)"
STRONG_PHENOTYPES_LAMBDA = 4 
WEAK_PHENOTYPES_LAMBDA = 4 

# UNIVERSAL DISTRACTORS
OBLIGATE_PHENOTYPES_LAMBDA = 4 
EXCLUDED_PHENOTYPES_LAMBDA = 4 
INTERSECT_PHENOTYPES_LAMBDA = 4 

# NON_SYNDROMIC PHENOTYPES
NON_SYNDROMIC_PHENOTYPES_LAMBDA = 4

# INSUFFICIENT EXPLANATORY PHENOTYES
INSUFFICIENT_EXPLAINER_PHENOTYPE_LAMBDA = 4

# TISSUE DISTRACTORS
# used to define the number of candidate tissue distractors for a given gene
# we only sample from the top N_TISSUE_GENES genes given the heavy tail of gtex expression cosine similarity
N_TISSUE_GENES = 100

# COMMON FALSE POSITIVES
# we only sample from the top common FP genes 
N_COMMON_FP_GENES = 100

# NOISY PHENOTYPES
PROB_NOISY_POSITIVE = 0.5
NOISY_POS_PHEN_SAMPLES_LAMBDA = 5 
NOISY_NEG_PHEN_SAMPLES_LAMBDA = 3

# PROBABILITY OF DIFFERENT DISTRACTOR GENE MODULES
NON_SYNDROM_PHEN_PROB = 0.09
COMMON_FP_PROB = 0.03
TISSUE_DIST_PROB = 0.08
PATH_PHEN_PROB = 0.42
INSUFF_EXPLAIN_PROB = 0.05 
UNIVERSAL_DIST_PROB = 0.03
PHENO_DIST_PROB = 0.30

# allowed % overlap between disease & phenotypes caused by insufficient explainer gene
INSUFF_EXPLAINER_HPO_OVERLAP = 0.20

# max number of times to try to use a gene module to add gene  
MAX_ADD_GENE_ATTEMPTS = 100 
MAX_ADD_ANY_DISTRACTOR_ATTEMPTS = 1000 

