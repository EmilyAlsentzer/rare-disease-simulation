# Simulation of undiagnosed patients with novel genetic conditions

[![made-with-python](https://img.shields.io/badge/Made%20with-Python-1f425f.svg)](https://www.python.org/)
[![GitHub Repo stars](https://img.shields.io/github/stars/EmilyAlsentzer/rare-disease-simulation)](https://github.com/EmilyAlsentzer/rare-disease-simulation/stargazers)
[![GitHub Repo forks](https://img.shields.io/github/forks/EmilyAlsentzer/rare-disease-simulation)](https://github.com/EmilyAlsentzer/rare-disease-simulation/network/members)
 [![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](https://opensource.org/licenses/MIT)

## Overview
We present a computational pipeline to simulate realistic undiagnosed rare disease patients that can be used to evaluate gene prioritization tools. Each simulated patient is represented by sets of candidate disease-causing genes and standardized phenotype terms. 

### Pipeline Featues
- :high_brightness:  **Realistic simulation of patient phenotypes and candidate genes:** We provide a taxonomy of categories of “distractor” genes that do not cause the patient’s presenting syndrome yet would be considered plausible candidates during the clinical process. We then introduce a simulation framework that jointly samples genes and phenotypes according to these categories to simulate nontrivial and realistic patients. 
- :high_brightness: **Modelling novel genetic conditions:** To simulate patients with novel genetic conditions, we curate a knowledge graph (KG) of known gene–disease and gene–phenotype annotations that is time-stamped to 2015. This enables us to define post-2015 medical genetics discoveries as novel with respect to our KG. We manually time-stamp each disease and disease–gene association according to the date of the Pubmed article that reported the discovery and use these time-stamps to annotate each patient according to each of the disease-gene novelty categories below. 




### Ways to use this simulation pipeline
- Probing Gene Prioritization Models: We provide our framework and simulated patients as a public resource to enable developers to internally validate and improve their tools by separately evaluating on simulated patients across novelty categories and distractor gene categories. 
- Training Machine-Learning based Gene Prioritization Models: We also suspect that our dataset of simulated patients that span diverse genetic disorders and reflect realistic clinical processes and imprecision can serve as invaluable training data for machine learning models for rare disease diagnosis.



## Environment Setup
This codebase leverages Python and many associated packages. To create a conda environment with all of the required packages, ensure that [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) is installed and run the following:

```
conda env create -f environment.yml
conda activate rare-dx-env
```

Alternatively, to install via pip, run:
```
pip install -r requirements.txt
```



## Installation
After the conda environment is created and activated, install the github repo with the following:

```
git clone https://github.com/EmilyAlsentzer/rare-disease-simulation
cd rare-disease-simulation
pip install -e .
```

## Set Up
Before you run any code in the repository, you will need to (1) download the required data and (2) specify the path to the project directory.

### Download Data
First, download the data from Harvard Dataverse

### Setting up config file
Go to `config.py` and set the project directory (`PROJECT_ROOT`) to be the path to the data folder downloaded in the previous step. 



## Usage

### :one: Simulate Realistic Rare Disease Patients
To run the simulation pipeline with default parameters, run the following code. This will run the simulation pipeline and output the simulation patients to a jsonl file named `simulated_patients.jsonl` by default in the directory specified by `config.SIMULATED_DATA_PATH`.

```
cd simulate_patients
python simulate_patients.py
```

The parameters of the simulation pipeline and their associated descriptions can be found in the `config.py` file. Modify these parameters to change the number of patients produced, the frequency of sampling each distractor gene module, etc.

### :two: Annotate patients with their disease-gene novelty category
To annotate the simulated patients from the prior step with their disease-gene novelty category, run:

```
python label_patients_with_novelty_categories.py --input simulated_patients.jsonl
```

### :three: Evaluate gene prioritization model performance on the patients

### :four: Create simulated patients with various pipeline components removed


## License
This project is covered under the MIT License.

## Questions
For questions, please leave a Github issue or contact Emily Alsentzer at emilya@mit.edu.

