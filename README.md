# Simulation of undiagnosed patients with novel genetic conditions

[![made-with-python](https://img.shields.io/badge/Made%20with-Python-1f425f.svg)](https://www.python.org/)
[![GitHub Repo stars](https://img.shields.io/github/stars/EmilyAlsentzer/rare-disease-simulation)](https://github.com/EmilyAlsentzer/rare-disease-simulation/stargazers)
[![GitHub Repo forks](https://img.shields.io/github/forks/EmilyAlsentzer/rare-disease-simulation)](https://github.com/EmilyAlsentzer/rare-disease-simulation/network/members)
 [![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](https://opensource.org/licenses/MIT)

## Overview

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



## Example Usage

### :boom: Simulating Patients
To run the simulation pipeline with default parameters, run the following code. This will run the simulation pipeline and output the simulation patients to a jsonl file in the directory specified by `config.SIMULATED_DATA_PATH`.

```
cd simulate_patients
python simulate_patients.py
```

The parameters of the simulation pipeline and their associated descriptions can be found in the `config.py` file. Adjusting these parameters will change the simulated patients output from the `simulate_patients.py` script.



## License
This project is covered under the MIT License.

## Questions
For questions, please leave a Github issue or contact Emily Alsentzer at emilya@mit.edu.


