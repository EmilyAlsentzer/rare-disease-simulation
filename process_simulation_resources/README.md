# Generating the 2015 Timestamped KG

## To combine the HPO-A files with the rest of the 2015 phenolyzer KG:
- Download the 2015 phenolyzer KG from [here](https://github.com/WGLab/phenolyzer/tree/ecec7410729276859b9023a00f20e75c2ce58862).
- Move the  [`phenotype_annotation.tab`](github.com/drseb/HPO-archive/tree/master/2014-2015/2015_week_4/annotations/artefacts) and [`ALL_SOURCES_ALL_FREQUENCIES_phenotype_to_genes.txt`](github.com/drseb/HPO-archive/tree/master/hpo.annotations.monthly/2015-01-28_14-15-03/archive/annotation) files into the folder with the KG (`lib/compiled_database`).
- Run `./generate_compiled_hpo.sh`.

## To convert the genes in the KG into ensembl IDs and older HPO terms to the 2019 HPO:
- Run `convert KG & simulation sources to ensembl & 2019 HPO.ipynb`.
