#!/bin/bash
perl integrate_hpo_annotation.pl
cat ALL_SOURCES_ALL_FREQUENCIES_phenotype_to_genes.txt | awk  '{a[$1]=$2; b[$1]=b[$1]?b[$1]","$4:$4;} END{for (i in a) print i, a[i], b[i];}' OFS='\t' FS='\t' | sort -k2,2  > DB_COMPILED_HPO_PHENOTYPE_GENE
