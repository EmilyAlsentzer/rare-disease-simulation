import sys
import argparse
import pandas as pd
sys.path.insert(0, '../..') # add config to path
sys.path.insert(0, '..') # add util to path

import config
from simulation_pipeline.utils.util import read_simulated_patients

'''
usage: java -jar exomizer-cli-6.0.0.jar [...]
    --batch-file <file>                  Path to batch file. This should
                                         contain a list of fully qualified
                                         path names for the settings files
                                         you wish to process. There should
                                         be one file name on each line.
    --candidate-gene <arg>               Gene symbol of known or suspected
                                         gene association e.g. FGFR2
 -D,--disease-id <arg>                   OMIM ID for disease being
                                         sequenced. e.g. OMIM:101600
 -E,--hiphive-params <type>              Comma separated list of optional
                                         parameters for hiphive: human,
                                         mouse, fish, ppi. e.g.
                                         --hiphive-params=human or
                                         --hiphive-params=human,mouse,ppi
 -F,--max-freq <arg>                     Maximum frequency threshold for
                                         variants to be retained. e.g.
                                         100.00 will retain all variants.
                                         Default: 100.00
 -f,--out-format <type>                  Comma separated list of format
                                         options: HTML, VCF, TAB-GENE or
                                         TAB-VARIANT,. Defaults to HTML if
                                         not specified. e.g.
                                         --out-format=TAB-VARIANT or
                                         --out-format=TAB-GENE,TAB-VARIANT
                                         ,HTML,VCF
    --full-analysis <true/false>         Run the analysis such that all
                                         variants are run through all
                                         filters. This will take longer,
                                         but give more complete results.
                                         Default is false
    --genes-to-keep <Entrez geneId>      Comma separated list of seed
                                         genes (Entrez gene IDs) for
                                         filtering
 -h,--help                               Shows this help
 -H,--help                               Shows this help
    --hpo-ids <HPO ID>                   Comma separated list of HPO IDs
                                         for the sample being sequenced
                                         e.g.
                                         HP:0000407,HP:0009830,HP:0002858
 -I,--inheritance-mode <arg>             Filter variants for inheritance
                                         pattern (AR, AD, X)
    --num-genes <arg>                    Number of genes to show in output
 -o,--out-file <arg>                     name of out file. Will default to
                                         vcf-filename-exomiser-results.htm
                                         l
 -p,--ped <file>                         Path to pedigree (ped) file.
                                         Required if the vcf file is for a
                                         family.
 -P,--keep-non-pathogenic <true/false>   Keep the predicted non-pathogenic
                                         variants that are normally
                                         removed by default. These are
                                         defined as syonymous, intergenic,
                                         intronic, upstream, downstream or
                                         intronic ncRNA variants. This
                                         setting can optionally take a
                                         true/false argument. Not
                                         including the argument is
                                         equivalent to specifying 'true'.
    --prioritiser <name>                 Name of the prioritiser used to
                                         score the genes. Can be one of:
                                         hiphive, exomewalker, phenix,
                                         phive, uber-pheno or none. e.g.
                                         --prioritiser=none
 -Q,--min-qual <arg>                     Mimimum quality threshold for
                                         variants as specifed in VCF
                                         'QUAL' column.  Default: 0
 -R,--restrict-interval <arg>            Restrict to region/interval
                                         (e.g., chr2:12345-67890)
    --remove-dbsnp <true/false>          Filter out all variants with an
                                         entry in dbSNP/ESP (regardless of
                                         frequency).
 -S,--seed-genes <Entrez geneId>         Comma separated list of seed
                                         genes (Entrez gene IDs) for
                                         random walk
    --settings-file <file>               Path to settings file. Any
                                         settings specified in the file
                                         will be overidden by parameters
                                         added on the command-line.
 -T,--keep-off-target <true/false>       Keep the off-target variants that
                                         are normally removed by default.
                                         These are defined as intergenic,
                                         intronic, upstream, downstream or
                                         intronic ncRNA variants. This
                                         setting can optionally take a
                                         true/false argument. Not
                                         including the argument is
                                         equivalent to specifying 'true'.
 -v,--vcf <file>                         Path to VCF file with mutations
                                         to be analyzed. Can be either for
                                         an individual or a family.
                                         '''

def create_settings(args, patients, filepath, patient_type, prioritiser):
    print('creating settings....')
    print(f'There are {len(patients)} patients')
    for p in patients:
        positive_hpo = p['positive_phenotypes']
        pid = p['id']

        with open(str(filepath / f'{pid}_prioritizer={prioritiser}.settings' ), 'w') as f:
            f.write(f'#settings for {pid}\n')
            f.write(f"vcf={config.GENERATED_VCF_PATH}/{patient_type}/{patient_type}_{pid}_hg19.vcf\n") #NOTE: assumes hg19
            f.write(f"prioritiser={prioritiser}\n")
            f.write("ped=\n")
            f.write("max-freq=100.00\n")
            f.write("restrict-interval=\n")
            f.write("min-qual=0\n")
            f.write("keep-non-pathogenic=true\n")
            f.write("remove-dbsnp=false\n")
            f.write("keep-off-target=true\n")
            f.write("full-analysis=false\n")
            f.write("candidate-gene=\n")
            f.write("seed-genes=\n")
            f.write("disease-id=\n")
            f.write("inheritance-mode=\n")
            f.write(f"hpo-ids={','.join(positive_hpo)}\n")
            f.write("num-genes=0\n")
            f.write(f"out-file=results/{patient_type}/{pid}_exomiser_results_prioritizer={prioritiser}\n")
            f.write(f"out-format=TSV-GENE,TSV-VARIANT,HTML\n")
            
                
def main():
    parser = argparse.ArgumentParser(description='Create exomiser settings files')
    parser.add_argument('--prioritiser', type=str, default='phenix')
    parser.add_argument('--patients', type=str, default='simulated_patients_formatted.jsonl')
    parser.add_argument('--output_filepath', type=str, default=config.PROJECT_ROOT / 'exomiser_6.0.0/exomiser-cli-6.0.0/settings/simulated')
    args = parser.parse_args()
    
    sim_fname = args.sim_input.split('.jsonl')[0]
    patients = read_simulated_patients(config.SIMULATED_DATA_PATH / args.patients)
    
    create_settings(args, patients, args.output_filepath, 'simulated', args.prioritiser)


if __name__ == "__main__":
    main()
