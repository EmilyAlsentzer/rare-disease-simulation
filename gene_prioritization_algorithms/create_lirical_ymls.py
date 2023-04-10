import sys
import argparse
sys.path.insert(0, '../') 
sys.path.insert(0, '../../') # add config to path

import config
from simulation_pipeline.utils.util import read_simulated_patients


def create_yml(args, patients, filepath, exomiser_filepath, patient_type, orphanet='true', use_negated_hpo=True, reference='hg19', transcriptdb='ucsc'):
    print('creating ymls....')
    print(f'Reference: {reference}')
    print(f'There are {len(patients)} patients')
    for p in patients:
        positive_hpo = p['positive_phenotypes']
        negative_hpo = p['negative_phenotypes']
        pid = p['id']
        use_global = 'true' if args.use_global else 'false'
        with open(str(filepath / f'{patient_type}_{pid}_orphanet={orphanet}_negatedhpo={use_negated_hpo}_vcf={args.use_vcf}_global={use_global}_reference={reference}_transcriptdb={transcriptdb}.yml' ), 'w') as f:
            f.write('---\n')
            f.write('analysis:\n')
            if args.use_vcf: 
                f.write(f'  genomeAssembly: {reference}\n')
                f.write(f'  vcf: {config.GENERATED_VCF_PATH}/{patient_type}/{patient_type}_{pid}_{reference}.vcf\n')
                f.write(f'  exomiser: {exomiser_filepath}/2209_{reference}/\n')
                f.write(f'  global: {use_global}\n') #If the YAML file contains the line global: true then it will not discard candidate diseases with no known disease gene or candidates for which no predicted pathogenic variant was found in the VCF.
                f.write(f'  transcriptdb: {transcriptdb}\n') 
            f.write('  threshold: 0.0\n')
            f.write('  tsv: true\n')
            f.write(f'  orphanet: {orphanet}\n')
            f.write(f'outdir: /home/ema30/zaklab/udn_data/LIRICAL/results/{patient_type}\n')
            f.write(f'prefix: patient_{pid}_{patient_type}_orphanet={orphanet}_negatedhpo={use_negated_hpo}_vcf={args.use_vcf}_global={use_global}_{reference}\n')
            f.write(f'hpoIds: {positive_hpo}\n')
            if use_negated_hpo: f.write(f'negatedHpoIds: {negative_hpo}')


def main():
    parser = argparse.ArgumentParser(description='Create LIRICAL ymls')
    parser.add_argument('--simulated', action='store_true')
    parser.add_argument('--reference', type=str, default='hg19')
    parser.add_argument('--use_global', action='store_true')
    parser.add_argument('--use_vcf', action='store_true')
    parser.add_argument('--patients', type=str, default='simulated_patients_formatted.jsonl')
    parser.add_argument('--output_filepath', type=str, default=config.PROJECT_ROOT / 'LIRICAL/yml/simulated')
    parser.add_argument('--exomiser_filepath', type=str, default= config.PROJECT_ROOT / 'Exomiser/data/')
    args = parser.parse_args()

    sim_fname = args.sim_input.split('.jsonl')[0]
    patients = read_simulated_patients(config.SIMULATED_DATA_PATH / args.patients)
    create_yml(args, patients, args.output_filepath, args.exomiser_filepath, 'simulated', reference=args.reference)


if __name__ == "__main__":
    main()
