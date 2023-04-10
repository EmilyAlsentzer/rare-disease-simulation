
import sys
import argparse
from pathlib import Path
sys.path.insert(0, '../') # add config to path
import config
from simulation_pipeline.utils.util import read_simulated_patients


def vcf_header(reference):
    """
    :return: header for the VCF file (same everytime)
    """
    if reference == 'hg19': ref_str = '##reference=GRCh37/hg19'
    elif reference == 'hg38': ref_str = '##reference=GRCh38/hg38'
    else: raise Exception
    vcf_lines = ['##fileformat=VCFv4.2', ref_str] + \
                ['##contig=<ID='+chrom+'>' for chrom in [str(i) for i in range(1, 23)] + ['X', 'Y', 'MT']] + \
                ['##INFO=<ID=ensembl_id,Number=.,Type=String,Description="Stable Ensembl Gene Identifier v106">',
                 '##INFO=<ID=name,Number=.,Type=String,Description="HGNC Gene Name">',
                 '##INFO=<ID=alias,Number=.,Type=String,Description="former HGNC symbol">',
                 '##INFO=<ID=deprecated_ensembl_id,Number=.,Type=String,Description="Ensembl ID that no longer maps to version 108">',
                 '##INFO=<ID=typo,Number=.,Type=String,Description="Mistyped gene name in manually curated clinical report">',
                 '##INFO=<ID=CADD,Number=1,Type=Float,Description="Phred-scaled CADD score">',
                 '##INFO=<ID=PolyPhen2,Number=1,Type=Float,Description="VEP-assigned PolyPhen2 score">',
                 '##INFO=<ID=MAF,Number=1,Type=Float,Description="Maximum gnomAD subpopulation AF or TOPMed AF">',
                 '##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequences assigned by Ensembl VEP">',
                 '##INFO=<ID=biotype,Number=.,Type=String,Description="Ensembl Gene Biotype">',
                 '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
                 '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">',
                 '\t'.join(['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'proband'])
                 ]
    return '\n'.join(vcf_lines)+'\n'


########################################################################################################

def gene_variant_map(gene_to_variant_mapping):
    """
    :param gene_to_variant_mapping: file name with gene name as the first column followed by a complete VCF line
    :return: dictionary of gene_name -> variant
    """

    gene_to_variant = {}  # gene_name -> variant line (no newline) in VCF format
    with open(gene_to_variant_mapping) as gnmap_handle:
        for gline in gnmap_handle:
            gene = gline.strip().split('\t')[0]
            vcf_line = '\t'.join(gline.strip().split('\t')[1:])
            gene_to_variant[gene] = vcf_line
    return gene_to_variant


########################################################################################################

def create_vcf(patient, candidate_gene_list, output_vcf, gn_to_var, reference, genes_with_no_variants):
    """
    :param candidate_gene_list: genes we are interested in
    :param output_vcf: write out properly formatted VCF file
    :return: running list of genes that don't map to a variant
    """

    sorted_vars = []
    for cgene in candidate_gene_list:
        try:
            var = gn_to_var[cgene]
            chrom, pos, _, ref, alt = var.split('\t')[:5]
            sorted_vars.append(((23 if chrom == 'X' else (24 if chrom == 'Y' else (25 if chrom == 'MT' else int(chrom))),
                                 int(pos),
                                 ref,
                                 alt), var))
        except:
            sys.stderr.write('Could not find a variant for '+cgene+'\n')
            genes_with_no_variants.append(cgene)

    outhandle = open(output_vcf, 'w')
    outhandle.write(vcf_header(reference) + '\n'.join([a[1]+'\tGT:DP\t0/1:40' for a in sorted(sorted_vars)])+'\n')
    outhandle.close()
    sys.stderr.write('VCF output in ' + output_vcf + '\n')
    return genes_with_no_variants

def create_vcfs(patients, patient_type, reference):
    if reference == 'hg38':
        GENE_TO_VARIANT_MAP_LOC = 'gene_to_variant_hg38.tsv'
    elif reference == 'hg19':
        GENE_TO_VARIANT_MAP_LOC = 'gene_to_variant_hg19.tsv'
    else: raise Exception
    gn_to_var = gene_variant_map(GENE_TO_VARIANT_MAP_LOC)

    genes_with_no_variants = []
    for patient in patients:
        cand_genes = patient['all_candidate_genes']
        patient_id = patient['id']
        print(f'creating VCF for patient {patient_id}...')
        output_vcf = str(config.GENERATED_VCF_PATH / patient_type /  f'{patient_type}_{patient_id}_{reference}.vcf' )
        genes_with_no_variants = create_vcf(patient, cand_genes, output_vcf, gn_to_var, reference, genes_with_no_variants)
    print(f'There are {len(set(genes_with_no_variants))} unique genes with no variants: {set(genes_with_no_variants)}')

def main():
    parser = argparse.ArgumentParser(description='Create VCFs')
    parser.add_argument('--reference', type=str, default='hg19')
    parser.add_argument('--patients', type=str, default='simulated_patients_formatted.jsonl')
    args = parser.parse_args()

    patients = read_simulated_patients(config.SIMULATED_DATA_PATH / args.patients)
    create_vcfs(patients, 'simulated', args.reference)

if __name__ == "__main__":
    main()