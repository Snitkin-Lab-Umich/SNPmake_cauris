import os
import argparse
import pandas as pd
import gffutils as gff


def create_gene_family_dict(pandf,gff_db,reference_name):
    # create a dictionary of gene families to gene names + coordinates
    # only include gene families that are actually present in the reference genome
    gf_dict = {}
    for gene in pandf['Gene']:
        gene_data = pandf.loc[pandf['Gene'] == gene, reference_name]
        if gene_data.iloc[0] != '':
            gene_data_list = gene_data.iloc[0].split(';')
            # remove any entries with 'refound', and ensure they either start with 'cds-' or end with '.cds'
            gene_data_list = [x for x in gene_data_list if 'refound' not in x and (x.startswith('cds-') or x.endswith('.cds'))]
            if len(gene_data_list) > 0:
                gene_name = gene_data_list[0]
                gene_info = gff_db[gene_name]
                gf_dict[gene] = (gene_name, gene_info.seqid, gene_info.start, gene_info.end)
    return gf_dict


def gene_coverage(gene_cov_df, start, end):
    """start/end must be 0-based, half-open coordinates."""
    if end <= start:
        raise ValueError("Gene interval must have positive length")

    overlap_lengths = (
        gene_cov_df['end'].clip(upper=end)
        - gene_cov_df['start'].clip(lower=start)
    )

    covered_depth = (overlap_lengths * gene_cov_df['coverage']).sum()
    return covered_depth / (end - start)


def make_coverage_matrix(input_dir,panaroo_matrix,reference_gff,reference_name,output_file):
    # read in the panaroo matrix
    pandf = pd.read_csv(panaroo_matrix,sep=',',keep_default_na=False)
    # create a new dataframe to hold the coverage data, with the same rows and columns as the panaroo matrix
    summary_matrix = pd.DataFrame(float('nan'), index=pandf.index, columns=pandf.columns)
    summary_matrix['Gene'] = pandf['Gene']
    # use the reference gff, the reference name, and the panaroo matrix to create a dictionary of gene families to gene names + coordinates
    gff_db = gff.create_db(reference_gff,dbfn=":memory:",force=True,keep_order=False,merge_strategy="create_unique",sort_attribute_values=True,from_string=False)
    gf_dict = create_gene_family_dict(pandf,gff_db,reference_name)
    # drop any gene families that are not present in gf_dict from the summary matrix
    summary_matrix = summary_matrix[summary_matrix['Gene'].isin(gf_dict.keys())]
    # add a reference_gene column to the summary matrix
    summary_matrix['reference_gene'] = summary_matrix['Gene'].apply(lambda x: gf_dict[x][0] if x in gf_dict else 'NA')
    # iterate through the input directory and read in each bedcov file
    for sample in summary_matrix.columns:
        if sample == 'Gene' or sample == 'reference_gene':
            continue
        bedcov_file = os.path.join(input_dir,'filtered_vcf',sample,f'{sample}_GenomeCoverage.bedgraph.gz')
        if not os.path.exists(bedcov_file):
            # drop this sample from the summary matrix
            summary_matrix.drop(columns=[sample],inplace=True)
            continue
        cov_df = pd.read_csv(bedcov_file,sep='\t',header=None,names=['chrom','start','end','coverage'],compression='gzip')
        for gene in summary_matrix['Gene']:
            gene_name, gene_chrom, gene_start, gene_end = gf_dict[gene]
            # correct for 0-based coordinates in bedcov file
            gene_start -= 1
            gene_cov_df = cov_df[(cov_df['chrom'] == gene_chrom) & (cov_df['start'] < gene_end) & (cov_df['end'] > gene_start)]
            if gene_cov_df.empty:
                # set coverage to -1 in the case of an error
                summary_matrix.loc[summary_matrix['Gene'] == gene, sample] = -1
            else:
                cov = gene_coverage(gene_cov_df,gene_start,gene_end)
                summary_matrix.loc[summary_matrix['Gene'] == gene, sample] = cov
    # write the coverage dataframe to a tsv file
    summary_matrix.to_csv(output_file,sep='\t',index=True)


def main():
    # define all args
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--input','-i',type=str,
        help='''Provide the path to input directory containing bedtools genomecov files.''',
        required=True
        )
    parser.add_argument(
        '--panaroo','-p',type=str,
        help='''Provide the path to the panaroo matrix. This should have gene names rather than numbers for the cells.''',
        required=True
        )
    parser.add_argument(
        '--reference_gff','-g',type=str,
        help='''Provide the path to the reference GFF file. This should match the reference genome used during variant calling.''',
        required=True
        )
    parser.add_argument(
        '--reference_name','-n',type=str,
        help='''Provide the name of the reference genome.''',
        required=True
        )
    parser.add_argument(
        '--output','-o',type=str,
        help='''Provide the path to the output summary .tsv file.''',
        required=True
        )
    args = parser.parse_args()
    make_coverage_matrix(args.input,args.panaroo,args.reference_gff,args.reference_name,args.output)
    

if __name__ == "__main__":
    main()




