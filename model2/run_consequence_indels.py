import pandas as pd
from multiprocessing import Pool, cpu_count
import numpy as np

import pysam

num_partitions = cpu_count() - 2

num_cores = cpu_count() - 2

#TODO: Make sure this code is varsatile; in that we can easily change the data
chalenge_data = pd.read_csv('../../data/Clinical_Data/sc1_Training_ClinAnnotations.csv')
data_path = '../../data/Genomic_Data/MMRF_IA9_CelgeneProcessed/Strelka_SnpSift_Annotated_vcfs/indels/'

def parallelize_dataframe(df, func):
    df_split = np.array_split(df, num_partitions)
    pool = Pool(num_cores)
    df = pd.concat(pool.map(func, df_split))
    pool.close()
    pool.join()
    return df


def run_consequence(bed_df):
    """
    
    """
    consequence = bed_df.apply(lambda row: count_consequence(row[0]),axis=1)
    
    return consequence


def count_consequence(vcf_r):
    vcf_r = data_path+vcf_r
    count = 0
    try:
        test = pysam.VariantFile(vcf_r, mode='r')
        for i in test.fetch():
            for t in i.info.get('ANN'):
                if variant in t:
                    count = count + 1
    except TypeError:
        count = None
    
    return count

variants = ['upstream_gene_variant',
 'intron_variant',
 'downstream_gene_variant',
 'non_coding_transcript_exon_variant',
 '3_prime_UTR_variant',
 '5_prime_UTR_variant',
 'intergenic_region',
 'frameshift_variant',
 'sequence_feature',
 'splice_region_variant&intron_variant',
 'TF_binding_site_variant',
 'splice_region_variant',
 'disruptive_inframe_deletion',
 'splice_region_variant&disruptive_inframe_deletion',
 'non_coding_transcript_variant']

consequence_data = pd.DataFrame()

for variant in variants:
    consequence_data[variant] = parallelize_dataframe(pd.DataFrame(chalenge_data['WES_mutationFileStrelkaIndel'].head(736)), run_consequence)

consequence_data = ((consequence_data.T/consequence_data.T.sum())*100).T

consequence_data.insert(0,'Patient',chalenge_data['Patient'])
consequence_data.to_csv('consequence_data_percentage_indels.csv', index=False)