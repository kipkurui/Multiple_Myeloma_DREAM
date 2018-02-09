import pandas as pd
from multiprocessing import Pool, cpu_count
import numpy as np

import pysam

num_partitions = cpu_count()

num_cores = cpu_count()


chalenge_data = pd.read_csv('/mnt/lustre/users/ckibet/MM_DREAM/Data/Clinical_Data/sc1_Training_ClinAnnotations.csv')


data_path = '/mnt/lustre/users/ckibet/MM_DREAM/Data/Genomic_Data/MMRF_IA9_CelgeneProcessed/MuTect2_SnpSift_Annotated_vcfs/'


def parallelize_dataframe(df, func):
    df_split = np.array_split(df, num_partitions)
    pool = Pool(num_cores)
    df = pd.concat(pool.map(func, df_split))
    pool.close()
    pool.join()
    return df

def run_nonsense(bed_df):
    
    return bed_df.apply(lambda row: count_nonsense(row))


def count_nonsense(vcf_r):
    vcf_r = data_path+vcf_r
    count = 0
    try:
        test = pysam.VariantFile(vcf_r, mode='r')
        for i in test.fetch():
            for t in i.info.get('ANN'):
                if 'nonsense_mediated_decay' in t:
                    count = count + 1
    except TypeError:
        count = None
    
    return count

nonsense = parallelize_dataframe(chalenge_data['WES_mutationFileMutect'], run_nonsense)

nonsense.to_csv('/mnt/lustre/users/ckibet/MM_DREAM/nonesense.csv',index=False)