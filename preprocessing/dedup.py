'''
Remove duplicate genes and non-coding sequences from each tissue dataset. 

While each row in the datasets has a unique ENSEMBL ID, multiple ENSEMBL IDs 
may correlate to the same gene. Duplicates for a given gene are removed by 
retaining the single entry with the highest variance in expression across all 
samples. 

Non-coding sequences are removed by cross-checking all Entrez symbols with 
a list of coding gene symbols from NCBI's Gene Database. 
'''

import numpy as np
import pandas as pd 

filelist = [
    'AdiposeTissue_counts.tsv',
    'AdrenalGland_counts.tsv',
    'Bladder_counts.tsv',
    'Blood_counts.tsv',
    'BloodVessel_counts.tsv',
    'BoneMarrow_counts.tsv',
    'Brain_counts.tsv',
    'Breast_counts.tsv',
    'CervixUteri_counts.tsv',
    'Colon_counts.tsv',
    'Esophagus_counts.tsv',
    'FallopianTube_counts.tsv',
    'Heart_counts.tsv',
    'Kidney_counts.tsv',
    'Liver_counts.tsv',
    'Lung_counts.tsv',
    'Muscle_counts.tsv',
    'Nerve_counts.tsv',
    'Ovary_counts.tsv',
    'Pancreas_counts.tsv',
    'Pituitary_counts.tsv',
    'Prostate_counts.tsv',
    'SalivaryGland_counts.tsv',
    'Skin_counts.tsv',
    'SmallIntestine_counts.tsv',
    'Spleen_counts.tsv',
    'Stomach_counts.tsv',
    'Testis_counts.tsv',
    'Thyroid_counts.tsv',
    'Uterus_counts.tsv',
    'Vagina_counts.tsv'
]

# load list of human protein-coding genes from NCBI's Gene Database
gene_info = pd.read_csv(
    'https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz', 
    compression='gzip', sep='\t')
symbols = gene_info['Symbol'].tolist()

# iterate through each tissue dataset
for f in filelist:
    df = pd.read_csv(f, sep='\t')

    # retain only rows for genes whose symbols are included in NCBI Gene
    df = df[df['Description'].isin(symbols)]

    # double check that number of rows is less than total number of genes
    if (df.shape[0] <= len(symbols)):
        print("Successfully removed non-coding genes!")
        print("Removing duplicates from " + f + "...")

        # identify genes with more than one entry
        vc = df['Description'].value_counts()
        dups = vc[vc > 1]
        dup_list = list(dups.index)

        # store all entries for each duplicated gene in a dictionary
        dup_dict = {}
        for d in dup_list:
            dup_dict[d] = []
        for row in df.iterrows():
            if row[1]['Description'] in dup_list:
                vals = row[1].drop(['Unnamed: 0', 'Name', 'Description'])
                dup_dict[row[1]['Description']].append((row[1]['Name'], vals))

        # initialize a list to store ENSEMBL IDs of entries to be deleted
        remove_ids = []

        # build a dictionary to store highest variance entry for each dup gene
        max_var_dict = {}
        for d in dup_list:
            max_var_dict[d] = None

        # iterate through all duplicated genes and their corresponding entries
        for key in dup_dict:
            for item in dup_dict[key]:
                var = item[1].var()
                # if no entries yet, set the value
                if max_var_dict[key] == None:
                    max_var_dict[key] = (item[0], var)
                # if value already set, compare with variance of new value
                else:
                    if var > max_var_dict[key][1]:
                        # make sure to identify IDs to be removed
                        remove_ids.append(max_var_dict[key][0])
                        max_var_dict[key] = (item[0], var)
                    else:
                        # make sure to identify IDs to be removed
                        remove_ids.append(item[0])

        # remove all low-variance duplicated entries by ENSEMBL ID
        no_dup_df = df[~df['Name'].isin(remove_ids)]

        # remove unneeded index and ENSEMBL ID columns
        # for purposes of analysis, only Entrez gene symbols are required
        no_dup_df = no_dup_df.drop(labels=['Unnamed: 0', 'Name'], axis=1)

        # ensure that the correct number of gene entries remain, once all 
        # duplicates have been removed
        correct_num = (len(no_dup_df.index) == 
            len(df.index) - dups.sum() + len(dups.index))

        # if no errors, create new tsv file for deduped dataset
        if correct_num:
            print("Final gene count correct")
            tis_str = f.replace('.tsv','')
            no_dup_df.to_csv(tis_str + '_deduped.tsv', sep='\t', index=False)
        else:
            print("Final gene count incorrect")

    # do not proceed with deduplication if failed to remove all non-coding seq
    else:
        print("Number of genes is NOT less than 19907, check again")