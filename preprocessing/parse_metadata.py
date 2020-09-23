'''
    Rename and organize metadata categories for each sample.
'''

import numpy as np
import pandas as pd

# list of deduplicated tissue rna-seq datasets
filelist = [
    'AdiposeTissue_counts_deduped.tsv',
    'AdrenalGland_counts_deduped.tsv',
    'Bladder_counts_deduped.tsv',
    'Blood_counts_deduped.tsv',
    'BloodVessel_counts_deduped.tsv',
    'BoneMarrow_counts_deduped.tsv',
    'Brain_counts_deduped.tsv',
    'Breast_counts_deduped.tsv',
    'CervixUteri_counts_deduped.tsv',
    'Colon_counts_deduped.tsv',
    'Esophagus_counts_deduped.tsv',
    'FallopianTube_counts_deduped.tsv',
    'Heart_counts_deduped.tsv',
    'Kidney_counts_deduped.tsv',
    'Liver_counts_deduped.tsv',
    'Lung_counts_deduped.tsv',
    'Muscle_counts_deduped.tsv',
    'Nerve_counts_deduped.tsv',
    'Ovary_counts_deduped.tsv',
    'Pancreas_counts_deduped.tsv',
    'Pituitary_counts_deduped.tsv',
    'Prostate_counts_deduped.tsv',
    'SalivaryGland_counts_deduped.tsv',
    'Skin_counts_deduped.tsv',
    'SmallIntestine_counts_deduped.tsv',
    'Spleen_counts_deduped.tsv',
    'Stomach_counts_deduped.tsv',
    'Testis_counts_deduped.tsv',
    'Thyroid_counts_deduped.tsv',
    'Uterus_counts_deduped.tsv',
    'Vagina_counts_deduped.tsv'
]

# new headers for sample attributes
headers = {
    'SAMPID': 'sample_id',
    'SMATSSCR': 'autolysis_score',
    'SMNABTCH': 'iso_batch_id',
    'SMNABTCHT': 'iso_batch_type',
    'SMNABTCHD': 'iso_batch_date',
    'SMGEBTCH': 'gtype_batch_id',
    'SMGEBTCHD': 'gtype_batch_date',
    'SMGEBTCHT': 'gtype_batch_type',
    'SMCENTER': 'bss_site',
    'SMPTHNTS': 'path_notes',
    'SMRIN': 'rin_num',
    'SMTS': 'tissue_type',
    'SMTSD': 'tissue_detail',
    'SMUBRID': 'uberon_id',
    'SMTSISCH': 'total_isch_time',
    'SMTSPAX': 'paxgene_time',
    'SMAFRZE': 'analysis_freeze',
    'SMGTC': 'gtype_file',
    'SME2MPRT': 'end2_map_rate',
    'SMCHMPRS': 'chimeric_pairs',
    'SMNTRART': 'intragenic_rate',
    'SMNUMGPS': 'num_gaps',
    'SMMAPRT': 'map_rate',
    'SMEXNCRT': 'exon_rate',
    'SM550NRM': 'norm_5end',
    'SMGNSDTC': 'genes_detec',
    'SMUNMPRT': 'uniq_map_rate',
    'SM350NRM': 'norm_3end',
    'SMRDLGTH': 'max_read_len',
    'SMMNCPB': 'mean_cov_pb',
    'SME1MMRT': 'end1_mismatch',
    'SMSFLGTH': 'frag_len_std',
    'SMESTLBS': 'est_lib_size',
    'SMMPPD': 'mapped',
    'SMNTERRT': 'intergenic_rate',
    'SMRRNANM': 'num_rrna',
    'SMRDTTL': 'total_reads',
    'SMVQCFL': 'failed_qc',
    'SMMNCV': 'mean_coef_var',
    'SMTRSCPT': 'trans_detec',
    'SMMPPDPR': 'mapped_pairs',
    'SMCGLGTH': 'cum_gap_len',
    'SMGAPPCT': 'gap_percent',
    'SMUNPDRD': 'unpaired_reads',
    'SMNTRNRT': 'intron_rate',
    'SMMPUNRT': 'map_uniq_rate',
    'SMEXPEFF': 'exp_prof_eff',
    'SMMPPDUN': 'mapped_uniq',
    'SME2MMRT': 'end2_mismatch',
    'SME2ANTI': 'end2_antisense',
    'SMALTALG': 'alt_align',
    'SME2SNSE': 'end2_sense',
    'SMMFLGTH': 'frag_len_mean',
    'SMSPLTRD': 'split_reads',
    'SME1ANTI': 'end1_antisense',
    'SMBSMMRT': 'base_mismatch',
    'SME1SNSE': 'end1_sense',
    'SME1PCTS': 'end1_pc_sense',
    'SMRRNART': 'rate_rrna',
    'SME1MPRT': 'end1_map_rate',
    'SMNUM5CD': 'num_cov_5end',
    'SMDPMPRT': 'dup_rate',
    'SME2PCTS': 'end2_pc_sense'
}

# load phenotype metadata
df_ptype = pd.read_csv('https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt', 
sep='\t',index_col=0).sort_index()

# load sample metadata
meta_df = pd.read_csv('https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt', sep='\t')

# re-format sex and death hardy scale attributes for readability
df_ptype["SEX"].replace({1: 'M', 2: 'F'}, inplace=True)
df_ptype["DTHHRDY"].replace({0.0: 'Ventilator', 
                             1.0: 'Fast violent', 
                             2.0: 'Fast natural', 
                             3.0: 'Intermediate',
                             4.0: 'Slow'}, inplace=True)

# iterate through each tissue to remove all non-relevant sample metadata
for f in filelist:
    print("Now parsing " + f + "...")

    # obtain all relevant sample IDs in a tissue dataset
    gene_df = pd.read_csv(f, sep='\t')
    samps = gene_df.columns

    # create new metadata table with only relevant sample IDs
    rel_meta_df = meta_df[meta_df['SAMPID'].isin(samps)]
    
    # check that correct number of samples remain in new metadata table
    if rel_meta_df.shape[0] == len(samps):
        print("Metadata parsed for " + f)
    else:
        print(rel_meta_df.shape[0], len(samps))
        raise Exception("Sample number mismatch")

    # combine phenotype and sample metadata
    print("Now formatting " + f + "...")

    # obtain current metadata column headers
    curr_head = rel_meta_df.columns.tolist()

    # match and assign new column headers in order 
    new_head = []
    for h in curr_head:
        new_head.append(headers[h])
    rel_meta_df.columns = new_head

    # initialize lists to store sample phenotype attributes in order
    samp_sex = []
    samp_age = []
    samp_death = []

    # build ordered list of phenotype data matched to samples by donor ID
    ptype_dict = df_ptype.to_dict('index')
    for entry in rel_meta_df['sample_id']:
        donor = '-'.join(entry.split('-')[:2])
        samp_sex.append(ptype_dict[donor]['SEX'])
        samp_age.append(ptype_dict[donor]['AGE'])
        samp_death.append(ptype_dict[donor]['DTHHRDY'])

    # prepend phenotype attribute columns to sample attributes
    rel_meta_df.insert(1, column='dthhrdy', value=samp_death)
    rel_meta_df.insert(1, column='age', value=samp_age)
    rel_meta_df.insert(1, column='sex', value=samp_sex)

    # save reformatted and combined metadata table
    print("Done parsing " + f + "!")
    rel_meta_df.to_csv('datasets/metadata/formatted/' + f, 
        sep='\t', index=False)
