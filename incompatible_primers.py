import pandas as pd
import os
import itertools
import numpy as np
import subprocess


def create_df(sam_file):
    # read the input file into a pandas dataframe
    primers_df = pd.read_csv(sam_file, sep='\t', header=None)

    # extract columns 1, 3, 4, 6 and 20, dropping the rest of the columns
    primers_df = primers_df[[0,1, 2, 3, 5 ,19]]

    # name the columns of the new dataset
    primers_df.columns = ['primer_name','strand', 'chromosome', 'pos', 'cigar' ,'extra_hits']

    primers_df['extra_hits'] = primers_df['extra_hits'].str.replace('XA:Z:', '')
    primers_df['strand'] = primers_df['strand'].replace(16, '-')
    primers_df['strand'] = primers_df['strand'].replace(0, '+')

    primers_df.fillna("nan",inplace=True)

    df_list = pd.DataFrame()

    # df_list.columns = ['primer_name', 'chromosome', 'pos', 'cigar', 'extra_hits']

    for index, row in primers_df.iterrows():
        if row['extra_hits'] != 'nan':
            hits = row['extra_hits'].split(';')

            for hit in hits:
                if hit != '':
                    if len(hits) < 50:
                        hit_chrom, hit_strand_pos, hit_cigar, _ = hit.split(',')
                        new_row = row.copy()
                        new_row['extra_hits'] = 'nan'
                        new_row['chromosome'] = hit_chrom
                        hit_strand = hit_strand_pos[:1]
                        hit_pos = hit_strand_pos[1:]
                        new_row['pos'] = hit_pos
                        new_row['cigar'] = hit_cigar
                        new_row['strand'] = hit_strand

                        df_list = pd.concat([df_list, new_row.to_frame(1).T])
                    else:
                        pass


    primers_df = pd.concat([primers_df, df_list], ignore_index=True)
    primers_df.drop('extra_hits', axis=1, inplace=True)

    primers_df.to_csv('all_primers_mapping.csv', sep='\t', index=False)

    sub_dfs = {}
    for chromosome in primers_df['chromosome'].unique():
        sub_dfs[chromosome] = primers_df[primers_df['chromosome'] == chromosome]


    # write the resulting dataframe to a new

    try:
        for chr in sub_dfs.keys():
            file_name = 'chr_dfs/'+ str(chr) + '_df.csv'
            sub_dfs[chr].to_csv(file_name, sep='\t', index=False)
    except OSError:
        subprocess.run(['mkdir', 'chr_dfs'])
        for chr in sub_dfs.keys():
            file_name = 'chr_dfs/'+ str(chr) + '_df.csv'
            sub_dfs[chr].to_csv(file_name, sep='\t', index=False)


####### using dataframes for

def find_incompatibilities():
    print('Buscando incompatibilidades entre los primers...')

    targets_bed = pd.read_csv('example/targets-CFTR.bed', sep='\t', header=None)
    targets_bed = targets_bed.loc[:,1:2]

    targets_coord = [(row[1][1], row[1][2]) for row in targets_bed.iterrows()]

    incompatible_set = set()

    for file in os.listdir('chr_dfs'):
        chr_df = pd.read_csv('chr_dfs/' + file, sep='\t')

        # Filtering out on target matchs
        for row in chr_df.iterrows():
            hit_pos = row[1]['pos']
            for target in targets_coord:
                if target[0] < hit_pos < target[1]:
                    chr_df.drop(row[0], inplace=True)

        plus_pos = chr_df[chr_df["strand"] == "+"]["pos"]
        plus_primer_name = chr_df[chr_df["strand"] == "+"]["primer_name"]
        minus_pos = chr_df[chr_df["strand"] == "-"]["pos"]
        minus_primer_name = chr_df[chr_df["strand"] == "-"]["primer_name"]

        plus_pos = np.array(plus_pos)[:, None]
        a, b = np.broadcast_arrays(plus_pos, minus_pos)
        df = pd.DataFrame(np.abs(a - b))
        df.columns = minus_primer_name
        df.index = plus_primer_name


        for index, row in df.iterrows():
            for column, value in row.items():
                if value < 2000:
                    primers = [index, column]
                    primers.sort()
                    pair = (primers[0], primers[1])
                    incompatible_set.add(pair)


    incompatiblity_dict = {}

    for (primer_1, primer_2) in incompatible_set:
        incompatiblity_dict[primer_1] = []
        incompatiblity_dict[primer_2] = []

    for (primer_1, primer_2) in incompatible_set:
        incompatiblity_dict[primer_1].append(primer_2)
        incompatiblity_dict[primer_2].append(primer_1)

    total_incomp = 0
    for value in incompatiblity_dict.values():
        total_incomp += len(value)

    if total_incomp > 0:
        print('Se encontaron', total_incomp, 'incompatiblidades entre los primers.')
        print('Estos excluiran mutuamente durante la fase de optimizacion.')
    else:
        print('No se encontaron incompatiblidades entre los primers.')

    return incompatiblity_dict
