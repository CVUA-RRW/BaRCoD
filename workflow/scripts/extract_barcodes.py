#! /usr/bin/env python3
# -*- coding: utf-8 -*-


import pandas as pd


def find_amplicon_pos(df, min_length, max_length):
    """
    Find amplicon position from a list of matches.
    Returns a tuple (start, stop).
    
    Will yield occurences of matches on the + strand followed by a match on the minus strand.
    """
    dfs = df.sort_values(by='start')
    start = None
    for index, row in dfs.iterrows():
        if row["strand"] == 'plus' and not start: 
            start = row['start']
        elif row['strand'] == 'minus' and start:
            end = row['start'] # blast reports end as the 3' position of the primer on the reverse strand
            length = int(end) - int(start)
            
            if length >= min_length and length <= max_length:
                yield (start, end, length)
            
            start = None  # find next matching seq


def main(blastfile, reportout, min_length, max_length):
    df = pd.read_csv(blastfile, 
                        sep="\t", 
                        header=0)

    # empty df to store amplicon informations
    dfout = pd.DataFrame(columns = ['seqid', 'taxid', 'start', 'end', 'length'])

    for s in set(df.seqid):
        sub_df = df.loc[df['seqid'] == s]
        for pos in find_amplicon_pos(sub_df, min_length, max_length):
            df_s = pd.DataFrame([{'seqid' : s,
                                    'taxid' : list(sub_df['taxid'])[0],
                                    'start' : pos[0],
                                    'end' : pos[1],
                                    'length' : pos[2],
                                    }])
            dfout = pd.concat([dfout, df_s])
    dfout.to_csv(reportout, sep='\t', header=False, index=False) 


if __name__ == '__main__':
    main(snakemake.input['blast'], snakemake.output['positions'],
         snakemake.params['min_length'], snakemake.params['max_length'])
