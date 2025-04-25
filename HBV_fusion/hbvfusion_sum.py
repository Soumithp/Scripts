import sys
import argparse
import os
import pandas as pd
import numpy as np

def main():
    args = parse_args()
    ################ define output files ################
    outexon = os.path.abspath(args.outprefix + "_hgene_exon.csv")
    # outmat = os.path.abspath(args.outprefix + "_hgene_matrix.csv")
    outvirus = os.path.abspath(args.outprefix + "_vgene.csv")
    outmerge = os.path.abspath(args.outprefix + "_hgene_vgene.csv")

    ################ annoatate events with human genes ################
    # read input table
    df = pd.read_csv(args.inputhuman, sep="\t", header=None, dtype= 'str', index_col=False)
    df.columns = ['chr','start','end','sample','read_count','gene','exon']

    # table with exon
    agg = df.groupby(by=['chr','start','end','sample','read_count','gene']) \
        .agg( lambda x: ';'.join(x.dropna().apply(str)) ).reset_index() \
        .replace(r'^\s*$', np.nan, regex=True)
    agg.to_csv(outexon, index=False, na_rep='Unknown')

    # # table with gene-by-sample matrix with values as event count
    # mat = agg[['gene','sample']] \
    #     .pivot_table(index='gene', columns='sample', aggfunc= lambda x: len(x)) \
    #     .reset_index().fillna(0).set_index('gene').rename_axis(None, axis=1)
    # # total for each sample
    # mat.loc['total'] = mat.sum()
    # # total for each gene
    # mat['total'] = mat.sum(axis=1)
    # mat.astype(int).to_csv(outmat)

    ################ annoatate events with hbv genes ################
    inputhbv = pd.read_csv(args.inputhbv, sep="\t", header=None, dtype=str, index_col=False)
    inputhbv.columns = ['chr','start','end','sample','event','ID','Name','product']
    inputhbv['ID'] = inputhbv['ID'].replace(to_replace = pattern, regex=True) \
        .replace(to_replace = simple_pattern, regex=True)
    inputhbv['Name'] = inputhbv['Name'].replace(to_replace = pattern, regex=True) \
        .replace(to_replace = simple_pattern, regex=True)
    inputhbv['product'] = inputhbv['product'].replace(to_replace = product_pattern, regex=True) \
        .replace(to_replace = simple_pattern, regex=True) 
    inputhbv['end'] = pd.to_numeric(inputhbv['end'], errors='coerce')

    # dataframe for each read (Note that dupcate reads were collapsed in the previous step)
    # aggread = df.groupby(by=['chr','start','end','sample','event']) \
    #     .agg( lambda x: ';'.join(x.dropna().apply(str).drop_duplicates()) ).reset_index() \
    #     .replace(r'^\s*$', "Unknown", regex=True)
    
    # dataframe for each event to aggregate virus position, 0-based coordinate
    aggviruspos = inputhbv[['sample','event','end']] \
        .groupby(by=['sample','event']) \
        .agg( func='median' ).reset_index() \
        .replace(r'^\s*$', np.nan, regex=True)
    aggviruspos['viral_pos'] = aggviruspos['end'].astype(int)

    # dataframe for each event
    aggvirus = inputhbv[['sample','event','ID','Name','product']] \
        .groupby(by=['sample','event']) \
        .agg( lambda x: ';'.join(x.dropna().apply(str).drop_duplicates().sort_values()) ).reset_index() \
        .replace(r'^\s*$', np.nan, regex=True)
    # take the union of ID, Name, and product
    aggvirus['consensus'] = aggvirus \
        .apply(lambda x: ';'.join([ i for i in x[['ID','Name','product']] if str(i) != 'nan' ]).split(';'), axis=1) \
        .apply(lambda x: ';'.join(sorted(set(x)))) \
        .replace(r'^\s*$', np.nan, regex=True)
    aggvirus['id_name_consensus'] = aggvirus \
        .apply(lambda x: ';'.join([ i for i in x[['ID','Name']] if str(i) != 'nan' ]).split(';'), axis=1) \
        .apply(lambda x: ';'.join(sorted(set(x)))) \
        .replace(r'^\s*$', np.nan, regex=True)
    aggvirus = pd.merge(aggvirus, aggviruspos[['sample','event','viral_pos']], on=['sample','event'], how='outer')
    aggvirus.to_csv(outvirus, index=False, na_rep='Unknown')

    ################ Combine human gene and virus gene annotation ################
    agghuman = df.drop(columns=['exon']) \
        .groupby(by=['chr','start','end','sample','read_count']) \
        .agg( lambda x: ';'.join(x.dropna().apply(str).drop_duplicates().sort_values()) ).reset_index() \
        .replace(r'^\s*$', np.nan, regex=True)
    agghuman['event'] = agghuman['chr'] + "_" + agghuman['start'].apply(str)
    merged = pd.merge(agghuman, aggvirus, on=['sample','event'], how='outer')
    merged = merged[['sample','event','read_count','gene','id_name_consensus','consensus','viral_pos']]
    merged.columns = ['sample','event','read_count','human','virus_id','virus_product','viral_pos']
    merged.to_csv(outmerge, index=False, na_rep='Unknown')

    
    for virus_annotation in ['virus_id','virus_product']:
        outmatevent = os.path.abspath(args.outprefix + "_human-" + virus_annotation + "_eventc_matrix.csv")
        outmatread = os.path.abspath(args.outprefix + "_human-" + virus_annotation + "_readc_matrix.csv")

        # hvgene-by-sample matrix with values as event count
        mergedmat = merged.replace(np.nan, 'Unknown')
        mergedmat['annotation'] = mergedmat['human'] + "-" + mergedmat[ virus_annotation ]
        mergedmat_event = mergedmat[['annotation','sample']] \
            .pivot_table(index='annotation', columns='sample', aggfunc= lambda x: len(x)) \
            .reset_index().fillna(0).set_index('annotation').rename_axis(None, axis=1)
        # total for each sample
        mergedmat_event.loc['total'] = mergedmat_event.sum()
        # total for each gene
        mergedmat_event['total'] = mergedmat_event.sum(axis=1)
        mergedmat_event.astype(int).to_csv(outmatevent)

        # hvgene-by-sample matrix with values as aggregated read count
        mergedmat_read = mergedmat[['annotation','sample','read_count']]
        mergedmat_read['read_count_num'] = pd.to_numeric(mergedmat['read_count'], errors='coerce')
        mergedmat_read = mergedmat_read[['annotation','sample','read_count_num']] \
            .pivot_table(values='read_count_num', index='annotation', columns='sample', aggfunc='sum') \
            .reset_index().fillna(0).set_index('annotation').rename_axis(None, axis=1)
        mergedmat_read.astype(int).to_csv(outmatread)



# dict to reformatting ID column
pattern = {
    r'^hbv_ref' : np.nan,
    r'^Gen[Bb]ank' : np.nan,
    r'^[A-Z]{3}[0-9]+': np.nan,
    r'""': np.nan,
    r'^HBAust.*$': np.nan, r'^pgRNA.*$': np.nan,
    r'^[P]?C.*$': "Pre-C/C", r'^[Pp]re[-]?\s?[Cc].*$': "Pre-C/C", r'^HBVgp4.*$': "Pre-C/C", r'^(HBcAg)|(HBe/c.*$)': "Pre-C/C",
    r'^(Large\s)?S.*$': "Pre-S1/Pre-S2/S", r'^[Pp]re[-]?\s?[Ss].*$': "Pre-S1/Pre-S2/S", r'^HBVgp2.*$': "Pre-S1/Pre-S2/S", r'^(LHBs.*$)|(MHBs.*$)': "Pre-S1/Pre-S2/S",
    r'^P[ol]?.*$': "P", r'^pol.*$': "P", r'^HBVgp1.*$': "P",
    r'^[xX].*$': "X", r'^HBVgp3.*$': "X", r'^HBx.*$': "X"
}
# dict to reformatting product column
product_pattern = {
    r'^HBAust.*$': np.nan, r'^pgRNA.*$': np.nan,
    r'^[P]?C.*$': "Pre-C/C", r'^[Pp]re[-]?\s?[Cc].*$': "Pre-C/C", r'^HBVgp4.*$': "Pre-C/C", r'^(HBcAg)|(HBe/c.*$)': "Pre-C/C", r'^.*capsid.*$': "Pre-C/C",
    r'^(Large\s)?S.*$': "Pre-S1/Pre-S2/S", r'^[Pp]re[-]?\s?[Ss].*$': "Pre-S1/Pre-S2/S", r'^HBVgp2.*$': "Pre-S1/Pre-S2/S", r'^(LHB[S].*$)|(MHB[S].*$)': "Pre-S1/Pre-S2/S",
    r'^P[ol]?.*$': "P", r'^pol.*$': "P", r'^HBVgp1.*$': "P",
    r'^[xX].*$': "X", r'^HBVgp3.*$': "X", r'^HBx.*$': "X",
    r'^.*core.*$': "Pre-C/C", r'^HBeAg.*$': "Pre-C/C",
    r'^.*surface.*$': "Pre-S1/Pre-S2/S", r'^.*S\sprotein.*$': "Pre-S1/Pre-S2/S", r'HBsAg': "Pre-S1/Pre-S2/S", r'^.*envelop.*$': "Pre-S1/Pre-S2/S",
    r'^.*polymerase.*$': "P"
}
# simplyfied labels
simple_pattern = {
    'Pre-S1/Pre-S2/S': 'S',
    'Pre-C/C': 'C',
}

def parse_args():
    parser = argparse.ArgumentParser(prog='hbvfusion_sum.py')

    ## required ##
    essential_args = parser.add_argument_group("Required")
    essential_args.add_argument("--inputhuman", type=str, help="Path to input file", required=True)
    essential_args.add_argument("--inputhbv", type=str, help="Path to input file", required=True)
    essential_args.add_argument("--outprefix", type=str, help="output file prefix", required=True)

    args = parser.parse_args()

    ## Required ##
    args.inputhuman = os.path.abspath(args.inputhuman)
    args.inputhbv = os.path.abspath(args.inputhbv)
    args.outprefix = str(args.outprefix)

    return args


if __name__ == "__main__":
    main()
