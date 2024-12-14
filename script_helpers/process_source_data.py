import argparse
import gzip

import pandas as pd

def parse_args():
    parser = argparse.ArgumentParser(description='Exctracts DNAm fraction')
    parser.add_argument('-i', '--intensity', required=True, nargs='?',
                        help='path to input signal intensities file')
    parser.add_argument('-d', '--description', required=True, nargs='?',
                        help='path to input signal intensities file')
    parser.add_argument('-o', '--out',  required=True, nargs='?',
                        help='path to output file (pandas .pkl)')
    parser.add_argument('-od', '--outdescr',  required=True, nargs='?',
                        help='path to output file (pandas .pkl)')
    parser.add_argument('-f', '--filter', default=None, type=str, nargs='?',
                        help='filter the particular DNAm sites from the .txt file with IDs of DNAm sites')
    return parser.parse_args()


def read_sample(path: str) -> pd.DataFrame:
    """
    extracts information about samples (ages)
    """
    sample_info = {}

    with gzip.open(path, 'rt') as out_f:
        for line in out_f:
            if line.startswith('!Sample'):
                dat = line.strip().split('\t')
                name = dat[0].replace('!', '')
                if name == 'Sample_characteristics_ch1':
                    name = dat[1][1:].split(' ')[0].replace(':', '')
                dat_list = [i.replace('"', '') for i in dat[1:]]
                sample_info[name] = dat_list
    samples_descr = pd.DataFrame.from_dict(sample_info, orient='columns')
    samples_descr.Sample_title = samples_descr.Sample_title.str.split('_').str[0]
    samples_descr.age = samples_descr.age.str.split(' ').str[-1].astype(int)

    # we are only interested in aging information
    return samples_descr[['Sample_title', 'age']]


def read_intensities(path: str, filter_file: str | None = None) -> pd.DataFrame:
    """
    read intensities into a pandas DataFrame object
    """
    signal_intensities = pd.read_csv(path,
                                    compression='gzip', sep='\t', 
                                    index_col='ID_REF')
    if filter_file is not None:
        ids = pd.read_pickle(filter_file)
        signal_intensities = filter_dat(signal_intensities, ids)
    return signal_intensities


def filter_dat(df: pd.DataFrame, ids: pd.Series) -> pd.DataFrame:
    """
    filter only target DNAm sites based on pd.Series of ids
    """
    sub_df = df.loc[ids]

    assert sub_df.shape[0] == ids.shape[0]  # all sites found
    assert sub_df.index.value_counts().max() == 1  # all sites unique
    
    return sub_df


def calculate_DNAm_frac(source_df: pd.DataFrame) -> pd.DataFrame:
    """
    for each DNAm calculates fraction of methylated sites  
    uses information from columns: "X.methylated" and "X.unmethylated"
    """
    df_T = source_df.T
    # get info
    df_T['sample'] = df_T.index.str.split('.').str[0]
    df_T['tag'] = df_T.index.str.split('.').str[1]
    # drop p-value colums
    df_T = df_T.loc[df_T['tag'] != 'detectionPval']
    # calc sum
    df_T_sum = df_T.groupby('sample').sum().drop(columns=['tag'])
    # get DNAm counts
    df_T_me = df_T.loc[df_T['tag'] == 'methylated'].reset_index(drop=True).set_index('sample').drop(columns=['tag'])

    # divide one and another
    df_T_frac = df_T_me / df_T_sum

    return df_T_frac.T


if __name__ == '__main__':
    intensity_file = parse_args().intensity
    description_file = parse_args().description
    out_file = parse_args().out
    out_descr_file = parse_args().outdescr
    filter_file = parse_args().filter

    # read and process intensity files 
    df = read_intensities(path=intensity_file, filter_file=filter_file)
    df = calculate_DNAm_frac(source_df=df)
    df.to_pickle(out_file)

    # read and process description file 
    descr_df = read_sample(description_file)
    descr_df.to_pickle(out_descr_file)
