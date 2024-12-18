import argparse
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(description='Exctracts DNAm fractions')
    parser.add_argument('-i', '--infile', required=True, nargs='?',
                        help='path to input signal intensities file')
    parser.add_argument('-o', '--outfile', required=True, nargs='?',
                        help='path to output file (pandas .pkl)')
    return parser.parse_args()


def parse_hv_s3_file(path: str) -> pd.DataFrame:
    """
    extracts from predictors table the IDs of correpsonding CpGs
    """ 
    df = pd.read_csv(path, skiprows=2).iloc[1:]
    df = df.astype({'Gene_ID': 'int', 'MapInfo': 'int', 'Chr': 'int', 'GenomeBuild': 'int'})
    return df


if __name__ == '__main__':
    in_file = parse_args().infile
    out_file = parse_args().outfile
    df = parse_hv_s3_file(in_file)
    df['CpGmarker'].to_pickle(out_file)