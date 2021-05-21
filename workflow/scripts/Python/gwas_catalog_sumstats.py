"""Code for downloading GWAS catalog summary statistics using the API endpoint

Based on documentation at https://www.ebi.ac.uk/gwas/summary-statistics/docs/

The summary statistics are formatted as .csv with the following columns:

    SNP: variant_id
    CHR: chromosome
    STUDY: study_accession
    P: p_value
    A1: effect_allele
    A2: other_allele
    OR: odds_ratio
    BETA: beta
    FRQ: effect_allele_frequency

Note that it is assumed that the user has checked that the summary statistics for the 
given study id supports API access (you can check this at the GWAS catalog website)
"""

import requests
import argparse
import pandas as pd

import re
import os
import ftplib
import ntpath

import sys
import time

import gzip
import shutil

from progressbar import AnimatedMarker, Bar, BouncingBar, Counter, ETA, \
    AdaptiveETA, FileTransferSpeed, FormatLabel, Percentage, \
    ProgressBar, ReverseBar, RotatingMarker, \
    SimpleProgress, Timer, UnknownLength


parser = argparse.ArgumentParser()

parser.add_argument('--study-id', default=None, type=str,
                    help="Study ID from the GWAS catalog")
parser.add_argument('--out', default=None, type=str,
                    help="Output filename prefix")


def get_request(url):
    """Make API requests and return the data in JSON format

    Args:
        url: The url for the GET request
    """
    try:
        response = requests.get(url,timeout=10)
        response.raise_for_status()
    except requests.exceptions.RequestException as err:
        print ("Error:",err)
    except requests.exceptions.HTTPError as errh:
        print ("Http Error:",errh)
    except requests.exceptions.ConnectionError as errc:
        print ("Error Connecting:",errc)
    except requests.exceptions.Timeout as errt:
        print ("Timeout Error:",errt)   

    return response.json()  


def get_study_details(study_id):
    """Returns detailed information about the study given in the GWAS catalog

    Args:
        study_id: The study id from the GWAS catalog 
    """
    url = f'https://www.ebi.ac.uk/gwas/rest/api/studies/{study_id}'

    study_data = get_request(url)

    return study_data


def get_sample_size(study_details):
    """Returns the sample size information about the study given in the GWAS catalog

    Args:
        study_details: The study information obtained from the API endpoint
    """
    # note there's no field for sample size in the GWAS catalog, so this is calculated here by summing individual ancestry populations
    # (therefore it's advised to specify the sample size in the input studies.tsv file)
    N = sum([_['numberOfIndividuals'] for _ in study_details['ancestries'] if _['type']=='initial'])

    return N


def cleanup_summary_statistics(df, study_id, n_gwas):
    """Reformat the summary statistics download from the GWAS catalog
    into a format that can be used by the QC script

    Args:
        df: The dataframe of summary statistics downloaded from the GWAS catalog
        study_id: The study id from the GWAS catalog 
        n_gwas: The sample size of the GWAS study
    """
    # rename the columns to those used by https://github.com/bulik/ldsc/blob/master/munge_sumstats.py
    rename_cols = {'variant_id':'SNP', 
                    'chromosome':'CHR',
                    'p_value':'P', 
                    'effect_allele':'A1', 
                    'other_allele':'A2', 
                    'odds_ratio':'OR', 
                    'beta':'BETA', 
                    'effect_allele_frequency':'FRQ'}

    df = df.rename(columns=rename_cols)

    df['STUDY'] = study_id

    # ensure A1,A2 are upper case - otherwise QC script disregards them
    df['A1'] = df['A1'].str.upper()
    df['A2'] = df['A2'].str.upper()

    # drop irrelevant columns
    # also drop anything with all NaNs (otherwise the QC script will remove all rows)
    df = df[rename_cols.values()]
    df = df.dropna(axis=1,how='all')

    df['N'] = n_gwas

    return df


def file_write(data, file, pbar):
   file.write(data) 
   pbar += len(data)


def download_ftp_sumstats(local_path, remote_path):
    """
    Download GWAS summary statistics using FTP download at given remote path
    Stores the unzipped file at the given local path

    Args: 
        local_path: Filepath for storing the downloaded (both raw and unzipped) files
        remote_path: FTP path for the file download
    """
    ftp = ftplib.FTP('ftp.ebi.ac.uk')
    ftp.login()

    size = ftp.size(remote_path)

    widgets = ['Downloading: ', Percentage(), ' ',
                        Bar(marker='#',left='[',right=']'),
                        ' ', ETA(), ' ', FileTransferSpeed()]

    pbar = ProgressBar(widgets=widgets, maxval=size)
    pbar.start()

    with open(local_path, 'wb') as outfile:
        ftp.retrbinary("RETR " + remote_path, lambda block: file_write(block, outfile, pbar))

    # unzip the file
    #with gzip.open(local_path, 'rb') as f_in:
    #    with open(local_path[:-3], 'wb') as f_out:
    #        shutil.copyfileobj(f_in, f_out)

    sumstats_df = pd.read_csv(local_path, sep='\t')
    return sumstats_df


def download_sumstats(args, p=True):
    if args.study_id is None:
        raise ValueError('The --study-id flag is required.')
    if args.out is None:
        raise ValueError('The --out flag is required.')

    studies = pd.read_csv('config/studies.tsv', dtype=str, sep='\t').set_index(["study_id"], drop=False)
    ftp_path = studies.loc[args.study_id]['ftp_address']
    local_path = '{}/{}'.format('/'.join(args.out.split('/')[:-1]), ftp_path.split('/')[-1])

    print(f'downloading summary statistics for {args.study_id}')
    sumstats = download_ftp_sumstats(local_path, ftp_path)
    n_gwas = int(studies.loc[args.study_id]['n_cases']) + int(studies.loc[args.study_id]['n_controls'])
    cleaned_sumstats = cleanup_summary_statistics(sumstats, args.study_id, n_gwas)
    cleaned_sumstats.to_csv(args.out, index=None)
    print(f'>> output saved at {args.out}')


if __name__ == '__main__':
    download_sumstats(parser.parse_args(), p=True)
