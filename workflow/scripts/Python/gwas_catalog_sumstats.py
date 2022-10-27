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
from pandas.errors import ParserError
import numpy as np

import re
import os
import ftplib
import ntpath

import sys
import time

import gzip
import shutil

parser = argparse.ArgumentParser()

parser.add_argument('--study-id', required=True, type=str,
                    help="Study ID from the GWAS catalog")
parser.add_argument('--out', required=True, type=str,
                    help="Output filename prefix")
parser.add_argument('--studies', required=True, type=str, help='Path to studies.tsv containing study metadata')


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
    # this won't do anything if the columns already have the correct names
    rename_cols = {'variant_id':'SNP',
                   'chromosome':'CHR',
                   'base_pair_location':'BP',
                   'p_value':'P',
                   'effect_allele':'A1', 
                   'other_allele':'A2', 
                   'odds_ratio':'OR', 
                   'beta':'BETA', 
                   'effect_allele_frequency':'FRQ'}
    
    df = df.rename(columns=rename_cols)
    
    # below some efforts are made to recover columns from summary statistics if they don't follow the recommended naming convention.

    if 'CHR' not in df.columns:
        chr_cols = list(c for c in df.columns if str(c).lower().startswith('chrom'))
        chr_cols += list(c for c in df.columns if str(c).lower().startswith('chr'))
        if len(chr_cols):
            print('Assuming "{}" is the CHR column'.format(chr_cols[0]))
            df.rename(columns={chr_cols[0]:'CHR'}, inplace=True)
    
    if 'BP' not in df.columns:
        pos_cols = list(c for c in df.columns if str(c).lower().startswith('pos'))
        pos_cols += list(c for c in df.columns if str(c).lower() == 'bp')
        if len(pos_cols):
            print('Assuming "{}" is the BP column'.format(pos_cols[0]))
            df.rename(columns={pos_cols[0]:'BP'}, inplace=True)
        else:
            print('Warning: could not find the BP (position) column.')
    
    if 'SNP' not in df.columns:
        snp_cols = list(c for c in df.columns if str(c).lower().startswith('rsid'))
        snp_cols += list(c for c in df.columns if str(c).lower().startswith('snp'))
        if len(snp_cols):
            print('Assuming "{}" is the SNP column'.format(snp_cols[0]))
            df.rename(columns={snp_cols[0]:'SNP'}, inplace=True)
        else:
            print('Unable to determine the SNP column. Using dummy values.')
            df['SNP'] = np.array(['snp_{}'.format(i) for i in range(len(df))])
                
    if 'A1' not in df.columns:
        a1_cols = []
        if 'ALLELE1' in df.columns:
            a1_cols = ['ALLELE1']
        a1_cols += list(c for c in df.columns if str(c).lower().startswith('a1') and ('freq' not in  str(c).lower()) and ('frq' not in str(c).lower()))
        a1_cols += list(c for c in df.columns if str(c).lower().startswith('allele') and str(c).endswith('1'))
        a1_cols += list(c for c in df.columns if str(c).lower().startswith('effect') and str(c).lower().endswith('allele'))
        a1_cols += list(c for c in df.columns if str(c).lower().startswith('alter') and str(c).lower().endswith('allele'))
        a1_cols += list(c for c in df.columns if str(c).lower().startswith('minor') and str(c).lower().endswith('allele'))
        if len(a1_cols):
            if len(a1_cols) > 1:
                print('Warning: found multiple possible columns containing A1 : {}'.format(a1_cols))
            print('Assuming "{}" is the A1 column'.format(a1_cols[0]))
            df.rename(columns={a1_cols[0]:'A1'}, inplace=True)
        else:
            print('Unable to determine A1 column.')            

    if 'A2' not in df.columns:
        a2_cols = []
        if ('ALLELE0' in df.columns):
            a2_cols = ['ALLELE0']
        a2_cols += list(c for c in df.columns if str(c).lower().startswith('a2') and ('freq' not in  str(c).lower()) and ('frq' not in str(c).lower()))
        a2_cols += list(c for c in df.columns if str(c).lower().startswith('allele') and str(c).endswith('2'))
        a2_cols += list(c for c in df.columns if str(c).lower().startswith('other') and str(c).lower().endswith('allele'))
        a2_cols += list(c for c in df.columns if str(c).lower().startswith('ref') and str(c).lower().endswith('allele'))
        if len(a2_cols):
            if len(a2_cols) > 1:
                print('Warning: found multiple possible columns containing A2 : {}'.format(a1_cols))
            print('Assuming "{}" is the A2 column'.format(a2_cols[0]))
            df.rename(columns={a2_cols[0]:'A2'}, inplace=True)
        else:
            print('Unable to determine A2 column.')
            
    if 'BETA' not in df.columns:
        beta_cols = list(c for c in df.columns if 'beta' in str(c).lower())
        beta_cols += list(c for c in df.columns if 'effect' in str(c).lower() and not 'allele' in str(c).lower())
        if len(beta_cols):
            print('Assuming "{}" is the BETA column'.format(beta_cols[0]))
            df.rename(columns={beta_cols[0]: 'BETA'}, inplace=True)
            
            
    if 'FRQ' not in df.columns:
        frq_cols = list(c for c in df.columns if ('freq' in str(c).lower()) or ('frq' in str(c).lower()) )
        if len(frq_cols):
            print('Assuming "{}" is the FRQ column'.format(frq_cols[0]))
            df.rename(columns={frq_cols[0]:'FRQ'}, inplace=True)
            
    if 'P' not in df.columns:
        pv_cols = []
        pv_cols += list(c for c in df.columns if str(c).lower() == 'p')
        pv_cols += list(c for c in df.columns if (str(c).lower().startswith('p')) and ('val' in c))
        if 'P_BOLT_LMM_INF' in df.columns:
            pv_cols += ['P_BOLT_LMM_INF']
        if 'P_LINREG' in df.columns:
            pv_cols +=  ['P_LINREG']
        if len(pv_cols):
            if len(pv_cols) > 1:
                print('Warning: found multiple possible columns containing p-values: {}'.format(pv_cols))
            print('Assuming "{}" is the P column'.format(pv_cols[0]))
            df.rename(columns={pv_cols[0]:'P'}, inplace=True)
    
    df['STUDY'] = study_id

    # ensure A1,A2 are upper case - otherwise QC script disregards them
    df['A1'] = df['A1'].str.upper()
    df['A2'] = df['A2'].str.upper()

    # drop irrelevant columns
    # also drop anything with all NaNs (otherwise the QC script will remove all rows)
    
    found_cols = list(col for col in df.columns if col in rename_cols.values())
    
    if 'N' in df.columns:
        # "N" is not a standard column in the GWAS catalog, but may be included in munged sumstats
        df = df[ found_cols + ['N']]
    else:
        df = df[ found_cols ]
        df['N'] = n_gwas
        
    df = df.dropna(axis=1,how='all')

    return df


def file_write(data, file):
    file.write(data)

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

    with open(local_path, 'wb') as outfile:
        ftp.retrbinary("RETR " + remote_path, lambda block: file_write(block, outfile))
        
    if local_path.endswith('gz'):
        with gzip.open(local_path, 'rt') as infile:
            l1 = infile.read()
    else:
        with open(local_path, 'r') as infile:
            l1 = infile.read()
            
    tabsep = len(l1.split('\t')) > len(l1.split(' '))
    if tabsep:
        print("Assuming summary statistics are tab-separated based on first line")
        sep = '\t'
    else:
        if len(l1.split(',')) > len(l1.split(' ')):
            print("Assuming summary statistics are comma-separated based on first line")
            sep = ','
        else:
            print("Assuming summary statistics are whitespace-separated based on first line")
            sep = ' '

    try:
        sumstats_df = pd.read_csv(local_path, sep=sep)
    except ParserError as e:
        print("ParserError while trying to read sumstats file. Trying delim_whitespace=True")
        sumstats_df = pd.read_csv(local_path, delim_whitespace=True)
        
    return sumstats_df


def download_sumstats(args, p=True):

    studies = pd.read_csv(args.studies, dtype=str, sep='\t').set_index(["study_id"], drop=False)

    ftp_path = studies.loc[args.study_id]['ftp_address']
    
    if pd.isna(ftp_path):
        # if ftp address is not given, try local path 
        path = studies.loc[args.study_id]['local_path']
        sumstats = pd.read_csv(path, sep='\t')
    else:
        local_path = '{}/{}'.format('/'.join(args.out.split('/')[:-1]), ftp_path.split('/')[-1])
        print(f'downloading summary statistics for {args.study_id}')
        sumstats = download_ftp_sumstats(local_path, ftp_path)
    
    try:
        n_cases = int(studies.loc[args.study_id]['n_cases'])
    except ValueError:
        n_cases = 0
    try:
        n_controls = int(studies.loc[args.study_id]['n_controls'])
    except ValueError:
        n_controls = 0
    
    n_gwas = n_cases + n_controls
    
    assert n_gwas > 0, 'Error: GWAS sample size specified is 0'
    
    cleaned_sumstats = cleanup_summary_statistics(sumstats, args.study_id, n_gwas)
    cleaned_sumstats.to_csv(args.out, index=None, sep='\t')
    print(f'>> output saved at {args.out}')


if __name__ == '__main__':
    download_sumstats(parser.parse_args(), p=True)
