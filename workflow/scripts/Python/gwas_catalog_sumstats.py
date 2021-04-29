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
    # TODO check I am calculating sample size correctly - there was no field for this so I'm adding all individuals in the initial ancestries (i.e. not replicates)
    N = sum([_['numberOfIndividuals'] for _ in study_details['ancestries'] if _['type']=='initial'])

    return N


def get_number_snps(study_details):
    """Returns the number of snps in the study given in the GWAS catalog

    Args:
        study_details: The study information obtained from the API endpoint
    """
    n_snps = study_details['snpCount']

    return n_snps


def get_summary_statistics(study_id, study_details):
    """Get summary statistics from the GWAS catalog for a specified study

    Args:
        study_id: The study id from the GWAS catalog 
        study_details: The study information obtained from the API endpoint
    """
    size = 500 # the API returns 20 items by default and you increase this by setting the size parameter (up to a limit of 500)
    n_snps = get_number_snps(study_details)

    url = f'http://www.ebi.ac.uk/gwas/summary-statistics/api/studies/{study_id}/associations?size={size}'

    df = pd.DataFrame()

    for _ in range(0,n_snps,size):
        counter = int(_/size)+1
        print(f'>> downloading page {counter}')
        data = get_request(url)
        url = data['_links']['next']['href'] # get url for next page
        df = pd.concat([df, pd.DataFrame(data['_embedded']['associations'].values())])
        if _==1000: break # TODO remove this later - this is here now to speed up execution for testing purposes

    return df


def cleanup_summary_statistics(df, study_id, study_details):
    """Reformat the summary statistics download from the GWAS catalog
    into a format that can be used by the QC script

    Args:
        df: The dataframe of summary statistics downloaded from the GWAS catalog
        study_id: The study id from the GWAS catalog 
        study_details: The study information obtained from the API endpoint
    """
    # rename the columns to those used by https://github.com/bulik/ldsc/blob/master/munge_sumstats.py
    rename_cols = {'variant_id':'SNP', 
                    'chromosome':'CHR',
                    'study_accession':'STUDY',
                    'p_value':'P', 
                    'effect_allele':'A1', 
                    'other_allele':'A2', 
                    'odds_ratio':'OR', 
                    'beta':'BETA', 
                    'effect_allele_frequency':'FRQ'}

    df = df.rename(columns=rename_cols)

    # drop irrelevant columns
    # also drop anything with all NaNs (otherwise the QC script will remove all rows)
    df = df[rename_cols.values()]
    df = df.dropna(axis=1,how='all')

    df['N'] = get_sample_size(study_details)

    return df


def download_sumstats(args, p=True):
    if args.study_id is None:
        raise ValueError('The --study-id flag is required.')
    if args.out is None:
        raise ValueError('The --out flag is required.')

    print(f'downloading summary statistics for {args.study_id}')
    study_details = get_study_details(args.study_id)
    sumstats = get_summary_statistics(args.study_id, study_details)
    cleaned_sumstats = cleanup_summary_statistics(sumstats, args.study_id, study_details)
    cleaned_sumstats.to_csv(args.out, index=None)
    print(f'>> output saved at {args.out}')


if __name__ == '__main__':
    download_sumstats(parser.parse_args(), p=True)