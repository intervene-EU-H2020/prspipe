
import pandas as pd
import argparse

def get_args():
    p = argparse.ArgumentParser()
    p.add_argument('-a','--annot')
    p.add_argument('-b','--bim')
    p.add_argument('-o','--out')
    args = p.parse_args()
    return args

def main():
    
    """
    Reads annot-file and intersects with bim.
    Will crash if not all variants in bim are in annot.
    Ldsc will crash if the variants are not ordered by position. This is not checked here.
    """
    
    args = get_args()
    
    var_bim = pd.read_csv(args.bim, header=None, sep='\t', usecols=[1], dtype=str)
    var_bim.columns = ['SNP']
    
    annot = pd.read_csv(args.annot, sep='\t', header=0)
    
    cols = annot.columns.tolist()
    
    annot['SNP'] = annot.SNP.astype(str)
    annot.set_index('SNP', inplace=True, drop=False)
    
    annot = annot.loc[var_bim.SNP.values]
    
    annot.to_csv(args.out, sep='\t', header=True, index=False)

if __name__ == '__main__':
    main()