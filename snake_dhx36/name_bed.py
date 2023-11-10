import pandas as pd
import sys
def name_bed(bedpath, outpath):
    bed = pd.read_csv(bedpath, 
                    sep = '\t',
                    header = None, 
                    names = ['chrom', 'start', 'end', 'pval', 'fc', 'strand'])
    bed['name'] = bed['chrom']+'_'+bed['start'].astype(str)+'_'+bed['end'].astype(str)+'_'+bed['strand']

    bed[['chrom', 'start', 'end', 'name', 'fc', 'strand']].to_csv(
        outpath,
    index = False,
    sep = '\t',
    header = False)

if __name__=='__main__':
    name_bed(sys.argv[1], sys.argv[2])