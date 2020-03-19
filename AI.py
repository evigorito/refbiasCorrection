""" Combine pre and post remapping AI across samples. Use 99%CI to decide if fSNPs with AI estimate >0.5 are used or discarded"""

import argparse
import sys
import pandas as pd
from functools import reduce


def parse_options():

    parser = argparse.ArgumentParser(description="Add reads mapping REF and ALT "
                                     "alleles before and after remapping")

    parser.add_argument("--initial-AI", "-i",  nargs='+', default=[],
                        help="Full name for files with initial AI per sample")
    parser.add_argument("--post-remap_AI", nargs='+', default=[],
                        help="Full name for files with post-remapping AI per sample")
    parser.add_argument("--output_file", "-o", action='store',
                        help="Full name for output file")

    options = parser.parse_args()

    return options


def count_reads_snp(files):
    
    dfs = [pd.read_csv(f, sep = " ") for f in files]
    for df in dfs:
        df.set_index(["CHROM", "POS", "REF", "ALT"], inplace = True)

    comb = reduce(lambda x, y: x.add(y), dfs)

    comb['Total'] = comb['NREF'] + comb['NALT']
    
    comb['AI'] = comb['NALT']/comb['Total']

    return comb


def binom_test(x,y):
    """ Computes AI null (0.5) plus 99% CI for binomial test of AI_post == 0.5, with x=AI_post and y=Total_post reads. Binomial test based on 99% CI"""
    if y == 0:
        return float('NaN')
    else:
        z = 0.5 + 2.576*(x*(1-x)/y)**0.5
        return z



def main(AI_before, AI_after, output):

    ai1 = count_reads_snp(AI_before)
    ai2 = count_reads_snp(AI_after)

    df = pd.merge(ai1, ai2, left_index=True, right_index=True, suffixes=["_pre", "_post"])
    
    df['Keep'] = df.apply(lambda y: "yes" if y['AI_post'] < binom_test(y['AI_post'], y['Total_post']) else "no", axis=1)

    df.to_csv(output, sep=" ", index=True, na_rep="NA")


if __name__ == '__main__':

    #sys.stderr.write("command line: %s\n" % " ".join(sys.argv))
    options = parse_options()

    main(options.initial_AI,
         options.post_remap_AI,
         options.output_file)
    
