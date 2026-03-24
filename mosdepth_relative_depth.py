import sys
import pandas as pd

def main():
    bed_file = sys.argv[1] # .regions.bed.gz from mosdepth
    out_file = sys.argv[2] # Output filename

    # output format of mosdepth(chrom, start, end, mean_depth)
    df = pd.read_csv(bed_file, sep='\t', header=None, names=['chrom', 'start', 'end', 'mean_depth'])

    # median of autosomes (filter chrX and chrY)
    autosomes = df[~df['chrom'].isin(['chrX', 'chrY'])]
    median_bindepth = autosomes['mean_depth'].median()

    # relative_depth
    denominator = (median_bindepth / 2) if median_bindepth > 0 else 1e-9
    df['relative_depth'] = df['mean_depth'] / denominator

    #  Convert bin_start to 1-based (mosdepth is 0-based), and generate bin_number
    df['bin_start'] = df['start'] + 1
    df['bin_number'] = df.index + 1

    # Output:chrom, bin_start, relative_depth, bin_number
    df[['chrom', 'bin_start', 'relative_depth', 'bin_number']].to_csv(
        out_file, sep='\t', index=False, header=False, quoting=3)

if __name__ == "__main__":
    main()
