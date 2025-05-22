#!/usr/bin/env python3

from itertools import islice
import sys
import pandas as pd 
import time

############

infi = sys.argv[1]  #input file specified in command line
outfi = sys.argv[2]  #output file specified in command line
hits = [] # list for storing hits for dataframe

############


start = time.perf_counter()

with open(infi, 'r') as fi, open(outfi, 'w') as ou:
    writethis = False
    for ln in fi:
        if ln.startswith("Query:"):
            query_id = ln.strip().split("Query:",1)[-1].split()[0].rsplit(".hmm", 1)[0] #pulls QueryID
            lines = list(islice(fi, 6)) #reads next 6 lines
            if not any("[No hits detected that satisfy reporting thresholds]" in line.lower() for line in lines):
                writethis = True
                continue

#if hits found then:
        if writethis:
            if ln.startswith("Domain"): #stops writing at "Domain"
                writethis = False
                continue 
            parts = ln.strip().split()
            if parts:
                try:
                    e_value = float(parts[0])
                    gene_id = (parts[8]) #don't need to specify variable type as it's already a str, doesnt need coverting
                    seq_score = (parts[1])
                    gene_start = (parts[10])
                    gene_end = (parts[12])

                    if e_value < 1e-15: #filter by evalue------ maybe i move the filtering to within the dataframe?
                        hits.append([gene_id, gene_start, gene_end, query_id]) #store hits in hit list
                    continue
                except ValueError: #skips if non-numerical (error handling)
                    pass


    if hits: #creates df with parts stored in hits list and sorts by gene id, first by the contig and then by the gene number
        df = pd.DataFrame(hits, columns=['gene_id', 'gene_start', 'gene_end', 'query_id'])

        # Extract sorting values with error handling
        df['sort_1'] = pd.to_numeric(
            df['gene_id'].str.extract(r'\|[^|]+\|(\d+)')[0], 
            errors='coerce'
        ).fillna(-1).astype(int)  # Convert NaN to -1

        # Extract the number after the last underscore, with error handling
        df['sort_2'] = pd.to_numeric(
            df['gene_id'].str.extract(r'_(\d+)(?!.*_\d+)')[0], 
            errors='coerce'
        ).fillna(-1).astype(int)  # Convert NaN to -1

        # Sort the dataframe
        df = df.sort_values(by=['sort_1', 'sort_2'])

        df = df.drop(columns=['sort_1', 'sort_2'])
        print(df)
        df.to_csv(outfi, mode='w', index=False, sep="\t")

############    

end = time.perf_counter()
elapsed = end - start
print(f'Hmmsearch file filtered, time taken: {elapsed:.6f} seconds')


