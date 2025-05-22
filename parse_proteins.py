import sys
from itertools import islice
import pandas as pd
import time 

############

infi_prot = sys.argv[1]  #input protein identification file (.out)
outfi = sys.argv[2]  #output file 
hits = []
start = time.perf_counter()

############


with open(infi_prot, 'r') as fiprot:
    writethis = False
    query_id = None
    query_desc = None

    lines = iter(fiprot)  #convert to an iterator to allow peeking ahead

    for ln in lines:
        if ln.startswith("Query:"):
            query_id = ln.strip().split("Query:", 1)[-1].strip()
        
        if ln.startswith("Accession:"):
            accession = ln.strip().split("Accession:", 1)[-1].strip()

            #get the next two lines
            next(lines)
            query_desc_line = next(lines).strip()

            #extract the description text
            if query_desc_line.startswith("Description:"):
                query_desc = query_desc_line.split("Description:", 1)[-1].strip()
            else:
                query_desc = None  #in case there's no description
            
            #read next 6 lines to check for "[No hits detected]"
            check_lines = list(islice(lines, 6))
            if not any("[No hits detected that satisfy reporting thresholds]" in line.lower() for line in check_lines):
                writethis = True
                continue

        if writethis:
            if ln.startswith("Domain"):  #stops writing at "Domain"
                writethis = False
                continue 

            parts = ln.strip().split()

            if len(parts) > 8: #extracts parts of the hmm file
                try:
                    gene_id = parts[8] 
                    seq_score = parts[1]
                    gene_start = parts[10]
                    gene_end = parts[12]
                    hits.append([gene_id, accession, gene_start, gene_end, query_id, query_desc])
                    continue
                except (ValueError, IndexError):
                    pass

    if hits: #creates df with parts stored in hits list and sorts by gene id, first by the contig and then by the gene number
        df = pd.DataFrame(hits, columns=['gene_id', 'accession', 'gene_start', 'gene_end', 'query_id', 'query_desc'])

        # Extract sorting values with error handling
        df['sort_key_1'] = pd.to_numeric(
            df['gene_id'].str.extract(r'\|[^|]+\|(\d+)')[0], 
            errors='coerce'
        ).fillna(-1).astype(int)  # Convert NaN to -1

        # Extract the number after the last underscore, with error handling
        df['sort_key_2'] = pd.to_numeric(
            df['gene_id'].str.extract(r'_(\d+)(?!.*_\d+)')[0], 
            errors='coerce'
        ).fillna(-1).astype(int)  # Convert NaN to -1

        # Sort the dataframe
        df = df.sort_values(by=['sort_key_1', 'sort_key_2'])

        df = df.drop(columns=['sort_key_1', 'sort_key_2'])
        print(df)
        df.to_csv(outfi, mode='w', index=False, sep="\t")

############

end = time.perf_counter()
elapsed = end - start
print(f'Annotated all proteins, time taken: {elapsed:.6f} seconds')
