import sys
import time
import os
from collections import defaultdict

#this script will take the listed genes in each PUL from the output of pul_mapper and will locate the genes in the .faa files, pull them out and generate .faa files for each PUL for each genome (and generate folders for each genome to keep it organised)



def log(msg):
    now = time.strftime('%H:%M:%S')
    print(f"[{now}] {msg}")


def create_genome_directory(genome_name, base_dir="pul_fasta_sequences"):
    """create's the output file for puls fasta files for each genome ran"""
    genome_dir = os.path.join(base_dir, genome_name)
    os.makedirs(genome_dir, exist_ok=True)
    return genome_dir

def parse_pul_file(pul_file):
    """parses pul mapper output into dictionary"""
    puls = defaultdict(list)
    current_pul = None

    with open(pul_file) as f:
        for line in f:
            
            if "PUL_" in line:  
                try:
                    current_pul = line.split("PUL_")[1].split()[0].rstrip(",")
                    log(f"Found PUL: {current_pul}")
                except IndexError:
                    log(f"Warning: Could not parse PUL from line: {line.strip()}")
                    continue
            elif line.strip() and not line.startswith(("Upstream", "Downstream")):
                parts = [p.strip() for p in line.strip().split(",")]
                try:
                    if len(parts) >= 3:
                        gene_id = parts[0]
                        try:
                            start = int(parts[1])
                            end = int(parts[2])
                        except ValueError:
                            continue
                        
                        if current_pul:  
                            puls[current_pul].append((gene_id, start, end))
                            log(f"Added gene {gene_id} to PUL_{current_pul}")
                except Exception as e:
                    log(f"Error processing line: {line.strip()}\nError: {str(e)}")
                    continue
    
    for pul_id, genes in puls.items():
        log(f"PUL_{pul_id} contains {len(genes)} genes")
    
    return puls

def extract_sequence(genome_file, target_contig, start, end):
    """Extract sequence from genome file using contig ID and coordinates"""
    current_contig = ""
    sequence = ""
    found_contig = False
    
    log(f"Searching for contig: {target_contig}")
    
    with open(genome_file) as f:
        for line in f:
            if line.startswith(">"):
                header = line.strip()  # Save full header for debugging
                # Get the full contig ID up to the #
                current_contig = header[1:].split('#')[0].strip()
                log(f"Found header: {header}")
                log(f"Extracted contig: {current_contig}")
                # More flexible matching
                found_contig = target_contig in current_contig
                if found_contig:
                    log(f"Match found for contig: {current_contig}")
            elif found_contig:
                sequence += line.strip()
    
    if not found_contig:
        log(f"Warning: Contig {target_contig} not found in genome file")
    
    if sequence:
        try:
            extracted = sequence[start-1:end]
            log(f"Extracted sequence length: {len(extracted)}")
            return extracted
        except IndexError:
            log(f"Warning: Coordinates {start}-{end} out of range for contig {target_contig}")
            return None
    return None

if __name__ == "__main__":
    in_fasta, in_puls, outfi = sys.argv[1:4]
    start = time.perf_counter()

    genome_name = os.path.basename(in_fasta).split(".fa")[0]


    genome_dir = create_genome_directory(genome_name)


    puls = parse_pul_file(in_puls)


    for pul_id, genes in puls.items():
        output_file = os.path.join(genome_dir, f"PUL_{pul_id}.fasta")
        log(f"Processing PUL_{pul_id}")

        with open(output_file, 'w') as out:
            for gene_id, start, end in genes:
                log(f"Processing gene: {gene_id}")
                contig_id = gene_id.split("|")[0]
                log(f"Extracted contig ID: {contig_id}")
                sequence = extract_sequence(in_fasta, contig_id, start, end)
                if sequence:
                    out.write(f">{gene_id}, {start}-{end}\n")
                    for i in range(0, len(sequence), 60):
                        out.write(f"{sequence[i:i+60]}\n")
                else:
                    log(f"Warning: Could not extract sequence for {gene_id}")


    end = time.perf_counter()
    elapsed = end - start
    log(f"Completed pul sequence extraction in {elapsed:.2f} seconds")






# store genome name
# parse pul mapper