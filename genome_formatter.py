import os
import sys 

############
# This script processes genome files in FASTA format and adds a unique identifier to the header lines.
############

input_path = sys.argv[1]
output_path = sys.argv[2]

input_directory = "./genomes"  
output_directory = "./genomes_formatted"  

os.makedirs(output_directory, exist_ok=True)

for file in os.listdir(input_directory):
    if file.endswith(".fa"):
        genome_name = os.path.splitext(file)[0] 
        input_path = os.path.join(input_directory, file)
        output_path = os.path.join(output_directory, f"{file}")

        with open(input_path, "r") as infile, open(output_path, "w") as outfile:
            contig = ""
            contig_counter = 1
            for line in infile:
                if line.startswith(">"):
                    if contig:
                        outfile.write(f"{contig}\n")
                    contig_name = line.strip()[1:].split()[0]
                    new_header = f">{genome_name}|{contig_name}|{contig_counter}"
                    outfile.write(new_header + "\n")
                    contig = ""
                    contig_counter += 1
                else:
                    contig += line.strip()
            if contig:
                outfile.write(f"{contig}\n")

print("Formatting complete. Formatted files are saved in the output directory.")