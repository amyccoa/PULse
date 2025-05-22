
genome, = glob_wildcards("genomes/{genome}.fa")

#~#~#~#~#~#~#~#~#

rule all:
    input:
        expand("genomes_formatted/formatted_{genome}.fa", genome=genome),
        expand("proteins/{genome}.faa", genome=genome), 
        expand("proteins/{genome}.gff", genome=genome), 
        expand("annotations/{genome}.out", genome=genome), 
        expand("cazymes/{genome}_cazyme.out", genome=genome), 
        expand("cazymes_filtered/{genome}_cazymes_f.out", genome=genome), 
        expand("proteins_filtered/{genome}_proteins_f.out", genome=genome),
        expand("markers/{genome}_markers.out", genome=genome)

#~#~#~#~#~#~#~#~#

rule reformat_genome:
    input:
        "genomes/{genome}.fa"
    output:
        formatted = 'genomes_formatted/formatted_{genome}.fa'
    shell:
        'python scripts/genome_formatter.py {input} {output.formatted}' 

rule protein_identify:
    input:
        "genomes_formatted/formatted_{genome}.fa"
    output:
        faa='proteins/{genome}.faa',
        gff='proteins/{genome}.gff'
    shell:
        'pyrodigal -i {input} -a {output.faa} -c -m -g 11 -p single -f gff > {output.gff}'

rule protein_annotate:
    input:
        "proteins/{genome}.faa"
    output:
        annotations='annotations/{genome}.out'
    shell: 
        'hmmsearch --cut_ga -o {output} databases/Pfam-A.hmm.gz {input}'

rule cazyme_annotate:
    input:
        "proteins/{genome}.faa"
    output:
        cazymes='cazymes/{genome}_cazyme.out'
    shell:
        'hmmsearch -o {output} databases/dbCAN-HMMdb-V13.txt {input}'

rule cazyme_filter:
    input:
        "cazymes/{genome}_cazyme.out"
    output:
        cazymes_f='cazymes_filtered/{genome}_cazymes_f.out'
    shell:
        'python scripts/parse_cazymes.py {input} {output}'

rule protein_filter:
    input:
        "annotations/{genome}.out"
    output:
        proteins_f='proteins_filtered/{genome}_proteins_f.out'
    shell:
        'python scripts/parse_proteins.py {input} {output}'

rule marker_spotter:
    input:
        "proteins_filtered/{genome}_proteins_f.out"
    output:
        markers='markers/{genome}_markers.out'
    shell:
        'python scripts/marker_spotter.py {input} {output}'

#rule PUL_spotter:
 #   input:
  #      "markers/{genome}_markers.out"
   # output:
    #    puls='puls/{genome}_puls.out'
    #shell:
     #   'python scripts/pul_mapper.py {input} {output}'
