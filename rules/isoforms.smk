rule filter_bam:
    '''Rule adapted from ProteomeGenerator'''
    input: "data/{sra}.sorted.bam" if check_sra() else "data/{fq}.sorted.bam",
    output: "data/{sra}.sorted.filtered.bam" if check_sra() else "data/{fq}.sorted.filtered.bam"
    benchmark: "data/{sra}.sorted.filtered.benchmark" if check_sra() else "data/{fq}.sorted.filtered.benchmark"
    log: "data/{sra}.sorted.filtered.log" if check_sra() else "data/{fq}.sorted.filtered.log"
    threads: 1
    shell:
        "(samtools view -b -h -F 4 -F 256 -F 512 -q 30 {input} > {output} && "
        "samtools index {output}) 2> {log}"

rule assemble_transcripts:
    '''Rule adapted from ProteomeGenerator'''
    input:
        bam="data/{sra}.sorted.bam" if check_sra() else "data/{fq}.sorted.bam",
        gff="data/ensembl/Homo_sapiens." + GENEMODEL_VERSION + ".gff3"
    output: "data/{sra}.sorted.gtf" if check_sra() else "data/{fq}.sorted.gtf"
    threads: 6
    log: "data/{sra}.sorted.gtf.log" if check_sra() else "data/{fq}.sorted.gtf"
    shell:
        "stringtie {input.bam} -p {threads} -G {input.gff} -o {output} -c 2.5 -m 300 -f .01 2> {log}" # strandedness: --fr for forwared or --rf for reverse

rule merge_transcripts:
    '''Rule adapted from ProteomeGenerator'''
    input:
        custom_gtfs=expand("data/{sra}.sorted.bam", sra=config["sra"]) if check_sra() is True else expand("data/{fq}.sorted.bam", fq=config["fq"]),
        gff="data/ensembl/Homo_sapiens." + GENEMODEL_VERSION + ".gff3"
    output: "data/combined.gtf"
    threads: 12
    log: "data/combined.gtf.log"
    shell:
        "strintie --merge -o {output} -c 2.5 -m 300 -T 1 -f .01 -p {threads} -i {input.custom_gtfs} 2> {log}"

rule convert2ucsc:
    '''Rule adapted from ProteomeGenerator'''
    input: "data/combined.gtf"
    output: "data/combined_ucsc.gtf"
    shell: "python script/convert_ensembl2ucsc.py {input} {output}"

rule gtf_file_to_cDNA_seqs:
    '''Rule adapted from ProteomeGenerator'''
    input: "data/combined_ucsc.gtf"
    output:
        fasta="data/combined.transcripts.fasta",
        gtf="data/combined.transcripts.gtf"
    benchmark: "data/combined.gtf_file_to_cDNA_seqs.benchmark"
    log: "data/combined.gtf_file_to_cDNA_seqs.log"
    threads: 1
    shell:
        "(gffread {input} -T -o {output.gtf} --no-pseudo --force-exons -M -Q && "
        "gffread -w {output.fasta} -g {FASTA} {output.gtf}) 2> {log}"

rule LongOrfs:
    '''Rule adapted from ProteomeGenerator'''
    input: "data/combined.transcripts.fasta"
    output: "data/combined.transcripts.fasta.transdecoder_dir/longest_orfs.pep"
    benchmark: "data/combined.LongOrfs.benchmark"
    log: "data/combined.LongOrfs.log"
    threads: 1
    shell: "TransDecoder.LongOrfs -t {input} -m 100 2> {log}"

rule makeblastdb:
    '''Rule adapted from ProteomeGenerator'''
    input: UNIPROTFASTA
    output: [UNIPROTFASTA+'.pin', UNIPROTFASTA+'.phr', UNIPROTFASTA+'.psq']
    benchmark: "data/makeblastdb.benchmark"
    log: "data/makeblastdb.log"
    threads: 1
    shell: "makeblastdb -in {UNIPROTFASTA} -dbtype prot 2> {log}"

rule blastp:
    '''Rule adapted from ProteomeGenerator'''
    input:
        fasta=UNIPROTFASTA,
        pep="data/combined.transcripts.fasta.transdecoder_dir/longest_orfs.pep",
        blastdb=[UNIPROTFASTA+'.pin', UNIPROTFASTA+'.phr', UNIPROTFASTA+'.psq']
    output: "data/combined.blastp.outfmt6"
    benchmark: "data/combined.blastp.benchmark"
    log: "data/combined.blastp.log"
    threads: 24
    shell: "blastp \
        -num_threads {threads} \
        -query {input.pep}  \
        -db {input.fasta}  \
        -max_target_seqs 1 \
        -outfmt 6 \
        -evalue 1e-5 \
        > {output} 2> {log}"

rule Predict:
    '''Rule adapted from ProteomeGenerator'''
    input:
        orfs="data/combined.transcripts.fasta.transdecoder_dir/longest_orfs.pep",
        fasta="data/combined.transcripts.fasta",
        blastp="data/combined.blastp.outfmt6"
    output:
        "data/combined.transcripts.fasta.transdecoder.pep",
        gff3="data/combined.transcripts.fasta.transdecoder.gff3"
    benchmark: "data/combined.Predict.benchmark"
    log: "data/combined.Predict.log"
    threads: 1
    shell:
        "TransDecoder.Predict -t {input.fasta} --single_best_only "
        "--retain_blastp_hits {input.blastp} 2> {log}"

rule gtf_to_alignment_gff3:
    '''Rule adapted from ProteomeGenerator'''
    input: "data/combined.transcripts.gtf"
    output: "data/combined/transcripts.gff3"
    benchmark: "data/combined.gtf_to_alignment_gff3.benchmark"
    log: "data/combined.gtf_to_alignment_gff3.log"
    threads: 1
    shell: "gtf_to_alignment_gff3.pl {input} > {output} 2> {log}"

rule cdna_alignment_orf_to_genome_orf:
    '''Rule adapted from ProteomeGenerator'''
    input:
        gff3="data/combined.transcripts.gff3",
        fasta_td="data/combined.transcripts.fasta",
        gff3_td="data/combined.transcripts.fasta.transdecoder.gff3"
    output: "data/combined.transcripts.genome.gff3"
    benchmark: "data/combined.cdna_alignment_orf_to_genome_orf.benchmark"
    log: "data/combined.cdna_alignment_orf_to_genome_orf.log"
    threads: 1
    shell: "cdna_alignment_orf_to_genome_orf.pl {input.gff3_td} {input.gff3} {input.fasta_td} > {output} 2> {log}"

rule gff3_file_to_bed:
    '''Rule adapted from ProteomeGenerator'''
    input: "data/combined.transcripts.genome.gff3"
    output: "data/combined.proteome.bed"
    benchmark: "data/combined.gff3_file_to_bed.benchmark"
    log: "data/combined.gff3_file_to_bed.log"
    threads: 1
    shell:
        "cat {input} | grep -P \"\tCDS\t\" | gffread --force-exons - -o- | gff3_file_to_bed.pl /dev/stdin | tail -n +2 > {output} 2> {log}"

rule gff3_file_to_proteins:
    '''Rule adapted from ProteomeGenerator'''
    input: "data/combined.transcripts.genome.gff3"
    output: "data/combined.proteome.fasta"
    benchmark: "data/combined.gff3_file_to_proteins.benchmark"
    log: "data/combined.gff3_file_to_proteins.log"
    threads: 1
    shell: "cat {input} | grep -P \"\tCDS\t\" | gffread --force-exons - -o- | gff3_file_to_proteins.pl --gff3 /dev/stdin --fasta {FA} | egrep -o '^[^*]+' > {output} 2> {log}"

rule reorderFASTA:
    '''Rule adapted from ProteomeGenerator
    Had to do this in Ubuntu 18.06: sudo ln -s /lib/x86_64-linux-gnu/libreadline.so.7 /lib/x86_64-linux-gnu/libreadline.so.6 '''
    input: "data/combined.proteome.fasta"
    output: "data/combined.proteome.unique.fasta"
    benchmark: "data/combined.reorderFASTA.benchmark"
    log: "data/combined.reorderFASTA.log"
    threads: 1
    script: "scripts/reorderFASTA.R"

rule generate_snpeff_database:
    input:
        jar="SnpEff/snpEff.jar",
        gtf="data/combined.transcripts.genome.gff3",
        pfa="data/ensembl/Homo_sapiens." + GENOME_VERSION + ".pep.all.fa",
        gfa="data/ensembl/Homo_sapiens." + GENOME_VERSION + ".dna.primary_assembly.karyotypic.fa"
    output:
        gtf="SnpEff/data/combined.transcripts.genome.gff3/genes.gff3",
        pfa="SnpEff/data/combined.transcripts.genome.gff3/protein.fa",
        gfa="SnpEff/data/genomes/combined.transcripts.genome.gff3.fa",
    params:
        ref="combined.transcripts.genome.gff3"
    resources:
        mem_mb=16000
    log:
        "data/combined.transcripts.genome.gff3.snpeffdatabase.log"
    shell:
        "cp {input.gtf} {output.gtf} && "
        "cp {input.pfa} {output.pfa} && "
        "cp {input.gfa} {output.gfa} && "
        "echo \"\n# {params.ref}\" >> SnpEff/snpEff.config && "
        "echo \"{params.ref}.genome : Human genome " + GENOME_VERSION + " using RefSeq transcripts\" >> SnpEff/snpEff.config && "
        "echo \"{params.ref}.reference : ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/\" >> SnpEff/snpEff.config && "
        "echo \"\t{params.ref}.M.codonTable : Vertebrate_Mitochondrial\" >> SnpEff/snpEff.config && "
        "echo \"\t{params.ref}.MT.codonTable : Vertebrate_Mitochondrial\" >> SnpEff/snpEff.config && "
        "(java -Xmx{resources.mem_mb}M -jar {input.jar} build -gff3 -v {params.ref}) 2> {log}"
