import os

rule directories:
    output: directory("data/ensembl/Homo_sapiens.GRCh38.dna.primary_assembly.karyotypic/")
    shell: "mkdir -p data/ensembl/Homo_sapiens.GRCh38.dna.primary_assembly.karyotypic/"

rule star_setup:
    output: "STAR-2.6.0c/bin/Linux_x86_64/STAR"
    shell:
        "wget https://github.com/alexdobin/STAR/archive/2.6.0c.tar.gz"
        "tar xvf 2.6.0c.tar.gz"
        "rm 2.6.0c.tar.gz"

rule star_genome_generate:
    input:
        star="STAR-2.7.0e/bin/Linux_x86_64/STAR",
        genomeDir=directory("data/ensembl/Homo_sapiens.GRCh38.dna.primary_assembly.karyotypic"),
        fa="data/ensembl/Homo_sapiens.GRCh38.dna.primary_assembly.karyotypic.fa",
        gff="data/ensembl/Homo_sapiens.GRCh38.81.gff3"
    output:
        "data/ensembl/Homo_sapiens.GRCh38.dna.primary_assembly.karyotypic/SA"
    shell:
        "{input.star} --runMode genomeGenerate --runThreadN {threads} --genomeDir {input.genomeDir} "
        "--genomeFastaFiles {input.fa} --sjdbGTFfile {input.gff} --sjdbGTFtagExonParentTranscript Parent --sjdbOverhang 100"

rule hisat_genome:
    input:
        fa="data/ensembl/Homo_sapiens.GRCh38.dna.primary_assembly.karyotypic.fa",
        gtf="data/ensembl/Homo_sapiens.GRCh38.81.gff3"
    output: "data/ensembl/Homo_sapiens.GRCh38.dna.primary_assembly.karyotypic.1.ht2"
    shell: "hisat2-build data/ensembl/Homo_sapiens.GRCh38.dna.primary_assembly.karyotypic.fa data/ensembl/Homo_sapiens.GRCh38.dna.primary_assembly.karyotypic"

rule hisat2_splice_sites:
    input: "data/ensembl/Homo_sapiens.GRCh38.81.gff3"
    output: "data/ensembl/Homo_sapiens.GRCh38.81.splicesites.txt"
    shell: "hisat2_extract_splice_sites.py {input} > {output}"

def input_fq_args(fastqs):
    fqs=fastqs.split()
    if len(fqs) == 1:
        return f"-U {fqs[0]}"
    else:
        return f"-1 {fqs[0]} -2 {fqs[1]}"

def check_sra():
    if 'sra' in config and config["sra"] is not None:
        if len(config["sra"]) > 0:
            return True
    return False

rule hisat2_align_bam:
    input:
        "data/ensembl/Homo_sapiens.GRCh38.dna.primary_assembly.karyotypic.1.ht2",
        fq1="{dir}/{sra}_1.fastq" if check_sra() is True else expand("{{dir}}/{fq1}_1.fastq", fq1=config["fq1"]),
        fq2="{dir}/{sra}_2.fastq" if check_sra() is True else expand("{{dir}}/{fq2}_2.fastq", fq2=config["fq2"]),
        ss="data/ensembl/Homo_sapiens.GRCh38.81.splicesites.txt"
    output:
        sorted="{dir}/{sra}.sorted.bam" if check_sra() is True else "{dir}/{fq1}.sorted.bam",
    threads: 12
    params:
        compression="9",
        tempprefix="{dir}/{sra}.sorted" if check_sra() is True else "{dir}/{fq1}.sorted",
    log: "{dir}/{sra}.hisat2.log" if check_sra() is True else "{dir}/{fq1}.hisat2.log"
    shell:
        "(hisat2 -p {threads} -x data/ensembl/Homo_sapiens.GRCh38.dna.primary_assembly.karyotypic -1 {input.fq1} -2 {input.fq2} --known-splicesite-infile {input.ss} | " # align the suckers
        "samtools view -h -F4 - | " # get mapped reads only
        "samtools sort -l {params.compression} -T {params.tempprefix} -o {output.sorted} -) 2> {log} && " # sort them
        "samtools index {output}"

rule hisat2_merge_bams:
    input:
        bams=expand("{{dir}}/{sra}.sorted.bam", sra=config["sra"]) if check_sra() is True else expand("{{dir}}/{fq1}.sorted.bam", fq1=config["fq1"])
    output:
        sorted="{dir}/combined.sorted.bam"
    params:
        compression="9",
        tempprefix="{dir}/combined.sorted"
    log: "{dir}/combined.sorted.log"
    shell:
        "(ls {input.bams} | "
        "{{ read firstbam; "
        "samtools view -h ""$firstbam""; "
        "while read bam; do samtools view ""$bam""; done; }} | "
        "samtools view -ubS - | "
        "samtools sort -l {params.compression} -T {params.tempprefix} -o {output.sorted} -) 2> {log} && "
        "samtools index {output.sorted}"
