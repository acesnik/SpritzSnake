import os

rule directories:
    output: directory("ensembl/202122/")
    shell: "mkdir -p ensembl/202122/"

rule star_setup:
    output: "STAR-2.6.0c/bin/Linux_x86_64/STAR"
    shell:
        "wget https://github.com/alexdobin/STAR/archive/2.6.0c.tar.gz"
        "tar xvf 2.6.0c.tar.gz"
        "rm 2.6.0c.tar.gz"

rule star_genome_generate:
    input:
        star="STAR-2.7.0e/bin/Linux_x86_64/STAR",
        genomeDir=directory("ensembl/Homo_sapiens.GRCh38.dna.primary_assembly.karyotypic"),
        fa="ensembl/Homo_sapiens.GRCh38.dna.primary_assembly.karyotypic.fa",
        gff="ensembl/202122.gff3"
    output:
        "ensembl/202122/SA"
    shell:
        "{input.star} --runMode genomeGenerate --runThreadN {threads} --genomeDir {input.genomeDir} "
        "--genomeFastaFiles {input.fa} --sjdbGTFfile {input.gff} --sjdbGTFtagExonParentTranscript Parent --sjdbOverhang 100"

rule hisat_genome:
    input:
        fa="ensembl/202122.fa",
        gtf="ensembl/202122.gff3"
    output: "ensembl/202122.1.ht2"
    shell: "hisat2-build ensembl/202122.fa ensembl/202122"

rule hisat2_splice_sites:
    input: "ensembl/202122.gff3"
    output: "ensembl/202122.splicesites.txt"
    shell: "hisat2_extract_splice_sites.py {input} > {output}"

def input_fqs():
    if "sra" in config:
        if os.path.exists(expand("TestData/{sample}_1.fastq", sample=config["sra"])[0]):
            return expand(["TestData/{sample}_1.fastq", "TestData/{sample}_2.fastq"], sample=config["sra"])
        else:
            return expand(["TestData/{sample}.fastq"], sample=config["sra"])
    elif "fastq" in config:
        print("Not implemented: ""fastq"" in config not implemented, yet.")
        exit(1) # todo: implement this option for fastqs that are specified directly
    else:
        print("Error: FASTQs must be specified in the config file using ""sra"" or ""fastq.""")
        exit(1)

def input_fq_args():
    fqs=input_fqs()
    if len(fqs) == 1:
        return f"-U {fqs[0]}"
    else:
        return f"-1 {fqs[0]} -2 {fqs[1]}"

rule hisat2_align_bam:
    input:
        "ensembl/202122.1.ht2",
        fq=lambda wildcards: input_fqs(),
        ss="ensembl/202122.splicesites.txt"
    output:
        sorted="TestData/{sample}.sorted.bam",
    threads: 12
    params:
        compression="9",
        tempprefix="TestData/{sample}.sorted"
    log: "TestData/{sample}.hisat2.log"
    run:
        infq = input_fq_args()
        shell(
            "(hisat2 -p {threads} -x ensembl/202122 {infq} --known-splicesite-infile {input.ss} | " # align the suckers
            "samtools view -h -F4 - | " # get mapped reads only
            "samtools sort -l {params.compression} -T {params.tempprefix} -o {output.sorted} -) 2> {log} && " # sort them
            "samtools index {output}")
