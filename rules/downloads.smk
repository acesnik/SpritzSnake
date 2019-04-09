rule download_genome_fasta:
    output: "ensembl/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
    shell:
        "wget -O - ftp://ftp.ensembl.org/pub/release-81//fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz | "
        "gunzip -c > ensembl/Homo_sapiens.GRCh38.dna.primary_assembly.fa"

rule download_gene_model:
    output: "ensembl/Homo_sapiens.GRCh38.81.gff3"
    shell:
        "wget -O - ftp://ftp.ensembl.org/pub/release-81/gff3/homo_sapiens/Homo_sapiens.GRCh38.81.gff3.gz | "
        "gunzip -c > ensembl/Homo_sapiens.GRCh38.81.gff3"

rule download_protein_fasta:
    output: "ensembl/Homo_sapiens.GRCh38.pep.all.fa"
    shell:
        "wget -O - ftp://ftp.ensembl.org/pub/release-81//fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz | "
        "gunzip -c > ensembl/Homo_sapiens.GRCh38.pep.all.fa"

rule download_common_known_variants:
    output:
        "ensembl/common_all_20170710.vcf",
        "ensembl/common_all_20170710.vcf.idx"
    shell:
        "wget -O - ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b150_GRCh38p7/VCF/GATK/common_all_20170710.vcf.gz | "
        "gunzip -c > ensembl/common_all_20170710.vcf;"
        "gatk IndexFeatureFile -F ensembl/common_all_20170710.vcf"

rule download_chromosome_mappings:
    output: "ChromosomeMappings/GRCh38_UCSC2ensembl.txt"
    shell: "git clone https://github.com/dpryan79/ChromosomeMappings.git"

rule reorder_genome_fasta:
    input: "ensembl/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
    output: "ensembl/Homo_sapiens.GRCh38.dna.primary_assembly.karyotypic.fa"
    script: "../scripts/karyotypic_order.py"

rule convert_ucsc2ensembl:
    input:
        "ensembl/common_all_20170710.vcf",
        "ChromosomeMappings/GRCh38_UCSC2ensembl.txt"
    output:
        "ensembl/common_all_20170710.ensembl.vcf",
    script:
        "../scripts/convert_ucsc2ensembl.py"

rule index_ucsc2ensembl:
    input: "ensembl/common_all_20170710.ensembl.vcf"
    output: "ensembl/common_all_20170710.ensembl.vcf.idx"
    shell: "gatk IndexFeatureFile -F {input}"

rule filter_gff3:
    input: "ensembl/Homo_sapiens.GRCh38.81.gff3"
    output: "ensembl/202122.gff3"
    shell: "grep \"^#\|20\|^21\|^22\" \"ensembl/Homo_sapiens.GRCh38.81.gff3\" > \"ensembl/202122.gff3\""

rule filter_fa:
    input: "ensembl/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
    output: "ensembl/202122.fa"
    script: "scripts/filter_fasta.py"

rule download_sras:
    # input:
    #     sample = lambda wildcards: config["sra"][wildcards.sample]
    output:
        "TestData/ERR315327_1.fastq"
        #"TestData/ERR315327_2.fastq"
    log:
        "TestData/ERR315327.log"
    threads: 4
    shell: # --outdir directory SRAaccession
        "fasterq-dump --progress --threads {threads} --split-files --outdir TestData ERR315327 2> {log}" # shell for downloading sras

# rule download_sras:
#     input:
#         sample = lambda wildcards: config["sra"][wildcards.sample]
#     output:
#         "TestData/{sample}_1.fastq"
#         "TestData/{sample}_2.fastq"
#     log:
#         "TestData/{sample}.log"
#     threads: 4
#     shell: # --outdir directory SRAaccession
#         "fasterq-dump --progress --threads {threads} --split-files --outdir TestData {wildcards.sample} 2> {log}" # shell for downloading sras
