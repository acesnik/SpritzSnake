GENOME_VERSION = "GRCh38"
ENSEMBL_VERSION = "96"
GENEMODEL_VERSION = GENOME_VERSION + "." + ENSEMBL_VERSION
GENOME_FA = f"ensembl/Homo_sapiens.{GENOME_VERSION}.dna.primary_assembly.fa"
ENSEMBL_GFF = f"ensembl/Homo_sapiens.{GENEMODEL_VERSION}.gff3"
TEST_GENOME_FA = f"ensembl/202122.fa"
TEST_ENSEMBL_GFF = f"ensembl/202122.gff3"
FA=GENOME_FA # for analysis; can also be TEST_GENOME_FA
GFF3=ENSEMBL_GFF # for analysis; can also be TEST_ENSEMBL_GFF
REFSTAR_PREFIX = f"ensembl/Homo_sapiens.{GENEMODEL_VERSION}RsemStar/RsemStarReference"
REFSTAR_FOLDER = f"ensembl/Homo_sapiens.{GENEMODEL_VERSION}RsemStar/"
REF_PREFIX = f"ensembl/Homo_sapiens.{GENEMODEL_VERSION}Rsem/RsemReference"
REF_FOLDER = f"ensembl/Homo_sapiens.{GENEMODEL_VERSION}Rsem/"

configfile: "config.yaml"

def check_dir():
    if ('analysisDirectory' in config and config["analysisDirectory"] is not None) and len(config["analysisDirectory"]) > 0:
        return True
    return False

# def output(wildcards):
#      if check_dir():
#          return expand(["{dir}/combined.spritz.snpeff.protein.withmods.xml.gz", "{dir}/combined.spritz.noindels.snpeff.protein.withmods.xml.gz", "{dir}/combined.spritz.isoformvariants.protein.withmods.xml.gz", "{dir}/combined.spritz.noindels.isoformvariants.protein.withmods.xml.gz", "{dir}/combined.spritz.isoform.protein.withmods.xml.gz", "{dir}/GRCh38.86.protein.withmods.xml.gz"], dir=config["analysisDirectory"])
#      return expand(["output/combined.spritz.snpeff.protein.withmods.xml.gz", "output/combined.spritz.noindels.snpeff.protein.withmods.xml.gz", "output/combined.spritz.isoformvariants.protein.withmods.xml.gz", "output/combined.spritz.noindels.isoformvariants.protein.withmods.xml.gz", "output/combined.spritz.isoform.protein.withmods.xml.gz", "output/GRCh38.86.protein.withmods.xml.gz"])

def output(wildcards):
     if check_dir():
         return expand(["{dir}/combined.spritz.snpeff.protein.withmods.xml.gz", "{dir}/combined.spritz.isoformvariants.protein.withmods.xml.gz", "{dir}/combined.spritz.isoform.protein.withmods.xml.gz", "{dir}/GRCh38.86.protein.withmods.xml.gz"], dir=config["analysisDirectory"])
     return expand(["output/combined.spritz.snpeff.protein.withmods.xml.gz", "output/combined.spritz.isoformvariants.protein.withmods.xml.gz", "output/combined.spritz.isoform.protein.withmods.xml.gz", "output/GRCh38.86.protein.withmods.xml.gz"])

rule all:
    input:
        "data/combined.spritz.snpeff.protein.withmods.xml.gz",
        "data/combined.spritz.noindels.snpeff.protein.withmods.xml.gz",
        "data/combined.spritz.isoformvariants.protein.withmods.xml.gz",
        "data/combined.spritz.noindels.isoformvariants.protein.withmods.xml.gz",
        "data/combined.spritz.isoform.protein.withmods.xml.gz", # no variants
        "data/GRCh38.86.protein.withmods.xml.gz", # no variants
        # "clean_snpeff"

# expand(["{dir}/GRCh38.86.protein.xml.gz", "{dir}/GRCh38.86.protein.withmods.xml.gz"], dir=config["analysisDirectory"])

rule clean:
    shell:
        "rm -rf data/ ChromosomeMappings/ SnpEff/ tmp/ fast.tmp/ && "
        "cd GtfSharp && dotnet clean && cd .. && "
        "cd TransferUniProtModifications && dotnet clean && cd .."

include: "rules/downloads.smk"
include: "rules/align.smk"
include: "rules/variants.smk"
include: "rules/isoforms.smk"
include: "rules/proteogenomics.smk"
include: "rules/qc.smk"
include: "rules/quant.smk"
include: "rules/fusion.smk"
