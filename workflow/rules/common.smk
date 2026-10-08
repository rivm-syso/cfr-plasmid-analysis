import re
from itertools import combinations

import pandas as pd


def rivm_abricate(wildcards):
    with open("resources/metadata/rivm_contigs.txt", "r") as file:
        samples = [line.strip() for line in file.readlines()]
    return expand("results/abricate/{db}/{sample}.tab", sample=samples, db="resfinder")

def gather_fasta_files(wildcards):
    rivm_samples = pd.read_csv(checkpoints.resistome_figure.get(**wildcards).output.fig_2c)["accession"].to_list()
    rivm_fasta = [f"resources/fasta/{sample}.fasta" for sample in rivm_samples]

    public_samples = pd.read_csv(checkpoints.get_ncbi_data.get(**wildcards).output.metadata)["accession"].to_list()
    public_fasta = [f"results/fasta/{sample}.fasta" for sample in public_samples]
    return rivm_fasta + public_fasta

def file_of_fasta_files(wildcards):
    rivm_samples = pd.read_csv(checkpoints.resistome_figure.get(**wildcards).output.fig_2c)["accession"].to_list()
    rivm_fasta = [f"resources/fasta/{sample}.fasta" for sample in rivm_samples]

    public_samples = pd.read_csv(checkpoints.get_ncbi_data.get(**wildcards).output.metadata)["accession"].to_list()
    public_fasta = [f"results/fasta/{sample}.fasta" for sample in public_samples]
    
    fasta_files = rivm_fasta + public_fasta
    
    fofn = "results/metadata/file_of_fasta_files.txt"
    with open(fofn, "w") as file:
        file.write("\n".join(fasta_files))
    return fofn

def fasta_file(wildcards):
    fasta_file = "{source}/fasta/{sample}.fasta".format(
        sample=wildcards.sample,
        source=("resources" if "RIVM" in wildcards.sample else "results")
        )
    return fasta_file

def gather_mobtyper(wildcards):
    rivm_samples = pd.read_csv(checkpoints.resistome_figure.get(**wildcards).output.fig_2c)["accession"].to_list()
    public_samples = pd.read_csv(checkpoints.get_ncbi_data.get(**wildcards).output.metadata)["accession"].to_list()

    samples = rivm_samples + public_samples
    return expand("results/mobtyper/typing/{sample}.tab", sample=samples)

def get_species(wildcards):
    if "RIVM" in wildcards.sample:
        rivm_metadata = pd.read_csv("resources/metadata/metadata.csv").set_index("sample_name")
        sample_name = re.search(r"(RIVM_[C|M]\d+)", wildcards.sample).group(0)
        species = rivm_metadata.loc[sample_name, "organism"]
    else:
        public_metadata = pd.read_csv(checkpoints.get_ncbi_data.get(**wildcards).output.metadata).set_index("accession")
        species = public_metadata.loc[wildcards.sample, "organism"]
    return species

def rivm_synteny_input(wildcards):
    rivm_samples = pd.read_csv(checkpoints.resistome_figure.get(**wildcards).output.fig_2c)["accession"].to_list()
    rivm_samples = [sample for sample in rivm_samples if sample.startswith("pRIVM_M")]
    return expand("resources/fasta/{sample}.fasta", sample=rivm_samples)


def rivm_links(wildcards):
    rivm_samples = pd.read_csv(checkpoints.resistome_figure.get(**wildcards).output.fig_2c)["accession"].to_list()
    rivm_samples = [sample for sample in rivm_samples if sample.startswith("pRIVM_M")]
    return ["results/synteny/blast/{sample_a}-{sample_b}.csv".format(sample_a=a, sample_b=b) for a, b in combinations(rivm_samples, 2)]

def rivm_genbank(wildcards):
    rivm_samples = pd.read_csv(checkpoints.resistome_figure.get(**wildcards).output.fig_2c)["accession"].to_list()
    rivm_samples = [sample for sample in rivm_samples if sample.startswith("pRIVM_M")]
    return expand("results/rotated/genbank/{sample}.gbff", sample=rivm_samples)

def rivm_biomarkers(wildcards):
    rivm_samples = pd.read_csv(checkpoints.resistome_figure.get(**wildcards).output.fig_2c)["accession"].to_list()
    rivm_samples = [sample for sample in rivm_samples if sample.startswith("pRIVM_M")]
    return expand("results/mobtyper/biomarkers/{sample}.tab", sample=rivm_samples)



