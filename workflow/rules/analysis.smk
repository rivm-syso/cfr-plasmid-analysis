checkpoint resistome_figure:
    input:
        resfinder = rivm_abricate,
        contig_metadata = "resources/metadata/contig_metadata.csv",
        isolate_metadata = "resources/metadata/metadata.csv"
    output: 
        fig_2a = "results/figures/fig2a_sau-resistome.png",
        fig_2b = "results/figures/fig2b_pmi-resistome.png",
        fig_2c = "results/metadata/rivm_plasmid_metadata.csv"
    log: notebook = "log/resistome_figure.ipynb"
    notebook: "../notebooks/resistome_figure.ipynb"

rule filter_plsdb_records:
    input:
        "resources/plsdb/plsdb_blastn.tsv",
        "resources/plsdb/plsdb_tblastn.tsv"
    output: "results/metadata/plsdb_records.txt"
    log: notebook = "log/filter_plsdb_records.ipynb"
    notebook: "../notebooks/filter_plsdb_records.ipynb"

checkpoint get_ncbi_data:
    input: "results/metadata/plsdb_records.txt"
    output:
        metadata = "results/metadata/public_plasmid_metadata.csv",
        fasta = directory("results/fasta")
    log: notebook = "log/get_ncbi_data.ipynb"
    notebook: "../notebooks/get_ncbi_data.ipynb"

rule mash_sensitivity:
    input: "results/mash/matrix.tab"
    output: "results/figures/sup-fig4_mash-sensitivity.png"
    log: notebook = "log/mash_sensitivity.ipynb"
    notebook: "../notebooks/mash_sensitivity.ipynb"

rule network_analysis:
    input:
        rivm_metadata = "results/metadata/rivm_plasmid_metadata.csv",
        public_metadata = "results/metadata/public_plasmid_metadata.csv",
        mash = "results/mash/edgelist.tab",
        pling = "results/pling"
    output:
        figure = "results/figures/fig5_network.png",
        node_table = "results/network_analysis/node_table.csv",
        ari = "results/network_analysis/adjusted_rand_indices.csv"
    params:
        mash_cluster_cutoff = config["mash"]["cluster_cutoff"],
        pling_containment_cutoff = config["pling"]["containment_cutoff"],
        pling_dcj_cutoff = config["pling"]["dcj_cutoff"]
    conda: "../envs/ggnetwork.yaml"
    script: "../scripts/network_analysis.R"

checkpoint mash_vs_pling:
    input: "results/network_analysis/node_table.csv",
    output: "results/network_analysis/mash-vs-pling.yaml"
    log: notebook = "log/mash-vs-pling.ipynb"
    notebook: "../notebooks/mash-vs-pling.ipynb"

rule distance_vs_size:
    input: 
        mash = "results/mash/edgelist.tab",
        pling = "results/pling",
        mobtyper = gather_mobtyper
    output: "results/figures/sup-fig8_distance-vs-size.png"
    log: notebook = "log/distance_vs_size.ipynb"
    notebook: "../notebooks/distance_vs_size.ipynb"

rule parse_genbank:
    input: "results/bakta/{sample}/{sample}.gbff"
    output: 
        genbank = "results/rotated/genbank/{sample}.gbff",
        fasta = "results/rotated/fasta/{sample}.fasta"
    log: notebook = "log/parse_genbank/{sample}/parse_genbank.ipynb"
    notebook: "../notebooks/parse_genbank.ipynb"

rule synteny_order:
    input: "results/mash/synteny_matrix.tab"
    output: "results/synteny/rivm/order.txt"
    log: notebook = "log/synteny_order_rivm.ipynb"
    notebook: "../notebooks/synteny_order_rivm.ipynb"
    
rule synteny_rivm_plasmids:
    input:
        genbank = rivm_genbank,
        links = rivm_links,
        biomarkers = rivm_biomarkers,
        order = "results/synteny/rivm/order.txt"
    output: "results/figures/fig3_synteny-rivm-plasmids.png"
    conda: "../envs/gggenomes.yaml"
    script: "../scripts/synteny_rivm_plasmids.R"

rule synteny_mash_vs_pling:
    input:
        genbank_a = lambda wildcards: gather_genbank(wildcards, "A"),
        links_a = lambda wildcards: gather_blastn(wildcards, "A"),
        biomarkers_a = lambda wildcards: gather_biomarkers(wildcards, "A"),
        genbank_b = lambda wildcards: gather_genbank(wildcards, "B"),
        links_b = lambda wildcards: gather_blastn(wildcards, "B"),
        biomarkers_b = lambda wildcards: gather_biomarkers(wildcards, "B"),
        genbank_c = lambda wildcards: gather_genbank(wildcards, "C"),
        links_c = lambda wildcards: gather_blastn(wildcards, "C"),
        biomarkers_c = lambda wildcards: gather_biomarkers(wildcards, "C")
    output: "results/figures/fig4_synteny-mash-vs-pling.png"
    conda: "../envs/gggenomes.yaml"
    script: "../scripts/synteny_mash-vs-pling.R"

rule generate_tables:
    input:
        node_table = "results/network_analysis/node_table.csv",
        public_metadata = "results/metadata/public_plasmid_metadata.csv",
        rivm_metadata = "results/metadata/rivm_plasmid_metadata.csv",
        resfinder = gather_resfinder,
        mobtyper = gather_mobtyper,
        tblastn = gather_tblastn
    output:
        metadata = "results/tables/merged_metadata.csv",
        resistance = "results/tables/resistance_profiles.csv",
        mobility = "results/tables/mobility.csv",
        cfr_count = "results/tables/cfr_count_table.csv",
        mash_clusters = "results/tables/mash_cluster_description.csv",
        pling_subcommunities = "results/tables/pling_subcommunity_description.csv"
    log: notebook = "log/generate_tables.ipynb"
    notebook: "../notebooks/generate_tables.ipynb"