rule run_abricate:
    input: fasta_file
    output: "results/abricate/{db}/{sample}.tab"
    log: "log/abricate/{db}/{sample}.txt"
    params:
        identity = config['abricate']['identity'],
        coverage = config['abricate']['coverage']
    conda: "../envs/abricate.yaml"
    shell:
        """
        abricate \
            --db {wildcards.db} \
            --minid {params.identity} \
            --mincov {params.coverage} \
            --nopath \
            {input} \
        > {output} \
        2> {log} 
        """

rule run_mash_matrix:
    input: gather_fasta_files
    output: "results/mash/matrix.tab"
    params:
        s = config["mash"]["s"],
        k = config["mash"]["k"]
    conda: "../envs/mash.yaml"
    log: "log/mash/log.txt"
    shell:
        """
        mash triangle \
            -k {params.k} \
            -s {params.s} \
            {input} \
        > {output} \
        2> {log}
        """

rule run_mash_edgelist:
    input: gather_fasta_files
    output: "results/mash/edgelist.tab"
    params:
        s = config["mash"]["s"],
        k = config["mash"]["k"]
    conda: "../envs/mash.yaml"
    log: "log/mash/log.txt"
    shell:
        """
        mash triangle \
            -E \
            -k {params.k} \
            -s {params.s} \
            {input} \
        > {output} \
        2> {log}
        """

rule run_mash_synteny:
    input: rivm_synteny_input
    output: "results/mash/synteny_matrix.tab"
    params:
        s = config["mash"]["s"],
        k = config["mash"]["k"]
    conda: "../envs/mash.yaml"
    log: "log/mash/log.txt"
    shell:
        """
        mash triangle \
            -k {params.k} \
            -s {params.s} \
            {input} \
        > {output} \
        2> {log}
        """

rule run_pling:
    localrule: True
    input: file_of_fasta_files
    output: directory("results/pling")
    conda: "../envs/pling.yaml"
    log: "log/pling/pling.txt"
    params: profile = "workflow/profile/pling_lsf"
    shell:
        """
        pling \
            {input} \
            {output} \
            align \
            --profile {params.profile} \
        2>&1> {log}
        """

rule run_tblastn:
    input:
        subject = fasta_file,
        query = "resources/card/ARO_3000202-protein.fasta"
    output: "results/blast/tblastn/{sample}.csv"
    params: 
        evalue = config["blast"]["tblastn"]["evalue"],
        qcov_hsp = config["blast"]["tblastn"]["qcov_hsp"],
        max_hsps = config["blast"]["tblastn"]["max_hsps"]
    conda: "../envs/blast.yaml"
    log: "log/blast/tblastn/{sample}.txt"
    shell:
        """
        tblastn \
            -task tblastn \
            -query {input.query} \
            -subject {input.subject} \
            -evalue {params.evalue} \
            -qcov_hsp_perc {params.qcov_hsp} \
            -max_hsps {params.max_hsps} \
            -out {output} \
            -outfmt '20 qseqid sseqid pident qcovhsp evalue' \
        2>&1> {log}
        """
rule run_blastn_pairwise:
    input:
        query = lambda wildcards: "results/rotated/fasta/{sample}.fasta".format(sample=wildcards.sample_a),
        subject = lambda wildcards: "results/rotated/fasta/{sample}.fasta".format(sample=wildcards.sample_b)
    output: "results/synteny/blast/{sample_a}-{sample_b}.csv"
    params: 
        evalue = config["blast"]["blastn"]["evalue"],
        perc_identity = config["blast"]["blastn"]["perc_identity"]
    conda: "../envs/blast.yaml"
    log: "log/blast/blastn/{sample_a}-{sample_b}.txt"
    shell:
        """
        blastn \
            -task blastn \
            -query {input.query}\
            -subject {input.subject}\
            -evalue {params.evalue} \
            -perc_identity {params.perc_identity} \
            -out {output} \
            -outfmt '20 qseqid qlen qstart qend qstrand sseqid slen sstart send sstrand' \
        2>&1> {log}
        """

rule run_mobtyper:
    input: fasta_file
    output: "results/mobtyper/typing/{sample}.tab"
    conda: "../envs/mob-suite.yaml"
    log: "log/mobtyper/typing/{sample}.txt"
    resources:
        runtime_min = 59
    shell:
        """
        mob_typer \
            --infile {input} \
            --out_file {output} \
        2>&1> {log}
        """
rule run_mobtyper_biomarkers:
    input: "results/rotated/fasta/{sample}.fasta"
    output:
        biomarkers = "results/mobtyper/biomarkers/{sample}.tab",
        out = "results/mobtyper/biomarkers/typing/{sample}.tab",
    conda: "../envs/mob-suite.yaml"
    log: "log/mobtyper/biomarkers/{sample}.txt"
    shell:
        """
        mob_typer \
            --infile {input} \
            --biomarker_report_file {output.biomarkers} \
            --out_file {output.out} \
        2>&1> {log}
        """

rule run_bakta:
    input: fasta_file
    output: "results/bakta/{sample}/{sample}.gbff"
    params:
        db = config["bakta"]["db"],
        species = get_species
    threads: config["bakta"]["threads"]
    resources:
        mem_mb = config["bakta"]["mem_mb"],
        runtime_min = config["bakta"]["runtime_min"]
    log: "log/bakta/{sample}.txt"
    conda: "../envs/bakta.yaml"
    shell:
        """
        bakta \
            {input} \
            --db {params.db} \
            --threads {threads} \
            --output results/bakta/{wildcards.sample} \
            --force \
            --species {params.species:q} \
            --complete \
            --keep-contig-headers \
        2>&1> {log}
        """ 