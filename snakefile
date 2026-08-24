rule all:
    input:
        "results/origin_ligands.yaml",
        "results/blast.tsv",
        "results/top_10_hits.fasta",
        "results/boltz_inputs",
        "results/boltz_results",
        "results/bKeyCatalyticDescriptors.csv",
        "results/ranked_candidates.csv"

rule download_swissprot:
    output:
        touch("db/.swissprot_downloaded")
    conda:
        "env/blast.yaml"
    shell:
        """
        mkdir -p db
        curl -L https://ftp.ncbi.nlm.nih.gov/blast/db/swissprot.tar.gz -o db/swissprot.tar.gz
        tar -xzf db/swissprot.tar.gz -C db
        rm db/swissprot.tar.gz
        """

rule get_ligands:
    input:
        "data/input.fasta"
    output:
        "results/origin_ligands.yaml"
    conda:
        "env/ligands.yaml"
    script:
        "scripts/get_ligands.py"        

rule blast:
    input:
        query = "data/input.fasta",
        db_check = "db/.swissprot_downloaded"
    output:
        "results/blast.tsv"
    conda:
        "env/blast.yaml"
    shell:
        """
        mkdir -p results
        blastp -query {input.query} -db db/swissprot -max_target_seqs 10 -outfmt '6 qseqid sseqid pident length evalue bitscore stitle' -out {output}
        """

rule extract_top_hits_fasta:
    input:
        tsv="results/blast.tsv",
        db_check="db/.swissprot_downloaded"
    output:
        ids="results/top_10_hit_ids.txt",
        fasta="results/top_10_hits.fasta"
    conda:
        "env/blast.yaml"
    shell:
        """
        cut -f2 {input.tsv} > {output.ids}

        blastdbcmd \
            -db db/swissprot \
            -entry_batch {output.ids} \
            -out {output.fasta}
        """

rule prepare_boltz_inputs:
    input:
        hits="results/top_10_hits.fasta",
        ligands="results/origin_ligands.yaml"
    output:
        directory("results/boltz_inputs")
    conda:
        "env/ligands.yaml"
    script:
        "scripts/prepare_boltz_inputs.py"
        
rule run_boltz_all:
    input:
        "results/boltz_inputs"
    output:
        directory("results/boltz_results")
    conda:
        "env/boltz.yaml"
    threads: 16
    shell:
        """
        mkdir -p {output}

        for yaml in {input}/*.yaml; do
            name=$(basename "$yaml" .yaml)

            boltz predict "$yaml" \
                --out_dir {output}/"$name" \
                --model boltz2 \
                --accelerator cpu \
                --devices 1 \
                --num_workers {threads} \
                --preprocessing-threads {threads} \
                --use_msa_server \
                --no_kernels
        done
        """
        
rule extract_boltz_key_descriptors:
    input:
        boltz_results="results/boltz_results"
    output:
        "results/bKeyCatalyticDescriptors.csv"
    conda:
        "env/descriptors.yaml"
    script:
        "scripts/get_catalytic_descriptors.py"

rule rank_candidates:
    input:
        "results/bKeyCatalyticDescriptors.csv"
    output:
        "results/ranked_candidates.csv"
    conda:
        "env/descriptors.yaml"
    script:
        "scripts/rank_candidates.py"
