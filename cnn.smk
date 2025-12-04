from pathlib import Path
GENOME = config.get("genome", "sample")
DIAMOND_DB = config.get("diamond_db", "eggnog_proteins.dmnd")
EGGNOG_DB = config.get("eggnog_db", "eggnog.db")
EGGNOG_DATA_DIR = config.get("eggnog_data_dir", "./")
EGGNOG_2_CNN = config.get("eggnog_2_cnn", "eggnog2cnnPredictions.py")
GENEFAMILIES = config.get("gene_families", "genefamilies.txt")
MODEL_FILES_DIR = config.get("model_files_dir", "./")
ANTIBIOTICS = [
    "amikacin", "amoxicillin", "ampicillin", "aztreonam", "cefepime", "ceftriaxone",
    "chloramphenicol", "ciprofloxacin", "clindamycin", "colistin", "doripenem",
    "ertapenem", "erythromycin", "fosfomycin", "gentamicin", "imipenem", "levofloxacin",
    "meropenem", "moxifloxacin", "nitrofurantoin", "tetracycline", "tigecycline", "tobramycin"
]

rule all:
    input:
        f"{GENOME}_results.tar.gz"


rule preprocess:
    input:
        fasta = f"{GENOME}.fna"
    output:
        faa = f"{GENOME}.faa"
    singularity:
        "/home/lucydillon.linux/ldillon_amrwebsite_amrmlpipeline_3_arm64.sif"
    shell:
        """
        prodigal -i {input.fasta} -a {output.faa} -f gff
        """

rule diamond:
    input:
        faa = f"{GENOME}.faa",
        db = DIAMOND_DB
    output:
        diamond = f"{GENOME}_Diamond.tsv"
    singularity:
        "/home/lucydillon.linux/ldillon_amrwebsite_amrmlpipeline_3_arm64.sif"
    resources:
        mem_mb = 8000
    threads: 1
    shell:
        """
        diamond blastp \
            --db {input.db} \
            --query {input.faa} \
            --out {output.diamond} \
            --outfmt 6 \
            --threads {threads}
        """

rule edit_diamond:
    input:
        f"{GENOME}_Diamond.tsv"
    output:
        f"{GENOME}_DiamondReduced.tsv"
    shell:
        """
        awk 'BEGIN{{prev = ""}} {{if ($1 != prev) {{ prev=$1; print $1"\\t"$2"\\t"$11"\\t"$12 }}}}' {input} > {output}
        """

rule eggnog:
    input:
        faa = f"{GENOME}_DiamondReduced.tsv"
    output:
        f"{GENOME}.emapper.annotations"
    params:
        prefix = GENOME,
        data_dir = EGGNOG_DATA_DIR
    resources:
        mem_mb = 6000
    singularity:
        "/home/lucydillon.linux/ldillon_amrwebsite_amrmlpipeline_3_arm64.sif"
    shell:
        """
        emapper.py --data_dir {params.data_dir} \
                   --annotate_hits_table {input.faa} \
                   --no_file_comments \
                   --override \
                   --output {params.prefix}
        """

rule eggnog2cnn:
    input:
        eggnog_output = f"{GENOME}.emapper.annotations",
        script = EGGNOG_2_CNN,
        file = GENEFAMILIES
    output:
        expand(f"CNN_predictions_{{ab}}_{GENOME}.txt", ab=ANTIBIOTICS)
    params:
        db_dir = config.get("original_model_files_dir", MODEL_FILES_DIR),
        antibiotics = ANTIBIOTICS
    singularity: 
        "/home/lucydillon.linux/lucyd_machinelearning_mlpackages.sif"
    shell:
        r"""
        for ab in {params.antibiotics}; do
            db="{params.db_dir}/${{ab}}_egg_model_K_fold_bin.h5"
            out="CNN_predictions_${{ab}}_{GENOME}.txt"
            echo "Running eggnog2cnn for $ab using $db"
            python3 {input.script} -file {input.eggnog_output} -m "$db" -o "$out" || echo "Failed: $ab" >&2
        done
        """

rule make_signature:
    input:
        fasta = f"{GENOME}.fna"
    output:
        sig = f"{GENOME}.sig"
    singularity:
        "/home/lucydillon.linux/ldillon_amrwebsite_amrmlpipeline_3_arm64.sif"
    shell:
        """
        sourmash sketch dna -p k=31 {input.fasta} -o {output.sig}
        """

rule sourmash_search:
    input:
        sig = f"{GENOME}.sig"
    output:
        expand(f"{GENOME}_{{ab}}.sm.txt", ab=ANTIBIOTICS)
    params:
        db_dir = config.get("original_model_files_dir", MODEL_FILES_DIR),
        antibiotics = ANTIBIOTICS
    singularity:
        "/home/lucydillon.linux/ldillon_amrwebsite_amrmlpipeline_3_arm64.sif"
    shell:
        r"""
        for ab in {params.antibiotics}; do
            db="{params.db_dir}/${{ab}}_training.sbt.zip"
            out="{GENOME}_${{ab}}.sm.txt"
            echo "Running sourmash search for $ab using $db"

            if [ -f "$db" ]; then
                sourmash search {input.sig} "$db" --threshold 0.0 > "$out" 2>&1 || \
                echo "Search failed" > "$out"
            else
                echo "Database not found: $db" > "$out"
            fi
        done
        """


rule package_results:
    input:
        test = f"{GENOME}.emapper.annotations",
        sourmash = expand(f"{GENOME}_{{ab}}.sm.txt", ab=ANTIBIOTICS),
        antibiotics = expand(f"CNN_predictions_{{ab}}_{GENOME}.txt", ab=ANTIBIOTICS)
    output:
        archive = f"{GENOME}_results.tar.gz"
    shell:
        """
        tar -czf {output.archive} {input}
        echo "Results packaged in: {output.archive}"
        """
