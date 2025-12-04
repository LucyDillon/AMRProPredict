from pathlib import Path

GENOME = config.get("genome", "sample")
DIAMOND_DB = config.get("diamond_db")
EGGNOG_DB = config.get("eggnog_db")
EGGNOG_DATA_DIR = config.get("eggnog_data_dir", "./")
EGGNOG_2_ARFF = config.get("eggnog_2_arff")
TOOL_DIR = config.get("tool_dir", "./")
MODEL_FILES_DIR = config.get("model_files_dir", TOOL_DIR + "/model_files")
EGGNOG_2_CNN = config.get("eggnog_2_cnn", "eggnog2cnnPredictions.py")
GENEFAMILIES = config.get("gene_families", "genefamilies.txt")
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

rule eggnog_2_arff:
    input:
        annotations = f"{GENOME}.emapper.annotations",
        script = EGGNOG_2_ARFF
    output:
        f"{GENOME}.arff"
    shell:
        """
        python {input.script} -file {input.annotations} -o {output}
        """

rule genome_names:
    input:
        f"{GENOME}.fna"
    output:
        f"{GENOME}_names.txt"
    shell:
        """
        ls {input} > {output}
        echo >> {output}
        """

rule run_decision_tree:
    input:
        arff_file = f"{GENOME}.arff",
        genome_names = f"{GENOME}_names.txt",
        dot_file = lambda wildcards: f"{MODEL_FILES_DIR}/{wildcards.ab}_EGGNOG_NEW.arff.dot"
    output:
        predictions = f"{GENOME}_{{ab}}.predictions.txt"
    resources:
        mem_mb = 4000
    params:
        ab = lambda wildcards: wildcards.ab
    singularity:
        "/home/lucydillon.linux/ldillon_amrwebsite_amrmlpipeline_3_arm64.sif"
    shell:
        """
        set -x
        echo "Starting decision tree analysis for {params.ab}..."
        echo "Dot file location: {input.dot_file}"

        # Increase stack size for C programs
        ulimit -s unlimited 2>/dev/null || ulimit -s 65536

        # Check if dot file exists
        if [ ! -f {input.dot_file} ]; then
            echo "ERROR: Dot file not found: {input.dot_file}" >&2
            exit 1
        fi

        # Copy to local file
        cp {input.dot_file} ./model.dot

        # Normalize line endings
        sed -i 's/\r$//' ./model.dot 2>/dev/null || true

        echo "Running apply_decision_tree..."
        apply_decision_tree ./model.dot {input.arff_file} {input.genome_names}
        exit_code=$?

        echo "apply_decision_tree exit code: $exit_code"

        # Handle output
        if [ -f model.dot.predictions.txt ]; then
            mv model.dot.predictions.txt {output.predictions}
        elif [ -f {input.genome_names}.predictions.txt ]; then
            mv {input.genome_names}.predictions.txt {output.predictions}
        else
            echo "ERROR: Could not find predictions output file" >&2
            touch {output.predictions}
            exit 1
        fi

        rm -f ./model.dot
        exit $exit_code
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

rule package_results:
    input:
        predictions = expand(f"{GENOME}_{{ab}}.predictions.txt", ab=ANTIBIOTICS),
        cnn = expand(f"CNN_predictions_{{ab}}_{GENOME}.txt", ab=ANTIBIOTICS),
        sourmash = expand(f"{GENOME}_{{ab}}.sm.txt", ab=ANTIBIOTICS),
        eggnog = f"{GENOME}.emapper.annotations"
    output:
        archive = f"{GENOME}_results.tar.gz"
    shell:
        """
        tar -czf {output.archive} {input}
        echo "Results packaged in: {output.archive}"
        """
