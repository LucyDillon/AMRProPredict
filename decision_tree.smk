singularity: "./ldillon_amrwebsite.sif"

GENOMES, = glob_wildcards("static/files/{genome}.fna")

rule all:
    input:
        expand("static/files/results/{genome}.txt", genome=GENOMES)


rule prodigal:
    input:
        "static/files/{genome}.fna"
    output:
        temp("static/files/prodigal_output/{genome}.faa")
    singularity:
        "ldillon_amrwebsite.sif"
    shell:
        '''
        prodigal -i {input} -a {output} -f gff
        '''

rule diamond:
    input:
        faa="static/files/prodigal_output/{genome}.faa",
        db="data/eggnog_proteins.dmnd"
    output:
        temp("static/files/diamond_annotation/{genome}_Diamond.faa")
    shell:
        '''
        diamond blastp -d {input.db} -q {input.faa} --more-sensitive --threads 8 -e 0.001000 -o {output} --top 3 --query-cover 0 --subject-cover 0
        '''

rule Edit_diamond:
    input:
        "static/files/diamond_annotation/{genome}_Diamond.faa"
    output:
        temp("static/files/diamond_annotation/{genome}_DiamondReduced.faa")
    run:
        shell("""cat {input} | awk 'BEGIN{{prev = ""}};{{if ($1 != prev) {{ prev=$1; print $1"\t"$2"\t"$11"\t"$12 }} }}' > {output}""")


rule eggnog:
    input:
        faa="static/files/diamond_annotation/{genome}_DiamondReduced.faa",
        db="data/"
    output:
        temp("static/files/{genome}.emapper.annotations")
    params:
        prefix="{genome}"
    shell:
        """
        emapper.py --data_dir {input.db} --annotate_hits_table {input.faa}  --no_file_comments --override --output {params.prefix} --cpu 8;
        mv {params.prefix}.emapper.annotations static/files/
        """

rule eggnog_2_arff:
    input:
        "static/files/{genome}.emapper.annotations"
    output:
        temp("static/files/decision_tree/{genome}.arff")
    shell:
        """
        python eggnog_output_to_arff.py -file {input} -o {output} -c_attr_name phenotype -catagory1 Susceptible -catagory2 Resistant
        """

rule genome_names:
    input:
        "static/files/{genome}.fna"
    output:
        temp("static/files/decision_tree/{genome}_names.txt")
    run:
        shell("""ls {input} > {output}; echo >> {output}""")


rule run_c_program: # This rule only works on Linux based systems :(
    input:
        dot_cip="dot_files/ciprofloxacin_EGGNOG_NEW.arff.dot",
        dot_gen="dot_files/gentamicin_EGGNOG_NEW.arff.dot",
        dot_col="dot_files/colistin_EGGNOG_NEW.arff.dot",
        dot_mer="dot_files/meropenem_EGGNOG_NEW.arff.dot",
        dot_az="dot_files/aztreonam_EGGNOG_NEW.arff.dot",
        genome_names="static/files/decision_tree/{genome}_names.txt",
        arff_file="static/files/decision_tree/{genome}.arff",
    output:
        paths=temp("static/files/decision_tree/{genome}_names.txt.paths.txt"),
        predictions=temp("static/files/decision_tree/{genome}_names.txt.predicitions.txt"),
        predictions_cip=temp("static/files/decision_tree/{genome}_cip.txt.predicitions.txt"),
        predictions_gen=temp("static/files/decision_tree/{genome}_gen.txt.predicitions.txt"),
        predictions_col=temp("static/files/decision_tree/{genome}_col.txt.predicitions.txt"),
        predictions_mer=temp("static/files/decision_tree/{genome}_mer.txt.predicitions.txt"),
        predictions_az=temp("static/files/decision_tree/{genome}_az.txt.predicitions.txt"),
    shell:
        """
        my_apply_decision_tree {input.dot_cip} {input.arff_file} {input.genome_names};
        tail -1 {output.predictions} > {output.predictions_cip};
        my_apply_decision_tree {input.dot_gen} {input.arff_file} {input.genome_names};
        tail -1 {output.predictions} > {output.predictions_gen};
        my_apply_decision_tree {input.dot_col} {input.arff_file} {input.genome_names};
       	tail -1 {output.predictions} > {output.predictions_col};
        my_apply_decision_tree {input.dot_mer} {input.arff_file} {input.genome_names};
      	tail -1  {output.predictions} > {output.predictions_mer};
        my_apply_decision_tree {input.dot_az} {input.arff_file} {input.genome_names};
      	tail -1 {output.predictions} > {output.predictions_az};
        """

rule combine_files:
    input:
        predictions_cip="static/files/decision_tree/{genome}_cip.txt.predicitions.txt",
        predictions_gen="static/files/decision_tree/{genome}_gen.txt.predicitions.txt",
        predictions_col="static/files/decision_tree/{genome}_col.txt.predicitions.txt",
        predictions_mer="static/files/decision_tree/{genome}_mer.txt.predicitions.txt",
        predictions_az="static/files/decision_tree/{genome}_az.txt.predicitions.txt"
    output:
        final=temp("static/files/decision_tree/{genome}_FINAL.txt.predicitions.txt")
    run:
        shell("""echo "Antibiotic	Genome_name	Genome_number	Node_number	Node_label" > {output.final}; awk '{{print "Ciprofloxacin", $0}}' {input.predictions_cip} >> {output.final}; awk '{{print "Gentamicin", $0}}' {input.predictions_gen} >> {output.final}; awk '{{print "Colistin", $0}}' {input.predictions_col} >> {output.final}; awk '{{print "Meropenem", $0}}' {input.predictions_mer} >> {output.final}; awk '{{print "Aztreonam", $0}}' {input.predictions_az} >> {output.final}""")


rule send_results:
    input:
        file="static/files/decision_tree/{genome}_FINAL.txt.predicitions.txt"
    params:
        name="static/files/decision_tree/{genome}.txt"
    output:
        result="static/files/results/{genome}.txt"
    shell:
        """
        mv {input.file} {params.name};
        echo "DT model smk" >> {params.name};
        mv {params.name} {output.result}
        """
