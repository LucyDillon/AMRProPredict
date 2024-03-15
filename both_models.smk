singularity: "./amrmlpipeline_6.sif"
GENOMES, = glob_wildcards("static/files/{genome}.fna")

rule all:
    input:
        expand("static/files/results/{genome}.txt", genome=GENOMES)

rule prodigal:
    input:
        "static/files/{genome}.fna"
    output:
        temp("static/files/prodigal_output/{genome}.faa")
    shell:
        '''
        prodigal -i {input} -a {output} -f gff
        '''

rule diamond:
    input:
        faa=temp("static/files/prodigal_output/{genome}.faa"),  # Mark as temporary
        db="data/eggnog_proteins.dmnd"
    output:
        temp("static/files/diamond_annotation/{genome}_Diamond.faa")
    shell:
        '''
        diamond blastp -d {input.db} -q {input.faa} --more-sensitive --threads 2 -e 0.001000 -o {output} --top 3 --query-cover 0 --subject-cover 0
        '''

rule Edit_diamond:
    input:
        "static/files/diamond_annotation/{genome}_Diamond.faa"
    output:
        "static/files/diamond_annotation/{genome}_DiamondReduced.faa"
    run:
        shell("""cat {input} | awk 'BEGIN{{prev = ""}};{{if ($1 != prev) {{ prev=$1; print $1"\t"$2"\t"$11"\t"$12 }} }}' > {output}""")

rule eggnog:
    input:
	faa="static/files/diamond_annotation/{genome}_DiamondReduced.faa",
        db="data/"
    output:
	"{genome}.emapper.annotations"  # Correct the output file name here
    params:
	prefix="{genome}"
    shell:
	"""
	emapper.py --data_dir {input.db} --annotate_hits_table {input.faa}  --no_file_comments --override --output {params.prefix} --cpu 8
        """

rule eggnog_2_arff:
    input:
	"static/files/{genomes}.emapper.annotations"
    output:
	"static/files/decision_tree/{genome}.arff"
    shell:
	"""
	python eggnog_output_to_arff.py -file {input} -o {output} -c_attr_name phenotype -catagory1 Susceptible -catagory2 Resistant
        """

rule genome_names:
    input:
        "static/files/{genome}.fna"
    output:
        "static/files/decision_tree/{genome}_names.txt"
    run:
        shell("""ls {input} > {output}; echo >> {output}""")

rule run_c_program: # this rule only works on linux based systems :(
    input:
        dot_file="dot_files/ciprofloxacin_EGGNOG_NEW.arff.dot",
        genome_names="static/files/decision_tree/{genome}_names.txt",
        arff_file="static/files/decision_tree/{genome}.arff"
    output:
        predictions="static/files/decision_tree/{genome}_names.predictions.txt"
    shell:
        """
        my_apply_decision_tree {input.dot_file} {input.arff_file} {input.genome_names}
        """

rule send_results: # Figure out a way of combining the results of DT and CNN into one file :D
    input:
        file="static/files/decision_tree/{genome}_names.txt"  
    params:
        name="static/files/decision_tree/{genome}.txt"
    output:
        result="static/files/results/{genome}.txt"
    shell:
        """
        mv {input.file} {params.name};
        echo "Both models smk" >> {params.name};
        mv {params.name} {output.result}
        """
