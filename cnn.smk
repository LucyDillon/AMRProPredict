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


# Now convert the eggNOG output to an appropriate format for python

# load CNN models
# test genome on CNN 

rule send_results:
    input:
        file="static/files/diamond_annotation/{genome}_DiamondReduced.faa"  
    params:
        name="static/files/diamond_annotation/{genome}.txt"
    output:
        result="static/files/results/{genome}.txt"
    shell:
        """
        mv {input.file} {params.name};
        echo "CNN model smk" >> {params.name};
        mv {params.name} {output.result}
        """
