singularity: "./ldillon_amrwebsite.sif"

GENOMES, = glob_wildcards("static/files/{genome}.fna")

rule all:
    input:
        expand("static/files/results/{genome}.tar.gz", genome=GENOMES)


rule prodigal:
    input:
        "static/files/{genome}.fna"
    output:
        "static/files/decision_tree/{genome}.faa"
    singularity:
        "ldillon_amrwebsite.sif"
    shell:
        '''
        prodigal -i {input} -a {output} -f gff
        '''

rule diamond:
    input:
        faa="static/files/decision_tree/{genome}.faa",
        db="data/eggnog_proteins.dmnd"
    output:
        "static/files/decision_tree/{genome}_Diamond.faa"
    shell:
        '''
        diamond blastp -d {input.db} -q {input.faa} --more-sensitive --threads 8 -e 0.001000 -o {output} --top 3 --query-cover 0 --subject-cover 0
        '''

rule Edit_diamond:
    input:
        "static/files/decision_tree/{genome}_Diamond.faa"
    output:
        "static/files/decision_tree/{genome}_DiamondReduced.faa"
    run:
        shell("""cat {input} | awk 'BEGIN{{prev = ""}};{{if ($1 != prev) {{ prev=$1; print $1"\t"$2"\t"$11"\t"$12 }} }}' > {output}""")


rule eggnog:
    input:
        faa="static/files/decision_tree/{genome}_DiamondReduced.faa",
        db="data/"
    output:
        "static/files/decision_tree/{genome}.emapper.annotations"
    params:
        prefix="{genome}"
    shell:
        """
        emapper.py --data_dir {input.db} --annotate_hits_table {input.faa}  --no_file_comments --override --output {params.prefix} --cpu 8;
        mv {params.prefix}.emapper.annotations static/files/decision_tree/
        """

rule eggnog_2_arff:
    input:
        "static/files/decision_tree/{genome}.emapper.annotations"
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
        temp("static/files/decision_tree/{genome}_names.txt")
    run:
        shell("""ls {input} > {output}; echo >> {output}""")

rule run_c_program:
    input:
        dot_ami="dot_files/amikacin_EGGNOG_NEW.arff.dot",
        dot_amo="dot_files/amoxicillin_EGGNOG_NEW.arff.dot",
        dot_amp="dot_files/ampicillin_EGGNOG_NEW.arff.dot",
        dot_az="dot_files/aztreonam_EGGNOG_NEW.arff.dot",
        dot_cip="dot_files/ciprofloxacin_EGGNOG_NEW.arff.dot",
        dot_clin="dot_files/clindamycin_EGGNOG_NEW.arff.dot",
        dot_col="dot_files/colistin_EGGNOG_NEW.arff.dot",
        dot_dor="dot_files/doripenem_EGGNOG_NEW.arff.dot",
        dot_ert="dot_files/ertapenem_EGGNOG_NEW.arff.dot",
        dot_gen="dot_files/gentamicin_EGGNOG_NEW.arff.dot",
        dot_mer="dot_files/meropenem_EGGNOG_NEW.arff.dot",
        dot_nit="dot_files/nitrofurantoin_EGGNOG_NEW.arff.dot",
        dot_tet="dot_files/tetracycline_EGGNOG_NEW.arff.dot",
        dot_tig="dot_files/tigecycline_EGGNOG_NEW.arff.dot",
        dot_tob="dot_files/tobramycin_EGGNOG_NEW.arff.dot",
        genome_names="static/files/decision_tree/{genome}_names.txt",
        arff_file="static/files/decision_tree/{genome}.arff",
    output:
        paths=temp("static/files/decision_tree/{genome}_names.txt.paths.txt"),
        predictions=temp("static/files/decision_tree/{genome}_names.txt.predicitions.txt"),
        predictions_ami=temp("static/files/decision_tree/{genome}_ami.txt.predicitions.txt"),
	predictions_amo=temp("static/files/decision_tree/{genome}_amo.txt.predicitions.txt"),
	predictions_amp=temp("static/files/decision_tree/{genome}_amp.txt.predicitions.txt"),
        predictions_az=temp("static/files/decision_tree/{genome}_az.txt.predicitions.txt"),
	predictions_cip=temp("static/files/decision_tree/{genome}_cip.txt.predicitions.txt"),
	predictions_clin=temp("static/files/decision_tree/{genome}_clin.txt.predicitions.txt"),
        predictions_col=temp("static/files/decision_tree/{genome}_col.txt.predicitions.txt"),
        predictions_dor=temp("static/files/decision_tree/{genome}_dor.txt.predicitions.txt"),
        predictions_ert=temp("static/files/decision_tree/{genome}_ert.txt.predicitions.txt"),
        predictions_gen=temp("static/files/decision_tree/{genome}_gen.txt.predicitions.txt"),
        predictions_mer=temp("static/files/decision_tree/{genome}_mer.txt.predicitions.txt"),
	predictions_nit=temp("static/files/decision_tree/{genome}_nit.txt.predicitions.txt"),
        predictions_tet=temp("static/files/decision_tree/{genome}_tet.txt.predicitions.txt"),
        predictions_tig=temp("static/files/decision_tree/{genome}_tig.txt.predicitions.txt"),
        predictions_tob=temp("static/files/decision_tree/{genome}_tob.txt.predicitions.txt"),
    shell:
        """
        my_apply_decision_tree {input.dot_ami} {input.arff_file} {input.genome_names};
        tail -1 {output.predictions} > {output.predictions_ami};
        my_apply_decision_tree {input.dot_amo} {input.arff_file} {input.genome_names};
        tail -1 {output.predictions} > {output.predictions_amo};
        my_apply_decision_tree {input.dot_amp} {input.arff_file} {input.genome_names};
        tail -1 {output.predictions} > {output.predictions_amp};
        my_apply_decision_tree {input.dot_az} {input.arff_file} {input.genome_names};
        tail -1 {output.predictions} > {output.predictions_az};
        my_apply_decision_tree {input.dot_cip} {input.arff_file} {input.genome_names};
        tail -1 {output.predictions} > {output.predictions_cip};
        my_apply_decision_tree {input.dot_clin} {input.arff_file} {input.genome_names};
        tail -1 {output.predictions} > {output.predictions_clin};
        my_apply_decision_tree {input.dot_col} {input.arff_file} {input.genome_names};
        tail -1 {output.predictions} > {output.predictions_col};
        my_apply_decision_tree {input.dot_dor} {input.arff_file} {input.genome_names};
        tail -1 {output.predictions} > {output.predictions_dor};
        my_apply_decision_tree {input.dot_ert} {input.arff_file} {input.genome_names};
        tail -1 {output.predictions} > {output.predictions_ert};
        my_apply_decision_tree {input.dot_gen} {input.arff_file} {input.genome_names};
        tail -1 {output.predictions} > {output.predictions_gen};
        my_apply_decision_tree {input.dot_mer} {input.arff_file} {input.genome_names};
        tail -1  {output.predictions} > {output.predictions_mer};
        my_apply_decision_tree {input.dot_nit} {input.arff_file} {input.genome_names};
        tail -1 {output.predictions} > {output.predictions_nit};
        my_apply_decision_tree {input.dot_tet} {input.arff_file} {input.genome_names};
        tail -1 {output.predictions} > {output.predictions_tet};
        my_apply_decision_tree {input.dot_tig} {input.arff_file} {input.genome_names};
        tail -1 {output.predictions} > {output.predictions_tig};
        my_apply_decision_tree {input.dot_tob} {input.arff_file} {input.genome_names};
        tail -1 {output.predictions} > {output.predictions_tob}
        """


rule combine_files:
    input:
        predictions_ami="static/files/decision_tree/{genome}_ami.txt.predicitions.txt",
	predictions_amo="static/files/decision_tree/{genome}_amo.txt.predicitions.txt",
	predictions_amp="static/files/decision_tree/{genome}_amp.txt.predicitions.txt",
        predictions_az="static/files/decision_tree/{genome}_az.txt.predicitions.txt",
        predictions_cip="static/files/decision_tree/{genome}_cip.txt.predicitions.txt",
	predictions_clin="static/files/decision_tree/{genome}_clin.txt.predicitions.txt",
        predictions_col="static/files/decision_tree/{genome}_col.txt.predicitions.txt",
        predictions_dor="static/files/decision_tree/{genome}_dor.txt.predicitions.txt",
        predictions_ert="static/files/decision_tree/{genome}_ert.txt.predicitions.txt",
        predictions_gen="static/files/decision_tree/{genome}_gen.txt.predicitions.txt",
        predictions_mer="static/files/decision_tree/{genome}_mer.txt.predicitions.txt",
	predictions_nit="static/files/decision_tree/{genome}_nit.txt.predicitions.txt",
        predictions_tet="static/files/decision_tree/{genome}_tet.txt.predicitions.txt",
        predictions_tig="static/files/decision_tree/{genome}_tig.txt.predicitions.txt",
        predictions_tob="static/files/decision_tree/{genome}_tob.txt.predicitions.txt"
    output:
        final="static/files/decision_tree/{genome}_FINAL.txt.predicitions.txt",
    run:
        shell("""echo "Antibiotic   Genome_name Genome_number   Node_number Node_label" > {output.final}; awk '{{print "Amikacin", $0}}' {input.predictions_ami} >> {output.final}; awk '{{print "Amoxicillin", $0}}' {input.predictions_amo} >> {output.final}; awk '{{print "Ampicillin", $0}}' {input.predictions_amp} >> {output.final}; awk '{{print "Aztreonam", $0}}' {input.predictions_az} >> {output.final}; awk '{{print "Ciprofloxacin", $0}}' {input.predictions_cip} >> {output.final}; awk '{{print "Clindamycin", $0}}' {input.predictions_clin} >> {output.final}; awk '{{print "Colistin", $0}}' {input.predictions_col} >> {output.final}; awk '{{print "Doripenem", $0}}' {input.predictions_dor} >> {output.final}; awk '{{print "Ertapenem", $0}}' {input.predictions_ert} >> {output.final}; awk '{{print "Gentamicin", $0}}' {input.predictions_gen} >> {output.final}; awk '{{print "Meropenem", $0}}' {input.predictions_mer} >> {output.final}; awk '{{print "Nitrofurantoin", $0}}' {input.predictions_nit} >> {output.final}; awk '{{print "Tetracycline", $0}}' {input.predictions_tet} >> {output.final}; awk '{{print "Tigecycline", $0}}' {input.predictions_tig} >> {output.final}; awk '{{print "Tobramycin", $0}}' {input.predictions_tob} >> {output.final}""")


rule zip_files:
    input:
        "static/files/decision_tree/{genome}_Diamond.faa",
	"static/files/decision_tree/{genome}_FINAL.txt.predicitions.txt",
        "static/files/decision_tree/{genome}.arff",
        "static/files/decision_tree/{genome}.emapper.annotations",
        "static/files/decision_tree/{genome}_DiamondReduced.faa",
    params:
        path="static/files/decision_tree"
    output:
        "static/files/{genome}.tar.gz"
    shell:
        """
        tar -czvf {output} {input};
        rm -r {params.path}
        """


rule send_results:
    input:
        file="static/files/{genome}.tar.gz"
    output:
        result="static/files/results/{genome}.tar.gz"
    shell:
        """
        mv {input.file} {output.result}
        """
