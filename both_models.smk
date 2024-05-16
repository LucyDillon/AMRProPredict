GENOMES, = glob_wildcards("static/files/{genome}.fna")

rule all:
    input:
        expand("static/files/results/{genome}.tar.gz", genome=GENOMES)

rule prodigal:
    input:
        "static/files/{genome}.fna"
    output:
        "static/files/BOTH/{genome}.faa"
    singularity: "./ldillon_amrwebsite.sif"
    shell:
        '''
        prodigal -i {input} -a {output} -f gff
        '''

rule diamond:
    input:
        faa="static/files/BOTH/{genome}.faa",  # Mark as temporary
        db="data/eggnog_proteins.dmnd"
    output:
        "static/files/BOTH/{genome}_Diamond.faa"
    singularity: "./ldillon_amrwebsite.sif"
    shell:
        '''
        diamond blastp -d {input.db} -q {input.faa} --more-sensitive --threads 2 -e 0.001000 -o {output} --top 3 --query-cover 0 --subject-cover 0
        '''

rule Edit_diamond:
    input:
        "static/files/BOTH/{genome}_Diamond.faa"
    output:
        "static/files/BOTH/{genome}_DiamondReduced.faa"
    run:
        shell("""cat {input} | awk 'BEGIN{{prev = ""}};{{if ($1 != prev) {{ prev=$1; print $1"\t"$2"\t"$11"\t"$12 }} }}' > {output}""")

rule eggnog:
    input:
        faa="static/files/BOTH/{genome}_DiamondReduced.faa",
        db="data/"
    output:
        "static/files/BOTH/{genome}.emapper.annotations"
    params:
        prefix="{genome}"
    singularity: "./ldillon_amrwebsite.sif"
    shell:
        """
        emapper.py --data_dir {input.db} --annotate_hits_table {input.faa}  --no_file_comments --override --output {params.prefix} --cpu 8;
        mv {params.prefix}.emapper.annotations static/files/BOTH/
        """

rule eggnog2cnn:
    input:
        eggnog_output="static/files/BOTH/{genome}.emapper.annotations",
        cnn_ami="Training_models/amikacin_egg_model_K_fold_bin.h5",
        cnn_amo="Training_models/amoxicillin_egg_model_K_fold_bin.h5",
        cnn_amp="Training_models/ampicillin_egg_model_K_fold_bin.h5",
        cnn_az="Training_models/aztreonam_egg_model_K_fold_bin.h5",
        cnn_cefe="Training_models/cefepime_egg_model_K_fold_bin.h5",
        cnn_ceft="Training_models/ceftriaxone_egg_model_K_fold_bin.h5",
        cnn_chlo="Training_models/chloramphenicol_egg_model_K_fold_bin.h5",
        cnn_cip="Training_models/ciprofloxacin_egg_model_K_fold_bin.h5",
        cnn_clin="Training_models/clindamycin_egg_model_K_fold_bin.h5",
        cnn_col="Training_models/colistin_egg_model_K_fold_bin.h5",
        cnn_dor="Training_models/doripenem_egg_model_K_fold_bin.h5",
        cnn_ert="Training_models/ertapenem_egg_model_K_fold_bin.h5",
        cnn_eryt="Training_models/erythromycin_egg_model_K_fold_bin.h5",
        cnn_fos="Training_models/fosfomycin_egg_model_K_fold_bin.h5",
        cnn_gen="Training_models/gentamicin_egg_model_K_fold_bin.h5",
        cnn_imi="Training_models/imipenem_egg_model_K_fold_bin.h5",
        cnn_lev="Training_models/levofloxacin_egg_model_K_fold_bin.h5",
        cnn_mer="Training_models/meropenem_egg_model_K_fold_bin.h5",
        cnn_mox="Training_models/moxifloxacin_egg_model_K_fold_bin.h5",
        cnn_nit="Training_models/nitrofurantoin_egg_model_K_fold_bin.h5",
        cnn_tet="Training_models/tetracycline_egg_model_K_fold_bin.h5",
        cnn_tig="Training_models/tigecycline_egg_model_K_fold_bin.h5",
        cnn_tob="Training_models/tobramycin_egg_model_K_fold_bin.h5",
    output:
        ami="static/files/BOTH/CNN_predictions_ami_{genome}.txt",
        amo="static/files/BOTH/CNN_predictions_amo_{genome}.txt",
        amp="static/files/BOTH/CNN_predictions_amp_{genome}.txt",
        az="static/files/BOTH/CNN_predictions_az_{genome}.txt",
        cefe="static/files/BOTH/CNN_predictions_cefe_{genome}.txt",
        ceft="static/files/BOTH/CNN_predictions_ceft_{genome}.txt",
        chlo="static/files/BOTH/CNN_predictions_chlo_{genome}.txt",
        cip="static/files/BOTH/CNN_predictions_cip_{genome}.txt",
        clin="static/files/BOTH/CNN_predictions_clin_{genome}.txt",
        col="static/files/BOTH/CNN_predictions_col_{genome}.txt",
        dor="static/files/BOTH/CNN_predictions_dor_{genome}.txt",
        ert="static/files/BOTH/CNN_predictions_ert_{genome}.txt",
        eryt="static/files/BOTH/CNN_predictions_eryt_{genome}.txt",
        fos="static/files/BOTH/CNN_predictions_fos_{genome}.txt",
        gen="static/files/BOTH/CNN_predictions_gen_{genome}.txt",
        imi="static/files/BOTH/CNN_predictions_imi_{genome}.txt",
        lev="static/files/BOTH/CNN_predictions_lev_{genome}.txt",
        mer="static/files/BOTH/CNN_predictions_mer_{genome}.txt",
        mox="static/files/BOTH/CNN_predictions_mox_{genome}.txt",
        nit="static/files/BOTH/CNN_predictions_nit_{genome}.txt",
        tet="static/files/BOTH/CNN_predictions_tet_{genome}.txt",
        tig="static/files/BOTH/CNN_predictions_tig_{genome}.txt",
        tob="static/files/BOTH/CNN_predictions_tob_{genome}.txt",
    singularity: "./lucyd_machinelearning.sif"
    shell:
        """
        python3 eggnog2cnnPredictions.py -file {input.eggnog_output} -model {input.cnn_ami} -o {output.ami};
        python3 eggnog2cnnPredictions.py -file {input.eggnog_output} -model {input.cnn_amo} -o {output.amo};
        python3 eggnog2cnnPredictions.py -file {input.eggnog_output} -model {input.cnn_amp} -o {output.amp};
        python3 eggnog2cnnPredictions.py -file {input.eggnog_output} -model {input.cnn_az} -o {output.az};
        python3 eggnog2cnnPredictions.py -file {input.eggnog_output} -model {input.cnn_cefe} -o {output.cefe};
        python3 eggnog2cnnPredictions.py -file {input.eggnog_output} -model {input.cnn_ceft} -o {output.ceft};
        python3 eggnog2cnnPredictions.py -file {input.eggnog_output} -model {input.cnn_chlo} -o {output.chlo};
        python3 eggnog2cnnPredictions.py -file {input.eggnog_output} -model {input.cnn_cip} -o {output.cip};
        python3 eggnog2cnnPredictions.py -file {input.eggnog_output} -model {input.cnn_clin} -o {output.clin};
        python3 eggnog2cnnPredictions.py -file {input.eggnog_output} -model {input.cnn_col} -o {output.col};
        python3 eggnog2cnnPredictions.py -file {input.eggnog_output} -model {input.cnn_dor} -o {output.dor};
        python3 eggnog2cnnPredictions.py -file {input.eggnog_output} -model {input.cnn_ert} -o {output.ert};
        python3 eggnog2cnnPredictions.py -file {input.eggnog_output} -model {input.cnn_eryt} -o {output.eryt};
        python3 eggnog2cnnPredictions.py -file {input.eggnog_output} -model {input.cnn_fos} -o {output.fos};
        python3 eggnog2cnnPredictions.py -file {input.eggnog_output} -model {input.cnn_gen} -o {output.gen};
        python3 eggnog2cnnPredictions.py -file {input.eggnog_output} -model {input.cnn_imi} -o {output.imi};
        python3 eggnog2cnnPredictions.py -file {input.eggnog_output} -model {input.cnn_lev} -o {output.lev};
        python3 eggnog2cnnPredictions.py -file {input.eggnog_output} -model {input.cnn_mer} -o {output.mer};
        python3 eggnog2cnnPredictions.py -file {input.eggnog_output} -model {input.cnn_mox} -o {output.mox};
        python3 eggnog2cnnPredictions.py -file {input.eggnog_output} -model {input.cnn_nit} -o {output.nit};
        python3 eggnog2cnnPredictions.py -file {input.eggnog_output} -model {input.cnn_tet} -o {output.tet};
        python3 eggnog2cnnPredictions.py -file {input.eggnog_output} -model {input.cnn_tig} -o {output.tig};
        python3 eggnog2cnnPredictions.py -file {input.eggnog_output} -model {input.cnn_tob} -o {output.tob}
        """

rule combine_files_CNN:
    input:
        predictions_ami="static/files/BOTH/CNN_predictions_ami_{genome}.txt",
        predictions_amo="static/files/BOTH/CNN_predictions_amo_{genome}.txt",
        predictions_amp="static/files/BOTH/CNN_predictions_amp_{genome}.txt",
        predictions_az="static/files/BOTH/CNN_predictions_az_{genome}.txt",
        predictions_cefe="static/files/BOTH/CNN_predictions_cefe_{genome}.txt",
        predictions_ceft="static/files/BOTH/CNN_predictions_ceft_{genome}.txt",
        predictions_chlo="static/files/BOTH/CNN_predictions_chlo_{genome}.txt",
        predictions_cip="static/files/BOTH/CNN_predictions_cip_{genome}.txt",
        predictions_clin="static/files/BOTH/CNN_predictions_clin_{genome}.txt",
        predictions_col="static/files/BOTH/CNN_predictions_col_{genome}.txt",
        predictions_dor="static/files/BOTH/CNN_predictions_dor_{genome}.txt",
        predictions_ert="static/files/BOTH/CNN_predictions_ert_{genome}.txt",
        predictions_eryt="static/files/BOTH/CNN_predictions_eryt_{genome}.txt",
        predictions_fos="static/files/BOTH/CNN_predictions_fos_{genome}.txt",
        predictions_gen="static/files/BOTH/CNN_predictions_gen_{genome}.txt",
        predictions_imi="static/files/BOTH/CNN_predictions_imi_{genome}.txt",
        predictions_lev="static/files/BOTH/CNN_predictions_lev_{genome}.txt",
        predictions_mer="static/files/BOTH/CNN_predictions_mer_{genome}.txt",
        predictions_mox="static/files/BOTH/CNN_predictions_mox_{genome}.txt",
        predictions_nit="static/files/BOTH/CNN_predictions_nit_{genome}.txt",
        predictions_tet="static/files/BOTH/CNN_predictions_tet_{genome}.txt",
        predictions_tig="static/files/BOTH/CNN_predictions_tig_{genome}.txt",
        predictions_tob="static/files/BOTH/CNN_predictions_tob_{genome}.txt",
    output:
        final="static/files/BOTH/{genome}_FINAL.txt_CNN_predicitions.txt"
    run:
        shell("""echo "Antibiotic predicted_pheno" > {output.final}; awk '{{print "Amikacin", $0}}' {input.predictions_ami} >> {output.final}; awk '{{print "Amoxicillin", $0}}' {input.predictions_amo} >> {output.final}; awk '{{print "Ampicillin", $0}}' {input.predictions_amp} >> {output.final}; awk '{{print "Aztreonam", $0}}' {input.predictions_az} >> {output.final}; awk '{{print "Cefepime", $0}}' {input.predictions_cefe} >> {output.final}; awk '{{print "Ceftriaxone", $0}}' {input.predictions_ceft} >> {output.final}; awk '{{print "Chloramphenicol", $0}}' {input.predictions_chlo} >> {output.final}; awk '{{print "Ciprofloxacin", $0}}' {input.predictions_cip} >> {output.final}; awk '{{print "Clindamycin", $0}}' {input.predictions_clin} >> {output.final}; awk '{{print "Colistin", $0}}' {input.predictions_col} >> {output.final}; awk '{{print "Doripenem", $0}}' {input.predictions_dor} >> {output.final}; awk '{{print "Ertapenem", $0}}' {input.predictions_ert} >> {output.final}; awk '{{print "Erythromcyin", $0}}' {input.predictions_eryt} >> {output.final}; awk '{{print "Fosfomycin", $0}}' {input.predictions_fos} >> {output.final}; awk '{{print "Gentamicin", $0}}' {input.predictions_gen} >> {output.final}; awk '{{print "Imipenem", $0}}' {input.predictions_imi} >> {output.final}; awk '{{print "Levofloxacin", $0}}' {input.predictions_lev} >> {output.final}; awk '{{print "Meropenem", $0}}' {input.predictions_mer} >> {output.final}; awk '{{print "Moxifloxacin", $0}}' {input.predictions_mox} >> {output.final}; awk '{{print "Nitrofurantoin", $0}}' {input.predictions_nit} >> {output.final}; awk '{{print "Tetracycline", $0}}' {input.predictions_tet} >> {output.final}; awk '{{print "Tigecycline", $0}}' {input.predictions_tig} >> {output.final}; awk '{{print "Tobramycin", $0}}' {input.predictions_tob} >> {output.final}""")


rule eggnog_2_arff:
    input:
        "static/files/BOTH/{genome}.emapper.annotations"
    output:
        "static/files/BOTH/{genome}.arff"
    singularity: "./ldillon_amrwebsite.sif"
    shell:
        """
        python eggnog_output_to_arff.py -file {input} -o {output} -c_attr_name phenotype -catagory1 Susceptible -catagory2 Resistant
        """

rule genome_names:
    input:
        "static/files/{genome}.fna"
    output:
        temp("static/files/BOTH/{genome}_names.txt")
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
        genome_names="static/files/BOTH/{genome}_names.txt",
        arff_file="static/files/BOTH/{genome}.arff",
    output:
        paths=temp("static/files/BOTH/{genome}_names.txt.paths.txt"),
        predictions=temp("static/files/BOTH/{genome}_names.txt.predicitions.txt"),
        predictions_ami=temp("static/files/BOTH/{genome}_ami.txt.DT_predicitions.txt"),
        predictions_amo=temp("static/files/BOTH/{genome}_amo.txt.DT_predicitions.txt"),
        predictions_amp=temp("static/files/BOTH/{genome}_amp.txt.DT_predicitions.txt"),
        predictions_az=temp("static/files/BOTH/{genome}_az.txt.DT_predicitions.txt"),
        predictions_cip=temp("static/files/BOTH/{genome}_cip.txt.DT_predicitions.txt"),
        predictions_clin=temp("static/files/BOTH/{genome}_clin.txt.DT_predicitions.txt"),
        predictions_col=temp("static/files/BOTH/{genome}_col.txt.DT_predicitions.txt"),
        predictions_dor=temp("static/files/BOTH/{genome}_dor.txt.DT_predicitions.txt"),
        predictions_ert=temp("static/files/BOTH/{genome}_ert.txt.DT_predicitions.txt"),
        predictions_gen=temp("static/files/BOTH/{genome}_gen.txt.DT_predicitions.txt"),
        predictions_mer=temp("static/files/BOTH/{genome}_mer.txt.DT_predicitions.txt"),
        predictions_nit=temp("static/files/BOTH/{genome}_nit.txt.DT_predicitions.txt"),
        predictions_tet=temp("static/files/BOTH/{genome}_tet.txt.DT_predicitions.txt"),
        predictions_tig=temp("static/files/BOTH/{genome}_tig.txt.DT_predicitions.txt"),
        predictions_tob=temp("static/files/BOTH/{genome}_tob.txt.DT_predicitions.txt"),
    singularity: "./ldillon_amrwebsite.sif"
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


rule combine_files_DT:
    input:
        predictions_ami="static/files/BOTH/{genome}_ami.txt.DT_predicitions.txt",
        predictions_amo="static/files/BOTH/{genome}_amo.txt.DT_predicitions.txt",
        predictions_amp="static/files/BOTH/{genome}_amp.txt.DT_predicitions.txt",
        predictions_az="static/files/BOTH/{genome}_az.txt.DT_predicitions.txt",
        predictions_cip="static/files/BOTH/{genome}_cip.txt.DT_predicitions.txt",
        predictions_clin="static/files/BOTH/{genome}_clin.txt.DT_predicitions.txt",
        predictions_col="static/files/BOTH/{genome}_col.txt.DT_predicitions.txt",
        predictions_dor="static/files/BOTH/{genome}_dor.txt.DT_predicitions.txt",
        predictions_ert="static/files/BOTH/{genome}_ert.txt.DT_predicitions.txt",
        predictions_gen="static/files/BOTH/{genome}_gen.txt.DT_predicitions.txt",
        predictions_mer="static/files/BOTH/{genome}_mer.txt.DT_predicitions.txt",
        predictions_nit="static/files/BOTH/{genome}_nit.txt.DT_predicitions.txt",
        predictions_tet="static/files/BOTH/{genome}_tet.txt.DT_predicitions.txt",
        predictions_tig="static/files/BOTH/{genome}_tig.txt.DT_predicitions.txt",
        predictions_tob="static/files/BOTH/{genome}_tob.txt.DT_predicitions.txt"
    output:
        final="static/files/BOTH/{genome}_FINAL.txt.DT_predicitions.txt",
    run:
        shell("""echo "Antibiotic   Genome_name Genome_number   Node_number Node_label" > {output.final}; awk '{{print "Amikacin", $0}}' {input.predictions_ami} >> {output.final}; awk '{{print "Amoxicillin", $0}}' {input.predictions_amo} >> {output.final}; awk '{{print "Ampicillin", $0}}' {input.predictions_amp} >> {output.final}; awk '{{print "Aztreonam", $0}}' {input.predictions_az} >> {output.final}; awk '{{print "Ciprofloxacin", $0}}' {input.predictions_cip} >> {output.final}; awk '{{print "Clindamycin", $0}}' {input.predictions_clin} >> {output.final}; awk '{{print "Colistin", $0}}' {input.predictions_col} >> {output.final}; awk '{{print "Doripenem", $0}}' {input.predictions_dor} >> {output.final}; awk '{{print "Ertapenem", $0}}' {input.predictions_ert} >> {output.final}; awk '{{print "Gentamicin", $0}}' {input.predictions_gen} >> {output.final}; awk '{{print "Meropenem", $0}}' {input.predictions_mer} >> {output.final}; awk '{{print "Nitrofurantoin", $0}}' {input.predictions_nit} >> {output.final}; awk '{{print "Tetracycline", $0}}' {input.predictions_tet} >> {output.final}; awk '{{print "Tigecycline", $0}}' {input.predictions_tig} >> {output.final}; awk '{{print "Tobramycin", $0}}' {input.predictions_tob} >> {output.final}""")

rule zip_files:
    input:
        "static/files/BOTH/{genome}.faa",
        "static/files/BOTH/{genome}_Diamond.faa",
        "static/files/BOTH/{genome}_DiamondReduced.faa",
        "static/files/BOTH/{genome}.emapper.annotations",
        "static/files/BOTH/{genome}.arff",
        "static/files/BOTH/{genome}_FINAL.txt_CNN_predicitions.txt",
        "static/files/BOTH/{genome}_FINAL.txt.DT_predicitions.txt",
    params:
        path="static/files/BOTH"
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
