# AMRcast 
Command line workflow to predict AMR phenotype.

## WARNING! This work is still under development.

## Quick Start
```bash
# Basic usage
python process.py \
    --input your_genome.fasta \
    --output results.tar.gz \
    --model decision_tree \
    --cores 4 \
    --tool-dir /path/to/this/directory

# Or use the wrapper
./run_analysis.sh --input genome.fasta --output results.tar.gz --model both_models
```

## Model Options
- `decision_tree` - more biological detail
- `cnn` - Deep learning (higher accuracy typically)
- `both_models` - Run both for comprehensive analysis


## Downloads and singularity 

Singularity definition files are in the repo for users to build (note: if you are using a Mac, you will need a VM for Singularity).

Due to large file sizes of the Sourmash DBs, please download the files from the following link on the Open Science Framework:
https://osf.io/qv3fs/overview?view_only=181294b8e47e4fd780e3b292158c8b43

Download eggNOG DBs using the instructions on their GitHub: https://github.com/eggnogdb/eggnog-mapper/ 

## Extra
I would advise checking the sourmash output file to ensure the genome is sufficiently similar to the training data to make an accurate prediction.
If you have any ideas to improve, please report a suggestion.



