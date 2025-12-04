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


## Edit paths:
You will need to add the respective paths and names of the singularity containers in the 'process.py' file:
```
container_sources = {
        'container': tool_path / 'amrmlpipeline_DT_arm64.sif', 
        'ml_container': tool_path / 'amrmlpipeline_CNN.sif'
    }
    
    # Also check in parent directory
    if not container_sources['container'].exists():
        container_sources['container'] = Path('/__YOUR__PATH___/amrmlpipeline_DT_arm64.sif')
    if not container_sources['ml_container'].exists():
        container_sources['ml_container'] = Path('/__YOUR__PATH___/amrmlpipeline_CNN.sif')
```

## Downloads and singularity 

Singularity definition files are in the repo for users to build (note: if you are using a Mac, you will need a VM for Singularity).

Due to large file sizes of the Sourmash DBs, please download the files from the following link on the Open Science Framework:
https://osf.io/qv3fs/overview?view_only=181294b8e47e4fd780e3b292158c8b43 \
Put these sourmash files into the same folder as the model files "data".

Download eggNOG DBs from: http://eggnog5.embl.de/download/emapperdb-5.0.2/ \
Put these files inside the AMRcast directory.

## Extra
I would advise checking the sourmash output file to ensure the genome is sufficiently similar to the training data to make an accurate prediction. \
Size of repo ~240MB (not including the extra downloads) \
If you have any ideas to improve, please report a suggestion.



