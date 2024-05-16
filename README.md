# AMRProPredict
Command line workflow to predict AMR phenotype

Two containers were needed due to memory-building restrictions on sylabs- we are working to integrate this into a single docker or singularity container.

To pull the singularity image used in the decision tree smk:
```
singularity pull --arch amd64 library://ldillon/amrwebsite/amrmlpipeline:6
```
To pull the singularity image used in the CNN smk:
```
singularity pull --arch amd64 library://lucyd/machinelearning/mlpackages:1
```

The eggNOG databases are required to run these workflows (please look at .smk files for the directory names). We are in the process of integrating a 'fast' option to allow users to use the eggNOG output (*.emapper.annotations).


