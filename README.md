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


