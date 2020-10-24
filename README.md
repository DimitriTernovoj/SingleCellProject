#Snakemake Workflow: single cell RNA-seq analysis
The chosen methods are based on the paper of Luecken and Theis, 2019 (https://doi.org/10.15252/msb.20188746), explaining the current best practice in single cell RNA-seq analysis.

##Installations
l. install Snakemake (workflow management system)
l. install conda (package management system)
l. install mamba (package management system) (optional)
l. install bamtofastq and cellranger by 10X Genomics (used for preprocessing of data)
l. install Sphetcher (downsampling algorithm) (optional)

##Usage
l. clone this repository recursively (because of the submodule Specter)
l. configure the workflow by editing the config.yaml-file (parameters described in the file)
l. start the execution in the folder the Snakefile is located by typing one of the following commands:
```
snakemake --use-conda --cores x
snakemake --use-conda --cores x --conda-frontend mamba
```
x is the number of cores used. The second command uses the package management system mamba instead of conda and should be used if the installation of the environments takes too long.

##Data
The workflow starts with bam-files, definied in two tsv-files, which are linked in the config.yaml. The Samples.tsv has two columns, of which the first one defines the sample and the second one the corresponding path to the bam-file. The Units.tsv has three to five columns and defines further information about the samples. The first column specifies the sample and the second one the alias used in the workflow. The third column defines regions for all samples (used in visualizations across the regions). The last two columns are optional and should be only used if differential testing is supposed to be performed. The fourth column is then named contrast and specifies the two groups between which differential testing is performed (definied by the letters A and B). The last column is optional and can be named by the user, it specifies another source of variability, which is accounted for in the differential testing.

example Samples.tsv:
sample | path
-------|------
sample1 | path_to_sample1.bam
sample2 | path_to_sample2.bam

example Units.tsv:
sample | sample_alias | region | contrast | variable_name
sample1 | s1 | reg1 | A | x
sample2 | s2 | reg2 | B | y
