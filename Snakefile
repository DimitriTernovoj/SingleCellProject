import os
import pandas as pd

configfile: "config.yaml"

samples = pd.read_table(config["data"]["samples"], index_col="sample")
units = pd.read_table(config["data"]["units"])

rule all:
        input:
                "figures/slingshot_data_vis.jpg","figures/monocle_trajectory.jpg", dynamic("figures/{cluster}_enrichment_vis.pdf") if len(units.columns) != 3 else [] 


include: "rules/preprocessing.smk"
include: "rules/read_in.smk"
include: "rules/quality_control.smk"
include: "rules/normalization.smk"
include: "rules/visualization.smk"
include: "rules/clustering.smk"
include: "rules/slingshot_R.smk"
include: "rules/monocle2.smk"
include: "rules/differential_testing.smk"
include: "rules/gene_set_analysis.smk"
