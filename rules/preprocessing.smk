rule bamtofastq:
        input:
                lambda wildcards: samples.at[wildcards.sample,"path"] if wildcards.sample in samples.index else ""
        output:
                directory("data/{sample}_results")
        params:
                path_to_exe = config["preprocessing"]["path_to_bamtofastq"],
        threads:
                workflow.cores * 0.75
        shell:
                """
                chmod 700 {params.path_to_exe}
                {params.path_to_exe} --cr11 --nthreads={threads} {input} data/{wildcards.sample}_results
                """

rule cellranger_count:
        input:
                "data/{sample}_results"
        output:
                directory("data/{sample}")
        params:
                path_to_fastq = "data/{sample}_results/fastqs/",
                path_to_cellranger = config["preprocessing"]["path_to_cellranger"],
                ref = config["preprocessing"]["path_to_ref"],
		cores = workflow.cores
        shell:
                """
                export PATH={params.path_to_cellranger}:$PATH
                mv data/{wildcards.sample}_results/* data/{wildcards.sample}_results/fastqs
                cellranger count --id={wildcards.sample} --fastqs={params.path_to_fastq} --transcriptome={params.ref} --localcores={params.cores}
		mv {wildcards.sample} data
                """

