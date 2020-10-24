rule normalization_1:
    input:
        "file_dir/filtered_data.h5ad"
    output:
        "file_dir/data_mat.csv",
        "file_dir/input_groups.csv"
    conda:
        "../env/oneforall.yaml"
    script:
        "../scripts/normalization.py"
        #"scripts/rpy2_tests.py"

rule normalization_R:
    input:
        "file_dir/data_mat.csv",
        "file_dir/input_groups.csv"
    output:
        "file_dir/size_factors.csv"
    conda:
        "../env/r_norm.yaml"
    script:
        "../scripts/scran.R"

rule normalization_2:
    input:
        "file_dir/filtered_data.h5ad",
        "file_dir/size_factors.csv"
    output:
        "file_dir/normed_data_normal.h5ad",
        report("figures/scatter_normed_ncounts_sizefactors.pdf", caption="../report/scatter_normed_ncounts_sizefactors.rst", category="Normalisierung"),
        report("figures/scatter_normed_ngenes_sizefactors.pdf", caption="../report/scatter_normed_ngenes_sizefactors.rst", category="Normalisierung"),
        report("figures/sizefactors_plot.png", caption="../report/sizefactors_plot.rst", category="Normalisierung"),
        report("figures/filter_genes_dispersion_high_variable_genes.pdf", caption="../report/gene_dispersion.rst", category="Visualisierung dimensionsreduzierter Daten", subcategory="stark variable Gene")
    conda:
        "../env/oneforall.yaml"
    params:
        downsampling = config["downsampling"]["downsampling_method"],
        amount_of_hvgs = config["general"]["amount_of_hvgs"]
    script:
        "../scripts/normalization_2.py"


####Sphetcher
rule downsampling_sphetcher_1:
	input:
		"file_dir/normed_data_normal.h5ad"
	output:
		"file_dir/subset.csv"
	params:
		path = config["downsampling"]["path_to_sphetcher"],
		amount_of_cells = config["downsampling"]["sketch_size"]
	shell:
		"{params.path}/sphetcher file_dir/sphetcher_pca.csv {params.amount_of_cells} file_dir/subset.csv"

rule downsampling_sphetcher_2:
	input:
		"file_dir/normed_data_normal.h5ad",
		"file_dir/subset.csv"
	output:
		"file_dir/normed_data_sphetcher.h5ad"
	conda:
		"../env/oneforall.yaml"
	script:
		"../scripts/downsampling.py"
