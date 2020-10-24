rule quality_control:
	input:
		"file_dir/data.h5ad"
	output:
		"file_dir/filtered_data.h5ad",
		report("figures/violin_ncounts_groupedby_sample.pdf", caption= "../report/violin_ncounts_groupedby_sample.rst", category="Qualitätskontrolle"),
		report("figures/violin_mt_frac_groupedby_sample.pdf", caption= "../report/violin_mt_frac_groupedby_sample.rst", category="Qualitätskontrolle"),
		report("figures/n_counts.png", caption= "../report/n_counts.rst", category="Qualitätskontrolle"),
		report("figures/lower_ncounts.png", caption= "../report/lower_ncounts.rst", category="Qualitätskontrolle"),
		report("figures/upper_ncounts.png", caption= "../report/upper_ncounts.rst", category="Qualitätskontrolle"),
		report("figures/ncounts_genes.png", caption= "../report/ncounts_genes.rst", category="Qualitätskontrolle"),
		report("figures/lower_ncounts_genes.png", caption= "../report/lower_ncounts_genes.rst", category="Qualitätskontrolle"),
                report("figures/scatter_ncounts_ngenes.pdf", caption= "../report/scatter_ncounts_ngenes.rst", category="Qualitätskontrolle")
	params:
		lower_quantile_genes = config["quality_control"]["cells"]["lower_quantile_genes"],
		lower_quantile_counts = config["quality_control"]["cells"]["lower_quantile_counts"],
		upper_quantile_counts = config["quality_control"]["cells"]["upper_quantile_counts"],
		min_counts = config["quality_control"]["cells"]["min_counts"],
		max_counts = config["quality_control"]["cells"]["max_counts"],
		mt_frac = config["quality_control"]["cells"]["mt_frac"],
		min_genes = config["quality_control"]["cells"]["min_genes"],
		min_cells = config["quality_control"]["genes"]["min_cells"]
	conda:
		"../env/oneforall.yaml"
	script:
		"../scripts/qc.py"
