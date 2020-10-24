rule differential_testing:
	input:
		"file_dir/adata_dt.h5ad"
	output:
		dynamic("results/{cluster}_diff_testing.csv")
	conda:
		"../env/differential_testing.yaml"
	script:
		"../scripts/differential_testing.R"
