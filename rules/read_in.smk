rule read_in:
	input:
		expand("data/{sample}", sample=list(samples.index))
	output:
		"file_dir/data.h5ad"
	conda:
		"../env/oneforall.yaml"
	params:
		units = config["data"]["units"]
	script:
		"../scripts/read_in.py"
