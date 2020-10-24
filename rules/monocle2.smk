rule monocle2:
	input:
		"file_dir/post_clustering_final.h5ad"
	output:
		report("figures/monocle_trajectory.jpg", caption ="../report/monocle_trajectory.rst", category = "Trajektorienbildung", subcategory="Monocle2"),
		report("figures/monocle_pseudotime.jpg", caption ="../report/monocle_pseudotime.rst", category = "Trajektorienbildung", subcategory="Monocle2")
	conda:
		"../env/monocle2.yaml"
	script:
		"../scripts/monocle2.R"
