rule visualization:
	input:
		lambda wildcards: expand("file_dir/normed_data_" + wildcards.downsampling + ".h5ad", downsampling=config["downsampling"]["downsampling_method"]) if wildcards.downsampling != "" else ""
	output:
		"file_dir/post_vis_data_{downsampling}.h5ad",
		"file_dir/specter_pca_data_{downsampling}.csv" if config["clustering"]["cluster_method"] == "specter" else [], 
		report("figures/pca_datavis_{downsampling}.pdf", caption="../report/pca_datavis.rst", category = "Visualisierung dimensionsreduzierter Daten"),
                report("figures/tsne_datavis_{downsampling}.pdf", caption="../report/tsne_datavis.rst", category = "Visualisierung dimensionsreduzierter Daten"),
                report("figures/umap_datavis_{downsampling}.pdf", caption="../report/umap_datavis.rst", category = "Visualisierung dimensionsreduzierter Daten"),
                report("figures/diffmap_datavis_{downsampling}.pdf", caption="../report/diffmap_datavis.rst", category = "Visualisierung dimensionsreduzierter Daten"),
                report("figures/draw_graph_fr_datavis_{downsampling}.pdf", caption="../report/draw_graph_fr_datavis.rst", category = "Visualisierung dimensionsreduzierter Daten"), #rausnehmen enventuell
                report("figures/umap_cellcycle_1_{downsampling}.pdf", caption="../report/umap_cellcycle_1.rst", category = "Visualisierung dimensionsreduzierter Daten", subcategory= "Zellzyklus Effekte") if config["cell_cycle_scoring"]["ref_genes"] != "" else [],
                report("figures/umap_cellcycle_2_{downsampling}.pdf", caption="../report/umap_cellcycle_2.rst", category = "Visualisierung dimensionsreduzierter Daten", subcategory= "Zellzyklus Effekte") if config["cell_cycle_scoring"]["ref_genes"] != "" else []
	conda:
		"../env/oneforall.yaml"
	params:
		ref_genes = config["cell_cycle_scoring"]["ref_genes"],
		clustering = config["clustering"]["cluster_method"]
	script:
		"../scripts/visualization.py"

