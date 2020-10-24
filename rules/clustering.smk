rule clustering_louvain:
	input:
		expand("file_dir/post_vis_data_{downsampling}.h5ad", downsampling=config["downsampling"]["downsampling_method"])
	output:
		"file_dir/post_clustering_louvain.h5ad",
                "file_dir/annotation_list_louvain.npy",
		report(expand("figures/umap_{gene}.pdf", gene = config["clustering"]["genes_to_vis"].split(",")), caption = "../report/gene_umap.rst", category = "Clustering"),
                report(expand("figures/violin_{gene}.pdf", gene = config["clustering"]["genes_to_vis"].split(",")), caption = "../report/gene_violin.rst", category = "Clustering"),
                report("figures/umap_clustering_result.pdf", caption= "../report/umap_clustering_result.rst", category="Clustering"),
                report("figures/umap_region_ncounts_oncluster.pdf", caption= "../report/umap_region_ncounts_oncluster.rst", category="Clustering" ),
                report("figures/umap_logcounts_mtfrac_oncluster.pdf",caption= "../report/umap_logcounts_mtfrac_oncluster.rst", category="Clustering"),
                report(f"figures/rank_genes_groups_louvain_r{config['clustering']['clustering_resolution']}_markergenes.pdf", caption= "../report/gene_ranking.rst", category="Clustering"), #Das ist problematisch, da 0.5 vom Nutzer besti$
                report("figures/markergenes_overlap.png", caption= "../report/markergenes_overlap.rst", category="Clustering"),
                report("figures/umap_annotated_clustering.pdf", caption= "../report/umap_annotated_clustering.rst", category="Clustering"),
	conda:
		"../env/oneforall.yaml"
	params:
		celltype_annotation = config["clustering"]["celltypes_markergenes"],
		genes = config["clustering"]["genes_to_vis"],
		subclustering = config["clustering"]["subclustering"],
		clustering_resolution = config["clustering"]["clustering_resolution"],
		#subclustering_resolution = config["clustering"]["subclustering_resolution"]
	script:
		"../scripts/clustering.py"

rule specter:
        input:
               expand("file_dir/specter_pca_data_{downsampling}.csv", downsampling=config["downsampling"]["downsampling_method"])
        output:
               "file_dir/specter_clustering.csv"
        params:
               number_of_clusters = config["specter"]["number_of_clusters"],
               ensemble_size = config["specter"]["ensemble_size"],
               mingamma = config["specter"]["mingamma"],
               downsampling = "0" if config["downsampling"]["downsampling_method"] =="normal" else "1"
        shell:
               """
               cd scripts/
               matlab -r 'n_clusters={params.number_of_clusters};ensemble_size={params.ensemble_size};mingamma={params.mingamma};downsampling=string({params.downsampling});' < Specter.m 
               """
	
rule clustering_specter:
        input:
                expand("file_dir/post_vis_data_{downsampling}.h5ad", downsampling=config["downsampling"]["downsampling_method"])
        output:
                "file_dir/post_clustering_specter.h5ad",
                "file_dir/annotation_list_specter.npy",
		report(expand("figures/umap_{gene}.pdf", gene = config["clustering"]["genes_to_vis"].split(",")), caption = "../report/gene_umap.rst", category = "Clustering"),
                report(expand("figures/violin_{gene}.pdf", gene = config["clustering"]["genes_to_vis"].split(",")), caption = "../report/gene_violin.rst", category = "Clustering"),
                report("figures/umap_clustering_result.pdf", caption= "../report/umap_clustering_result.rst", category="Clustering"),
                report("figures/umap_region_ncounts_oncluster.pdf", caption= "../report/umap_region_ncounts_oncluster.rst", category="Clustering" ),
                report("figures/umap_logcounts_mtfrac_oncluster.pdf",caption= "../report/umap_logcounts_mtfrac_oncluster.rst", category="Clustering"),
                report("figures/rank_genes_groups_specter_markergenes.pdf", caption= "../report/gene_ranking.rst", category="Clustering"), #Das ist problematisch, da 0.5 vom Nutzer besti$
                report("figures/markergenes_overlap.png", caption= "../report/markergenes_overlap.rst", category="Clustering"),
                report("figures/umap_annotated_clustering.pdf", caption= "../report/umap_annotated_clustering.rst", category="Clustering"),
        conda:
                "../env/oneforall.yaml"
        params:
                celltype_annotation = config["clustering"]["celltypes_markergenes"],
                genes = config["clustering"]["genes_to_vis"],
                subclustering = config["clustering"]["subclustering"],
                clustering_resolution = config["clustering"]["clustering_resolution"],
                #subclustering_resolution = config["clustering"]["subclustering_resolution"]
        script:
                "../scripts/clustering_specter.py"

rule subclustering:
	input:
		expand(["file_dir/post_clustering_{cluster}.h5ad","file_dir/annotation_list_{cluster}.npy"], cluster=config["clustering"]["cluster_method"])
	output:
		"file_dir/post_clustering_final.h5ad",
		"file_dir/adata_dt.h5ad",
		"file_dir/adata_ent.h5ad",
		report("figures/umap_named_annotated_clustering.pdf", caption = "../report/named_unannotated.rst", category = "Clustering", subcategory = "Subclustering") if config["subclustering"]["names_for_unannotated"] != "" else [],
		report("figures/umap_final_clustering.pdf", caption = "../report/umap_final_clustering.rst", category = "Clustering", subcategory = "Subclustering"),
		report("figures/umap_regions_on_finalclustering.pdf", caption = "../report/umap_regions_on_finalclustering.rst", category = "Clustering", subcategory = "Subclustering"),
		report("figures/pca_subset_for_ti.pdf", caption = "../report/pca_subset_for_ti.rst", category = "Vorbereitung für Trajektorienbildung"),
		report("figures/pca_variance_ratio_for_ti.pdf", caption = "../report/pca_variance_ratio_for_ti.rst", category = "Vorbereitung für Trajektorienbildung"),
		report("figures/pca_subset_for_ti_nbc.pdf", caption = "../report/pca_subset_for_ti_nbc.rst", category = "Vorbereitung für Trajektorienbildung"),
		report("figures/pca_variance_ratio_for_ti_nbc.pdf", caption = "../report/pca_variance_ratio_for_ti_nbc.rst", category = "Vorbereitung für Trajektorienbildung"),
		report("figures/diffmap_diffusion_pseudotime_1.pdf", caption = "../report/diffmap_diffusion_pseudotime_1.rst", category = "Vorbereitung für Trajektorienbildung"),
		report("figures/diffmap_diffusion_pseudotime_2.pdf", caption = "../report/diffmap_diffusion_pseudotime_2.rst", category = "Vorbereitung für Trajektorienbildung")
	conda:
		"../env/oneforall_subclustering.yaml"
	params:
		clustering_resolution = config["clustering"]["clustering_resolution"],
		subclustering_resolution = config["subclustering"]["subclustering_resolution"],
		names_for_unannotated = config["subclustering"]["names_for_unannotated"],
		further_subclusterings = config["subclustering"]["further_subclusterings"],
		clusters_to_include = config["trajectory_inference"]["clusters_to_include"],
		amount_of_hvgs = config["general"]["amount_of_hvgs"],
		normed_type = config["downsampling"]["downsampling_method"],
		cluster_method = config["clustering"]["cluster_method"]	
	script:
		"../scripts/subclustering.py"

