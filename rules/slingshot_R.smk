rule slingshot:
	input:
		"file_dir/adata_ent.h5ad"
	output:
		report("figures/slingshot_data_vis.jpg", caption = "../report/slingshot_data_vis.rst", category = "Trajektorienbildung", subcategory = "Slingshot"),
                report("figures/slingshot_pst_1.jpg", caption = "../report/slingshot_pst_1.rst", category = "Trajektorienbildung", subcategory = "Slingshot"),
                report("figures/slingshot_pst_2.jpg", caption = "../report/slingshot_pst_2.rst", category = "Trajektorienbildung", subcategory = "Slingshot"),
                report("figures/slingshot_data_vis_nbc.jpg", caption = "../report/slingshot_data_vis_nbc.rst", category = "Trajektorienbildung", subcategory = "Slingshot"),
                report("figures/slingshot_pst_1_nbc.jpg", caption = "../report/slingshot_pst_1_nbc.rst", category = "Trajektorienbildung", subcategory = "Slingshot"),
                report("figures/slingshot_pst_2_nbc.jpg", caption = "../report/slingshot_pst_2_nbc.rst", category = "Trajektorienbildung", subcategory = "Slingshot"),
		report("figures/slingshot_annotated.jpg", caption = "../report/annotated_trajectory.rst", category = "Trajektorienbildung", subcategory = "Slingshot"),
		report("figures/slingshot_annotated_nbc.jpg", caption = "../report/annotated_trajectory_nbc.rst", category = "Trajektorienbildung", subcategory = "Slingshot"),
		report("figures/Heatmap.pdf", caption = "../report/heatmap.rst", category = "Trajektorienbildung", subcategory = "Slingshot")
	conda:
		"../env/slingshot_R.yaml"
	script:
		"../scripts/slingshot.R"
