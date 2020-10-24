rule gprofiler:
	input:
		"results/{cluster}_diff_testing.csv"
	output:
		"results/{cluster}_enrich_results.csv"
	conda:
		"../env/gene_set_analysis.yaml"
	params:
		EN_threshold = config["gene_set_enrichment_analysis"]["enrichment_threshold"],
		DE_threshold = config["differential_testing"]["DE_threshold"],
		organism = config["gene_set_enrichment_analysis"]["organism"]
	script:
		"../scripts/gene_set_analysis.py"

rule enrichment_vis:
	input:
		"results/{cluster}_enrich_results.csv"
	output:
		report("figures/{cluster}_enrichment_vis.pdf", caption= "../report/enrichment_vis.rst", category="Gene Set Enrichment Analyse")
	conda:
		"../env/enrichment_vis.yaml"
	script:
		"../scripts/enrichment_vis.R"
