import pandas as pd
from gprofiler import GProfiler
import gprofiler

file = pd.read_csv("results/"+ snakemake.wildcards.cluster +"_diff_testing.csv")

gp = GProfiler(return_dataframe=True, user_agent='g:GOSt')

enrichment = gp.profile(organism= snakemake.params.organism, sources=['GO:BP'], user_threshold=float(snakemake.params.EN_threshold),
                              significance_threshold_method='fdr',
                              query=file['primerid'][file['FDR']<float(snakemake.params.DE_threshold)].tolist())

enrich_results = enrichment.set_index('native').sort_values('p_value').iloc[:,[2,5,7,10,1]]

enrich_results.to_csv("results/"+snakemake.wildcards.cluster+"_enrich_results.csv", sep=",")

