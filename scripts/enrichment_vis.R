library("dplyr")
library("ggplot2")

print(snakemake@wildcards[["cluster"]])

name <- snakemake@wildcards[["cluster"]]
data <- read.csv(file = snakemake@input[[1]])

file_name <- paste("figures/",name,"_enrichment_vis.pdf", sep="")
pdf(file_name)

data %>%
  top_n(10,wt=-p_value) %>%
  ggplot(aes(x=recall,
             y=name,
             colour=p_value,
             size=intersection_size)) +
         geom_point() +
         expand_limits(x=0) +
         labs(x="Gene Ratio", y="GO term", size="Count", colour="p value")
dev.off()
