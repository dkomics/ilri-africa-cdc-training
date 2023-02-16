# Fileformat is
# Gene Chr Pos Cov
# gene_A A01 40 0
# gene_A A01 41 2
# ...
# gene_D A01 508 41
# gene_D A01 509 42
# ...

geneCov <- read.table("./test_genes.txt", sep="\t", header=F)
colnames(geneCov) <- c("gene", "chr", "start", "cov")

library(ggplot2)
density<-ggplot(data=geneCov, aes(x=start, y=cov, group=gene)) + 
	geom_line(aes(linetype=gene)) + # linetype=gene is a bit unnecessary - 
	# styles don't distinguish anything here as every gene has own facet
	facet_wrap(~ gene, scales="free_x", ncol=1) +
	xlab("Position in the genome") + 
	ylab("Coverage density") + 
	theme_bw() 
png("test_cov_gene_density.png",500,750)
print(density)
dev.off()
