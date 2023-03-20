# Fileformat is
# gene Pos Cov
# gene_A 40 0
# gene_A 41 2
# ...
# gene_D 508 41
# gene_D 509 42
# ...

args = commandArgs(trailingOnly=TRUE)
coverageFile = args[1]

geneCov <- read.table(coverageFile, sep="\t", header=F)
colnames(geneCov) <- c("gene", "start", "cov")

library(ggplot2)
density<-ggplot(data=geneCov, aes(x=start, y=cov, group=gene)) + 
	geom_line(aes(linetype=gene)) + # linetype=gene is a bit unnecessary - 
	# styles don't distinguish anything here as every gene has own facet
	facet_wrap(~ gene, scales="free_x", ncol=1) +
	xlab("Position in the genome") + 
	ylab("Coverage density") + 
	theme_bw()
ggsave(plot = density, "test_cov_gene_density.png",
	width = 50, height = 75, units = "cm", res=300) 
#png("test_cov_gene_density.png",1000,1500,units="px",pointsize=12)
#print(density)
#dev.off()
