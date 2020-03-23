download.file("https://raw.githubusercontent.com/HelenaLC/muscat/master/vignettes/vignette.Rmd",destfile = "muscat.Rmd")

library("splatter")
library("muscat")
library("mitch")

rmarkdown::render("muscat.Rmd")

download.file("https://reactome.org/download/current/ReactomePathways.gmt.zip",destfile="ReactomePathways.gmt.zip")
unzip("ReactomePathways.gmt.zip")
genesets<-gmt_import("ReactomePathways.gmt")
y<-mitch_import(tbl,DEtype="muscat",geneID="gene",joinType = "full")
correl<-cor(y,use = "pairwise.complete.obs")

pdf("correl.pdf")
heatmap(correl,cexRow=1.5,cexCol = 1.5,scale="none",col = heat.colors(256),margins = c(12,12))
dev.off()

res<-mitch_calc(y,genesets,resrows = 25, priority="effect")
mitch_plots(res,outfile ="scrna_effect.pfd")

res<-mitch_calc(y,genesets,resrows = 25, priority="significance")
mitch_plots(res,outfile ="scrna_sig.pdf")

res<-mitch_calc(y,genesets,resrows = 25, priority="SD")
mitch_plots(res,outfile ="scrna_sd.pdf")



