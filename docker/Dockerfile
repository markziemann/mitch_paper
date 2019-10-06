# docker run -d -p 8787:8787 -e USER=rstudio -e PASSWORD=mitch -e ROOT=TRUE rocker/rstudio
FROM rocker/rstudio

RUN apt update && \
  apt install -y \
    zlib1g-dev \
    libxml2-dev \
    tcl8.6-dev \
    tk8.6-dev
    
RUN cat << EOF > installation.R
# cran install
install.packages("devtools")
install.packages("plyr")
install.packages("tidyverse")
install.packages("reshape2")
install.packages("ashr")
install.packages("parallel")
install.packages("grid")
install.packages("gridExtra")
install.packages("beeswarm")
install.packages("GGally")
install.packages("knitr")
install.packages("markdown")
install.packages("gplots")
install.packages("purrr")
install.packages("Seurat")
install.packages("Rmisc")
install.packages("UpSetR")
install.packages("vioplot")
install.packages("rmarkdown")
install.packages("xfun")

# github install
devtools::install_github("markziemann/dee2/getDEE2")
devtools::install_github("hrbrmstr/taucharts")
devtools::install_github("markziemann/Mitch")

# bioc install
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# bioc install
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(version='devel',ask = FALSE, update=TRUE)

BiocManager::install("ABSSeq", ask = FALSE, update=TRUE)
BiocManager::install("edgeR", ask = FALSE, update=TRUE)
BiocManager::install("DESeq2", ask = FALSE, update=TRUE)
BiocManager::install("fgsea", ask = FALSE, update=TRUE)
BiocManager::install("topconfects", ask = FALSE, update=TRUE)
BiocManager::install("biomaRt", ask = FALSE, update=TRUE)
BiocManager::install("ExperimentHub", ask = FALSE, update=TRUE)
BiocManager::install("muscat", ask = FALSE, update=TRUE)
BiocManager::install("scater",ask = FALSE , update=TRUE)
EOF

RUN Rscript installation.R
