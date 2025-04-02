# ICArus
Matrix Decomposition Based pipeline that extracts gene signatures from RNA-Seq data

Installation 

``` R

install.packages('devtools')
devtools::install_github('Zha0rong/ICARus',force = T)

```

For windows user:

Please install GO.db first from Biocmanager, then use devtools to install ICARus.

``` R
install.packages("BiocManager")
BiocManager::install("GO.db")
install.packages('devtools')
devtools::install_github('Zha0rong/ICARus',force = T)

```

Vignette in [R markdown format](https://github.com/Zha0rong/ICArus/blob/main/vignettes/ICARus.Rmd).

Pre-built Vignette in [HTML format](https://html-preview.github.io/?url=https://github.com/Zha0rong/ICArus/blob/main/vignettes/ICARus.htm).
