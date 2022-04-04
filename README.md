<p align="center">
<img src="vignettes/logo.png" align="center" class="center" alt="" width="250"/> 
</p> <br>
<p align="center">
<b>OGRE</b> - <b>O</b>verlapping annotated <b>G</b>enomic <b>RE</b>gions <br>
Calculate, visualize and analyse overlap between genomic regions
</p>


The OGRE package calculates overlap between user defined genomic region datasets. 
Any regions can be supplied i.e. genes, SNPs, or reads from sequencing experiments. 
Key numbers help analyse the extend of overlaps which can also be visualized at a genomic level.
 
[OGRE_Vignette.pdf](https://github.com/svenbioinf/OGRE/files/8412826/OGRE_Vignette.pdf)


https://user-images.githubusercontent.com/68910402/161344697-da98cb01-4104-41b2-b441-4ef0dbfecd51.mp4






## Installation via GitHub
For installing OGRE by GitHub use the install_github() function.

```{bash
if(!requireNamespace("remotes", quietly = TRUE)){
    install.packages("remotes")
}
remotes::install_github("svenbioinf/OGRE")

```

OGRE depends on the following packages:<br>

- IRanges<br>
- GenomicRanges<br>
- S4Vectors<br>
- methods<br>
- data.table<br>
- assertthat<br>
- ggplot2<br>
- Gviz<br>
- AnnotationHub<br>
- shiny<br>
- shinyFiles<br>
- DT<br>
- rtracklayer<br>
- shinydashboard<br>
- shinyBS<br>
- tidyr<br>
- GenomeInfoDb<br>

if not automatically installed, those can be installed with:
```{r}
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c(
"IRanges",
"GenomicRanges",
"S4Vectors",
"methods",
"data.table",
"assertthat",
"ggplot2",
"IRanges",
"Gviz",
"AnnotationHub",
"shiny",
"shinyFiles",
"DT",
"rtracklayer",
"shinydashboard",
"shinyBS",
"tidyr",
"GenomeInfoDb"
))
```    

The OGRE package itself can then be loaded with the following commands:
```{r}
library(OGRE) # load package
vignette("OGRE") #some information on how to use the package
```

## Installation via Bioconductor

(Once OGRE is available on Bioconductor, if not- use devel version below.)

Start R and enter:

```{r}
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("OGRE")
library(OGRE)    # load package
vignette("OGRE") # some information on how to use the package
```
## Installation via Bioconductor (develop version)

Start R and enter:

```{r}
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("OGRE",version = "devel")
library(OGRE)    # load package
vignette("OGRE") # some information on how to use the package
```



## Contact

This software was developed by Sven Berres svenbioinf(AT)gmail.com
