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



[Click here to download OGRE's Vignette as PDF](https://github.com/svenbioinf/OGRE/files/10057271/The.OGRE.user.guide.pdf)






https://user-images.githubusercontent.com/68910402/203071312-ff92a42d-05c6-4c83-9905-f47e477305e9.mp4







## Installation via GitHub
For installing OGRE by GitHub use the install_github() function.

```{bash
if(!requireNamespace("remotes", quietly = TRUE)){
    install.packages("remotes")
}
remotes::install_github("svenbioinf/OGRE",build_vignettes = TRUE)

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

## Installation via docker

Using docker you can easily install OGRE with the latest R and RStudio software
and all required packages already included. 
For this to work, your local computer needs to have docker installed. 
```{r}
sudo apt install docker.io #on linux via console
```
For windows just follow instruction on: [docker.com](https://docs.docker.com/desktop/install/windows-install/)
You can then pull OGRE's docker file from dockerHub and install by copying the
following to your shell (console)
```{r}
sudo docker pull svenbioinf/ogre:1 #on linux shell
sudo docker run -e PASSWORD=ogre -p 8787:8787 svenbioinf/ogre:1 #on linux shell

docker pull svenbioinf/ogre:1 #on windows console
docker run -e PASSWORD=ogre -p 8787:8787 svenbioinf/ogre:1 #on windows console

```
After installation, open RStudio using your favorite browser on 
[http://localhost:8787](http://localhost:8787)
username=rstudio, password=ogre. 
You can then load OGRE, browse the vignette or run the GUI with
```{r}
library("OGRE")
vignette("OGRE")
SHREC()
```

## Contact

This software was developed by Sven Berres svenbioinf(AT)gmail.com
