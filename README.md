# CACTA internship project
This project was completed during my internship with the [Catoni Lab](https://catonilab.com/), University of Birmingham, which focuses on plant epigenetics and genome plasticity. The project investigates a family of transposable elements (TEs) known as the CACTA family. This family is characterised by having an initial start sequence of CACTA, along with some other characteristics explained later. Here, I document my work looking for these TEs within 3 tree species: *Quercus robur* (English oak), *Corylus avelana* (Common hazel), and *Fraxinus excelsoir* (European hazel).

## CACTA transposons background
Transposable elements are sequences of DNA with the ability to translocate across a genome. Transposon families are categorised into being either class I retrotransposons, requiring reverse transcription of mRNA, or class II transposases, requiring transposase proteins to be encoded first. The target family for this research, CACTA, is a class II transposon, requiring the synthesis of a transpoase gene initially before translocation can occur. The paper by [Muñoz-López and García-Pérez (2010)][PMC2874221] provides a more detailed explanation of general transposons.

##Prerequisite R libaries
* [ggplot2](https://ggplot2.tidyverse.org/): used for all graphical plots
* [Bioconductor](https://www.bioconductor.org): used to install all bioinformatics-related R packaged
* [packFinder](http://www.bioconductor.org/packages/release/bioc/html/packFinder.html): used to return possible CACTA sequences given specific characteristics
* seqinr

## Prerequisite genome data
* Haploid *Q. robur* genome: [Oak Genome Sequencing FR](https://www.oakgenome.fr/?page_id=587)


## Key Findings




## References
Catoni, M., Jonesman, T., Cerruti, E. and Paszkowski, J., 2018. Mobilization of Pack-CACTA transposons in Arabidopsis suggests the mechanism of gene shuffling. *Nucleic Acids Research*, 47(3), pp.1311-1320. doi: <https://doi.org/10.1093/nar/gky1196>

Gisby, J. and Catoni, M., 2021. *packFinder: de novo Annotation of Pack-TYPE Transposable Elements*. R package version 1.4.0, <https://github.com/jackgisby/packFinder>.

Munoz-Lopez, M. and Garcia-Perez, J., 2010. DNA Transposons: Nature and Applications in Genomics. *Current Genomics*, 11(2), pp.115-128. PMCID: [PMC2874221]

[PMC2874221]: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2874221/






