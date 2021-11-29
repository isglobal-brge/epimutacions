# epimutacions

`epimutacions` is a R package to detect epimutations (rare alterations in the methylation pattern at specific loci) from DNA methylation data.


# Installation

```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("epimutacions")
```

The `epimutacions` package depends on `epimutacionsData`  package which can be installed  as following: 


```
install.packages("devtools")
devtools::install_github("LeireAbarrategui/epimutacionsData")

```

# Usage 

* See vignette:  https://github.com/isglobal-brge/epimutacions/blob/master/vignettes/epimutacions.html

# Acknowledgements

The authors would like to thank the team who collaborated in the initial design of the package in the European BioHackathon 2020:  Lordstrong Akano, James Baye, Alejandro Caceres, Pavlo Hrab, Raquel Manzano and Margherita Mutarelli. The authors also want to thank the organization of European BioHackathon 2020 for its support.


# Contact

* Leire Abarrategui: leire.abarrategui-martinez@newcastle.ac.uk
* Juan R. Gonzalez: juanr.gonzalez@isglobal.org
* Carlos Ruiz-Arenas: carlos.ruiza@upf.edu
* Carles Hernandez-Ferrer: carles.hernandez@cnag.crg.eu
