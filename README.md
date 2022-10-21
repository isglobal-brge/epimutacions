# epimutacions

`epimutacions` is an R/Biconductor package.  The package provides 
6 statistical methods for outlier detection in DNA methylation data. 
It also contains function for results visualization and annotation. 


# Installation

The `epimutacions` package can be installed using the following code:

```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("epimutacions")
```

The `epimutacions` package depends on `epimutacionsData`  
package which can be installed  as following: 


```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("epimutacionsData")

```

# Usage 

* See vignette: shorturl.at/gwFOQ

# Acknowledgements

The authors would like to thank the team who collaborated in the initial 
design of the package in the European BioHackathon 2020:  
Lordstrong Akano, James Baye, 
Alejandro Caceres, Pavlo Hrab, Raquel Manzano and Margherita Mutarelli. 
The authors also want to thank the 
organization of European BioHackathon 2020 for its support.


# Contact

* Dolors Pelegri-Siso: dolors.pelegri@isglobal.org
* Juan R. Gonzalez: juanr.gonzalez@isglobal.org
* Carlos Ruiz-Arenas: carlos.ruiza@upf.edu
* Carles Hernandez-Ferrer: carles.hernandez@cnag.crg.eu
