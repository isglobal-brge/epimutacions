# epimutacions v1.17.2

* Reduced the dependency footprint: removed the unused `ensembldb` dependency
  and moved the optional annotation/visualization packages (Homo.sapiens, Gviz,
  rtracklayer, AnnotationHub, ExperimentHub, the UCSC TxDb packages, the
  Illumina manifest/annotation packages, reshape2, purrr, ggrepel and
  gridExtra) from Imports to Suggests, as they are already used conditionally
  via `requireNamespace()`.

# epimutacions v1.17.1

* Fixed `plot_epimutations()` failing with a non-numeric `value` column.
* Made the vignette annotation robust to Ensembl/biomaRt network outages.
* Fixed the license DCF stub and an Rd `\itemize` NOTE.

# epimutacions v1.1.2

* fixed plots
* fixed annotations

# epimutacions v0.99.33

* Bugs and notes (if possible) fixed


# epimutacions v0.99.0

* Submitted to Bioconductor
