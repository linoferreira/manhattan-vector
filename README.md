# Make Manhattan plots as vector graphics in R with `ggplot2`

[Manhattan plots](https://en.wikipedia.org/wiki/Manhattan_plot) are a widely used tool in statistical genetics to visualise the results of genome-wide association studies (GWAS). While they are just basic scatter plots, since the number of points to display is often in the millions they become impractical to render as vector graphics in formats such as PDF or SVG. This is unfortunate as these formats are the gold standard for technical plots and are often requested by academic journals when submitting an article for publication.

This repo contains an R function to generate Manhattan plots with ggplot2 that can quickly be exported into a moderately-sized PDF or SVG file with `ggplot2::ggsave`. It is based on [Holtz Yan](https://github.com/holtzy/)'s excellent [Manhattan plot function](https://www.r-graph-gallery.com/101_Manhattan_plot.html) produced for [the R Graph Gallery](https://www.r-graph-gallery.com/index.html) which I extended by merging overlapping points into single shapes to simplify the resulting output.

The general idea is to use software for processing and plotting [geographic features](https://en.wikipedia.org/wiki/Simple_Features) (spefically, I use the [`sf` R package](https://r-spatial.github.io/sf/)) to convert each data point into a circle ([a polygon in simple features language](https://r-spatial.github.io/sf/articles/sf1.html)), merge overlapping circles into a single shape ([a union operation](https://r-spatial.github.io/sf/articles/sf3.html)) and [plot these simplified shapes](https://r-spatial.github.io/sf/articles/sf5.html).

The function (`ggmanh_vec`) is provided in the `fn-ggmanh_vec.R` script and I give a simple reproducible example in the script `make-manhattan-vec.R`. In this latter script, I download summary statistics of a GWAS of standing height in European-ancestry samples in the UK Biobank (includes 12 million variants) from [Watanabe et al. (2019)](https://doi.org/10.1038/s41588-019-0481-0) (this is available from the [GWAS Atlas](https://atlas.ctglab.nl/)) and make a basic Manhattan plot that is then exported into PDF format, producing the file `ukb-height-gwas.pdf` whose size is **751kB**.


## Runtime

Making the example Manhattan plot with 12 million points takes approximately 24min using four Intel Skylake 2.6GHz processors, each with 16GB of RAM. A simplified version in which we only plot points with p-value lower than or equal to 0.01 (approximately 1.2 million points) takes only 3min12s with the same resources (see the file `ukb-height-gwas-pv2.pdf`). For comparison, a standard rasterised PNG plot takes 2min14s to make with a slightly modified (and not parallelised) version of Holtz Yan's original function.


## Dependencies

The following R packages are required:

- [doParallel](https://cran.r-project.org/web/packages/doParallel/index.html)
- [dplyr](https://dplyr.tidyverse.org/)
- [foreach](https://cran.r-project.org/web/packages/foreach/index.html)
- [ggplot2](https://ggplot2.tidyverse.org/)
- [sf](https://r-spatial.github.io/sf/)

In addition, the `sf` R package requires that the [GDAL](https://gdal.org/) library be available on the system.

This code has been tested in R 3.6.2 (on CentOS 7.8.2003 with GDAL 3.0.2 installed) with `doParallel` 1.0.16, `dplyr` 1.0.6, `foreach` 1.5.1, `ggplot2` 3.3.3 and `sf` 1.0-0.
