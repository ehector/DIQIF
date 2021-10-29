# About

Joint integrative analysis of multiple data sources with correlated vector outcomes

This repository is an R package implementing methods from the following publication:
- Emily C. Hector and Peter X.-K. Song (2021). Joint integrative analysis of multiple data sources with correlated vector outcomes. arXiv, arXiv:2011.14996. To appear in The Annals of Applied Statistics.

Briefly, the method performs regression analysis of high-dimensional correlated responses. It divides data into blocks according to supplied indicators, 
analyses blocks using quadratic inference functions, and selectively combines blocks using a meta-estimator asymptotically equivalent to the optimal generalized method of moments equation.

Please email ehector@ncsu.edu with any questions or bug-reports.

# Installation

The DIQIF R package can be installed in one of two ways:
- from the downloaded gzipped tarball as R CMD INSTALL DIQIF_1.0-1.tar.gz
- from the downloaded and renamed DIQIF folder as R CMD build DIQIF and R CMD INSTALL DIQIF_1.0-1.tar.gz

Please make sure to have all packages listed in the DESCRIPTION file already installed.

# Citation

If you use the DIQIF R package, please consider citing the relevant manuscript: Hector and Song (2021).

# References

The posdef.matrix function was written by Ravi Varadhan: https://stat.ethz.ch/pipermail/r-help/2008-February/153708.

The GEE implementation is through the R package geepack: Højsgaard, S., Halekoh, U., Yan, J. (2006). The R package geepack for generalized estimating equations. Journal of Statistical Software, 15(2):1–11.

The QIF Rcpp implementation is courtesy of Lan Luo (https://github.com/luolsph). 
