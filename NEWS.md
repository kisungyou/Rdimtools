# Rdimtools 1.1.0

* Removed arbitrary regularization in `do.dne()`.
* Change of authorship.

# Rdimtools 1.0.9

* Change of some estimation methods to return implicit variables.
* `do.dppca()` added : dual probabilistic principal component analysis.

# Rdimtools 1.0.8

* Change of maintainer's contact.
* Some functions are converted to pure C++.

# Rdimtools 1.0.7

* Fixed an error in `do.olpp()` thanks to Nicholas Spyrison.

# Rdimtools 1.0.6

* Replaced nearest neighbor library.

# Rdimtools 1.0.5

* Methods added; `do.mmds()`, `do.phate()`.
* Bibliographic errors fixed.

# Rdimtools 1.0.4

* Methods for feature selection `do.wdfs()`, `do.uwdfs()`, `do.procrustes()`, and `do.mifs()` are added.
* `README` now shows numbers of currently available functions for DR and IDE.

# Rdimtools 1.0.3

* Fixed memory leaks in `do.sne()` and `do.tsne()`.
* `do.lsls()` added as a supervised feature selection method.

# Rdimtools 1.0.2

* README contains minimal examples for both dimension reduction and estimation.
* Porting to pure C++ implementations started, gaining computational efficiency. 
* `do.lmds` function is fixed for its discrepancy in nested structure.


# Rdimtools 1.0.1

* NEWS reformatted and [package website](https://www.kisungyou.com/Rdimtools/) is now available.
* Dependency on R package [ADMM](https://CRAN.R-project.org/package=ADMM) is removed.


# Rdimtools 1.0.0

## Major changes
* LDA solves trace ratio problem directly.
* Many of dependencies are removed.
* 133 dimension reduction methods available.
* 17 intrinsic dimension estimation methods available.

## Bug fixes
* Error fixed in `do.lscore` function (thanks to Jordan Lin).

