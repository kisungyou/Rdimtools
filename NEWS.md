# Rdimtools 0.4.3

 * MINOR CHANGE: LDA solves trace ratio problem directly.

 * FUNCTIONALITY: Added auxiliary functions:
   - `aux.traceratio()`: Solve trace ratio problem with 2012 Ngo's algorithm.
   
* DOCUMENTATION: Added NEWS.md to keep record of changes.

# Rdimtools 0.4.2 (2018-12-21)

 * DOCUMENTATION: Minor changes in documentation.

# Rdimtools 0.4.1 (2018-11-15)

 * FUNCTIONALITY: Added linear methods:
   - `do.disr()`: Diversity-Induced Self-Representation.
   - `do.lspe()`: Locality and Similarity Preserving Embedding.
   - `do.nrsr()`: Non-convex Regularized Self-Representation.
   - `do.rsr()`: Regularized Self-Representation.
   - `do.spc()`: Supervised Principal Component Analysis.
   - `do.spufs()`: Structure Preserving Unsupervised Feature Selection.
   - `do.udfs()`: Unsupervised Discriminative Features Selection.
   
 * FUNCTIONALITY: Added auxiliary functions:
   - `aux.which.mink()`: Returns index of smallest element.
   - `aux.which.maxk()`: Returns index of largest element.
   
# Rdimtools 0.4.0 (2018-09-28)

 * FUNCTIONALITY: Added linear methods:
   - `do.fastmap()`: FastMap.
   - `do.idmap()`: Interactive Document Map.
   - `do.lamp()`: Local Affine Multidimensional Scaling.
   - `do.llle()`: Local Linear Laplacian Eigenmaps.
   - `do.nnp()`: Nearest Neighbor Projection.
   - `do.spmds()`: Spectral Multidimensional Scaling.
     
 * FUNCTIONALITY: Added auxiliary functions:
   - `aux.findmaxidx()`: Find the row and column index of maximal elements.
   - `aux.randpartition()`: Given 1:n, divide it into K random partitions without replacement.
   
 * FUNCTIONALITY: Added dimensionality estimation methods:
   - `est.clustering()`: Clustering-based Estimation
   - `est.incisingball()`: Intrinsic Dimension Estimation with Incising Ball
   - `est.made()`: Manifold-Adaptive Dimension Estimation.
   - `est.mle1()`: MLE using Poisson Process.
   - `est.mle2()`: MLE using Poisson Process with Bias Correction.
   - `est.nearneighbor1()`: Near-Neighbor Information.
   - `est.nearneighbor2()`: Near-Neighbor Information with Bias Correction.
   - `est.packing()`: Estimation using Packing Numbers.
   - `est.pcathr()`: PCA Thresholding with Accumulated Variance.
   - `est.twonn()`: Minimal Neighborhood Information.
   - `est.Ustat()`: Convergence Rate of U-statistic.

# Rdimtools 0.3.2 (2018-03-06)

 * FUNCTIONALITY: Added linear methods:
   - `do.enet()`: Elastic Net Regularization.
     
 * FUNCTIONALITY: Added auxiliary functions:
   - `aux.oosprocess()`: Data processing for OOS prediction.
   
 * FUNCTIONALITY: Added out-of-sample prediction: 
   - `oos.linproj()`: Linear Projection.
   
# Rdimtools 0.3.1 (2018-02-06)

 * FUNCTIONALITY: Added auxiliary functions: 
   - `aux.pinv()`: Use SVD and NumPy scheme.
   - `aux.bicgstab()`: Rlinsolve can't be used.
   
 * MINOR CHANGE: `aux.geigen()` is now implemented using RcppArmadillo.

# Rdimtools 0.3.0 (2018-01-30)
 
 * FUNCTIONALITY: Added linear methods:
   - `do.adr()`: Adaptive Dimension Reduction.
   - `do.ammc()`: Adaptive Maximum Margin Criterion.
   - `do.asi()`: Adaptive Subspace Iteration.
   - `do.cnpe()`: Complete Neighborhood Preserving Embedding.
   - `do.crp()`: Collaborative Representation-based Projection.
   - `do.dagdne()`: Double-Adjacency Graphs-based Discriminant Neighborhood Embedding.
   - `do.dne()`: Discriminant Neighborhood Embedding.
   - `do.elde()`: Exponential Local Discriminant Embedding.
   - `do.elpp2()`: Enhanced Locality Preserving Projection (2013).
   - `do.fscore()`: Fisher Score.
   - `do.lasso()`: Least Absolute Shrinkage and Selection Operator.
   - `do.ldakm()`: Combination of LDA and K-means.
   - `do.ldp()`: Locally Discriminating Projection.
   - `do.lltsa()`: Linear Local Tangent Space Alignment.
   - `do.lpca()`: Locally Principal Component Analysis.
   - `do.lpe()`: Locality Pursuit Embedding.
   - `do.lpfda()`: Locality Preserving Fisher Discriminant Analysis.
   - `do.lpmip()`: Locality-Preserved Maximum Information Projection.
   - `do.lqmi()`: Linear Quadratic Mutual Information.
   - `do.lscore()`: Laplacian Score.
   - `do.lsda()`: Locality Sensitive Discriminant Analysis.
   - `do.lsdf()`: Locality Sensitive Discriminant Feature.
   - `do.lsir()`: Localized Sliced Inverse Regression.
   - `do.mcfs()`: Multi-Cluster Feature Selection.
   - `do.mfa()`: Marginal Fisher Analysis.
   - `do.mlie()`: Maximal Local Interclass Embedding.
   - `do.mmp()`: Maximum Margin Projection.
   - `do.mmsd()`: Multiple Maximum Scatter Difference.
   - `do.modp()`: Modified Orthogonal Discriminant Projection.
   - `do.msd()`: Maximum Scatter Difference.
   - `do.nolpp()`: Nonnegative Orthogonal Locality Preserving Projection.
   - `do.nonpp()`: Nonnegative Orthogonal Neighborhood Preserving Projections.
   - `do.npca()`: Nonnegative Principal Component Analysis.
   - `do.odp()`: Orthogonal Discriminant Projection.
     
 * FUNCTIONALITY: Added auxiliary functions:
   - `aux.adjqr()`: Adjust projection matrix using QR decomposition.
   - `aux.nbdlogical()`: Find homogeneous and heterogeneous neighborhood indexing.
   - `aux.geigen()`: geigen in my taste.
   - `aux.featureindicator()`: Generate (p-by-ndim) indicator matrix for projection.
   - `aux.traceratio.max()`: Compute trace ratio problem for maximal basis.

# Rdimtools 0.2.0 (2018-01-03)

 * FUNCTIONALITY: Added intrinsic dimension estimation methods for exploratory analysis.
 
 * DOCUMENTATION: Many fixes in documentation.
 
 * FUNCTIONALITY: Added linear methods:
   - `do.anmm()`: Average Neighborhood Margin Maximization.
   - `do.bpca()`: Bayesian Principal Component Analysis.
   - `do.dspp()`: Discriminative Sparsity Preserving Projection.
   - `do.eslpp()`: Extended Supervised Locality Preserving Projection.
   - `do.extlpp()`: Extended Locality Preserving Projection.
   - `do.isoproj()`: Isometric Projection.
   - `do.kmvp()`: Kernel-Weighted Maximum Variance Projection.
   - `do.kudp()`: Kernel-Weighted Unsupervised Discriminant Projection.
   - `do.lde()`: Local Discriminant Embedding.
   - `do.lea()`: Locally Linear Embedded Eigenspace Analysis.
   - `do.lfda()`: Local Fisher Discriminant Analysis-
   - `do.llp()`: Local Learning Projections.
   - `do.lspp()`: Local Similarity Preserving Projection.
   - `do.mmc()`: Maximum Margin Criterion.
   - `do.mvp()`: Maximum Variance Projection.
   - `do.ppca()`: Probabilistic Principal Component Analysis.
   - `do.sdlpp()`: Sample-Dependent Locality Preserving Projection.
   - `do.slpp()`: Supervised Locality Preserving Projection.
   - `do.spca()`: Sparse Principal Component Analysis.
   - `do.spp()`: Sparsity Preserving Projection.
   - `do.udp()`: Unsupervised Discriminant Projection.
 
 * FUNCTIONALITY: Added non-linear methods:
   - `do.crca()`: Curvilinear Component Analysis.
   - `do.crda()`: Curvilinear Distance Analysis.
   - `do.ispe()`: Isometric Stochastic Proximity Embedding.
   - `do.klde()`: Kernel Local Discriminant Embedding.
   - `do.klfda()`: Kernel Local Fisher Discriminant Analysis.
   - `do.kmmc()`: Kernel Maximium Margin Criterion.
     
 * FUNCTIONALITY: Added auxiliary functions:
   - `aux.graphnbdD()`: Construct nearest-neighborhood graph when distance matrix is already given.
   - `aux.kernelcentering()`: Centering the kernel/gram matrix.
   - `aux.kernelprojection()`: Given uncentered gram matrix, find the projected data.
   - `aux.adjprojection()`: Adjust projection matrix by simply normalizing each column.
   
 * FUNCTIONALITY: Added out-of-sample prediction: 
   - `oos.linear()`: Only for linear methods.

# Rdimtools 0.1.3 (2017-11-25)

 * FUNCTIONALITY: Added linear methods: CCA, OLPP, OPLS, PLS.
 
 * FUNCTIONALITY: Added auxiliary functions: 
   - `aux.pkgstat()`: Show the number of available functions.
 
# Rdimtools 0.1.2 (2017-11-14)
 
 * FUNCTIONALITY: Added linear method:
   - `do.lda()`: Linear Discriminant Analysis.
   
 * FUNCTIONALITY: Added non-linear method: 
   - `do.ree()`: Robust Euclidean Embedding.
 
 * DOCUMENTATION: Move references to separate bib file.

# Rdimtools 0.1.1 (2017-09-24)

 * DOCUMENTATION: Many fixes in roxygen documentation.
 
 * FUNCTIONALITY: Added non-linear methods: 
   - `do.keca()`: Kernel Entropy Component Analysis.
   - `do.ltsa()`: Local Tangent Space Alignment.
 
# Rdimtools 0.1.0 (2017-09-21)
Initial release of Rdimtools, a rich collection of linear and nonlinear
dimension reduction techniques implemented using 'RcppArmadillo'.

Linear dimensionality reduction methods:
 * `do.fa()`: Factor analysis.
 * `do.ica()`: Independent Component Analysis.
 * `do.lmds()`: Landmark Multidimensional Scaling.
 * `do.lpp()`: Locality Preserving Projection. 
 * `do.mds()`: (Metric) Multidimensional Scaling.
 * `do.npe()`: Neighborhood Preserving Embedding. 
 * `do.pca()`: Principal Component Analysis. 
 * `do.rndproj()`: Random Projection.
 
 Non-linear dimensionality reduction methods:
 * `do.cisomap()`: Conformal Isometric Feature Mapping.
 * `do.dm()`: Diffusion Maps.
 * `do.isomap()`: Isometric Feature Mapping.
 * `do.kpca()`: Kernel Principal Component Analysis.
 * `do.lapeig()`: Laplacian Eigenmaps.
 * `do.lisomap()`: Landmark Isometric Feature Mapping.
 * `do.lle()`: Locally Linear Embedding.
 * `do.mvu()`: Maximum Variance Unfolding / Semidefinite Embedding.
 * `do.plp()`: Piecewise Laplacian Projection.
 * `do.sammon()`: Sammon Mapping.
 * `do.sne()`: Stochastic Neighbor Embedding.
 * `do.tsne()`: t-distributed Stochastic Neighbor Embedding.
 
Dimensionality estimation methods:
 * `est.boxcount()`: Box-counting dimension, also known as Minkowski-Bouligand dimension.
 * `est.correlation()`: Correlation dimension.
 
Auxiliary functions:
 * `aux.typecheck()`: Check whether the data is poorly given.
 * `aux.preprocess()`: Preprocessing, Centering, decorrelation, or whitening.
 * `aux.gensamples()`: Generate a few popular data examples.
 * `aux.graphnbd()`: Compute neighborhood graph.
 * `aux.shortestpath()`: Compute shortest paths given neighborhood graph.
 * `aux.MaxMinLandmark()`: Choose a single landmark point.
 * `aux.kernekcov()`: Build K and centerd K matrix for kernel tricks.
 * `aux.eigendec()`: Use Armadillo in a descending order.

