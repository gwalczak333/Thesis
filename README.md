# Deep Dive into Gaussian Processes: Model Comparisons in Assessing Variable Importance

In statistical genetics, one may want to understand how a genetic variant affects a complex trait along with how the variant could disproportionately affect individuals from different subpopulations, and this can be achieved using the ”GlObal And Local Score” (GOALS) operator, which provides a way to assess local and global variable importance simultaneously in nonlinear models (Winn-Nu ̃nez et al.,2023).

In this paper, in addition to fitting a larger dataset using traditional linear models, a nonlinear, non-parametric Bayesian approach called Gaussian process (GP) regression will be utilized to implement the GOALS operator. To expand on this, I will use deep Gaussian process regression (DGP), which is essentially a multi-layer Gaussian process that can be implemented via Markov Chain Monte Carlo (Booth, 2023). Using a large dataset relating gene expression measurements to observed traits in a population of mice (Valdar et al.,
2006), I will use Gaussian process regression and deep Gaussian process regression to implement the GOALS operator. To compare the models, I will evaluate predicted values against the true values at the testing locations (Booth, 2023). 


## R-package for GOALS
[deepgp](https://cran.r-project.org/web/packages/deepgp/vignettes/deepgp.html)

[BGLR](https://cran.r-project.org/package=BGLR)

[LearnPCA](https://cran.r-project.org/package=LearnPCA)

[ggplot2](https://cran.r-project.org/package=ggplot2)

[dplyr](https://cran.r-project.org/package=dplyr)

[fields](https://cran.r-project.org/package=fields)


## Relevent Citations

Booth, A. S. (2023). deepgp: an R-package for Bayesian Deep Gaussian Processes.

Bryan A. Hanson, D. T. H. (2022). LearnPCA: Functions, Data Sets and Vignettes to Aid
in Learning Principal Components Analysis (PCA).

Gustavo de los Campos, P. P. R. (2023). BGLR: Bayesian Generalized Linear Regression.

Valdar, W., L. C. Solberg, D. Gauguier, S. Burnett, P. Klenerman, W. O. Cookson, M. S.
Taylor, J. N. Rawlins, R. Mott, and J. Flint (2006). Genome-wide genetic association of
complex traits in heterogeneous stock mice. Nature Genetics 38 (8), 879–887.

Winn-Nu ̃nez, E. T., M. Griffin, and L. Crawford (2023). A simple approach for local and
global variable importance in nonlinear regression models.
