# An Ising Similarity Regression Model

This repository code contains template code associated with the manuscript "[An Ising Similarity Regression Model for Modeling Multivariate Binary Data](https://aus01.safelinks.protection.outlook.com/?url=http%3A%2F%2Fwww3.stat.sinica.edu.tw%2Fss_newpaper%2FSS-2024-0021_na.pdf&data=05%7C02%7Czhiyang.tho%40anu.edu.au%7Cc183b7dbfd264377735308dd1eff913a%7Ce37d725cab5c46249ae5f0533e486437%7C0%7C0%7C638700805362576902%7CUnknown%7CTWFpbGZsb3d8eyJFbXB0eU1hcGkiOnRydWUsIlYiOiIwLjAuMDAwMCIsIlAiOiJXaW4zMiIsIkFOIjoiTWFpbCIsIldUIjoyfQ%3D%3D%7C0%7C%7C%7C&sdata=xJCgNK%2F2LDsF%2BBuL3XjbXOa17iHu2p%2BuIAavzrNRvow%3D&reserved=0)" by [Tho](https://rsfas.anu.edu.au/about/staff-directory/zhi-yang-tho), [Hui](https://francishui.netlify.app/), and [Zou](https://cbe.anu.edu.au/about/staff-directory/dr-tao-zou), which is accepted for publication in [Statistica Sinica](https://www3.stat.sinica.edu.tw/statistica/).

# Getting started

There are currently three directories in this repository:

-   `Code`, which contains `R` scripts implementing the proposed penalized pseudo-likelihood estimator along with other scripts for simulating data from the Ising similarity regression model and calculating empirical sandwich covariance matrix.

-   `Simulation`, which contains script to implement the simulation study in the manuscript. **Users are recommended to start here by examining `simulation_script.R` to understand how to obtain the penalized pseudo-likelihood estimator, among other estimators.**

-  `Application`, which contains script applying the Ising similarity regression model to the U.S. Senate roll call voting data of the 117-th congress. Also contained in this folder are a set of csv files that include the roll call voting data and senators' attributes from the [U.S. Senate's website](https://www.senate.gov/), additional senators' information from Wikipedia, and senators' Twitter data.

# If you find any bugs and issues...

If you find something that looks like a bug/issue, please use Github issues and post it up there. As much as possible, please include in the issue:

1.  A description of the bug/issue;
2.  Paste-able code along with some comments that reproduces the problem e.g., using the [reprex](https://cran.r-project.org/web/packages/reprex/index.html) package. If you also have an idea of how to fix the problem, then that is also much appreciated;
3.  Required data files etc...

Alternatively, please contact the corresponding author at [ZhiYang.Tho\@anu.edu.au](mailto:ZhiYang.Tho@anu.edu.au).
