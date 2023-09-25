# dimada <img src="man/figures/badge.png" align="right" alt="" width="155" />

[![PRs Welcome](https://img.shields.io/badge/PRs-welcome-brightgreen.svg?style=flat-square)](https://makeapullrequest.com)
![GitHub last commit](https://img.shields.io/github/last-commit/ccfang2/dimada?logo=GitHub)
![GitHub repo size](https://img.shields.io/github/repo-size/ccfang2/dimada?logo=GitHub)
![GitHub R package version](https://img.shields.io/github/r-package/v/ccfang2/dimada?logo=R)
![GitHub Repo stars](https://img.shields.io/github/stars/ccfang2/dimada?style=social)

> This package implements dimension adaptive estimation, which uses the sieve estimation and Lasso-type methods to obtain an estimator with the property of adapting its convergence rate to the unknown dimensionality of underlying true model. Specifically, if the underlying model is parametric, this estimator converges as fast as a parametric estimator. If the underlying model is non-parametric, its convergence rate might be slower but it still converges whilst the usual parametric estimator does not converge (in probability) to the true function. Hence, dimension adaptive estimator is weakly dominant in terms of approximation errors across different underlying models. The package 'dimada' allows users to easily perform such an estimator, together with other functions of plotting or summarizing the estimation results.

## Installation

You can install the development version of dimada from [GitHub](https://github.com/) with:
      
``` r
# install.packages("devtools")
devtools::install_github("ccfang2/dimada")
```

## Key Commands

The main command in this package is `dimada()`. Users are required to import response variable and a dataset of original independent variables. To construct sieves, users can choose from a variety of basis functions including power series, Legendre polynomials, B-splines and trigonometric polynomials, etc. Since a multivariate setup is considered, users also need to define the maximum number of interacting variables in a single term of sieves. If users already have the dataset of sieves, they could directly import the sieves instead of the original variables in this command. 

In terms of Lasso-type methods, users can choose from Lasso, [adaptive Lasso](https://www.tandfonline.com/doi/abs/10.1198/016214506000000735) and [twin adaptive Lasso](https://www.sciencedirect.com/science/article/abs/pii/S030440762100049X). Lasso is the baseline method, so it has to be included in the argument of `methods` in `dimada()`. An example of this command is given below.

```r
library(dimada)

# a data frame with 3 independent variables
set.seed(200)
df <- data.frame(a=runif(100), b=runif(100), c=runif(100))

# simulate a response variable with a non-parametric additive model and a normal error term
set.seed(200)
response2 <- sin(2*df[,1])+tan(df[,2])+log(df[,3])+rnorm(100,0.04)

# dimension adaptive estimation with B-splines and no interactions among variables
dimada2 <- dimada(y=response2, x=df, basis="bspline", max.interaction=1, methods=c("Lasso","adaLasso"))
```

To quickly plot out the main results of `dimada()`, one can use `dimada.plot()`. This command depicts the change of cross-validated approximation errors against lambdas used in Lasso-type methods or the consequent numbers of non zero coefficients. 

```r
plot <- dimada.plot(dimada2,x="nzeros")
```

The plot is as follows.

<img src="man/figures/plot.png" alt="Plot" width="600"/>    

In particular, post-selection OLS model is performed after each Lasso-type method to obtain unbiased coefficients because the estimated coefficients from Lasso-type methods are biased due to the penalty term in the loss function. Then, the coefficients from post-selection OLS model could be used for prediction in new dataset. The results of Lasso-type methods and the post-selection OLS models can be quickly printed out with the command `dimada.summary()`.

```r
dimada.summary(dimada2,verbose=TRUE)
```

The summary is as follows.

```
======================================== 
          Sieve 
======================================== 
Basis function: B-Splines
Number of basis in each term: 5
Maximum number of interactions in each term: 1
Time used: 0.2960498 secs
 
======================================== 
          LASSO
======================================== 
Selection rule of lambda: Value of lambda that gives minimum MSE
Selected lambda: 3.28e-02
Number of selected terms: 11
Cross-validated Empirical MSE: 1.22e+00
Time used: 1.135178 secs
 
Coefficients: 
| a.bs1   -0.17783| 
| b.bs1   -0.12352| 
| c.bs1   -3.94617| 
| b.bs2   -0.89202| 
| c.bs2   -0.48332| 
| a.bs3    0.27788| 
| a.bs4    0.77288| 
| b.bs4    0.11791| 
| c.bs4    0.94437| 
| b.bs5    0.69366| 
| c.bs5    0.30913| 

======================================== 
          Post LASSO
======================================== 
Method: OLS on selected terms
In-sample MSE: 9.01e-01
 
Coefficients: 
|(Intercept)   0.35890| 
|` a.bs1`     -0.15855| 
|` b.bs1`     -0.19250| 
|` c.bs1`     -4.09278| 
|` b.bs2`     -0.91587| 
|` c.bs2`     -0.54928| 
|` a.bs3`      0.39223| 
|` a.bs4`      0.89797| 
|` b.bs4`      0.14660| 
|` c.bs4`      0.95196| 
|` b.bs5`      0.81853| 
|` c.bs5`      0.34737| 

======================================== 
          Adaptive LASSO
======================================== 
Selection rule of lambda: Value of lambda that gives minimum MSE
Selected lambda: 7.68e-02
Number of selected terms: 8
Cross-validated Empirical MSE: 1.19e+00
Time used: 1.235119 secs
 
Coefficients: 
| c.bs1   -4.02574| 
| b.bs2   -1.10738| 
| c.bs2   -0.39484| 
| a.bs3    0.20018| 
| a.bs4    0.86268| 
| c.bs4    1.09413| 
| b.bs5    0.65998| 
| c.bs5    0.11182| 

======================================== 
          Post Adaptive LASSO
======================================== 
Method: OLS on selected terms
In-sample MSE: 9.06e-01
 
Coefficients: 
|(Intercept)   0.31862| 
|` c.bs1`     -4.06092| 
|` b.bs2`     -1.06783| 
|` c.bs2`     -0.54647| 
|` a.bs3`      0.49315| 
|` a.bs4`      0.93811| 
|` c.bs4`      0.98340| 
|` b.bs5`      0.82966| 
|` c.bs5`      0.36596| 
```

As is seen, the cross-validated (i.e., out-of-sample) MSEs for Lasso and adaptive Lasso are 1.22 and 1.19. However, the coefficients are biased due to the penalty term in loss function. To remove such an impact, post-Lasso and post-adaptive Lasso are performed, and MSEs for these post-selection methods are 0.901 and 0.906. Let's compare them to the MSE using parametric OLS estimator. The MSE from OLS `lm(response2~df[,1]+df[,2]+df[,3])` is 1.257. Hence, it is noted that in this non-parametric additive underlying model, our dimension adaptive estimator achieves a much smaller in-sample approximation error than parametric OLS estimator, about 2/3 of it. 

However, in practice, it is adviced to split the original dataset into train data and test data, then perform our dimension adaptive estimator in train data and finally use the coefficients from post-selection methods to make predictions in test data. The out-of-sample MSE of our estimator can thus be computed and compared with other estimators.

## Other Commands

In this package, there are also a series of functions that help to generate multivariate sieves from original dataset. Particularly, they consider the interactions among variables in the construction of sieves. These functions are suffixed with `.gen`, and the data frame in the output of these functions can be used as input of the argument `x.sieve` in command `dimada()`. Following is an example of generating sieves with power series by using `poly.gen()`.

```r
# a data frame with 3 variables and 100 observations
set.seed(200)
df <- data.frame(a=runif(100), b=runif(100), c=runif(100))
# sieves with power series, considering full interactions of these 3 variables
sieve1 <- poly.gen(data=df, test.data=NULL, n.basis=10, max.interaction=3, legendre=FALSE)
```

It is worthy of mentioning that users could also compute sieves for a test dataset, which can be used for evaluation of models estimated with sieves of the main dataset. 

## Log of Change

<details><summary>Version 1.0.2 (Date: 10.9.2023) </summary>
      <ol>
            <li> Fix an error in <code>trig.gen()</code> so it considers the interactions between cosine and sine components in trigonometric polynomials.
            <li> The function <code>create_index_matrix()</code> in package <code>Sieve</code> is now exportable, so I could call the function directly in my package.
            <li> In <code>dimada()</code>, outliers in the estimated coefficients from Lasso-type methods are excluded in post-selection methods.
            <li> In <code>dimada()</code>, the argument <code>seed</code> is removed.
      </ol>
</details>

<details><summary>Version 1.0.1 (Date: 30.8.2023) </summary>
      <ol>
            <li> Fix an error in the output of <code>trig.gen()</code>: the output now contains the data frames of sieves for both train and test datasets if original test dataset is given.
            <li> Add an argument <code>save.sieve</code> in the function <code>dimada()</code> so that one can choose if he hopes to save the data frame of generated or given sieves in the output of <code>dimada()</code>. The sieves can sometimes be very large in size, so if a user won't use it in future, it is advised not to save it in the output which may take up your machine memory.
      </ol>
</details>

## Note
- This package is a part of my [thesis](https://github.com/ccfang2/Masters_Thesis), which is supervised by Prof. Dr. Joachim Freyberger. This is an ongoing project, and future improvement on this package is expected.
- The picture on the hexagon badge of this package is drawn with an [AI image creator](https://openai.com/dall-e-2) of Dall·E, OpenAI.

## Contact

Chencheng Fang, Email: [ccfang[at]uni-bonn.de](mailto:ccfang@uni-bonn.de),
Bonn Graduate School of Economics, University of Bonn, Germany

