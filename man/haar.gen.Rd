% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/haar.gen.R
\name{haar.gen}
\alias{haar.gen}
\title{Generating Multivariate Sieves with Haar Wavelets}
\usage{
haar.gen(data, test.data = NULL, n.basis = NULL, max.interaction = 1)
}
\arguments{
\item{data}{a data frame or a matrix containing the original variables.}

\item{test.data}{an optional data frame or matrix containing the same variables as the argument \code{data}. If \code{NULL}, no test dataset is considered.}

\item{n.basis}{a non-negative integer that specifies the level of haar bases for each variable. If \code{NULL}, it will be calculated according to the formulas established in this \href{https://github.com/ccfang2/Masters_Thesis}{thesis}.}

\item{max.interaction}{a positive integer that specifies the largest number of interacting variables in a single term of the resulting sieve. For example, \code{max.interaction=1} indicates no interaction among different variables. It shouldn't
be larger than the number of original variables.}
}
\value{
\code{haar.gen} generates multivariate sieves with haar wavelets, considering the interactions between basis of different variables. If \code{test.data} is supplied,
the output is a list of two data frames. One contains the resulting sieves for main dataset, and the other includes sieves for test dataset. \code{n.basis} used to generate sieves for both datasets are the same, so the sieves of test dataset could be applied
for the evaluation of model estimated with sieves of main dataset. The data frame in the output of \code{haar.gen} could be directly used as an input of the argument \code{x.sieve} in command \code{\link{dimada}}.
}
\description{
The \code{haar.gen} command is used to generate multivariate sieves with bases of haar wavelets. It considers the interactions among these bases of different variables in a dataset.
In particular, if a user defines a test dataset which he hopes to use for his own model evaluation, the \code{haar.gen} command also generates the corresponding sieves of
the test dataset with the same number of basis used for each variable of the main dataset.
}
\note{
It is worthy of mentioning that \code{n.basis} is the level of haar wavelets bases, not the number of bases. In definition of haar wavelets, the number of wavelet bases on each level increases exponentially, so the resulting bases may be of extremely
large quantity if the defined level is high. Please read this \href{https://github.com/ccfang2/Masters_Thesis}{thesis} for detailed description of wavelet bases. \code{haar.gen} is not recommended if \code{n.basis} is too large because it will take much
longer time to compute the sieve than other types of basis functions like B-splines.

The command \code{\link[FlexCoDE]{haar_basis}} in package \pkg{FlexCoDE} is applied to create haar wavelets for each variable, and
the function \code{create_index_matrix} in package \pkg{Sieve} is used to help with the generation of multivariate sieves in this function \code{haar.gen}.
}
\examples{
# a data frame with 3 variables and 100 observations
set.seed(200)
df <- data.frame(a=runif(100), b=runif(100), c=runif(100))
# sieves with haar wavelets, considering partial interactions of these 3 variables and
# a self defined number of level of wavelet bases for each variable
sieve1 <- haar.gen(data=df, test.data=NULL, n.basis=NULL, max.interaction=2)

# adding a test dataset with the same variables as the above main dataset
set.seed(200)
df.test <- data.frame(a=runif(100), b=runif(100), c=runif(100))
sieve2 <- haar.gen(data=df, test.data=df.test, n.basis=NULL, max.interaction=2)
}
\seealso{
\link{poly.gen}; \link{bspline.gen}; \link{cosine.gen}; \link{sine.gen}; \link{trig.gen}; \link{daubechies.gen}; \link{dimada}.
}
\author{
Chencheng Fang, Bonn Graduate School of Economics, University of Bonn. Email: \email{ccfang@uni-bonn.de}
}
