#' Generating Multivariate Sieves with Trigonometric Polynomials
#' @description The \code{trig.gen} command is used to generate multivariate sieves with trigonometric polynomials. It considers the interactions among trigonometric polynomial bases of different variables in a dataset.
#' In particular, if a user defines a test dataset which he hopes to use for his own model evaluation, the \code{trig.gen} command also generates the corresponding sieves of
#' the test dataset with the same number of basis used for each variable of the main dataset.
#'
#' @param data a data frame or a matrix containing the original variables.
#' @param test.data an optional data frame or matrix containing the same variables as the argument \code{data}. If \code{NULL}, no test dataset is considered.
#' @param n.basis a non-negative integer that specifies the number of trigonometric polynomial basis for each variable. If \code{NULL}, it will be calculated according to the formulas established in this \href{https://github.com/ccfang2/Masters_Thesis}{thesis}.
#' @param max.interaction a positive integer that specifies the largest number of interacting variables in a single term of the resulting sieve. For example, \code{max.interaction=1} indicates no interaction among different variables. It shouldn't
#' be larger than the number of original variables.
#'
#' @return \code{trig.gen} generates multivariate sieves with trigonometric polynomials, considering the interactions between basis of different variables. If \code{test.data} is supplied,
#' the output is a list of two data frames. One contains the resulting sieves for main dataset, and the other includes sieves for test dataset. \code{n.basis} used to generate sieves for both datasets are the same, so the sieves of test dataset could be applied
#' for the evaluation of model estimated with sieves of main dataset. The data frame in the output of \code{trig.gen} could be directly used as an input of the argument \code{x.sieve} in command \code{\link{dimada}}.
#' @seealso \link{poly.gen}; \link{bspline.gen}; \link{cosine.gen}; \link{sine.gen}; \link{haar.gen}; \link{daubechies.gen}; \link{dimada}.
#' @export
#'
#' @note By definition, trigonometric polynomials are sum of cosine and sine polynomials. Hence, in the construction of \code{trig.gen} command, functions \code{cosine.gen} and \code{sine.gen} are called upon and the resulting sieves from them are then combined
#' to get the final sieve for \code{trig.gen}.
#'
#' @author Chencheng Fang, Bonn Graduate School of Economics, University of Bonn. Email: \email{ccfang@uni-bonn.de}
#'
#' @examples
#' # a data frame with 3 variables and 100 observations
#' set.seed(200)
#' df <- data.frame(a=runif(100), b=runif(100), c=runif(100))
#' # sieves with trigonometric polynomials, considering interactions of at most 2 variables
#' sieve1 <- trig.gen(data=df, test.data=NULL, n.basis=10, max.interaction=2)
#'
#' # adding a test dataset with the same variables as the above main dataset
#' set.seed(200)
#' df.test <- data.frame(a=runif(100), b=runif(100), c=runif(100))
#' sieve2 <- trig.gen(data=df, test.data=df.test, n.basis=10, max.interaction=2)
trig.gen <- function(data,
                     test.data=NULL,
                     n.basis=NULL,
                     max.interaction=1)
  {

  # --------------------------------
  # check arguments and preparation
  # --------------------------------

  options(warn=0)

  # check if 'data' is a matrix or data frame
  if (all(!is.data.frame(data), !is.matrix(data))) stop("The argument 'data' has to be a data frame or a matrix.")
  if (is.data.frame(data) && base::length(data)==0) stop("The argument 'data' cannot be an empty data frame.")
  if (is.matrix(data) && all(is.na(data))) stop("The argument 'data' cannot be an empty matrix.")

  # check if 'test.data' is a matrix or data frame, and if it has the same number of column(s) as 'data'
  if (all(!is.data.frame(test.data), !is.matrix(test.data), !is.null(test.data))) stop("The argument 'test.data' has to be a data frame, a matrix or NULL.")
  if (is.data.frame(test.data) && base::length(test.data)==0) stop("The argument 'test.data' cannot be an empty data frame.")
  if (is.matrix(test.data) && all(is.na(test.data))) stop("The argument 'test.data' cannot be an empty matrix.")
  if (!is.null(test.data) && ncol(data)!=ncol(test.data)) stop("The argument 'test.data' should have the same number of columns as 'data'.")

  # check if 'n.basis' is NULL or a single non-negative integer
  if (all(!is.null(n.basis), (n.basis<0 | n.basis%%1!=0 | base::length(n.basis)!=1))) stop("The argument 'n.basis' has to be either NULL or a single non-negative integer.")

  # check if 'max.interaction' is a single positive integer and smaller than the number of regressors
  if (max.interaction<=0 | max.interaction>ncol(data) | max.interaction%%1!=0 | base::length(max.interaction) !=1) stop("The argument 'max.interaction' has to be a single positive integer and smaller than the number of regressors.")

  # convert 'data' and 'test.data' to matrix if they are data frames
  if (is.data.frame(data)) data <- as.matrix(data)
  if (!is.null(test.data) && is.data.frame(test.data)) test.data <- as.matrix(test.data)

  # define variable names if not yet
  if (is.null(colnames(data))) {
    var.names <- colnames(data) <- base::paste("X", 1:ncol(data), sep="")
  } else {
    var.names <- colnames(data)
  }

  if (!is.null(test.data) && is.null(colnames(test.data))) var.names.test <- colnames(test.data) <- base::paste("X", 1:ncol(test.data), sep="")
  if (!is.null(test.data) && !is.null(colnames(test.data))) var.names.test <- colnames(test.data)

  # check if the set of variables in 'data' is the same as the set of variables in 'test.data'
  # if yes, sort variables in 'test.data' according to that in 'data'
  if (!is.null(test.data) && !identical(sort(var.names), sort(var.names.test))) {
    stop("The variable names in 'data' and 'test.data' are not the same. Please check.")
  } else {
    test.data <- test.data[,var.names]
  }

  # compute the lower bound of smoothness for cosine polynomials.
  # this set of basis functions theoretically includes the best approximation for all three cases of underlying model in my paper
  smoothness <- (ncol(data)+1)/2
  if (is.null(n.basis)) n.basis <- base::floor(log(nrow(data))*nrow(data)^(1/(2*smoothness+1)))

  # --------------------------------
  # use cosine.gen and sine.gen
  # --------------------------------

  cosine.object <- cosine.gen(data=data,
                              test.data=test.data,
                              n.basis=n.basis,
                              max.interaction=max.interaction)

  sine.object <- sine.gen(data=data,
                          test.data=test.data,
                          n.basis=n.basis,
                          max.interaction=max.interaction)

  # --------------------------------
  # output
  # --------------------------------

  output <- list(train=cbind(cosine.object$train, sine.object$train))
  if (!is.null(test.data)) output <- append(output, list(test=cbind(cosine.object$test, sine.object$test)))

  attr(output, "n.basis") <- n.basis
  return(output)

  options(warn=0)

}

