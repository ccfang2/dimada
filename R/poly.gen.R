#' Generating Multivariate Sieves with Power Series or Legendre Polynomials
#' @description The \code{poly.gen} command is used to generate multivariate sieves with power series or Legendre polynomials bases. It considers the interactions among these bases of different variables in a dataset.
#' In particular, if a user defines a test dataset which he hopes to use for his own model evaluation, the \code{poly.gen} command also generates the corresponding sieves of
#' the test dataset with the same number of basis used for each variable of the main dataset. The command \code{\link[plde]{legendre_polynomial}} in package \pkg{plde} is applied to create Legendre polynomials.
#'
#' @param data a data frame or a matrix containing the original variables.
#' @param test.data an optional data frame or matrix containing the same variables as the argument \code{data}. If \code{NULL}, no test dataset is considered.
#' @param n.basis a non-negative integer that specifies the number of basis for each variable. If \code{NULL}, it will be calculated according to the formulas established in this \href{https://github.com/ccfang2/Masters_Thesis}{thesis}.
#' @param max.interaction a positive integer that specifies the largest number of interacting variables in a single term of the resulting sieve. For example, \code{max.interaction=1} indicates no interaction among different variables. It shouldn't
#' be larger than the number of original variables.
#' @param legendre a logical value that indicates if power series or Legendre polynomials are used for construction of basis for each variable. If \code{TRUE}, Legendre polynomials are applied; otherwise, power series.
#'
#' @return \code{poly.gen} generates multivariate sieves with power series or Legendre polynomials, considering the interactions between basis of different variables. If \code{test.data} is supplied,
#' the output is a list of two data frames. One contains the resulting sieves for main dataset, and the other includes sieves for test dataset. \code{n.basis} used to generate sieves for both datasets are the same, so the sieves of test dataset could be applied
#' for the evaluation of model estimated with sieves of main dataset. The data frame in the output of \code{poly.gen} could be directly used as an input of the argument \code{x.sieve} in command \code{\link{dimada}}.
#' @seealso \link{bspline.gen}; \link{cosine.gen}; \link{sine.gen}; \link{trig.gen}; \link{haar.gen}; \link{daubechies.gen}; \link{dimada}.
#' @export
#'
#' @note An internal function \code{create_index_matrix} in package \pkg{Sieve} is applied to help with the generation of multivariate sieves in function \code{poly.gen}.
#'
#' @author Chencheng Fang, Bonn Graduate School of Economics, University of Bonn. Email: \email{ccfang@uni-bonn.de}
#'
#' @examples
#' # a data frame with 3 variables and 100 observations
#' set.seed(200)
#' df <- data.frame(a=runif(100), b=runif(100), c=runif(100))
#' # sieves with power series, considering full interactions of these 3 variables
#' sieve1 <- poly.gen(data=df, test.data=NULL, n.basis=10, max.interaction=3, legendre=FALSE)
#'
#' # sieves with Legendre polynomials, considering no interactions of variables
#' # and only 5 basis for each variable
#' sieve2 <- poly.gen(data=df, test.data=NULL, n.basis=5, max.interaction=1, legendre=TRUE)
poly.gen <- function(data,
                     test.data=NULL,
                     n.basis=NULL,
                     max.interaction=1,
                     legendre=TRUE)
  {

  # --------------------------------
  # check arguments and preparation
  # --------------------------------

  options(warn=-1)

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

  # check if 'legendre' is logical
  if (!is.logical(legendre) | base::length(legendre)!=1) stop("The argument 'legendre' has to be a single logical value.")

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

  # compute the lower bound of smoothness for polynomials.
  # this set of basis functions theoretically includes the best approximation for all three cases of underlying model in my paper
  smoothness <- (3*ncol(data)-1)/2
  if (is.null(n.basis)) n.basis <- base::floor(log(nrow(data))*nrow(data)^(1/(2*smoothness+1)))

  # --------------------------------
  # generate index tables
  # --------------------------------

  if (ncol(data) %in% 1:30) {
    index.table <- base::eval(base::parse(text = base::paste("index.table.xdim",ncol(data),sep="")))
  } else {
    #index.table <- Sieve:::create_index_matrix(xdim = ncol(data), basisN = 10e5)[,-1] - 1
    stop("The number of variables is so large that the consequent sieve is of huge dimension,
         and it will then be computationally hard to perform dimension adaptive estimation.")
  }

  keep.rows.degree <- base::apply(index.table, 1, function(row) sum(row) <= n.basis)
  index.table <- index.table[keep.rows.degree, ,drop=FALSE]
  keep.rows.interact <- base::apply(index.table, 1, function(row) length(row[row!=0]) <= max.interaction)
  index.table <- index.table[keep.rows.interact, ,drop=FALSE]
  #index.table <- index.table[-1, ,drop=FALSE]
  colnames(index.table) <- var.names

  # --------------------------------
  # generate names for polynomial terms
  # --------------------------------

  # define a function to generate name for each term
  terms.names.gen <- function(index.row) {
    if (sum(index.row)==0) {
      name <- "Intercept"
    } else {
      name <- character()
      for (i in 1:length(index.row)) {
        if (index.row[i] != 0){
          name.temp <- base::paste(var.names[i],index.row[i], sep="^")
          name <- base::paste(name,name.temp,sep=" * ")
        } else name <- name
      }
      name <- base::substr(name,3,base::nchar(name))
    }
    return(name)
  }

  # generate names for all terms
  terms.names <- base::apply(index.table,1,terms.names.gen)

  if (!legendre) {
    terms.names <- base::unlist(terms.names)
  } else {
    terms.names <- base::paste(base::unlist(terms.names), "(Lg)")
  }

  # --------------------------------
  # compute values for polynomial terms
  # --------------------------------

  if (!legendre | n.basis==1) {
    if (!is.null(test.data)) {
      terms.values <- as.data.frame(t(base::apply(data,1, function(x) base::apply(index.table,1,function(y) base::prod(x^y)))))
      terms.values.test <- as.data.frame(t(base::apply(test.data,1, function(x) base::apply(index.table,1,function(y) base::prod(x^y)))))
      base::colnames(terms.values.test) <- terms.names
      base::rownames(terms.values.test) <- NULL
    } else {
      terms.values <- as.data.frame(t(base::apply(data,1, function(x) base::apply(index.table,1,function(y) base::prod(x^y)))))
    }
  } else {
    if (!is.null(test.data)) {
      legendre.poly <- plyr::llply(1:ncol(data),function(x) base::cbind(1,plde::legendre_polynomial(data[,x],list(dimension = n.basis))))
      terms.values.int <- base::apply(index.table,1, function(x) plyr::llply(1:ncol(index.table), function(y) legendre.poly[[y]][,(x[y]+1)]))
      terms.values <- as.data.frame(t(plyr::ldply(terms.values.int, function(x) Reduce("*", x))))
      legendre.poly.test <- plyr::llply(1:ncol(test.data),function(x) base::cbind(1,plde::legendre_polynomial(test.data[,x],list(dimension = n.basis))))
      terms.values.int.test <- base::apply(index.table,1, function(x) plyr::llply(1:ncol(index.table), function(y) legendre.poly.test[[y]][,(x[y]+1)]))
      terms.values.test <- as.data.frame(t(plyr::ldply(terms.values.int.test, function(x) Reduce("*", x))))
      base::colnames(terms.values.test) <- terms.names
      base::rownames(terms.values.test) <- NULL
    } else {
      legendre.poly <- plyr::llply(1:ncol(data),function(x) base::cbind(1,plde::legendre_polynomial(data[,x],list(dimension = n.basis))))
      terms.values.int <- base::apply(index.table,1, function(x) plyr::llply(1:ncol(index.table), function(y) legendre.poly[[y]][,(x[y]+1)]))
      terms.values <- as.data.frame(t(plyr::ldply(terms.values.int, function(x) Reduce("*", x))))
    }
  }

  base::colnames(terms.values) <- terms.names
  base::rownames(terms.values) <- NULL

  # --------------------------------
  # output
  # --------------------------------

  output <- list(train=terms.values)
  if (!is.null(test.data)) output <- append(output, list(test=terms.values.test))

  attr(output, "n.basis") <- n.basis
  return(output)

  options(warn=0)

}

