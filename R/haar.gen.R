#' Generating Multivariate Sieves with Haar Wavelets
#' @description The \code{haar.gen} command is used to generate multivariate sieves with bases of haar wavelets. It considers the interactions among these bases of different variables in a dataset.
#' In particular, if a user defines a test dataset which he hopes to use for his own model evaluation, the \code{haar.gen} command also generates the corresponding sieves of
#' the test dataset with the same number of basis used for each variable of the main dataset.
#'
#' @param data a data frame or a matrix containing the original variables.
#' @param test.data an optional data frame or matrix containing the same variables as the argument \code{data}. If \code{NULL}, no test dataset is considered.
#' @param n.basis a non-negative integer that specifies the level of haar bases for each variable. If \code{NULL}, it will be calculated according to the formulas established in this \href{https://github.com/ccfang2/Masters_Thesis}{thesis}.
#' @param max.interaction a positive integer that specifies the largest number of interacting variables in a single term of the resulting sieve. For example, \code{max.interaction=1} indicates no interaction among different variables. It shouldn't
#' be larger than the number of original variables.
#'
#' @return \code{haar.gen} generates multivariate sieves with haar wavelets, considering the interactions between basis of different variables. If \code{test.data} is supplied,
#' the output is a list of two data frames. One contains the resulting sieves for main dataset, and the other includes sieves for test dataset. \code{n.basis} used to generate sieves for both datasets are the same, so the sieves of test dataset could be applied
#' for the evaluation of model estimated with sieves of main dataset. The data frame in the output of \code{haar.gen} could be directly used as an input of the argument \code{x.sieve} in command \code{\link{dimada}}.
#' @seealso \link{poly.gen}; \link{bspline.gen}; \link{cosine.gen}; \link{sine.gen}; \link{trig.gen}; \link{daubechies.gen}; \link{dimada}.
#' @export
#'
#' @note It is worthy of mentioning that \code{n.basis} is the level of haar wavelets bases, not the number of bases. In definition of haar wavelets, the number of wavelet bases on each level increases exponentially, so the resulting bases may be of extremely
#' large quantity if the defined level is high. Please read this \href{https://github.com/ccfang2/Masters_Thesis}{thesis} for detailed description of wavelet bases. \code{haar.gen} is not recommended if \code{n.basis} is too large because it will take much
#' longer time to compute the sieve than other types of basis functions like B-splines.
#'
#' The command \code{\link[FlexCoDE]{haar_basis}} in package \pkg{FlexCoDE} is applied to create haar wavelets for each variable, and
#' the function \code{create_index_matrix} in package \pkg{Sieve} is used to help with the generation of multivariate sieves in this function \code{haar.gen}.
#'
#' @author Chencheng Fang, Bonn Graduate School of Economics, University of Bonn. Email: \email{ccfang@uni-bonn.de}
#'
#' @examples
#' # a data frame with 3 variables and 100 observations
#' set.seed(200)
#' df <- data.frame(a=runif(100), b=runif(100), c=runif(100))
#' # sieves with haar wavelets, considering partial interactions of these 3 variables and
#' # a self defined number of level of wavelet bases for each variable
#' sieve1 <- haar.gen(data=df, test.data=NULL, n.basis=NULL, max.interaction=2)
#'
#' # adding a test dataset with the same variables as the above main dataset
#' set.seed(200)
#' df.test <- data.frame(a=runif(100), b=runif(100), c=runif(100))
#' sieve2 <- haar.gen(data=df, test.data=df.test, n.basis=NULL, max.interaction=2)
haar.gen <- function(data,
                     test.data=NULL,
                     n.basis=NULL,
                     max.interaction=1)
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

  # compute the lower bound of smoothness for haar wavelets
  smoothness <- (ncol(data)+1)/2
  if (is.null(n.basis)) n.basis <- base::floor(log(nrow(data), base=10)*nrow(data)^(1/(2*smoothness+1)))

  # rescale the columns of matrix to [0,1]
  data <- base::apply(data, 2, function(x) (x-min(x))/(max(x)-min(x)))
  if(!is.null(test.data)) test.data <- base::apply(test.data, 2, function(x) (x-min(x))/(max(x)-min(x)))

  # --------------------------------
  # generate index tables
  # --------------------------------

  # expand.grid may be an alternative, but it can easily run out of vector memory
  if (ncol(data) %in% 1:20) {
    index.table <- base::eval(base::parse(text = base::paste("index.table.xdim",ncol(data),sep="")))
  } else {
    index.table <- Sieve::create_index_matrix(xdim = ncol(data), basisN = 10e4)[,-1] - 1
  }

  keep.rows.degree <- base::apply(index.table, 1, function(row) sum(row) <= 2^(n.basis+1)-2)
  index.table <- index.table[keep.rows.degree, , drop=FALSE]

  keep.rows.interact <- base::apply(index.table, 1, function(row) base::length(row[row!=0]) <= max.interaction)
  index.table <- index.table[keep.rows.interact, ,drop=FALSE]
  #index.table <- index.table[-1, ,drop=FALSE]
  colnames(index.table) <- var.names

  # --------------------------------
  # generate names for haar wavelets
  # --------------------------------

  # define a function to generate name for each term
  terms.names.gen <- function(index.row) {
    if (sum(index.row)==0) {
      name <- "Intercept"
    } else {
      name <- character()
      for (i in 1:length(index.row)) {
        if (index.row[i] !=0){
          name.temp <- base::paste(var.names[i],index.row[i], sep=".haar")
          name <- base::paste(name,name.temp,sep=" * ")
        } else name <- name
      }
      name <- base::substr(name,3,base::nchar(name))
    }
    return(name)
  }

  # generate names for all terms
  terms.names <- base::apply(index.table,1,terms.names.gen)
  terms.names <- base::unlist(terms.names)

  # --------------------------------
  # define a function of haar_basis (from R package 'FlexCoDE': https://github.com/rizbicki/FlexCoDE)
  # --------------------------------

  .pow2seq <- function(n) {
    seq <- rep(NA, n)
    pow <- 1
    for (ii in 1:n) {
      if (ii >= (2 ^ pow)) {
        pow <- pow + 1
      }
      seq[ii] <- pow
    }
    return(seq)
  }

  .haar_phi <- function(x) {
    if (0 <= x && x < 0.5) {
      return(1)
    } else if (0.5 <= x && x < 1) {
      return(-1)
    } else {
      return(0)
    }
  }

  haar_basis <- function(z, n_basis) {
    basis <- matrix(NA, length(z), n_basis)
    basis[, 1] <- 1.0

    kk <- 0
    jj <- 0

    for (ii in 2:n_basis) {
      if (jj == 2 ^ kk - 1) {
        kk <- kk + 1
        jj <- 0
      } else {
        jj <- jj + 1
      }
      basis[, ii] <- 2 ^ (kk / 2) * sapply(2 ^ kk * z - jj, .haar_phi)
    }

    attr(basis, "levels") <- .pow2seq(n_basis)
    return(basis)
  }

  # --------------------------------
  # compute values for haar wavelets
  # --------------------------------

  if (!is.null(test.data)) {

    haars <- plyr::llply(1:ncol(data),function(x) haar_basis(z=data[,x],n_basis=2^(n.basis+1)-1))
    terms.values.int <- base::apply(index.table,1, function(x) plyr::llply(1:ncol(index.table), function(y) haars[[y]][,(x[y]+1)]))
    terms.values <- as.data.frame(t(plyr::ldply(terms.values.int, function(x) Reduce("*", x))))
    colnames(terms.values) <- terms.names
    rownames(terms.values) <- NULL

    haars.test <- plyr::llply(1:ncol(test.data),function(x) haar_basis(z=test.data[,x],n_basis=2^(n.basis+1)-1))
    terms.values.int.test <- base::apply(index.table,1, function(x) plyr::llply(1:ncol(index.table), function(y) haars.test[[y]][,(x[y]+1)]))
    terms.values.test <- as.data.frame(t(plyr::ldply(terms.values.int.test, function(x) Reduce("*", x))))
    colnames(terms.values.test) <- terms.names
    rownames(terms.values.test) <- NULL

  } else {
    haars <- plyr::llply(1:ncol(data),function(x) haar_basis(z=data[,x],n_basis=2^(n.basis+1)-1))
    terms.values.int <- base::apply(index.table,1, function(x) plyr::llply(1:ncol(index.table), function(y) haars[[y]][,(x[y]+1)]))
    terms.values <- as.data.frame(t(plyr::ldply(terms.values.int, function(x) Reduce("*", x))))
    colnames(terms.values) <- terms.names
    rownames(terms.values) <- NULL
  }

  # --------------------------------
  # output
  # --------------------------------

  output <- list(train=terms.values)
  if (!is.null(test.data)) output <- append(output, list(test=terms.values.test))

  attr(output, "n.basis") <- n.basis
  return(output)

  options(warn=0)

}
