#' Generating Multivariate Sieves with Daubechies Wavelets
#' @description The \code{daubechies.gen} command is used to generate multivariate sieves with bases of daubechies wavelets. It considers the interactions among these bases of different variables in a dataset.
#' In particular, if a user defines a test dataset which he hopes to use for his own model evaluation, the \code{daubechies.gen} command also generates the corresponding sieves of
#' the test dataset with the same number of basis used for each variable of the main dataset.
#'
#' @param data a data frame or a matrix containing the original variables.
#' @param test.data an optional data frame or matrix containing the same variables as the argument \code{data}. If \code{NULL}, no test dataset is considered.
#' @param n.basis a non-negative integer that specifies the level of daubechies bases for each variable. If \code{NULL}, it will be calculated according to the formulas established in this \href{https://github.com/ccfang2/Masters_Thesis}{thesis}.
#' @param max.interaction a positive integer that specifies the largest number of interacting variables in a single term of the resulting sieve. For example, \code{max.interaction=1} indicates no interaction among different variables. It shouldn't
#' be larger than the number of original variables.
#'
#' @return \code{daubechies.gen} generates multivariate sieves with daubechies wavelets, considering the interactions between basis of different variables. If \code{test.data} is supplied,
#' the output is a list of two data frames. One contains the resulting sieves for main dataset, and the other includes sieves for test dataset. \code{n.basis} used to generate sieves for both datasets are the same, so the sieves of test dataset could be applied
#' for the evaluation of model estimated with sieves of main dataset. The data frame in the output of \code{daubechies.gen} could be directly used as an input of the argument \code{x.sieve} in command \code{\link{dimada}}.
#' @seealso \link{poly.gen}; \link{bspline.gen}; \link{cosine.gen}; \link{sine.gen}; \link{trig.gen}; \link{haar.gen}; \link{dimada}.
#' @export
#'
#' @note It is worthy of mentioning that \code{n.basis} is the level of daubechies wavelets bases, not the number of bases. In definition of daubechies wavelets, the number of wavelet bases on each level increases exponentially, so the resulting bases may be of extremely
#' large quantity if the defined level is high. Please read this \href{https://github.com/ccfang2/Masters_Thesis}{thesis} for detailed description of wavelet bases. Here, \code{daubechies.gen} is not recommended if \code{n.basis} is too large because it will take much
#' longer time to compute the sieve than other types of basis functions like B-splines.
#'
#' The command \code{\link[FlexCoDE]{daubechies_basis}} in package \pkg{FlexCoDE} is applied to create daubechies wavelets for each variable, and
#' an internal function \code{create_index_matrix} in package \pkg{Sieve} is used to help with the generation of multivariate sieves in this function \code{daubechies.gen}.
#'
#' @author Chencheng Fang, Bonn Graduate School of Economics, University of Bonn. Email: \email{ccfang@uni-bonn.de}
#'
#' @examples
#' # a data frame with 3 variables and 100 observations
#' set.seed(200)
#' df <- data.frame(a=runif(100), b=runif(100), c=runif(100))
#' # sieves with daubechies wavelets, considering full interactions of these 3 variables
#' sieve1 <- daubechies.gen(data=df, test.data=NULL, n.basis=NULL, max.interaction=3)
#'
#' # sieves with daubechies wavelets, considering no interactions of variables and
#' # a self defined number of level of wavelet bases for each variable
#' sieve2 <- daubechies.gen(data=df, test.data=NULL, n.basis=5, max.interaction=1)
daubechies.gen <- function(data,
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

  # compute the lower bound of smoothness for daubechies wavelets
  smoothness <- (ncol(data)+1)/2
  if (is.null(n.basis)) n.basis <- base::floor(log(nrow(data), base=5)*nrow(data)^(1/(2*smoothness+1)))

  # rescale the columns of matrix to [0,1]
  data <- base::apply(data, 2, function(x) (x-min(x))/(max(x)-min(x)))
  if(!is.null(test.data)) test.data <- base::apply(test.data, 2, function(x) (x-min(x))/(max(x)-min(x)))

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

  keep.rows.degree <- base::apply(index.table, 1, function(row) sum(row) <= 2^(n.basis+1)-2)
  index.table <- index.table[keep.rows.degree, , drop=FALSE]

  keep.rows.interact <- base::apply(index.table, 1, function(row) base::length(row[row!=0]) <= max.interaction)
  index.table <- index.table[keep.rows.interact, ,drop=FALSE]
  #index.table <- index.table[-1, ,drop=FALSE]
  colnames(index.table) <- var.names

  # --------------------------------
  # generate names for daubechies wavelets
  # --------------------------------

  # define a function to generate name for each term
  terms.names.gen <- function(index.row) {
    if (sum(index.row)==0) {
      name <- "Intercept"
    } else {
      name <- character()
      for (i in 1:length(index.row)) {
        if (index.row[i] !=0){
          name.temp <- base::paste(var.names[i],index.row[i], sep=".db")
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
  # define a function of daubechies_basis ( in reference to R package 'FlexCoDE': https://github.com/rizbicki/FlexCoDE )
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

  daubechies_basis <- function(z, n_basis, n_aux_basis = max(n_basis, 2^12),
                               filter_number = 10, family = "DaubExPhase") {
    # Strategy is to project z onto a fine grid to use NN approximation
    aux <- wavethresh::GenW(n_aux_basis, filter.number = filter_number,
                            family = family)

    z_grid <- seq(0, 1, length.out = n_aux_basis)
    which_closest <- FNN::get.knnx(z_grid, as.matrix(z), k = 1)$nn.index

    aux[, 2:ncol(aux)] <- aux[, ncol(aux):2]
    aux <- aux[, 1:n_basis, drop = FALSE]
    basis <- aux[which_closest, ] / aux[1, 1]

    if (n_basis == 1) {
      basis <- as.matrix(basis)
    }

    attr(basis, "levels") <- .pow2seq(n_basis)
    return(basis)
  }

  # --------------------------------
  # compute values for daubechies wavelets
  # --------------------------------

  if (!is.null(test.data)) {

    daubechies <- plyr::llply(1:ncol(data),function(x) daubechies_basis(z=data[,x],n_basis=2^(n.basis+1)-1))
    terms.values.int <- base::apply(index.table,1, function(x) plyr::llply(1:ncol(index.table), function(y) daubechies[[y]][,(x[y]+1)]))
    terms.values <- as.data.frame(t(plyr::ldply(terms.values.int, function(x) Reduce("*", x))))
    colnames(terms.values) <- terms.names
    rownames(terms.values) <- NULL

    daubechies.test <- plyr::llply(1:ncol(test.data),function(x) daubechies_basis(z=test.data[,x],n_basis=2^(n.basis+1)-1))
    terms.values.int.test <- base::apply(index.table,1, function(x) plyr::llply(1:ncol(index.table), function(y) daubechies.test[[y]][,(x[y]+1)]))
    terms.values.test <- as.data.frame(t(plyr::ldply(terms.values.int.test, function(x) Reduce("*", x))))
    colnames(terms.values.test) <- terms.names
    rownames(terms.values.test) <- NULL

  } else {
    daubechies <- plyr::llply(1:ncol(data),function(x) daubechies_basis(z=data[,x],n_basis=2^(n.basis+1)-1))
    terms.values.int <- base::apply(index.table,1, function(x) plyr::llply(1:ncol(index.table), function(y) daubechies[[y]][,(x[y]+1)]))
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
