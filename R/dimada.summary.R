#' Summarizing the Output of Function \code{dimada}
#' @description The \code{dimada.summary} command aims to quickly summarize the output of function \code{dimada}. The summary includes an overview of the sieve that is used in dimension adaptive estimation and the information
#' of selected model for each Lasso-type method and the consequent terms with non-zero coefficients. Moreover, as the selected terms are also used to fit an OLS model, the results of these post selection models are summarized as well.
#'
#' @param object an object of S3 class \code{"dimada"}, which is an output of the function \code{dimada}.
#' @param verbose a logical value that specifies whether \code{dimada.summary} should suppress the summary to be directly printed on console. If \code{FALSE}, no summary will be printed out directly but users can save it to
#' an object and later print it out. See example below. Default is \code{TRUE}.
#'
#' @return \code{dimada.summary} uses \code{cat()} to print out the summary of \code{"dimada"} object on console, if not suppressed. Also, \code{dimada.summary} saves the same output as a character vector invisibly, so users could use
#' \code{cat()} to print it out later at their demand, as shown in examples. The summary is composed of important information in the output of function \code{dimada}.
#' @seealso \link{dimada}; \link{dimada.plot}.
#' @export
#'
#' @author Chencheng Fang, Bonn Graduate School of Economics, University of Bonn. Email: \email{ccfang@uni-bonn.de}
#'
#' @examples
#' # summarize the output 'dimada2' in the example of function 'dimada'
#' # and save the summary to an object called 'summa'
#' set.seed(200)
#' df <- data.frame(a=runif(100), b=runif(100), c=runif(100))
#' response2 <- sin(2*df[,1])+tan(df[,2])+log(df[,3])+rnorm(100,0.04)
#' dimada2 <- dimada(y=response2, x=df, basis="bspline", max.interaction=1,
#'                   methods=c("Lasso","adaLasso"))
#' summary <- dimada.summary(dimada2,verbose=TRUE)
#'
#' # print out the saved summary
#' cat(paste(summary,"\n"))
dimada.summary <- function(object,
                           verbose=TRUE)
  {

  # --------------------------------
  # checking arguments
  # --------------------------------

  # check if 'object' is of class 'dimada'
  if (!base::inherits(object,"dimada")) stop("The argument 'object' should be an output of function dimada.")

  # check if 'verbose' is a single logical value
  if (!is.logical(verbose) | base::length(verbose)!=1) stop("The argument 'verbose' has to be a single logical value.")

  # --------------------------------
  # data preparation
  # --------------------------------

  basis <- dplyr::case_when(attr(object, "basis")=="power" ~ "Power Series",
                            attr(object, "basis")=="legendre" ~ "Legendre Polynomials",
                            attr(object, "basis")=="bspline" ~ "B-Splines",
                            attr(object, "basis")=="haar" ~ "Haar Wavelets",
                            attr(object, "basis")=="daubechies" ~ "Daubechies Wavelets",
                            attr(object, "basis")=="cosine" ~ "Cosine Polynomials",
                            attr(object, "basis")=="sine" ~ "Sine Polynomials",
                            attr(object, "basis")=="trig" ~ "Trigonometric Polynomials")

  lambda.rule <- dplyr::case_when(attr(object, "s")=="lambda.min" ~ "Value of lambda that gives minimum MSE",
                                  attr(object, "s")=="lambda.1se" ~ "Largest value of lambda such that error is within 1 standard error of the minimum")

  methods <- base::names(object)

  data <- plyr::llply(1:length(methods), function(x) list(method=object[[x]]$method,
                                                          parameters.final=object[[x]]$parameters.final,
                                                          coefs.final=format(object[[x]]$coefs.final, digits=5),
                                                          post.lm=object[[x]]$post.lm,
                                                          elapse=object[[x]]$elapse))

  # coefficients
  tables.space <- plyr::llply(1:length(methods), function(x) max(nchar(data[[x]]$coefs.final$terms))-nchar(data[[x]]$coefs.final$terms)+4-(data[[x]]$coefs.final$coefs<0))
  tables <- plyr::llply(1:length(methods), function(x) plyr::llply(1:nrow(data[[x]]$coefs.final), function(y)  c("|",base::paste(data[[x]]$coefs.final[y,],collapse = paste0(rep('\40', tables.space[[x]][y]), collapse = '') ),"| \n") ))

  # coefficients for Post Methods
  tables.post.space <- plyr::llply(1:length(methods), function(x) max(nchar(names(data[[x]]$post.lm$coefficients)))-nchar(names(data[[x]]$post.lm$coefficients))+2) #-(data[[x]]$post.lm$coefficients<0)
  df.post <- plyr::llply(1:length(methods), function(x) format(data.frame(terms=names(data[[x]]$post.lm$coefficients), coefficients=data[[x]]$post.lm$coefficients), digits=5))
  tables.post <- plyr::llply(1:length(methods), function(x) plyr::llply(1:nrow(summary(data[[x]]$post.lm)$coefficients), function(y) c("|",base::paste(df.post[[x]][y,],collapse = paste0(rep('\40', tables.post.space[[x]][y]), collapse = '')),"| \n")  ))

  # --------------------------------
  # syntax
  # --------------------------------

  output <- plyr::llply(1:length(methods), function(x) c("======================================== \n",
                                                         "\40 \40 \40 \40 \40 ", data[[x]]$method, "\n",
                                                         "======================================== \n",
                                                         "Selection rule of lambda: ", lambda.rule,"\n",
                                                         "Selected lambda: ",base::sprintf("%.2e",data[[x]]$parameters.final$lambdas), "\n",
                                                         "Number of selected terms: ", data[[x]]$parameters.final$nzeros, "\n",
                                                         "Cross-validated Empirical MSE: ", base::sprintf("%.2e",data[[x]]$parameters.final$cvms), "\n",
                                                         "Time used: ", base::format(data[[x]]$elapse, units="auto"), "\n \n",
                                                         "Coefficients: \n",
                                                         tables[[x]],"\n",
                                                         "======================================== \n",
                                                         "\40 \40 \40 \40 \40 ","Post ", data[[x]]$method, "\n",
                                                         "======================================== \n",
                                                         "Method: OLS on selected terms\n",
                                                         "In-sample MSE: ", base::sprintf("%.2e",mean(data[[x]]$post.lm$residuals^2)), "\n \n",
                                                         "Coefficients: \n",
                                                         tables.post[[x]],"\n"))

  sieve <- c("======================================== \n",
             "\40 \40 \40 \40 \40 Sieve \n",
             "======================================== \n",
             base::ifelse(is.na(basis), "Sieve is directly given by the user, so the following is unknown.\n",""),
             "Basis function: ", basis,"\n",
             ifelse(attr(object, "basis") %in% c("haar", "daubechies"),"Number of levels for each original variable: ", "Number of basis for each original variable: "), attr(object, "n.basis"),"\n",
             "Maximum number of interactions in each term: ", attr(object, "max.interaction"),"\n",
             "Time used: ", attr(object, "sieve.time"), "\n \n")

  # --------------------------------
  # output
  # --------------------------------

  if (verbose) {
    cat(sieve, unlist(output), sep="")
  } else base::message("The output is suppressed by your command, but it is saved in an object, if assigned. See Examples to find out how to print out the saved output.")

  syntax <- utils::capture.output(cat(sieve, unlist(output), sep=""))

}

