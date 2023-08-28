#' Plotting the Output of Function \code{dimada}
#' @description The \code{dimada.plot} command aims to quickly plot the output of function \code{dimada}. The plot depicts the change of cross-validated approximation errors against lambdas used in Lasso-type methods
#' or the consequent numbers of non zero coefficients. If more than one method is defined in the argument \code{methods} of function \code{dimada}, all subplots will be combined into a single plot.
#'
#' @param object an object of S3 class \code{"dimada"}, which is an output of the function \code{dimada}.
#' @param x a character string that specifies the x-axis. Available options are \code{"nzeros"} and \code{"lambdas"}. The former represents the numbers of non zero coefficients and the latter indicates lambdas
#' used in Lasso-type methods. Default is \code{"nzeros"}.
#'
#' @return An object of \code{ggplot2} plot is returned. Specifically, cross-validated approximation errors (i.e., mean squared errors) are plotted against the numbers of non zero coefficients or lambdas. A red star is
#' used to mark the selected lambda or the consequently chosen number of non-zero coefficients, and the value of corresponding mean squared error is denoted on the red star.
#' @seealso \link{dimada}; \link{dimada.summary}.
#' @importFrom magrittr %>%
#' @export
#'
#' @author Chencheng Fang, Bonn Graduate School of Economics, University of Bonn. Email: \email{ccfang@uni-bonn.de}
#'
#' @examples
#' # plot the output 'dimada2' in the example of function 'dimada'
#' set.seed(200)
#' df <- data.frame(a=runif(100), b=runif(100), c=runif(100))
#' response2 <- sin(2*df[,1])+tan(df[,2])+log(df[,3])+rnorm(100,0.04)
#' dimada2 <- dimada(y=response2, x=df, basis="bspline", max.interaction=1,
#'                   methods=c("Lasso","adaLasso"))
#' (plot <- dimada.plot(dimada2,x="nzeros"))
dimada.plot <- function(object,
                        x="nzeros")
  {

  # --------------------------------
  # checking arguments
  # --------------------------------

  # check if 'object' is of class 'dimada'
  if (!base::inherits(object,"dimada")) stop("The argument 'object' should be an output of dimada function.")

  # check if 'x' is defined correctly
  if (any(base::length(x)!=1, !(x %in% c("nzeros","lambdas")))) stop("The argument 'x' has to be 'nzeros' or 'lambdas'.")

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

  xlabel <- dplyr::case_when(x=="nzeros" ~ "Number of Non Zero Coefficients",
                             x=="lambdas" ~ "Penalty Parameter: Lambda")

  methods <- base::names(object)

  data <- do.call(rbind, plyr::llply(1:length(methods), function(x) cbind(object[[x]]$parameters,method=object[[x]]$method))) %>%
    dplyr::mutate(method=factor(method, ordered=TRUE, levels=c("LASSO","Adaptive LASSO", "Twin Adaptive LASSO")))

  data.final <- do.call(rbind, plyr::llply(1:length(methods), function(x) cbind(object[[x]]$parameters.final, method=object[[x]]$method))) %>%
    dplyr::mutate(method=factor(method, ordered=TRUE, levels=c("LASSO","Adaptive LASSO", "Twin Adaptive LASSO")))

  # --------------------------------
  # plot
  # --------------------------------

  plots <- ggplot2::ggplot(data = data, ggplot2::aes(x=eval(parse(text=x)), y=cvms))+
    ggplot2::geom_point(shape=1)+
    ggplot2::geom_line()+
    ggplot2::geom_point(data = data.final, ggplot2::aes(x = eval(parse(text=x)), y = cvms), shape = 8, size = 6, color = "red") +
    ggplot2::geom_text(data = data.final, ggplot2::aes(x = eval(parse(text=x)), y = cvms, label = base::sprintf("%.2e",cvms)), vjust = -0.5, hjust=0.5) +
    ggplot2::labs(title="Cross-Validated Empirical MSE", subtitle = paste("Basis Functions: ",basis,sep=""))+
    ggplot2::xlab(xlabel)+
    ggplot2::ylab("MSE")+
    ggplot2::labs(caption=paste("The cross indicates the point selected by rule '", attr(object, "s"),"'.", sep=""))+
    ggplot2::facet_wrap(facets=ggplot2::vars(method), ncol=1,scales="free")+
    ggthemes::theme_base()

  return(plots)

}
