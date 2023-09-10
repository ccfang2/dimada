#' Dimension Adaptive Estimation with Sieves and Lasso-type Regularization
#' @description The \code{dimada} command implements dimension adaptive estimation with sieves and Lasso-type methods as proposed in this \href{https://github.com/ccfang2/Masters_Thesis}{thesis}.
#' The key advantage of dimension adaptive estimator is that its rate of convergence is adaptive to the dimension of underlying true model, which is usually unknown. Specifically, if the underlying
#' model is parametric, the convergence rate of this estimator is as fast as parametric estimator (e.g., OLS); if the underlying model is non-parametric, the convergence rate of this estimator is
#' slower but it still converges whilst the usual parametric method doesn't converge at all. To put it simply, the dimension adaptive estimator is weakly dominant in terms of approximation
#' error across different underlying models. One can implement this estimator easily with this \code{dimada} command.
#'
#' @param y a numeric vector of response variable
#' @param x a data frame or a matrix containing the original independent variables. The number of observations must be the same as the length of \code{y}. If the argument \code{x.sieve} is given directly, then \code{x} should be \code{NULL}.
#' @param basis a character string that specifies the basis with which the sieves of \code{x} are constructed. Available bases in \code{dimada} include \code{"power"}, \code{"legendre"}, \code{"bspline"}, \code{"haar"}, \code{"daubechies"},
#' \code{"cosine"}, \code{"sine"} and \code{"trig"}. Default is \code{"power"}, i.e., power series. If the argument \code{x.sieve} is given directly, there is no need to define \code{basis}, and one can simply leave it as default.
#' @param n.basis a non-negative integer that specifies the number of basis (or the level of basis in cases of haar and daubechies wavelets) for each variable. If \code{NULL}, it will be calculated according to the formulas established in this \href{https://github.com/ccfang2/Masters_Thesis}{thesis}.
#' Default is \code{NULL}. If the argument \code{x.sieve} is given, there is no need to define \code{n.basis}, and one can just leave it as default.
#' @param max.interaction a positive integer that specifies the largest number of interacting variables in a single term of the resulting sieve. For example, \code{max.interaction=1} indicates no interaction among different variables. It shouldn't
#' be larger than the number of original variables.
#' @param x.sieve a data frame or a matrix already containing the sieves of original dataset. If it is given, the arguments \code{x}, \code{basis}, \code{n.basis} and \code{max.interaction} are not needed. The data frame in the output of functions like \code{poly.gen} can be directly used as
#' an input of \code{x.sieve}. Default is \code{NULL}, when \code{x} must be provided.
#' @param methods a character vector that specifies the methods one hopes to use to select significant terms in the sieve. Available methods comprise \code{"Lasso"}, \code{"adaLasso"} and \code{"taLasso"}. Importantly, \code{"Lasso"} is the baseline method which must be included in this argument.
#' \code{"adaLasso"} stands for \href{https://www.tandfonline.com/doi/abs/10.1198/016214506000000735}{adaptive Lasso}, and \code{"taLasso"} represents \href{https://www.sciencedirect.com/science/article/abs/pii/S030440762100049X}{twin adaptive Lasso}, which is literally post-selection adaptive Lasso.
#' @param family a character string that specifies the quantitative family of response variable. Available families include \code{"gaussian"}, \code{"binomial"}, \code{"poisson"}, \code{"multinomial"}, \code{"cox"} and \code{"mgaussian"}. For more explanation of \code{family}, see function
#' \code{\link[stats]{glm}} in package \pkg{stats}. Default is \code{"gaussian"}.
#' @param nfolds an integer that specifies the number of folds in cross validation. It should be between 3 and the number of observations. Cross validation is needed to select the terms of sieves that give the best out-of-sample approximation error in Lasso-type methods. Default is \code{10}.
#' @param scale a non-negative number that regulates the weight given to the L1 penalization in adaptive Lasso (also twin adaptive Lasso). According to equation (4) in \href{https://www.tandfonline.com/doi/abs/10.1198/016214506000000735}{Zou (2006)}, the weight is \eqn{\frac{1}{|\hat{\beta}|^{\mathrm{scale}}}}. Hence, if \code{scale=0}, the weight for all \eqn{\beta} is 1,
#' which reduces adaptive Lasso to Lasso. If \code{scale=1}, the weight for higher \eqn{\beta} is smaller, and this is the key insight of adaptive Lasso. Thus, default value of \code{scale} is \code{1}.
#' @param parallel a logical value that indicates if a user hopes to use parallel \code{foreach} to fit each fold. If \code{TRUE}, the user must register parallel beforehand, such as \code{doMC}. See example below. Default is \code{FALSE}.
#' @param lambda.min.ratio a ratio that specifies the smallest value for \code{lambda} (i.e., penalty parameter in Lasso-type methods) as a fraction of maximal \code{lambda}.
#' @param s a character string that specifies the selection criterion of \code{lambda} in Lasso-type regularization. Available \code{"s"} includes \code{"lambda.min"} and \code{"lambda.1se"}. The former chooses the \code{lambda} that minimizes the cross-validated (i.e., out-of-sample) approximation error,
#' while the latter selects the largest value of \code{lambda} such that the cross-validated approximation error is within 1 standard error of the minimum. Default is \code{"lambda.min"}.
#' @param save.sieve a logical value that indicates if a user hopes to save the data frame of sieves in the output. Default is \code{FALSE}.
#'
#' @return The \code{dimada} command gives out an object of S3 class \code{"dimada"}, which is a list of results from all methods defined in the argument \code{methods} above. The result generally consists of the following items.
#' \itemize{
#'   \item \code{sieve}: a data frame containing the response variable and the sieves that are either generated from \code{x} or directly given by \code{x.sieve}. If the argument \code{save.sieve} is \code{FALSE}, then the data frame of sieves will not be included in the output.
#'   \item \code{method}: a character string that denotes the method, i.e., "Lasso", "adaLasso" or "taLasso".
#'   \item \code{parameters}: a data frame containing all lambdas (i.e., penalty parameters) that are used in above \code{method} and the consequent cross-validated approximation errors (i.e., mean squared errors) as well as corresponding numbers of non zero coefficients of terms in the sieve.
#'   \item \code{parameters.final}: a data frame containing the selected lambda and the corresponding cross-validated approximation error and number of non-zero coefficients. The selection rule is defined in argument \code{s}.
#'   \item \code{coefs.final}: a data frame containing the selected terms in the sieve and their coefficients.
#'   \item \code{lasso.all}, \code{adaLasso.all} or \code{taLasso.all}: an object of S3 class \code{"cv.glmnet"} containing all Lasso-type models estimated with a full list of lambdas as presented in \code{parameters}.
#'   \item \code{elapse}: time used to estimate all models with the above \code{method}. It is noted that "adaLasso" is based on the result of "Lasso", and "taLasso" is based on the result of "adaLasso", so time used for the methods of "Lasso", "adaLasso" and "taLasso" increases accordingly.
#'   \item \code{post.lm}: an object of S3 class \code{"lm"} containing the OLS model with all selected terms from Lasso-type methods.
#' }
#' Attributes are returned that correspond to the arguments \code{basis}, \code{n.basis}, \code{max.interaction}, \code{methods} and \code{s}. Also, if \code{x.sieve} is not given directly, the time used to compute sieves for \code{x} is then included as \code{sieve.time} in attributes. However, if \code{x.sieve} is given,
#' the attributes of \code{basis}, \code{n.basis}, \code{max.interaction} and \code{sieve.time} are \code{NA}.
#' @seealso  \link{dimada.plot}; \link{dimada.summary}; \link{poly.gen}; \link{bspline.gen}; \link{cosine.gen}; \link{sine.gen}; \link{trig.gen}; \link{haar.gen}; \link{daubechies.gen}.
#' @export
#'
#' @note The function \code{\link[glmnet]{cv.glmnet}} in package \pkg{glmnet} is applied to conduct Lasso-type methods in \code{dimada}. Cross-validation is demanded because out-of-sample approximation errors are needed to select the best model, as required in the dimension adaptive estimator.
#' Moreover, post-selection OLS model is performed to obtain unbiased coefficients because the estimated coefficients from Lasso-type methods are biased due to the penalty term in the loss function. So, the coefficients from post-selection OLS model could be used for prediction.
#'
#' @author Chencheng Fang, Bonn Graduate School of Economics, University of Bonn. Email: \email{ccfang@uni-bonn.de}
#'
#' @examples
#' # a data frame with 3 independent variables
#' set.seed(200)
#' df <- data.frame(a=runif(100), b=runif(100), c=runif(100))
#' # simulate a response variable with a parametric model with coefficients of c(1,2,1)
#' # and a normal error term
#' response1 <- as.matrix(df)%*%c(1,2,1)+rnorm(100,0.04)
#' # dimension adaptive estimation with power series and full interactions among independent
#' # variables as well as methods of Lasso and adaptive Lasso
#' dimada1 <- dimada(y=as.vector(response1), x=df, basis="power", max.interaction=3,
#'                   methods=c("Lasso","adaLasso"))
#'
#' # simulate a response variable with a non-parametric additive model and a normal error term
#' set.seed(200)
#' response2 <- sin(2*df[,1])+tan(df[,2])+log(df[,3])+rnorm(100,0.04)
#' # dimension adaptive estimation with B-splines and no interactions among variables
#' dimada2 <- dimada(y=response2, x=df, basis="bspline", max.interaction=1,
#'                   methods=c("Lasso","adaLasso"))
#'
#' # the same as above, but compute sieve with function 'bspline.gen' at first and then use
#' # the output as an input of 'x.sieve' in function 'dimada' directly
#' sieve <- bspline.gen(data=df, max.interaction=1)$train
#' dimada3 <- dimada(y=response2, x=NULL, x.sieve=sieve, methods=c("Lasso","adaLasso"))
#'
#' # Parallel
#' # install.packages("doMC")
#' require(doMC)
#' registerDoMC(cores = 2)
#' system.time(dimada(y=response2, x=df, basis="bspline", max.interaction=1,
#'             methods=c("Lasso","adaLasso")))
#' system.time(dimada(y=response2, x=df, basis="bspline", max.interaction=1,
#'             methods=c("Lasso","adaLasso"), parallel=TRUE))
dimada <- function(y,
                   x,
                   basis="power",
                   n.basis=NULL,
                   max.interaction=1,
                   x.sieve=NULL,
                   methods="Lasso",
                   family="gaussian",
                   nfolds=10,
                   scale=1,
                   parallel=FALSE,
                   lambda.min.ratio=1e-10,
                   s="lambda.min",
                   save.sieve=FALSE)
  {

  # --------------------------------
  # checking arguments
  # --------------------------------

  # check if 'x' is defined correctly
  if (all(!is.null(x), !is.matrix(x), !is.data.frame(x))) stop("The argument 'x' should be a matrix, a data frame or NULL")

  # check if 'y' is a vector, the length of which is the same as the number of rows in 'x'
  if (!is.vector(y)) stop("The argument 'y' should be a vector.")
  if (!is.null(x) && is.vector(y) && base::length(y)!= nrow(x)) stop("The length of 'y' should be the same as the number of observations in 'x'.")

  # check if 'basis' is defined correctly
  if (!is.null(x) && any(!is.character(basis), base::length(basis)!=1, !(basis %in% c("power","legendre","bspline","haar","daubechies","cosine", "sine", "trig")))) stop("Elements in the argument 'basis' have to be 'power','legendre', 'bspline', 'haar', 'daubechies', 'cosine', 'sine' or 'trig'.")

  # check if 'n.basis' is NULL or a single non-negative integer
  if (!is.null(x) && all(!is.null(n.basis), (n.basis<0 | n.basis%%1!=0 | base::length(n.basis) !=1))) stop("The argument 'n.basis' has to be either NULL or a single non-negative integer.")

  # check if 'max.interaction' is a single positive integer and smaller than the number of regressors
  if (!is.null(x) && (max.interaction<=0 | max.interaction>ncol(x) | max.interaction%%1!=0 | base::length(max.interaction) !=1)) stop("The argument 'max.interaction' has to be a single positive integer and smaller than the number of regressors.")

  # check if 'x.sieve' is defined correctly
  if (is.null(x) + is.null(x.sieve) !=1) stop("If 'x.sieve' is given, then the arguments 'x', 'basis', 'n.basis' and 'max.interaction' are not needed and 'x' must be NULL. Otherwise, if 'x.sieve' is NULL, then those arguments including 'x' must be given.")
  if (all(!is.null(x.sieve), !is.matrix(x.sieve), !is.data.frame(x.sieve))) stop("The argument 'x.sieve' has to be either NULL or a matrix / a data frame, but not a list.")
  if (!is.null(x.sieve) && base::length(y)!= base::nrow(x.sieve)) stop("If 'x.sieve' is a matrix, the length of 'y' should be the same as the number of rows in 'x.sieve'.")

  # check if 'methods' is defined correctly
  if (!("Lasso" %in% methods)) stop("The baseline method 'Lasso' has to be in the argument 'methods'.")
  if (any(!(methods %in% c("Lasso","adaLasso","taLasso")))) stop("Elements in the argument 'methods' have to be among 'Lasso','adaLasso' and 'taLasso'.")

  # check if 'family' is defined correctly
  if (any(base::length(family)!=1, !(family %in% c("gaussian", "binomial", "poisson", "multinomial", "cox", "mgaussian")))) stop("Elements in the argument 'family' have to be 'gaussian', 'binomial', 'poisson', 'multinomial', 'cox' or 'mgaussian'.")

  # check if 'nfolds' is a single integer between 3 and 'n'
  if (any(nfolds%%1!=0, nfolds>base::length(y), nfolds<3, base::length(nfolds)!=1)) stop("The argument 'nfolds' has to be a single integer between 3 and the number of observations.")

  # check if 'scale' is a single non-negative number
  if (any(!is.numeric(scale), scale <0, base::length(scale)!=1)) stop("The argument 'scale' has to be a single non-negative number.")

  # check if 'parallel' is a single logical value
  if (!is.logical(parallel) | base::length(parallel)!=1) stop("The argument 'parallel' has to be a single logical value.")

  # check if 'lambda.min.ratio' is defined correctly
  if (any(lambda.min.ratio>1, lambda.min.ratio<0, base::length(lambda.min.ratio)!=1)) stop("The argument 'lambda.min.ratio' has to be a single ratio number.")

  # check if 's' is defined correctly
  if (any(base::length(s)!=1, !(s %in% c("lambda.min","lambda.1se")))) stop("The argument 's' has to be 'lambda.min' or 'lambda.1se'.")

  # check if 'save.sieve' is a single logical value
  if (!is.logical(save.sieve) | base::length(save.sieve)!=1) stop("The argument 'save.sieve' has to be a single logical value.")

  # --------------------------------
  # data preparation
  # --------------------------------

  basis.input <- basis

  if (is.matrix(x.sieve) | is.data.frame(x.sieve)) {
    terms.values <- as.data.frame(x.sieve)
    basis <- NA
    n.basis <- NA
    max.interaction <- NA
  }

  if (is.matrix(x) | is.data.frame(x)) sieve.start <- base::Sys.time()
  if ((is.matrix(x) | is.data.frame(x)) & basis=="power") {terms.values.list <- poly.gen(data=x, test.data=NULL, n.basis=n.basis, max.interaction=max.interaction, legendre=FALSE); terms.values <- terms.values.list[[1]]}
  if ((is.matrix(x) | is.data.frame(x)) & basis=="legendre") {terms.values.list <- poly.gen(data=x, test.data=NULL, n.basis=n.basis, max.interaction=max.interaction, legendre=TRUE); terms.values <- terms.values.list[[1]]}
  if ((is.matrix(x) | is.data.frame(x)) & basis=="bspline") {terms.values.list <- bspline.gen(data=x, test.data=NULL, n.basis=n.basis, max.interaction=max.interaction, spline.degree=2); terms.values <- terms.values.list[[1]]}
  if ((is.matrix(x) | is.data.frame(x)) & basis=="haar") {terms.values.list <- haar.gen(data=x, test.data=NULL, n.basis=n.basis, max.interaction=max.interaction); terms.values <- terms.values.list[[1]]}
  if ((is.matrix(x) | is.data.frame(x)) & basis=="daubechies") {terms.values.list <- daubechies.gen(data=x, test.data=NULL, n.basis=n.basis, max.interaction=max.interaction); terms.values <- terms.values.list[[1]]}
  if ((is.matrix(x) | is.data.frame(x)) & basis=="cosine") {terms.values.list <- cosine.gen(data=x, test.data=NULL, n.basis=n.basis, max.interaction=max.interaction); terms.values <- terms.values.list[[1]]}
  if ((is.matrix(x) | is.data.frame(x)) & basis=="sine") {terms.values.list <- sine.gen(data=x, test.data=NULL, n.basis=n.basis, max.interaction=max.interaction); terms.values <- terms.values.list[[1]]}
  if ((is.matrix(x) | is.data.frame(x)) & basis=="trig") {terms.values.list <- trig.gen(data=x, test.data=NULL, n.basis=n.basis, max.interaction=max.interaction); terms.values <- terms.values.list[[1]]}
  if (is.matrix(x) | is.data.frame(x)) sieve.end <- base::Sys.time()

  response <- y

  # --------------------------------
  # LASSO
  # --------------------------------

  lasso.start <- base::Sys.time()
  lasso <- glmnet::cv.glmnet(x=as.matrix(terms.values), y=response, family=family, alpha =1, nfolds = nfolds, standardize=TRUE, parallel=parallel, relax = FALSE, nlambda=500, lambda.min.ratio=lambda.min.ratio, type.measure = "mse", intercept=TRUE)
  lasso.end <- base::Sys.time()
  lasso.lambda <- lasso$lambda
  lasso.final.coefs <- as.vector(as.matrix(lasso$glmnet.fit$beta, drop=FALSE)[,which(lasso$lambda==lasso[[s]])[1]]) #adding [1] to avoid repetition
  # kick out the outliers in coefficients
  if (basis.input=="bspline") lasso.final.coefs <- ifelse(lasso.final.coefs< stats::quantile(lasso.final.coefs[lasso.final.coefs!=0], prob=0.25, na.rm = TRUE)-20*stats::IQR(lasso.final.coefs[lasso.final.coefs!=0], na.rm = TRUE) | lasso.final.coefs> stats::quantile(lasso.final.coefs[lasso.final.coefs!=0], prob=0.75, na.rm = TRUE)+20*stats::IQR(lasso.final.coefs[lasso.final.coefs!=0], na.rm = TRUE),0, lasso.final.coefs)
  lasso.final.coefs.index <- lasso.final.coefs!=0

  # post lasso with ols
  # the term of intercept is already included in sieve. If it is important, it would remain in the selected term, so there is no need to have an additional intercept in ols.
  # it would be better to use glm because the family may not be 'gaussian'. This remains to be revised in future version.
  # post.lasso <- stats::glm(response ~ ., data=cbind(response=response, terms.values[,colnames(terms.values)[lasso.final.coefs.index],drop=FALSE]), family = base::eval(base::parse(text=family)))
  post.lasso <- stats::lm(response ~ ., data=cbind(response=response, terms.values[,colnames(terms.values)[lasso.final.coefs.index],drop=FALSE]))

  # --------------------------------
  # Adaptive LASSO
  # --------------------------------

  if (any(c("adaLasso","taLasso") %in% methods) & sum(lasso.final.coefs!=0)>=2) {

    # terms that are left for adaptive LASSO
    weight.adaLasso<- abs(lasso.final.coefs)^(-scale)
    index.w.adaLasso <- weight.adaLasso!=Inf
    weight.adaLasso <- weight.adaLasso[index.w.adaLasso]
    adaLasso.terms.values <- terms.values[,index.w.adaLasso]

    # Adaptive LASSO estimation
    adaLasso.start <- base::Sys.time()
    adaLasso <- glmnet::cv.glmnet(x=as.matrix(adaLasso.terms.values), y=response, family=family, alpha =1, nfolds = nfolds, penalty.factor=weight.adaLasso, standardize=TRUE, parallel=parallel, relax = FALSE, nlambda=500, lambda.min.ratio=lambda.min.ratio, type.measure = "mse", intercept=TRUE)
    adaLasso.end <- base::Sys.time()
    adaLasso.lambda <-  adaLasso$lambda
    adaLasso.final.coefs <- as.vector(as.matrix(adaLasso$glmnet.fit$beta, drop=FALSE)[,which(adaLasso$lambda==adaLasso[[s]])[1]])
    # kick out outliers in coefficients
    if (basis.input=="bspline") adaLasso.final.coefs <- ifelse(adaLasso.final.coefs< stats::quantile(adaLasso.final.coefs[adaLasso.final.coefs!=0], prob=0.25,na.rm = TRUE)-20*stats::IQR(adaLasso.final.coefs[adaLasso.final.coefs!=0], na.rm = TRUE) | adaLasso.final.coefs>stats::quantile(adaLasso.final.coefs[adaLasso.final.coefs!=0], prob=0.75,na.rm = TRUE)+20*stats::IQR(adaLasso.final.coefs[adaLasso.final.coefs!=0], na.rm = TRUE),0,adaLasso.final.coefs)
    adaLasso.final.coefs.index <- adaLasso.final.coefs!=0

    # post adaptive lasso
    # post.adaLasso <- stats::glm(response ~ ., data=cbind(response=response, adaLasso.terms.values[,colnames(adaLasso.terms.values)[adaLasso.final.coefs.index],drop=FALSE]), family = base::eval(base::parse(text=family)))
    post.adaLasso <- stats::lm(response ~ ., data=cbind(response=response, adaLasso.terms.values[,colnames(adaLasso.terms.values)[adaLasso.final.coefs.index],drop=FALSE]))

    # --------------------------------
    # Twin Adaptive LASSO
    # --------------------------------

    if (("taLasso" %in% methods) &(sum(adaLasso.final.coefs!=0)>=2)) {

      # terms that are left for the first step of Twin Adaptive LASSO
      index.taLasso <- adaLasso.final.coefs!=0
      taLasso.terms.values <- adaLasso.terms.values[,index.taLasso]

      # the first step of Twin Adaptive LASSO
      taLasso.s1.start <- base::Sys.time()
      taLasso.s1 <- glmnet::cv.glmnet(x=as.matrix(taLasso.terms.values), y=response, family=family, alpha =1, nfolds = nfolds, standardize=TRUE, parallel=parallel, relax = FALSE, nlambda=500, lambda.min.ratio=lambda.min.ratio, type.measure = "mse", intercept=TRUE)
      taLasso.s1.end <- base::Sys.time()
      taLasso.s1.lambda <-  taLasso.s1$lambda
      taLasso.s1.final.coefs <- as.vector(as.matrix(taLasso.s1$glmnet.fit$beta, drop=FALSE)[,which(taLasso.s1$lambda==taLasso.s1[[s]])[1]])
      # kick out outliers in coefficients
      if (basis.input=="bspline") taLasso.s1.final.coefs <- ifelse(taLasso.s1.final.coefs< stats::quantile(taLasso.s1.final.coefs[taLasso.s1.final.coefs!=0], prob=0.25, na.rm = TRUE)-20*stats::IQR(taLasso.s1.final.coefs[taLasso.s1.final.coefs!=0], na.rm = TRUE) | taLasso.s1.final.coefs>stats::quantile(taLasso.s1.final.coefs[taLasso.s1.final.coefs!=0], prob=0.75, na.rm = TRUE)+20*stats::IQR(taLasso.s1.final.coefs[taLasso.s1.final.coefs!=0], na.rm = TRUE),0,taLasso.s1.final.coefs)
      taLasso.s1.final.coefs.index <- taLasso.s1.final.coefs!=0

      if (sum(taLasso.s1.final.coefs!=0)>=2) {

        # terms that are left for the second step of Twin Adaptive LASSO
        weight.taLasso<- abs(taLasso.s1.final.coefs)^(-scale)
        index.w.taLasso <- weight.taLasso!=Inf
        weight.taLasso <- weight.taLasso[index.w.taLasso]
        taLasso.terms.values <- taLasso.terms.values[,index.w.taLasso]

        # the second step of Twin Adaptive LASSO
        taLasso.start <- base::Sys.time()
        taLasso <- glmnet::cv.glmnet(x=as.matrix(taLasso.terms.values), y=response, family=family, alpha =1, nfolds = nfolds, penalty.factor=weight.taLasso, standardize=TRUE, parallel=parallel, relax = FALSE,  nlambda=500, lambda.min.ratio=lambda.min.ratio, type.measure = "mse", intercept=TRUE)
        taLasso.end <- base::Sys.time()
        taLasso.lambda <-  taLasso$lambda
        taLasso.final.coefs <- as.vector(as.matrix(taLasso$glmnet.fit$beta, drop=FALSE)[,which(taLasso$lambda==taLasso[[s]])[1]])
        # kick out outliers in coefficients
        if (basis.input=="bspline") taLasso.final.coefs <- ifelse(taLasso.final.coefs< stats::quantile(taLasso.final.coefs[taLasso.final.coefs!=0], prob=0.25, na.rm = TRUE)-20*stats::IQR(taLasso.final.coefs[taLasso.final.coefs!=0], na.rm = TRUE) | taLasso.final.coefs>stats::quantile(taLasso.final.coefs[taLasso.final.coefs!=0], prob=0.75, na.rm = TRUE)+20*stats::IQR(taLasso.final.coefs[taLasso.final.coefs!=0], na.rm = TRUE),0,taLasso.final.coefs)
        taLasso.final.coefs.index <- taLasso.final.coefs!=0

        # post twin adaptive lasso
        # post.taLasso <- stats::glm(response ~ ., data=cbind(response=response, taLasso.terms.values[,colnames(taLasso.terms.values)[taLasso.final.coefs.index],drop=FALSE]), family = base::eval(base::parse(text=family)))
        post.taLasso <- stats::lm(response ~ ., data=cbind(response=response, taLasso.terms.values[,colnames(taLasso.terms.values)[taLasso.final.coefs.index],drop=FALSE]))

      }
    }
  }

  # --------------------------------
  # Output
  # --------------------------------

  if (save.sieve) {
    output <- list(Lasso=list(sieve=as.data.frame(cbind(response, terms.values)),
                              method="LASSO",
                              parameters=data.frame(cvms=lasso$cvm, lambdas=lasso.lambda, nzeros=lasso$nzero),
                              parameters.final=data.frame(cvms=lasso$cvm[which(lasso$lambda==lasso[[s]])[1]],lambdas=lasso[[s]],nzeros=lasso$nzero[which(lasso$lambda==lasso[[s]])[1]]),
                              coefs.final=data.frame(terms=colnames(terms.values)[lasso.final.coefs.index], coefs=lasso.final.coefs[lasso.final.coefs.index]),
                              lasso.all=lasso,
                              elapse=lasso.end-lasso.start,
                              post.lm=post.lasso
    ))

    if (("adaLasso" %in% methods) & (sum(lasso.final.coefs!=0)>=2)) {
      output <- append(output, list(adaLasso=list(sieve=as.data.frame(cbind(response, adaLasso.terms.values)),
                                                  method="Adaptive LASSO",
                                                  parameters=data.frame(cvms=adaLasso$cvm, lambdas=adaLasso.lambda, nzeros=adaLasso$nzero),
                                                  parameters.final=data.frame(cvms=adaLasso$cvm[which(adaLasso$lambda==adaLasso[[s]])[1]], lambdas=adaLasso[[s]], nzeros=adaLasso$nzero[which(adaLasso$lambda==adaLasso[[s]])[1]]),
                                                  terms.index = index.w.adaLasso,
                                                  coefs.final=data.frame(terms=colnames(adaLasso.terms.values)[adaLasso.final.coefs.index], coefs=adaLasso.final.coefs[adaLasso.final.coefs.index]),
                                                  adaLasso.all=adaLasso,
                                                  elapse=(adaLasso.end-adaLasso.start)+(lasso.end-lasso.start),
                                                  post.lm=post.adaLasso
      )))
      if (("taLasso" %in% methods) & (sum(adaLasso.final.coefs!=0)>=2)) {
        if(sum(taLasso.s1.final.coefs!=0)>=2) {
          output <- append(output, list(taLasso=list(sieve=as.data.frame(cbind(response, taLasso.terms.values)),
                                                     method="Twin Adaptive LASSO",
                                                     parameters=data.frame(cvms=taLasso$cvm, lambdas=taLasso.lambda, nzeros=taLasso$nzero),
                                                     parameters.final=data.frame(cvms=taLasso$cvm[which(taLasso$lambda==taLasso[[s]])[1]], lambdas=taLasso[[s]], nzeros=taLasso$nzero[which(taLasso$lambda==taLasso[[s]])[1]]),
                                                     terms.index=index.w.adaLasso[index.taLasso][index.w.taLasso],
                                                     coefs.final=data.frame(terms=colnames(taLasso.terms.values)[taLasso.final.coefs.index], coefs=taLasso.final.coefs[taLasso.final.coefs.index]),
                                                     taLasso.all=taLasso,
                                                     elapse=(taLasso.s1.end-taLasso.s1.start)+(taLasso.end-taLasso.start)+(adaLasso.end-adaLasso.start)+(lasso.end-lasso.start),
                                                     post.lm=post.taLasso
          )))
        }
      }
    }
  } else {
    output <- list(Lasso=list(method="LASSO",
                              parameters=data.frame(cvms=lasso$cvm, lambdas=lasso.lambda, nzeros=lasso$nzero),
                              parameters.final=data.frame(cvms=lasso$cvm[which(lasso$lambda==lasso[[s]])[1]],lambdas=lasso[[s]],nzeros=lasso$nzero[which(lasso$lambda==lasso[[s]])[1]]),
                              coefs.final=data.frame(terms=colnames(terms.values)[lasso.final.coefs.index], coefs=lasso.final.coefs[lasso.final.coefs.index]),
                              lasso.all=lasso,
                              elapse=lasso.end-lasso.start,
                              post.lm=post.lasso
    ))

    if (("adaLasso" %in% methods) & (sum(lasso.final.coefs!=0)>=2)) {
      output <- append(output, list(adaLasso=list(method="Adaptive LASSO",
                                                  parameters=data.frame(cvms=adaLasso$cvm, lambdas=adaLasso.lambda, nzeros=adaLasso$nzero),
                                                  parameters.final=data.frame(cvms=adaLasso$cvm[which(adaLasso$lambda==adaLasso[[s]])[1]], lambdas=adaLasso[[s]], nzeros=adaLasso$nzero[which(adaLasso$lambda==adaLasso[[s]])[1]]),
                                                  terms.index = index.w.adaLasso,
                                                  coefs.final=data.frame(terms=colnames(adaLasso.terms.values)[adaLasso.final.coefs.index], coefs=adaLasso.final.coefs[adaLasso.final.coefs.index]),
                                                  adaLasso.all=adaLasso,
                                                  elapse=(adaLasso.end-adaLasso.start)+(lasso.end-lasso.start),
                                                  post.lm=post.adaLasso
      )))
      if (("taLasso" %in% methods) & (sum(adaLasso.final.coefs!=0)>=2)) {
        if(sum(taLasso.s1.final.coefs!=0)>=2) {
          output <- append(output, list(taLasso=list(method="Twin Adaptive LASSO",
                                                     parameters=data.frame(cvms=taLasso$cvm, lambdas=taLasso.lambda, nzeros=taLasso$nzero),
                                                     parameters.final=data.frame(cvms=taLasso$cvm[which(taLasso$lambda==taLasso[[s]])[1]], lambdas=taLasso[[s]], nzeros=taLasso$nzero[which(taLasso$lambda==taLasso[[s]])[1]]),
                                                     terms.index=index.w.adaLasso[index.taLasso][index.w.taLasso],
                                                     coefs.final=data.frame(terms=colnames(taLasso.terms.values)[taLasso.final.coefs.index], coefs=taLasso.final.coefs[taLasso.final.coefs.index]),
                                                     taLasso.all=taLasso,
                                                     elapse=(taLasso.s1.end-taLasso.s1.start)+(taLasso.end-taLasso.start)+(adaLasso.end-adaLasso.start)+(lasso.end-lasso.start),
                                                     post.lm=post.taLasso
          )))
        }
      }
    }
  }

  # --------------------------------
  # Attributes
  # --------------------------------

  attr(output, "basis") <- basis
  attr(output, "n.basis") <- ifelse(is.null(x), NA, attr(terms.values.list, "n.basis"))
  attr(output, "max.interaction") <- max.interaction
  attr(output, "methods") <- methods
  attr(output, "s") <- s
  attr(output, "sieve.time") <- ifelse(is.null(x), NA, base::format(sieve.end-sieve.start, units="auto"))

  class(output) <- "dimada"

  return(output)

}
