#' Empirical sorted partial effects (SPE) and inference
#'
#' \code{SPE} conducts SPE estimation and inference at user-specifed quantile index. The bootstrap procedures
#' follows algorithm 2.1 as in Chernozhukov, Fernandez-Val and Luo (2018). All estimates are bias-corrected
#' and all confidence bands are monotonized. For graphical results, please use \code{\link{SPEplot}}.
#'
#' @return The output is a list with 4 components: (1) \code{spe} stores spe estimates and confidence bounds;
#' (2) \code{ape} stores ape estimates and confidence bounds; (3) \code{us} stores percentile index as in
#' \code{SPE} command; (4) \code{alpha} stores significance level as in \code{SPE} command.
#'
#' @param fm          Regression formula.
#' @param data        Data in use.
#' @param method      Models to be used for estimating partial effects. Four options: \code{"logit"} (binary response),
#'                    \code{"probit"} (binary response), \code{"ols"} (interactive linear with additive errors), \code{"QR"} (linear model
#'                    with non-additive errors). Default is \code{"ols"}.
#' @param var.type    The type of parameter in interest. Three options: \code{"binary"}, \code{"categorical"}, \code{"continuous"}. Default
#'                    is \code{"binary"}.
#' @param var.T       Variable T in interset. Should be a character type.
#' @param compare     If parameter in interest is categorical, then user needs to specify which two category to
#'                    compare with. Should be a 1 by 2 character vector. For example, if the two levels to compare
#'                    with is 1 and 3, then \code{c=("1", "3")}, which will calculate partial effect from 1 to 3. To use
#'                    this option, users first need to specify var.T as a factor variable.
#' @param subgroup    Subgroup in interest. Default is \code{NULL}. Specifcation should be a logical variable. For example, suppose data contains
#'                    indicator variable for women (female if 1, male if 0). If users are interested in women SPE, then users
#'                    should specify \code{subgroup = data[, "female"] == 1}.
#' @param samp_weight Sampling weight of data. If null then function implements empirical bootstrap.
#'                    If data specifies sampling weight, the function implements weighted bootstrap. Input
#'                    should be a n by 1 vector, where n denotes sample size. Default is \code{NULL}.
#' @param us          Percentile of interest for SPE. Should be a vector of values between 0 and 1. Default
#'                    is \code{c(1:9)/10}.
#' @param alpha       Size for confidence interval. Shoule be between 0 and 1. Default is 0.1
#' @param taus        Indexes for quantile regression. Default is \code{c(1:9)/10}.
#' @param B           Number of bootstrap draws. Default is set to be 10. For more accurate results, we recommend 500.
#' @param ncores      Number of cores for computation. Default is set to be 1. For large dataset, parallel computing
#'                    is highly recommended since bootstrap is time-consuming.
#' @param seed        Pseudo-number generation for reproduction. Default is 1.
#' @param bc          Whether want the estimate to be bias-corrected. Default is \code{TRUE}. If \code{FALSE} uncorrected
#'                    estimate and corresponding confidence bands will be reported.
#' @param boot.type   Type of bootstrap. Default is \code{boot.type = "nonpar"}, and the package implements nonparametric
#'                    bootstrap. An alternative is \code{boot.type = "weighted"}, and the package implements weighted
#'                    bootstrap.
#' @examples
#' data("mortgage")
#' fm <- deny ~ black + p_irat + hse_inc + ccred + mcred + pubrec + ltv_med +
#' ltv_high + denpmi + selfemp + single + hischl
#' test <- SPE(fm = fm, data = mortgage, var.T = "black", method = "logit",
#' us = c(1:9)/10)
#'
#' @importFrom Hmisc wtd.quantile
#' @importFrom boot boot
#' @importFrom stats quantile rexp qnorm
#' @export
SPE <- function(fm, data, method = "ols", var.type = "binary", var.T, compare,
                subgroup = NULL, samp_weight = NULL, us = c(1:9)/10, alpha = 0.1,
                taus = c(1:9)/10, B = 10, ncores = 1, seed = 1, bc = TRUE,
                boot.type = "nonpar"){
  # --------------------------- Stopping Condition -------------------------------------
  if(alpha >= 1 || alpha <= 0) stop("Please specify a correct size for hypothesis testing.")
  if(boot.type != "nonpar" && boot.type != "weighted") stop("Please specify a bootstrap type as a char: nonpar or weighted.")
  # --------------------------- Replace Null samp_weight Specification ---------------------
  if(is.null(samp_weight)) samp_weight <- rep(1, dim(data)[1])
  # --------------------------- 1. Call to estimate PE ---------------------------
  output <- PEestimate(fm, data, samp_weight, var.type, var.T, compare, method, subgroup, taus)
  PE.est <- output$PE.est
  # --------------------------- 2. Now get estimated SPE (full and subgroup samples)---------------------------
  if(method != "QR"){
    SPE.est <- wtd.quantile(PE.est, samp_weight, us)
  } else {
    SPE.est <- wtd.quantile(PE.est, matrix(samp_weight, ncol = 1, nrow = nrow(PE.est), byrow = FALSE), us)
  }
  if(!is.null(subgroup)){
    PEsub.est <- output$PEsub.est
    PEsub.w <- output$samp_weight_sub
    SPEsub.est <- wtd.quantile(PEsub.est, PEsub.w, us)
  }
  # --------------------------- 3. Bootstrap Samples ---------------------------
  #  statistic in one bootstrap for SPE
  #  Need to consider two scenarios: whether the regression is weighted. If not then use nonpar boot, if yes use
  #  exchangeably weighted boot. Turns out the statistic functions for boot are a bit different for these two
  #  scenarios, so I write two stastistics functions.
  stat.boot.noweight <- function(data, indices){
    data$.w <- samp_weight
    # create bootstrap sample, the "indices" is automatically provided by the boot package, this provides the
    # sampling with replacement part of bootstrap
    data.bs <- data[indices, ]
    output.bs <- PEestimate(fm, data.bs, samp_weight=data.bs$.w, var.type, var.T, compare, method, subgroup, taus)
    est.PE.bs <- output.bs$PE.est
    est.APE.bs <- mean(est.PE.bs) # scalar, length = 1
    est.SPE.bs <- wtd.quantile(est.PE.bs, data.bs$.w, us) # length = length(us)
    if(!is.null(subgroup)){
      est.PEsub.bs <- output.bs$PEsub.est
      PEsub.w.bs <- output.bs$samp_weight_sub
      est.APEsub.bs <- mean(est.PEsub.bs) # length = 1
      est.SPEsub.bs <- wtd.quantile(est.PEsub.bs, PEsub.w.bs, us) # length = length(us)
      return(c(est.APEsub.bs, est.SPEsub.bs)) # Each boot stat has length: 1+length(us)
    }else{
      return(c(est.APE.bs, est.SPE.bs)) # Each boot stat has length: 1+length(us)
    }
    data.bs$.w <- NULL
  }
  # The resampling is wrt sample weight, so the data.rg function does that
  data.rg <- function(data, mle){
    n <- dim(data)[1]
    # Exponential weights
    multipliers  <- rexp(n)
    # Sampling weight of data.bs
    weight <- (multipliers/sum(multipliers)) * samp_weight
    data$.w <- weight
    return(data)
  }
  # The following function implements nonparametric bootstrap for quantile regression
  data.non <- function(data, mle){
    n <- dim(data)[1]
    multipliers <- as.vector(table(factor(sample(n,n,replace=T), levels = c(1:n))))
    # Sampling weight of data.bs
    weight <- (multipliers/sum(multipliers)) * samp_weight
    data$.w <- weight
    return(data)
  }
  # The following is boot stat function for exchageably weighted boot
  stat.boot.weight <- function(data){
    output.bs <- PEestimate(fm, data, samp_weight=data$.w, var.type, var.T, compare, method, subgroup, taus)
    est.PE.bs <- output.bs$PE.est
    est.APE.bs <- mean(est.PE.bs) # scalar, length = 1
    if (method != "QR"){
      est.SPE.bs <- wtd.quantile(est.PE.bs, samp_weight, us)
    } else {
      est.SPE.bs <- wtd.quantile(est.PE.bs, matrix(samp_weight, ncol = 1, nrow = nrow(PE.est), byrow = FALSE), us)
    } # length = length(us)
    if(!is.null(subgroup)){
      est.PEsub.bs <- output.bs$PEsub.est
      PEsub.w.bs <- output.bs$samp_weight_sub
      est.APEsub.bs <- mean(est.PEsub.bs) # length = 1
      est.SPEsub.bs <- wtd.quantile(est.PEsub.bs, PEsub.w.bs, us) # length = length(us)
      return(c(est.APEsub.bs, est.SPEsub.bs)) # Each boot stat has length: 1+length(us)
    }else{
      return(c(est.APE.bs, est.SPE.bs)) # Each boot stat has length: 1+length(us)
    }
  }
  # Use boot command
  set.seed(seed)
  if(boot.type == "nonpar"){
    if(method != "QR"){
      result.boot <- boot(data = data, statistic = stat.boot.noweight, parallel="multicore", ncpus = ncores, R = B)
    } else{
      data$.w <- samp_weight
      result.boot <- boot(data = data, statistic = stat.boot.weight, sim = "parametric", ran.gen = data.non, mle = 0, parallel = "multicore", ncpus = ncores, R = B)
      data$.w <- NULL
    }
  } else if(boot.type == "weighted"){
    data$.w <- samp_weight
    result.boot <- boot(data = data, statistic = stat.boot.weight, sim = "parametric", ran.gen = data.rg, mle = 0, parallel = "multicore", ncpus = ncores, R = B)
    data$.w <- NULL
  }
  # --------------------------- 4. Inference ---------------------------
  ##############################
  ### Full Sample
  ##############################
  # SPE inference
  draws.SPE.bs <- result.boot$t[, 2:(length(us)+1)]
  # bundle function implements algorithm 2.1, remark 2.2 and 2.3
  inf.SPE <- bundle(draws.SPE.bs, SPE.est, alpha)
  # APE inference
  est.APE <- mean(PE.est)
  draws.APE.bs <- result.boot$t[, 1]
  inf.APE <- ape(draws.APE.bs, est.APE, alpha)

  ##############################
  ### Subgroup Sample (if subgroup is not NULL)
  ##############################
  if(!is.null(subgroup)){
    # SPE
    draws.SPEsub.bs <- result.boot$t[, 2:(length(us)+1)]
    inf.SPEsub <- bundle(draws.SPEsub.bs, SPEsub.est, alpha)
    # APE
    est.APEsub <- mean(PEsub.est)
    draws.APEsub.bs <- result.boot$t[, 1]
    inf.APEsub <- ape(draws.APEsub.bs, est.APEsub, alpha)
  }
  # --------------------------- 5. Return Results ---------------------------
  # Depends on whether subgroup is NULL and whether bias correction is wanted
  # The resulting outputs are APE(sub) & SPE(sub) with correpsonding confidence intervals. All are stored in an "inf" bundle list
  if(!is.null(subgroup)){
    if(bc == TRUE) {
      output <- list(spe = inf.SPEsub[1:3], ape = inf.APEsub[1:3], us = us, alpha = alpha)
    } else {
      output <- list(spe = inf.SPEsub[4:6], ape = inf.APEsub[4:6], us = us, alpha = alpha)
    }
  } else {
    if(bc == TRUE) {
      output <- list(spe=inf.SPE[1:3], ape=inf.APE[1:3], us = us, alpha = alpha)
    } else {
      output <- list(spe=inf.SPE[4:6], ape=inf.APE[4:6], us = us, alpha = alpha)
    }
  }
}

# ---------------------------------- Auxiliary Functions ------------------------------------------
# Implementing algorithm 2.1 and get (bias corrected) estimate and confidence bands for sorted effects
bundle <- function(bs, est, alpha){
  z.bs <- bs - matrix(est, nrow = nrow(bs), ncol = ncol(bs), byrow = TRUE)
  sigma <- (apply(z.bs, 2, quantile, .75, na.rm=TRUE) - apply(z.bs, 2, quantile, .25, na.rm=TRUE))/(qnorm(0.75) - qnorm(.25))
  t.hat <- apply(abs(z.bs /matrix(sigma, nrow = nrow(bs), ncol = ncol(bs), byrow=TRUE)), 1, max)
  crt <- quantile(t.hat, 1-alpha)
  # bias-correction
  est.bc <- sort(2*est - apply(bs, 2, mean))
  ubound.est.bc <-  sort(est.bc + crt * sigma)
  lbound.est.bc <-  sort(est.bc - crt * sigma)
  # uncorrected
  ubound.est <-  sort(est + crt * sigma)
  lbound.est <-  sort(est - crt * sigma)
  out <- list(est.bc = est.bc, ubound.est.bc = ubound.est.bc, lbound.est.bc = lbound.est.bc,
              est = est, ubound.est = ubound.est, lbound.est = lbound.est)
}

# Get (biased corrected) average partial effects and confidence intervals
ape <- function(bs, est, alpha){
  sigma <- (quantile(bs, .75, na.rm = TRUE) - quantile(bs, .25, na.rm = TRUE))/(qnorm(0.75) - qnorm(.25))
  mzs <- abs(bs - est) / sigma
  crt2 <- quantile(mzs, 1-alpha)
  # bias-correction
  est.bc    <- 2 * est - mean(bs)
  ubound.APE.bc <- est.bc + crt2 * sigma
  lbound.APE.bc <- est.bc - crt2 * sigma
  # uncorrected
  ubound.APE <- est + crt2 * sigma
  lbound.APE <- est - crt2 * sigma
  out <- list(est.bc = est.bc, ubound.est.bc = ubound.APE.bc, lbound.est.bc = lbound.APE.bc,
              est = est, ubound.APE = ubound.APE, lbound.APE = lbound.APE)
}

# ----------------------- Plotting (Ingredients: average effect, sorted effect, correspondent confidence bands) ----------------
#' Plot output of \code{\link{SPE}} command.
#' @param output      Output of \code{\link{SPE}} command.
#' @param xlim        x-axis limits. Default is range of percentile index.
#' @param ylim        y-axis limits. Default is NULL.
#' @param main        Main title of the plot. Defualt is NULL.
#' @param sub         Sub title of the plot. Default is NULL.
#' @param xlab        x-axis label. Default is "Percentile Index".
#' @param ylab        y-axis label. Default is "Sorted Effects".
#'
#' @examples
#' data("mortgage")
#' fm <- deny ~ black + p_irat + hse_inc + ccred + mcred + pubrec + ltv_med +
#' ltv_high + denpmi + selfemp + single + hischl
#' test <- SPE(fm = fm, data = mortgage, var.T = "black", method = "logit",
#' us = c(1:9)/10)
#'
#' SPEplot(output = test, main="APE and SPE of Being Black on the prob of Mortgage Denial",
#' sub="Logit Model", ylab="Change in Probability")
#'
#' @importFrom graphics plot polygon lines abline points legend
#' @export
SPEplot <- function(output, xlim = NULL, ylim = NULL, main = NULL, sub = NULL,
                    xlab = "Percentile Index", ylab = "Sorted Effects"){
  # SPE and confidence bands
  SE <- output$spe
  AE <- output$ape
  us <- output$us
  alpha <- output$alpha
  xlim <- range(us)
  plot(us, SE[[1]], type="l", xlim, ylim, log = "", main, sub, xlab, ylab, col=4, lwd=2)
  polygon(c(us, rev(us)),c(SE[[2]], rev(SE[[3]])), density = 60, border = F,
          col='light blue', lty = 1, lwd = 1);
  lines(us, SE[[1]], lwd = 2, col = 4 )
  # APE and CI
  abline(h = AE[[1]], col = 1);
  abline(h = AE[[2]], col = 1, lty = 2);
  abline(h = AE[[3]], col = 1, lty = 2);
  points(c(min(us),.2,.4,.6,.8,max(us)), rep(AE[[1]], 6), col = 1, pch = 15)
  legend(x = "topleft", legend = c("SPE","APE", paste0((1-alpha)*100,"% CB(SPE)"), paste0((1-alpha)*100,"% CB(APE)")), col = c(4, 1, "light blue", 1),
         lwd = c(1, 1, 5, 1), lty = c(1, 1, 1, 2), pch = c(NA, 15, NA, NA), pt.cex = c(2, 1), bty = 'n')
}
