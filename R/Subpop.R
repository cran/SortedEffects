#' Inference on Most and Least Affected Groups
#'
#' \code{Subpop} conducts set inference on the groups of most and least affected. When \code{subgroup = NULL}, output is for whole sample.
#' Otherwise the results are subgroup. The results can be visualized using the
#' \code{\link{Subpopplot}} command. The output of \code{Subpop} is a list containing four components: \code{most}, \code{least}, \code{u} and \code{sub}.
#' As the names indicate, \code{most} and \code{least} denote the confidence sets for the most and least affected units. \code{u} stores the u-th most
#' and least affected index and \code{sub} stores the indicators for subpopulations.
#'
#' @param fm          Regression formula
#' @param data        The data in use
#' @param method      Models to be used for estimating partial effects. Four options: \code{"logit"} (binary response),
#'                    \code{"probit"} (binary response), \code{"ols"} (interactive linear with additive errors), \code{"QR"} (linear model
#'                    with non-additive errors). Default is \code{"ols"}.
#' @param var.type    The type of parameter in interest. Three options: \code{"binary"}, \code{"categorical"}, \code{"continuous"}. Default
#'                    is \code{"binary"}.
#' @param var.T       Variable T in interset. Should be a character.
#' @param compare     If parameter in interest is categorical, then user needs to specify which two category to
#'                    compare with. Should be a 1 by 2 character vector. For example, if the two levels to compare
#'                    with is 1 and 3, then \code{c=("1", "3")}, which will calculate partial effect from 1 to 3. To use
#'                    this option, users first need to specify var.T as a factor variable.
#' @param subgroup    Subgroup in interest. Default is \code{NULL}. Specifcation should be a logical variable. For example, suppose data contains
#'                    indicator variable for women (female if 1, male if 0). If users are interested in women SPE, then users
#'                    should specify \code{subgroup = data[, "female"] == 1}.
#' @param samp_weight Sampling weight of data. If null then function implements empirical bootstrap.
#'                    If data specifies sampling weight, put that in and the function implements
#'                    weighted (i.i.d exponential weights) bootstrap.
#' @param taus        Indexes for quantile regression. Default is \code{c(1:9)/10}.
#' @param u           Percentile of most and least affected. Default is set to be 0.1.
#' @param alpha       Size for confidence interval. Shoule be between 0 and 1. Default is 0.1
#' @param B           Number of bootstrap draws. Default is set to be 10. For more accurate results, we recommend 500.
#' @param ncores      Number of cores for computation. Default is set to be 1. For large dataset, parallel computing
#'                    is highly recommended since bootstrap is time-consuming.
#' @param seed        Pseudo-number generation for reproduction. Default is 1.
#' @param boot.type   Type of bootstrap. Default is \code{boot.type = "nonpar"}, and the package implements nonparametric
#'                    bootstrap. An alternative is \code{boot.type = "weighted"}, and the package implements weighted
#'                    bootstrap.
#'
#' @examples
#' data("mortgage")
#' fm <- deny ~ black + p_irat + hse_inc
#' result <- Subpop(fm = fm, data = mortgage, var.T = "black", method = "logit")
#'
#' @importFrom Hmisc wtd.quantile
#' @importFrom boot boot
#' @importFrom stats quantile rexp qnorm
#' @export
Subpop <- function(fm, data, method = "ols", var.type = "binary", var.T, compare,
                   subgroup = NULL, samp_weight = NULL, taus = c(1:9)/10, u = 0.1,
                   alpha = 0.1, B = 10, ncores = 1, seed = 1, boot.type = "nonpar"){
  #---------------------- Stopping Conditions ---------------------------
  if(alpha >= 1 || alpha <= 0) stop("Please specify a correct size for hypothesis testing.")
  if(u >= 1 || u <= 0) stop("Please provide a meaningful group classification.")
  if(boot.type != "nonpar" && boot.type != "weighted") stop("Please specify type of bootstrap as a char: nonpar or weighted.")
  # --------------------------- Replace Null samp_weight specification --------------------
  if(is.null(samp_weight)) samp_weight <- rep(1, dim(data)[1])
  # --------------------------- 1. Call to estimate PE and PEsub ---------------------------
  output <- PEestimate(fm, data, samp_weight, var.type, var.T, compare, method, subgroup, taus)
  PE_est <- output$PE.est

  # --------------------------- 2. Threshold Values for u-most/least Affected ---------------------------
  # Full sample
  if(method == "QR"){
    effect.high <- wtd.quantile(PE_est, matrix(samp_weight, ncol = 1, nrow = nrow(PE_est), byrow = FALSE), 1-u)
    effect.low <- wtd.quantile(PE_est, matrix(samp_weight, ncol = 1, nrow = nrow(PE_est), byrow = FALSE), u)
  } else {
    effect.high <- wtd.quantile(PE_est, samp_weight, 1-u)
    effect.low <- wtd.quantile(PE_est, samp_weight, u)
  }
  # Subgroup sample
  if(!is.null(subgroup)){
    PEsub_est <- output$PEsub.est
    PEsub_w <- output$samp_weight_sub
    effect_sub.high <- wtd.quantile(PEsub_est, PEsub_w, 1-u)
    effect_sub.low <- wtd.quantile(PEsub_est, PEsub_w, u)
  }
  # ---------------------------- 3. Bootstrap Samples --------------------------
  # 1. Statistics in one boot if no weight is specifed
  boot.stat.noweight <- function(data, indices){
    data$.w <- samp_weight
    data <- data[indices, ]
    out.bs <- PEestimate(fm, data, samp_weight = data$.w, var.type, var.T, compare, method, subgroup, taus)
    PE_est.bs <- out.bs$PE.est
    effect.high.bs <- wtd.quantile(PE_est.bs, data$.w, 1-u) # scalar
    effect.low.bs <- wtd.quantile(PE_est.bs, data$.w, u) # scalar
    if(!is.null(subgroup)){
      PEsub_est.bs <- out.bs$PEsub.est
      PEsub_w <- out.bs$samp_weight_sub
      effect_sub.high.bs <- wtd.quantile(PEsub_est.bs, PEsub_w, 1-u)
      effect_sub.low.bs <- wtd.quantile(PEsub_est.bs, PEsub_w, u)
      return(c(effect_sub.high.bs, effect_sub.low.bs, PEsub_est.bs))
    }else{
      return(c(effect.high.bs, effect.low.bs, PE_est.bs))
    }
  }
  # 2. Stats in one boot if weight is specified
  data.rg <- function(data, mle){
    n <- dim(data)[1]
    # Exponential weights
    multipliers  <- rexp(n)
    # Sampling weight of data
    weight <- samp_weight * multipliers / sum(multipliers) * 20000
    data$.w <- weight
    return(data)
  }

  # The following function implements nonparametric bootstrap for quantile regression
  data.non <- function(data, mle){
    n <- dim(data)[1]
    multipliers <- as.vector(table(factor(sample(n,n,replace=T), levels = c(1:n))))
    # Sampling weight of data.bs
    weight <- (multipliers/sum(multipliers)) * samp_weight * 20000
    data$.w <- weight
    return(data)
  }

  boot.stat.weight <- function(data){
    out.bs <- PEestimate(fm, data, samp_weight=data$.w, var.type, var.T, compare, method, subgroup, taus)
    PE_est.bs <- out.bs$PE.est
    if(method == "QR"){
      effect.high.bs <- wtd.quantile(PE_est.bs, matrix(data$.w, ncol = 1, nrow = nrow(PE_est), byrow = FALSE), 1-u) # scalar
      effect.low.bs <- wtd.quantile(PE_est.bs, matrix(data$.w, ncol = 1, nrow = nrow(PE_est), byrow = FALSE), u) # scalar
    } else {
      effect.high.bs <- wtd.quantile(PE_est.bs, data$.w, 1-u) # scalar
      effect.low.bs <- wtd.quantile(PE_est.bs, data$.w, u) # scalar
    }
    if(!is.null(subgroup)){
      PEsub_est.bs <- out.bs$PEsub.est
      PEsub_w.bs <- out.bs$samp_weight_sub
      effect_sub.high.bs <- wtd.quantile(PEsub_est.bs, PEsub_w.bs, 1-u)
      effect_sub.low.bs <- wtd.quantile(PEsub_est.bs, PEsub_w.bs, u)
      return(c(effect_sub.high.bs, effect_sub.low.bs, PEsub_est.bs))
    }else{
      return(c(effect.high.bs, effect.low.bs, PE_est.bs))
    }
  }
  # --------------------------- 4. Conduct Bootstrap ---------------------------------------
  set.seed(seed)
  if(boot.type == "nonpar"){
    if(method != "QR"){
      result.boot <- boot(data = data, statistic = boot.stat.noweight, parallel="multicore", ncpus = ncores, R = B)
    } else{
      data$.w <- samp_weight
      result.boot <- boot(data = data, statistic = boot.stat.weight, sim = "parametric", ran.gen = data.non, mle = 0, parallel = "multicore", ncpus = ncores, R = B)
      data$.w <- NULL
    }
  } else if(boot.type == "weighted"){
    data$.w <- samp_weight
    result.boot <- boot(data = data, statistic = boot.stat.weight, sim = "parametric", ran.gen = data.rg, mle = 0, parallel = "multicore", ncpus = ncores, R = B)
    data$.w <- NULL
  }
  # --------------------------- 5. Analysis for the subpopulation ---------------------------
  if(is.null(subgroup)){
    # (a) PE
    # MOST Affected Group confidence set
    is.he <- submost(PE_est, result.boot$t[, 3:(length(PE_est)+2)], result.boot$t[, 1], effect.high, alpha, B)
    # LEAST Affected Group confidence set
    is.le <- subleast(PE_est, result.boot$t[, 3:(length(PE_est)+2)], result.boot$t[, 2], effect.low, alpha, B)
    output <- list(most = is.he, least = is.le, u = u, sub = subgroup)
  } else{
    # (b) PEsub
    is.he_sub <- submost(PEsub_est, result.boot$t[, 3:(length(PEsub_est)+2)], result.boot$t[, 1], effect_sub.high, alpha, B)
    is.le_sub <- subleast(PEsub_est, result.boot$t[, 3:(length(PEsub_est)+2)], result.boot$t[, 2], effect_sub.low, alpha, B)
    output <- list(most = is.he_sub, least = is.le_sub, u = u, sub = subgroup)
  }
}

# --------------------- Two Auxiliary Functions -------------------------
# Implementing algorithm: Output is confidence set indicator
# Most affected group
submost <- function(est_pe, bs_pe, bs_u, est_u, alpha, B){
  is.0 <- rank(abs(est_pe - est_u)) == min(rank(abs(est_pe - est_u))) # Find the min (to implement the sup condition as in the paper)
  draws <- bs_pe - matrix(bs_u, nrow = B, ncol = length(est_pe)) - matrix(est_pe - est_u, nrow = B, ncol = length(est_pe), byrow = TRUE)
  bse <- (apply(draws, 2, quantile, .75, na.rm = TRUE) - apply(draws, 2, quantile, .25, na.rm = TRUE))/(qnorm(0.75) - qnorm(.25))
  bm <- apply(draws, 2, mean) # bias estimator
  zs <- apply(-draws[, is.0] / matrix(bse[is.0], nrow = B, ncol = length(bse[is.0]), byrow = TRUE), 1, max, na.rm = TRUE)
  crt <- quantile(zs, 1 - alpha)  #critical value
  is.he <- -(est_pe - est_u) / bse <= crt
  return(is.he)
}
# Least affected group
subleast <- function(est_pe, bs_pe, bs_u, est_u, alpha, B){
  is.0 <- rank(abs(est_pe - est_u)) == min(rank(abs(est_pe - est_u))) # Find the min (to implement the sup condition as in the paper)
  draws <- bs_pe - matrix(bs_u, nrow = B, ncol = length(est_pe)) - matrix(est_pe - est_u, nrow = B, ncol = length(est_pe), byrow = TRUE)
  bse <- (apply(draws, 2, quantile, .75, na.rm = TRUE) - apply(draws, 2, quantile, .25, na.rm = TRUE))/(qnorm(0.75) - qnorm(.25))
  bm <- apply(draws, 2, mean) # bias estimator
  zs <- apply(draws[, is.0] / matrix(bse[is.0], nrow = B, ncol = length(bse[is.0]), byrow = TRUE), 1, max, na.rm = TRUE)
  crt <- quantile(zs, 1 - alpha)  #critical value
  is.le <- (est_pe - est_u) / bse <= crt
  return(is.le)
}

# ------------------------ Plotting (2-dimensional projection plots of two specified variables) --------------------
#'
#' Plot 2-dimensional projections of variables in interest.
#'
#' \code{Subpopplot} takes output from \code{\link{Subpop}} command as inputs and plots 2-dimensional projection plots of two specified variables.
#' A factor variable has to be put on the y-axis, otherwise the code breaks. The users need to specify the two variables for the projection
#' with \code{varx} and \code{vary}, and \code{output} should be specified as the output of \code{\link{Subpop}}.
#'
#' @param varx        Variable to be plotted on the x-axis.
#' @param vary        Variable to be plotted on the y-axis. If a variable in interest is of type factor, then user
#'                    must put it on the y-axis.
#' @param output      Output of \code{Subpop} command.
#' @param xlim        x-axis limits. Default is \code{range(data[, varx])}
#' @param ylim        y-axis limits. Default is \code{NULL}.
#' @param main        Main title of the plot. Default is \code{NULL}.
#' @param sub         Sub title of the plot. Default is NULL.
#' @param xlab        x-axis label. Default is \code{NULL}.
#' @param ylab        y-axis label. Default is \code{NULL}.
#'
#' @examples
#' data("mortgage")
#' fm <- deny ~ black + p_irat + hse_inc
#' result <- Subpop(fm = fm, data = mortgage, var.T = "black", method = "logit")
#' Subpopplot(varx = mortgage$p_irat, vary = mortgage$ccred, output = result)
#'
#' @importFrom graphics plot polygon lines legend
#' @export
Subpopplot <- function(varx, vary, output, xlim = NULL, ylim = NULL,
                       main = NULL, sub = NULL, xlab = NULL, ylab = NULL){
  if(is.factor(varx)) stop("Variable on the x-axis cannot be a factor.")
  is.he <- output$most
  is.le <- output$least
  subgroup <- output$sub
  u <- output$u
  xlim <- range(varx)
  if(!is.factor(vary)){
    plot(varx, vary, type = "n", xlim, ylim, log="", main, sub, xlab, ylab, col = 4, pch = 20, lwd = 2)
  } else {
    plot(varx, vary, type = "n", xlim, ylim, log="", main, sub, xlab, ylab, col = 4, pch = 20, lwd = 2, yaxt = 'n')
    axis(side = 2, at = c(1 : nlevels(vary)), labels=levels(vary))
  }
  if(is.null(subgroup)){
    # Full sample
    points(varx[is.he & !is.le], vary[is.he & !is.le], col = 'lightblue1', pch = 1, lwd = 2)
    points(varx[is.le & !is.he], vary[is.le & !is.he], col = 4, pch = 20, lwd = 2)
  }else{
    # Sub sample
    subx <- varx[subgroup]
    suby <- vary[subgroup]
    points(subx[is.he & !is.le], suby[is.he & !is.le], col = 'lightblue1', pch = 1, lwd = 2)
    points(subx[is.le & !is.he], suby[is.le & !is.he], col = 4, pch = 20, lwd = 2)
  }
  legend('topleft', c(paste0(u*100,"% Most"), paste0(u*100,"% Least")), col = c(4, 'lightblue1'), pch = c(20, 1), horiz = F, bty = 'n')
}
