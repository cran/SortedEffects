#' Empirical classification analysis (CA) and inference
#'
#' \code{CA} conducts CA estimation and inference on user-specified objects of interest: first (weighted)
#' moment or (weighted) distribution. Users can use \code{t} to specify variables in interest. When object
#' of interest is moment, use \code{cl} to specify linear combinations for hypothesis testing. All estimates
#' are bias-corrected and all confidence bands are monotonized. The bootstrap procedures follow algorithm 2.2
#' as in Chernozhukov, Fernandez-Val and Luo (2018).
#'
#' @return
#' If \code{subgroup = NULL}, all outputs are whole sample. Otherwise output are subgroup results. When
#' \code{interest = "moment"}, the output is a list showing
#' \itemize{
#'   \item \code{est} Estimates of variables in interest.
#'   \item \code{bse} Bootstrap standard errors.
#'   \item \code{joint_p} P-values that are adjusted for multiplicity to account for joint testing for all variables.
#' }
#' If users have further specified \code{cat} (e.g., \code{!is.null(cat)}), the output has a fourth component
#' \itemize{
#'   \item \code{p_cat} P-values that are adjusted for multiplicity to account for joint testing for all variables within a category.
#' }
#' When \code{interest = "dist"}, the output is a list of two components:
#' \itemize{
#'   \item \code{infresults} A list that stores estimates, upper and lower confidence bounds for all variables
#'                           in interest for least and most affected groups.
#'   \item \code{sortvar} A list that stores sorted and unique variables in interest.
#' }
#' We recommend using \code{\link{CAplot}} command for result visualization.
#'
#' @param fm          Regression formula
#' @param data        The data in use (full sample or subpopulation in interset)
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
#'                    weighted (i.i.d exponential weights) bootstrap. Default is \code{NULL}.
#' @param taus        Indexes for quantile regression. Default is \code{c(1:9)/10}.
#' @param u           Percentile of most and least affected. Default is set to be 0.1.
#' @param interest    Generic objects in the least and most affected subpopulations. Two options:
#'                    (1) \code{"moment"}: weighted mean of Z in the u-least/most affected subpopulation.
#'                    (2) \code{"dist"}: distribution of Z in the u-least/most affected subpopulation.
#'                    Default is \code{interest = "moment"}.
#' @param cl          A pre-specified linear combination. Should be a 2 by L matrix. Default is \code{matrix(c(1,0), nrow=2)}. L-th column denotes L-th hypothesis
#'                    For "moment" interest L means the number of hypotheses. \code{cl} must be specified as a matrix
#' @param t           An index for CA object. Should be a 1 by ncol(data) indicator vector. Users can either
#'                    (1) specify names of variables of interest directly, or
#'                    (2) use 1 to indicate the variable of interest. For example, total number of variables is 5 and interested in the 1st and 3rd vars, then specify
#'                    \code{t = c(1, 0, 1, 0, 0)}.
#' @param cat         P-values in classification analysis are adjusted for multiplicity to account for joint testing of
#'                    zero coefficients on for all variables within a category. Specify all variables in interest in a list using numbers
#'                    to denote relative positions. For example, if variables in interest are "educ", "male", "female", "low income", "middle income",
#'                    and "high income", cat should be specified as \code{cat = list(a=1, b=c(2,3), c=c(4,5,6))}. Default of cat is \code{NULL}.
#' @param alpha       Size for confidence interval. Shoule be between 0 and 1. Default is 0.1
#' @param B           Number of bootstrap draws. Default is 10. For more accurate results, we recommend 500.
#' @param ncores      Number of cores for computation. Default is set to be 1. For large dataset, parallel computing
#'                    is highly recommended since bootstrap is time-consuming.
#' @param seed        Pseudo-number generation for reproduction. Default is 1.
#' @param bc          Whether want the estimate to be bias-corrected. Default is \code{TRUE}. If \code{FALSE} uncorrected
#'                    estimate and corresponding confidence bands will be reported.
#' @param range.cb    When \code{interest = "dist"}, we sort and unique variables in interest to estimate weighted CDF. For large dataset there can be
#'                    memory problem storing very many of observations, and thus users can provide a Sort value and the package will sort and unique
#'                    based on the weighted quantile of Sort. If users don't want this feature, set \code{range.cb = NULL}. Default is \code{c(0.5:99.5)/100}.
#'                    To see how \code{range.cb} makes a difference in the plot, refer to the examples in the companion vignette.
#' @param boot.type   Type of bootstrap. Default is \code{boot.type = "nonpar"}, and the package implements nonparametric
#'                    bootstrap. An alternative is \code{boot.type = "weighted"}, and the package implements weighted
#'                    bootstrap.
#' @examples
#' data("mortgage")
#' fm <- deny ~ black + p_irat
#' t <- c(rep(1, 2), rep(0, 14)) # Specify variables in interest
#' cl <- matrix(c(1,0,0,1), nrow=2) # Meaning: show variables in interest for both groups
#' CA <- CA(fm = fm, data = mortgage, var.T = "black", method = "logit", cl = cl, t = t)
#'
#' # Tabulate the results
#' est <- matrix(CA$est, ncol=2)
#' se <- matrix(CA$bse, ncol=2)
#' Table <- matrix(0, ncol=4, nrow=2)
#' Table[, 1] <- est[, 1] # Least Affected Bias-corrected estimate
#' Table[, 2] <- se[, 1] # Corresponding SE
#' Table[, 3] <- est[, 2] # Most affected
#' Table[, 4] <- se[, 2] # Corresponding SE
#' rownames(Table) <- colnames(CA$est)[1:2] # assign names to each row
#' colnames(Table) <- rep(c("Estimate", "SE"), 2)
#'
#' @importFrom boot boot
#' @importFrom rlist list.cbind
#' @importFrom Hmisc wtd.quantile wtd.mean
#' @importFrom stats quantile rexp qnorm
#' @export
CA <- function(fm, data, method = "ols", var.type = "binary", var.T, compare,
               subgroup = NULL, samp_weight = NULL, taus = c(1:9)/10, u = 0.1,
               cl = matrix(c(1,0), nrow=2), t = c(1, 1, rep(0, dim(data)[2]-2)),
               interest = "moment", cat = NULL, alpha = 0.1, B = 10, ncores = 1,
               seed = 1, bc = TRUE, range.cb = c(0.5:99.5)/100, boot.type = "nonpar"){
  # --------------------------- Stopping Conditions ---------------------
  if(alpha >= 1 || alpha <= 0) stop("Please specify a correct size for hypothesis testing.")
  if(u >= 1 || u <= 0) stop("Please provide a meaningful group classification.")
  if(typeof(t) == "double" && length(t) != dim(data)[2]) stop("Number of variable indicators don't match.")
  if(!is.matrix(cl)) stop("cl must be specified as a matrix.")
  if(dim(cl)[1] != 2) stop("Please specify cl as a 2 by L matrix.")
  if(interest != "moment" && interest != "dist") stop("Please specify CA object of interest: moment or dist.")
  if(!is.null(cat) && !is.list(cat)) stop("Please specify cat as a list.")
  if(boot.type != "nonpar" && boot.type != "weighted") stop("Please specify a bootstrap type as a char: nonpar or weighted.")
  # --------------------------- Replace Null Weight Specification ------------------------
  if(is.null(samp_weight)) samp_weight <- rep(1, dim(data)[1])
  # --------------------------- 1. Call to estimate PE ---------------------------
  output <- PEestimate(fm, data, samp_weight, var.type, var.T, compare, method, subgroup, taus)
  PE_est <- output$PE.est
  # if(is.null(subgroup)) is equivalent to if(groupCA=="full")
  # if(!is.null(subgroup)) is equivalent to if(groupCA=="subgroup")

  # --------------------------- 2. u-most and least Affected Groups ---------------------------
  if(method != "QR"){
    if(is.null(subgroup)){
      # Threshold Values for u-most/least Affected for WHOLE sample
      effect.high <- wtd.quantile(PE_est, samp_weight, 1-u)
      effect.low <- wtd.quantile(PE_est, samp_weight, u)
      data$.w <- samp_weight # This is to get sampling weight for the two groups. If weight is NULL, then it does nothing.
      # Two affected groups
      high.affected <- data[PE_est >= effect.high, ]
      low.affected <- data[PE_est <= effect.low, ]
    } else {
      # SUB sample
      PEsub_est <- output$PEsub.est
      PEsub_w <- output$samp_weight_sub
      # Threshold Values for u-most/least Affected
      effect.high <- wtd.quantile(PEsub_est, PEsub_w, 1-u)
      effect.low <- wtd.quantile(PEsub_est, PEsub_w, u)
      data$.w <- samp_weight
      subdata <- data[subgroup, ]
      high.affected <- subdata[PEsub_est >= effect.high, ]
      low.affected <- subdata[PEsub_est <= effect.low, ]
    }
  } # QR requires special treatment due to stacking of quantile indices
  if (method == "QR"){
    if (is.null(subgroup)){
      effect.high <- wtd.quantile(PE_est, matrix(samp_weight, ncol = 1, nrow = nrow(PE_est), byrow = FALSE), 1-u) # This is where QR is different.
      effect.low <- wtd.quantile(PE_est, matrix(samp_weight, ncol = 1, nrow = nrow(PE_est), byrow = FALSE), u)
      # Obtain sampling weight of each group
      data$.w <- samp_weight
      mesh <- data[rep(1:nrow(data), times = length(taus)), ] #Kronecker function doesn't work for all data types, so I wrote this alternative
      # Two Affected Groups
      high.affected <- mesh[PE_est >= effect.high, ]
      low.affected <- mesh[PE_est <= effect.low, ]
    } else { # Subsample: use PEsub_est
      PEsub_est <- output$PEsub.est
      PEsub_w <- output$samp_weight_sub
      # Threshold Values for u-most/least Affected
      effect.high <- wtd.quantile(PEsub_est, PEsub_w, 1-u)
      effect.low <- wtd.quantile(PEsub_est, PEsub_w, u)
      data$.w <- samp_weight
      subdata <- data[subgroup, ]
      mesh <- subdata[rep(1:nrow(subdata), times = length(taus)), ]
      # Two Affected Groups
      high.affected <- mesh[PEsub_est >= effect.high, ]
      low.affected <- mesh[PEsub_est <= effect.low, ]
    }
  }
  # Obtain the sampling weight for each group (If samp_weight = NULL, then these two vars are NULL)
  weight.high <- high.affected$.w
  weight.low <- low.affected$.w
  # Remove the .w
  high.affected$.w <- NULL
  low.affected$.w <- NULL
  data$.w <- NULL

  ############################ Classicification Analysis Starts Here ############################

  # ----------------------------  3. Index Configuration -----------------------------
  # if t is variable names, convert to index
  if(typeof(t) == "character"){
    posit <- match(t, colnames(data))
    ind <- rep(0, dim(data)[2])
    ind[posit] <- 1
    t <- ind
  }
  # Now that t is an index, use it to get prod and vars
  n2 <- dim(data)[2]
  col.num <- c(1:n2)
  # prod is index vector what variables are selected
  prod <- col.num * t
  # The new data is data[, prod]
  # Number of characateristics in interest
  m <- dim(data[, prod])[2]
  if(is.null(m)) m <- 1
  # Extract Names of variables
  vars <- colnames(data)[which(prod != 0)]

  # ---------------------------- 4. Specify object of interest -----------------------
  if(interest == "moment"){
    # weighted mean of var in interest for the u-least/most affected groups
    interest.high <- apply(as.matrix(high.affected[, prod]), 2, wtd.mean, weights = weight.high, na.rm = TRUE)
    interest.low <- apply(as.matrix(low.affected[, prod]), 2, wtd.mean, weights = weight.low, na.rm = TRUE)
    # Average of var in interest for the two groups (m * 2), where m is number of variables in interest
    est.interest <- cbind(interest.low, interest.high)
    # Hypothesis in interest (m * L matrix)
    # Component (a, b) meaning: b-th hypothesis statistics for a-th variable in interest
    H <- est.interest %*% cl
    # Reshape: 1 * mL
    H <- matrix(H, nrow=1)
    # Assign names
    colnames(H) <- rep(vars, length(H)/m)
  } else if(interest == "dist"){
    # weighted CDF of var in interest for the u-least/most affected groups
    # Sort and eliminate repetitive obs. Note the list "temp" is ragged.
    # temp0 will be used at the end of the function
    # Depending on whether user provides Sort
    if(is.null(range.cb)){
      temp0 <- temp <- lapply(as.data.frame(data[, prod]), function(x) sort(unique(x)))
    } else {
      temp0 <- temp <- lapply(as.data.frame(data[, prod]), function(x) sort(unique(wtd.quantile(x, samp_weight, range.cb))))
    }
    names(temp0) <- names(temp) <- vars
    # Double each list and then split each list into a k by 2 matrix, where k denotes number of elements for each var
    temp <- lapply(temp, rep, 2)
    temp <- lapply(temp, matrix, ncol = 2)
    # I hate for loops in general, but this loop is short: it assigns names to each matrix within the list
    # In general, I find it easier to use for loop if the purpose is to modify part of dataframe
    for (x in 1:length(temp)) colnames(temp[[x]]) <- rep(names(temp)[[x]], 2)
    # Get the trim tails (for inference and plotting)
    trim <- function(x){
      varname <- colnames(x)[1]
      notails <- (x[, 1] >= wtd.quantile(unlist(data[varname]), samp_weight, 0.02) &
                    x[, 1] <= wtd.quantile(unlist(data[varname]), samp_weight, 0.98))
      notails <- rep(notails, 2)
    }
    trimtails <- lapply(temp, trim)
    # The following is a function to calculate weighted cdf
    fun <- function(a){ # a is a name component in the list
      if(is.null(samp_weight)){
        # Take care of no samp weight scenario
        weight.low <- rep(1, length(unlist(low.affected[, colnames(a)[1]])))
        weight.high <- rep(1, length(unlist(high.affected[, colnames(a)[1]])))
      }
      cdf.low <- lapply(a[, 1], wcdf, x = unlist(low.affected[, colnames(a)[1]]), w = weight.low) # low group cdf
      cdf.high <- lapply(a[, 2], wcdf, x = unlist(high.affected[, colnames(a)[[1]]]), w = weight.high) # high group cdf
      output <- cbind(cdf.low, cdf.high)
    }
    # Call 'fun' to calculate wcdfs for both groups for all vars in interest (least group 1st col, most group 2nd col)
    cdf_bundle <- lapply(temp, fun)
    index <- lengths(cdf_bundle)
    # Reshape
    H <- lapply(cdf_bundle, matrix, nrow = 1)
    H <- list.cbind(H) # estimates of weighted cdf
    trimtails <- lapply(trimtails, matrix, nrow=1)
    trimtails <- list.cbind(trimtails) # long vector of trimtail logics
  }
  #  ---------------------------- 6. The Bootstrap Statistics Function -----------------------
  # No weight bootstrap
  stat.boot.CA.noweight <- function(data, indices){
    data$.w <- samp_weight
    data <- data[indices, ]
    output.bs <- PEestimate(fm, data, samp_weight = data$.w, var.type, var.T, compare, method, subgroup, taus)
    PE_est <- output.bs$PE.est
    if(method != "QR"){
      if(is.null(subgroup)){
        # Threshold Values for u-most/least Affected
        effect.high <- wtd.quantile(PE_est, data$.w, 1-u)
        effect.low <- wtd.quantile(PE_est, data$.w, u)
        # Two affected groups
        high.affected <- data[PE_est >= effect.high, ]
        low.affected <- data[PE_est <= effect.low, ]
      } else {
        PEsub_est <- output.bs$PEsub.est
        PEsub_w <- output.bs$samp_weight_sub
        # Threshold Values for u-most/least Affected
        effect.high <- wtd.quantile(PEsub_est, PEsub_w, 1-u)
        effect.low <- wtd.quantile(PEsub_est, PEsub_w, u)
        subdata <- data[subgroup, ]
        high.affected <- subdata[PEsub_est >= effect.high, ]
        low.affected <- subdata[PEsub_est <= effect.low, ]
      }
    } # QR requires special treatment due to stacking of quantile indices
    if (method == "QR"){
      # Full sample: use PE_est
      if (is.null(subgroup)){
        # Threshold Values for u-most/least Affected
        effect.high <- wtd.quantile(PE_est, matrix(data$.w, ncol = 1, nrow = nrow(PE_est), byrow = FALSE), 1-u) # This is where QR is different.
        effect.low <- wtd.quantile(PE_est, matrix(data$.w, ncol = 1, nrow = nrow(PE_est), byrow = FALSE), u)
        mesh <- data[rep(1:nrow(data), times = length(taus)), ]
        # Two Affected Groups
        high.affected <- mesh[PE_est >= effect.high, ]
        low.affected <- mesh[PE_est <= effect.low, ]
      } else {
        PEsub_est <- output.bs$PEsub.est
        PEsub_w <- output.bs$samp_weight_sub
        # Subsample: use PEsub_est
        # Threshold Values for u-most/least Affected
        effect.high <- wtd.quantile(PEsub_est, PEsub_w, 1-u)
        effect.low <- wtd.quantile(PEsub_est, PEsub_w, u)
        subdata <- data[subgroup, ]
        mesh <- subdata[rep(1:nrow(subdata), times = length(taus)), ]
        # Two Affected Groups
        high.affected <- mesh[PEsub_est >= effect.high, ]
        low.affected <- mesh[PEsub_est <= effect.low, ]
      }
    }
    # Obtain the sampling weight for each group (If samp_weight = NULL, then these two vars are NULL)
    weight.high <- high.affected$.w
    weight.low <- low.affected$.w
    # Remove the .w
    high.affected$.w <- NULL
    low.affected$.w <- NULL
    data$.w <- NULL

    if(interest == "moment"){
      # weighted mean of var in interest for the u-least/most affected groups
      interest.high <- apply(as.matrix(high.affected[, prod]), 2, weighted.mean, weights = weight.high, na.rm = TRUE)
      interest.low <- apply(as.matrix(low.affected[, prod]), 2, weighted.mean, weights = weight.low, na.rm = TRUE)
      # Average of var in interest for the two groups (m * 2), where m is number of variables in interest
      est.interest <- cbind(interest.low, interest.high)
      # Hypothesis in interest (m * L matrix)
      # Component (a, b) meaning: b-th hypothesis statistics for a-th variable in interest
      H <- est.interest %*% cl
      # Reshape: 1 * mL
      H <- matrix(H, nrow=1)
      # Assign names
      colnames(H) <- rep(vars, length(H)/m)
    }
    if (interest == "dist"){
      fun <- function(a){ # a is a name component in the list
        if(is.null(samp_weight)){
          # Take care of no samp weight scenario
          weight.low <- rep(1, length(unlist(low.affected[, colnames(a)[1]])))
          weight.high <- rep(1, length(unlist(high.affected[, colnames(a)[1]])))
        }
        cdf.low <- lapply(a[, 1], wcdf, x = unlist(low.affected[, colnames(a)[1]]), w = weight.low) # low group cdf
        cdf.high <- lapply(a[, 2], wcdf, x = unlist(high.affected[, colnames(a)[[1]]]), w = weight.high) # high group cdf
        output <- cbind(cdf.low, cdf.high)
      }
      # Call 'fun' to calculate wcdfs for both groups for all vars in interest (least group 1st col, most group 2nd col)
      cdf_bundle <- lapply(temp, fun)
      index <- lengths(cdf_bundle)
      # Reshape
      H <- lapply(cdf_bundle, matrix, nrow = 1)
      H <- unlist(list.cbind(H)) # estimates of weighted cdf
    }
    out <- H
  }

  data.rg <- function(data, mle){
    n <- dim(data)[1]
    # Exponential weights
    multipliers  <- rexp(n)
    # Sampling weight of data.bs
    # Multiply 20000 to enlarge the weight (otherwise the effect.high and effect.low in stat.boot.CA.weight
    # function would be the same. Same strategy applies to u.Subpop. For SPE it doesn't matter since we plot
    # the entire quantile index. Computation is subtle...)
    weight <- (multipliers/sum(multipliers)) * samp_weight * 20000
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

  stat.boot.CA.weight <- function(data){
    output.bs <- PEestimate(fm, data = data, samp_weight = data$.w, var.type, var.T, compare, method, subgroup, taus)
    PE_est <- output.bs$PE.est
    if(method != "QR"){
      if(is.null(subgroup)){
        # Threshold Values for u-most/least Affected
        effect.high <- wtd.quantile(PE_est, data$.w, 1-u)
        effect.low <- wtd.quantile(PE_est, data$.w, u)
        # Two affected groups
        high <- PE_est >= effect.high
        low <- PE_est <= effect.low
        high.affected <- data[PE_est >= effect.high, ]
        low.affected <- data[PE_est <= effect.low, ]
      } else {
        PEsub_est <- output.bs$PEsub.est
        PEsub_w <- output.bs$samp_weight_sub
        # Threshold Values for u-most/least Affected
        effect.high <- wtd.quantile(PEsub_est, PEsub_w, 1-u)
        effect.low <- wtd.quantile(PEsub_est, PEsub_w, u)
        subdata <- data[subgroup, ]
        high.affected <- subdata[PEsub_est >= effect.high, ]
        low.affected <- subdata[PEsub_est <= effect.low, ]
      }
    } # QR requires special treatment due to stacking of quantile indices
    if (method == "QR"){
      # Full sample: use PE_est
      if (is.null(subgroup)){
        # Threshold Values for u-most/least Affected
        effect.high <- wtd.quantile(PE_est, matrix(data$.w, ncol = 1, nrow = nrow(PE_est), byrow = FALSE), 1-u)
        effect.low <- wtd.quantile(PE_est, matrix(data$.w, ncol = 1, nrow = nrow(PE_est), byrow = FALSE), u)
        mesh <- data[rep(1:nrow(data), times = length(taus)), ]
        # Two Affected Groups
        high.affected <- mesh[PE_est >= effect.high, ]
        low.affected <- mesh[PE_est <= effect.low, ]
      } else {
        PEsub_est <- output.bs$PEsub.est
        PEsub_w <- output.bs$samp_weight_sub
        # Subsample: use PEsub_est
        # Threshold Values for u-most/least Affected
        effect.high <- wtd.quantile(PEsub_est, PEsub_w, 1-u)
        effect.low <- wtd.quantile(PEsub_est, PEsub_w, u)
        subdata <- data[subgroup, ]
        mesh <- subdata[rep(1:nrow(subdata), times = length(taus)), ]
        # Two Affected Groups
        high.affected <- mesh[PEsub_est >= effect.high, ]
        low.affected <- mesh[PEsub_est <= effect.low, ]
      }
    }
    # Obtain the sampling weight for each group
    weight.high <- high.affected$.w
    weight.low <- low.affected$.w
    # Remove the .w
    high.affected$.w <- NULL
    low.affected$.w <- NULL
    data$.w <- NULL

    if(interest == "moment"){
      # weighted mean of var in interest for the u-least/most affected groups
      interest.high <- apply(as.matrix(high.affected[, prod]), 2, wtd.mean, weights = weight.high, na.rm = TRUE)
      interest.low <- apply(as.matrix(low.affected[, prod]), 2, wtd.mean, weights = weight.low, na.rm = TRUE)
      # Average of var in interest for the two groups (m * 2), where m is number of variables in interest
      est.interest <- cbind(interest.low, interest.high)
      # Hypothesis in interest (m * L matrix)
      # Component (a, b) meaning: b-th hypothesis statistics for a-th variable in interest
      H <- est.interest %*% cl
      # Reshape: 1 * mL
      H <- matrix(H, nrow=1)
      # Assign names
      colnames(H) <- rep(vars, length(H)/m)
    }
    if (interest == "dist"){
      fun2 <- function(a, weight.high = weight.high, weight.low = weight.low){ # a is a name component in the list
        cdf.low <- lapply(a[, 1], wcdf, x = unlist(low.affected[, colnames(a)[1]]), w = weight.low) # low group cdf
        cdf.high <- lapply(a[, 2], wcdf, x = unlist(high.affected[, colnames(a)[1]]), w = weight.high) # high group cdf
        output <- cbind(cdf.low, cdf.high)
      }
      # Call 'fun' to calculate wcdfs for both groups for all vars in interest (least group 1st col, most group 2nd col)
      cdf_bundle <- lapply(temp, fun2)
      index <- lengths(cdf_bundle)
      # Reshape
      H <- lapply(cdf_bundle, matrix, nrow = 1)
      H <- unlist(list.cbind(H)) # estimates of weighted cdf
    }
    out <- H
  }

  # -----------------------------------  7. Call boot for bootstrap ----------------------------------------
  set.seed(seed)
  if(boot.type == "nonpar"){
    if(method != "QR"){
      result.boot <- boot(data = data, statistic = stat.boot.CA.noweight, parallel="multicore", ncpus = ncores, R = B)
    } else{
      data$.w <- samp_weight
      result.boot <- boot(data = data, statistic = stat.boot.CA.weight, sim = "parametric", ran.gen = data.non, mle = 0, parallel = "multicore", ncpus = ncores, R = B)
      data$.w <- NULL
    }
  } else if(boot.type == "weighted"){
    data$.w <- samp_weight
    result.boot <- boot(data = data, statistic = stat.boot.CA.weight, sim = "parametric", ran.gen = data.rg, mle = 0, parallel = "multicore", ncpus = ncores, R = B)
    data$.w <- NULL
  }

  # -----------------------------------  8. Inference ----------------------------------------
  if(interest == "moment"){
    # draws.H is B * mL
    draws.H <- result.boot$t
    colnames(draws.H) <- rep(vars, dim(cl)[2])
    # Zc is B * mL
    Zc <- draws.H - matrix(H, nrow = B, ncol = length(H), byrow = TRUE)
    # sig is 1 * mL (bootstrap standard error)
    sig <- (apply(Zc, 2, quantile, 0.75, na.rm = TRUE) - apply(Zc, 2, quantile, 0.25, na.rm = TRUE))/(qnorm(0.75) - qnorm(0.25))
    # t.tilde is B * mL
    t.tilde <- apply(abs(Zc)/matrix(sig, nrow = B, ncol = length(sig), byrow = TRUE), 1, max, na.rm = TRUE)
    # Bias correction
    H.bc <- 2 * H - apply(draws.H, 2, mean, na.rm = TRUE) # 1 * mL
    # Joint p-values
    stat <- abs(H.bc)/sig
    # Proportions of the B draws of t.tilde that are greater than s
    j_pvals <- sapply(stat, epvals, zs = t.tilde) # 1 * mL. First m columns represent joint p-vals for the m vars for the first hypothesis.
    # Not bias-corrected joint p-values
    stat.un <- abs(H)/sig
    j_pvals.un <- sapply(stat.un, epvals, zs = t.tilde)
    # Now work on p-vals accounting for testing all vars within one category
    if(!is.null(cat)){
      # The aux function subset variables within one category among all hypotheses. (which function breaks for number of vars >= 4,
      # so I need to come up this function...)
      aux <- function(x, target){
        M <- rep(x, dim(cl)[2]) # Recall dim(cl)[2] denotes number of hypotheses L
        i <- rep(m, dim(cl)[2])
        i <- cumsum(i)
        i <- i - m
        index <- M + i
        mat <- target[, index]
      }
      Zc.cat <- lapply(cat, aux, target = Zc)
      H.bc.cat <- lapply(cat, aux, target = H.bc)
      H.cat <- lapply(cat, aux, target = H)
      # The following function takes each list component as input and calculates bootstrap se
      joint_se <- function(x){
        if(!is.null(dim(x))){
          sig <- (apply(x, 2, quantile, .75, na.rm=TRUE) - apply(x, 2, quantile, .25, na.rm=TRUE))/(qnorm(0.75) - qnorm(.25))
        } else { # Special case: one variable one hypothesis
          sig <- (quantile(x, 0.75, na.rm = TRUE) - quantile(x, 0.25, na.rm = TRUE))/(qnorm(0.75) - qnorm(0.25))
        }
        return(sig)
      }
      # The following function takes each list component as input and calculates t.tilde
      joint_t <- function(x){
        if(!is.null(dim(x))){
          sig <- (apply(x, 2, quantile, .75, na.rm=TRUE) - apply(x, 2, quantile, .25, na.rm=TRUE))/(qnorm(0.75) - qnorm(.25))
          t.tilde.cat <- apply(abs(x)/matrix(sig, nrow = B, ncol = length(sig), byrow=TRUE), 1, max, na.rm = TRUE)
        } else {
          sig <- (quantile(x, 0.75, na.rm = TRUE) - quantile(x, 0.25, na.rm = TRUE))/(qnorm(0.75) - qnorm(0.25))
          t.tilde.cat <- abs(x)/matrix(sig, nrow = B, ncol = length(sig), byrow = TRUE)
        }
        return(t.tilde.cat)
      }
      bse <- lapply(Zc.cat, joint_se) # k lists (k denotes # of categories)
      t.tilde.cat <- lapply(Zc.cat, joint_t)
      H.bc.cat <- unlist(H.bc.cat)
      H.cat <- unlist(H.cat)
      bse <- unlist(bse)
      # a vector of stat for p-value evaluation.
      stat.cat <- abs(H.bc.cat)/bse # bias-corrected
      stat.cat.un <- abs(H.cat)/bse # For uncorrected estimate
      t.tilde.cat <- unlist(t.tilde.cat) # B * k elements
      t.tilde.cat <- matrix(t.tilde.cat, ncol = B, byrow = TRUE) # convert to a k by B matrix
      # Create a temp matrix container
      A <- matrix(0, nrow = length(cat), ncol = length(stat.cat))
      B <- matrix(0, nrow = length(cat), ncol = length(stat.cat.un))
      # For all rows of A, calculate "p-values". Of course many values will be nonsensical.
      for (i in 1:nrow(A)) A[i, ] <- sapply(stat.cat, epvals, zs = t.tilde.cat[i, ])
      for (i in 1:nrow(B)) B[i, ] <- sapply(stat.cat.un, epvals, zs = t.tilde.cat[i, ])

      # Now extract the p-values that we want
      # splitting criterion (how to distinguish the stat.cat variables)
      sp_crit <- lengths(cat)*dim(cl)[2]
      CumSum <- cumsum(sp_crit)
      for(i in 2:length(cat)) A[1, (CumSum[i-1]+1) : CumSum[i]] <- A[i, (CumSum[i-1]+1) : CumSum[i]]
      for(i in 2:length(cat)) B[1, (CumSum[i-1]+1) : CumSum[i]] <- B[i, (CumSum[i-1]+1) : CumSum[i]]
      # This is the p-values accounting for testing within a category
      p_cat <- A[1, ]
      p_cat.un <- B[1, ]
      # Return output
      if(bc == TRUE){
        output <- list(est = H.bc, bse = sig, joint_p = j_pvals, p_cat = p_cat)
      } else {
        output <- list(est = H, bse = sig, joint_p = j_pvals.un, p_cat = p_cat.un)
      }
    } else { # if cat is null
      if(bc == TRUE){
        output <- list(est = H.bc, bse = sig, joint_p = j_pvals)
      } else {
        output <- list(est = H, bse = sig, joint_p = j_pvals.un)
      }
    }
  } else if(interest == "dist"){
    draws.H <- result.boot$t # nrow is B
    # Zc <- mapply("-", draws.H, lapply(H, rep, B))
    Zc <- draws.H - matrix(unlist(H), nrow = B, ncol = length(H), byrow = TRUE)
    bse <- (apply(Zc, 2, quantile, .75, na.rm = TRUE) - apply(Zc, 2, quantile, .25, na.rm = TRUE))/ (qnorm(0.75) - qnorm(.25))
    # Bias corrected estimation
    bm <- apply(draws.H, 2, mean, na.rm = TRUE)
    H.bc <- 2*unlist(H) - bm
    # Extract the var in interest by Assigning names to H
    # The following chunk of code is to get A: column of variable names
    indexname <- names(index)
    A <- NULL
    A[indexname] <- list(NULL)
    for (i in 1:length(A)) A[i] <- names(A[i])
    # A small function to assign names
    assign <- function(x, y=index){
      num <- which(names(y)==x)
      x <- rep(names(y)[num], y[num])
    }
    A <- lapply(A, assign)
    A <- lapply(A, matrix, nrow=1)
    A <- list.cbind(A)
    colnames(Zc) <- A
    # The following function outputs a bundle of bias-corrected estimate, upper and lower bound for both groups for a variable in interest
    CDFinf <- function(x){
      Zc.sub <- x[1:B, ]
      bse.sub <- bse.sub1 <- x[B+1, ]
      bse.sub1[bse.sub1==0] <- NA
      est.sub <- x[B+2, ]
      trim.sub <- as.logical(x[B+3, ])
      tempo <- abs(Zc.sub)/matrix(bse.sub1, nrow = B, ncol = length(bse.sub1), byrow = TRUE)
      t.sub <- apply(tempo[, trim.sub], 1, max, na.rm = TRUE) # B by 1
      crt.sub <- quantile(t.sub, 1 - alpha) # scalar
      # For Least-Affected Group (bias-corrected)
      least_upper <- est.sub[1:(length(est.sub)/2)] + crt.sub * bse.sub[1:(length(bse.sub)/2)]
      least_lower <- est.sub[1:(length(est.sub)/2)] - crt.sub * bse.sub[1:(length(bse.sub)/2)]
      # For Most-Affected Group (bias-corrected)
      most_upper <- est.sub[-(1:(length(est.sub)/2))] + crt.sub * bse.sub[-(1:(length(bse.sub)/2))]
      most_lower <- est.sub[-(1:(length(est.sub)/2))] - crt.sub * bse.sub[-(1:(length(bse.sub)/2))]
      bundle <- list(least = est.sub[1:(length(est.sub)/2)], least_up = least_upper, least_low = least_lower,
                     most = est.sub[-(1:(length(est.sub)/2))], most_up = most_upper, most_low = most_lower)
      # Impose shape restriction
      bundle <- lapply(bundle, function(y) sort(y * (y > 0 & y < 1) + (y >= 1)))
    }
    if(bc == TRUE){
      # Returns a list storing statistics in interest
      info <- rbind(Zc, bse, H.bc, trimtails) # Zc (row 1:B), bse (row B+1), H.bc (row B+2), trimtails (row B+3)
      info.sub <- lapply(split.data.frame(t(info), colnames(info)), t)
      # Now call function "CDFinf" to get the bundled result as a list.
      infresults <- lapply(info.sub, CDFinf)
      # Return output
      output <- list(infresults = infresults, sortvar = temp0, alpha = alpha)
    }else if(bc == FALSE){
      # For non bias-correction case
      info_n <- rbind(Zc, bse, unlist(H), trimtails) # For non bias-corrected
      info.sub_n <- lapply(split.data.frame(t(info_n), colnames(info_n)), t)
      infresults_n <- lapply(info.sub_n, CDFinf)
      output <- list(infresults = infresults_n, sortvar = temp0, alpha = alpha)
    }
  }
}

# ---------------------------- Two Auxiliary Functions -----------------------
# 1. Function to obtain empirical p-values
epvals <- function(z, zs){
  return(mean(zs > z))
}

# 2. Weighted distribution function estimation
#' @importFrom stats weighted.mean
wcdf <- function(y, x, w){
  Ff <- weighted.mean((x <= y), w)
  return(Ff)
}

# ----------------- Plotting (CDFs of a continuous variable in interest for least/most affected groups) ------------
#' Distribution plotting
#'
#' \code{CAplot} plots distributions and joint uniform confidence bands of variables in interest from \code{\link{CA}} command.
#'
#' @param var         The variable in interest after having been sorted. Must be specified as a character.
#' @param output      Output of \code{CA} command.
#' @param xlim        x-axis limits. Default is range of percentile index.
#' @param ylim        y-axis limits. Default is NULL.
#' @param main        Main title of the plot. Defualt is NULL.
#' @param sub         Sub title of the plot. Default is NULL.
#' @param xlab        x-axis label. Default is "Percentile Index".
#' @param ylab        y-axis label. Default is "Sorted Effects".
#'
#' @examples
#' data("mortgage")
#' fm <- deny ~ black + p_irat
#' t <- c(rep(1, 2), rep(0, 14)) # Specify variables in interest
#' cl <- matrix(c(1,0,0,1), nrow=2) # Meaning: show variables in interest for both groups
#' # Drawing Distributions
#' CAgraph <- CA(fm=fm, data = mortgage, var.T = "black", method = "logit", cl = cl,
#'               t = t, interest = "dist")
#' CAplot(var = "p_irat", output = CAgraph)
#'
#' @importFrom graphics axis legend lines plot polygon
#' @export
CAplot <- function(var, output, xlim = NULL, ylim = NULL, main = NULL,
                   sub = NULL, xlab = NULL, ylab = NULL){
  alpha = output$alpha
  if(is.na(match(var, names(output$sortvar)))) stop("The variable must be consistent with interest specification in u.CA command.")
  xvar <- unlist(output$sortvar[var])
  yvars <- unlist(output$infresults[var])
  q <- length(xvar)
  plot(xvar, yvars[(3*q+1):(4*q)], type = "n", xlim, ylim, log = "", main, sub, xlab, ylab, col=4, lwd=2) # frame
  polygon(c(xvar, rev(xvar)),c(yvars[(q+1):(2*q)], rev(yvars[(2*q+1):(3*q)])), density = 60, border = F,
          col = 'darkcyan', lty = 1, lwd = 1, type = "l") # CB of least affected
  polygon(c(xvar, rev(xvar)), c(yvars[(4*q+1):(5*q)], rev(yvars[(5*q+1):(6*q)])), density = 60, border = F,
          col = 'tomato', lty = 1, lwd = 1, type = "l") # CB of most affected
  lines(xvar, yvars[1:q], lwd = 2, col = 4, type = "l") # Est of least affected
  lines(xvar, yvars[(3*q+1):(4*q)], lwd=2, col = 2, type = "l") # Est of most affected
  legend(x = "topleft", legend = c("Least Affected","Most Affected", paste0((1-alpha)*100,"% CB"), paste0((1-alpha)*100,"% CB")),
         col = c(4, 2, "darkcyan","tomato"), lwd = c(1, 1, 5, 5), lty = c(1, 1, 1, 1), bty = 'n')
}
