#' @title Sensitivity analysis for effects of multiple exposures in the presence of unmeasured confounding based on an alternative assumption.
#'
#' @description \code{GSAMU.alt()} is an alternative sensitivity analysis method using assumption for the conditional distribution of U | (L, X) rather than for the joint distribution of (U, L, X).
#'
#' @param data A data frame in which contains outcome, measured confounders, and exposures.
#' @param outcome The name of variable for outcome. For continuous, count, and binary outcomes, a single character for outcome is specified. For time-to-event outcome, two characters corresponding to survival time and status are specified (i.e., outcome=c("survivaltime", "status")).
#' @param outcome.type The type of outcome variable. Possible values are "continuous" (for continuous outcome), "count" (for count outcome), "binary" (for binary outcome), or "timetoevent" (for time-to-event outcome).
#' @param link The specification for the model link function. Possible values are "identity" (for continuous outcome), "log" (for count outcome), "probit" (for binary outcome), or "logit" (for binary outcome). Not used for outcome.type="timetoevent".
#' @param hazard.model The specification for the hazard model. Possible values are "coxph" (for Cox PH model) or "ah" (for additive hazard model). Not used for outcome.type="continuous", "binary", or "count".
#' @param confounder The vector of variable names for confounders. See Examples.
#' @param exposure The vector of variable names for exposures. See Examples.
#' @param delta The values of \eqn{\delta}s. The length of delta must be less than 6. See Examples.
#' @param bound The range of \eqn{(\rho_{1}, \ldots, \rho_{p})}. The order of inputs is (\eqn{\rho_{1,min}, \ldots, \rho_{p,min}, \rho_{1,max}, \ldots, } \eqn{\rho_{p,max}}). See Examples.
#' @param bound.sigma2 The range of \eqn{\sigma_{c}^{2}}. The length of bound.sigma2 must be 1 (for constant variance) or 2 (for lower and upper bounds of variance). See Examples.
#' @param bootsCI Logical value: if TRUE, the process to obtain the bootstrap confidence interval for the population sensitivity interval is carried out. Default: TRUE.
#' @param B The number of bootstrap replicates. Not used for bootsCI=FALSE. Default: 1000
#' @param seed The value of seed number when bootstrap works. Not used for bootsCI=FALSE. Default: 231111.
#' @param verbose Logical value: if TRUE, bootstrap process is output every 100 times. Not used for bootsCI=FALSE. Default: TRUE.
#'
#' @return An object of class \code{GSAMU}. The object is a data.frame with the following components:
#' \item{label}{Exposures}
#' \item{delta}{The values of \eqn{\delta}}
#' \item{coef}{The coefficients of exposures}
#' \item{Lower.SI}{The lower sensitivity interval of conditional single- or joint-exposure effect}
#' \item{Upper.SI}{The upper sensitivity interval of conditional single- or joint-exposure effect}
#' \item{Lower.CI}{The lower confindece interval of sensitivity interval. Results are displayed only when bootstrap is performed (i.e., bootsCI=TRUE).}
#' \item{Upper.CI}{The upper confindece interval of sensitivity interval. Results are displayed only when bootstrap is performed (i.e., bootsCI=TRUE).}
#'
#' The results for the \code{GSAMU.alt} are printed with the \code{\link{print.GSAMU}} function.
#' To generate the plot of results for the \code{GSAMU.alt}, use the \code{\link{autoplot.GSAMU}} function.
#'
#' @details See Lee et al. (2024+) for details.
#'
#' @examples
#' ############################################
#' ## Real data analysis for binary outcome ##
#' ############################################
#' ### Log transformed of exposures.
#' library(readr)
#' ul=paste0("https://static-content.springer.com/esm/art%3A10.1186%2Fs12940-020-00644-4/",
#'           "MediaObjects/12940_2020_644_MOESM2_ESM.csv")
#' data_orig <- as.data.frame(read_csv(url(ul))); data_orig <- na.omit(data_orig)
#'
#' ### Data cleaning
#' ## Make binary outcome
#' data_orig$triglyceride_b <- ifelse(data_orig$triglyceride >= 150, 1, 0)
#' data_orig$triglyceride <- NULL
#'
#' ### Environmental factors
#' ## All exposure variables, except for retinol, were log-transformed,
#' # Retinyl-palmitate
#' data_orig$a3.Retinyl.palmitate <- log(data_orig$a3.Retinyl.palmitate)
#' # Retinol
#' summary(data_orig$a5.Retinol)
#' # trans.b.carotene
#' data_orig$a1.trans.b.carotene <- log(data_orig$a1.trans.b.carotene)
#' # a-Tocopherol
#' data_orig$a7.a.Tocopherol <- log(data_orig$a7.a.Tocopherol)
#' # g-tocopherol
#' data_orig$a6.g.tocopherol <- log(data_orig$a6.g.tocopherol)
#'
#' ### Covariates
#' # Age
#' summary(data_orig$age)
#' # Sex
#' data_orig$sex <- ifelse(data_orig$sex == 2, 1, 0)
#' # Race [Race/Ethnicity (1: Non-Hispanic white; 2: Non-Hispanic black; 3: Mexican American;
#' #                       4: Other race - Including multi-racial; 5: Other Hispanic)]
#' data_orig$race <- ifelse(data_orig$rac == 1, 1, 0)
#'
#' ## using only needed variables
#' usedata <- data_orig[,c("triglyceride_b", "sex", "race", "age",
#'                         "a3.Retinyl.palmitate", "a5.Retinol", "a1.trans.b.carotene",
#'                         "a7.a.Tocopherol", "a6.g.tocopherol")]
#' ## Rename the varoables
#' names(usedata) <- c("triglyceride_b", "sex", "race", "age",
#'                     "Retinyl palmitate", "Retinol", "trans-beta-carotene",
#'                     "alpha-Tocopherol", "gamma-Tocopherol")
#' data_r4 <- usedata
#'
#' ## glm fitting (working regression model)
#' fit.glm <- glm(triglyceride_b~., family=binomial(link="logit"), data=data_r4)
#' summary(fit.glm)
#' round(cbind("Est"=exp(coef(fit.glm)), exp(confint(fit.glm)),
#'       "P value"=summary(fit.glm)$coefficients[,4]), 3)
#'
#' ##
#' bound <- c(## lower bound
#'            # Retinyl palmitate    Retinol  trans-beta-carotene
#'            rep(-0.19, 3),
#'            # alpha-Tocopherol  gamma-Tocopherol
#'            rep(-0.19, 2),
#'
#'            ## upper bound
#'            # Retinyl palmitate    Retinol  trans-beta-carotene
#'            rep(0.15, 3),
#'            # alpha-Tocopherol  gamma-Tocopherol
#'            rep(0.10, 2))
#'
#' ## Sensitivity analysis
#' binary.re0 <- GSAMU.alt(data=data_r4, outcome="triglyceride_b", outcome.type="binary",
#'                         link="logit", hazard.model=NULL,
#'                         confounder=c("sex", "race", "age"),
#'                         exposure=c("Retinyl palmitate", "Retinol", "trans-beta-carotene",
#'                                    "alpha-Tocopherol", "gamma-Tocopherol"),
#'                         delta=c(0.11, 0.22, 0.33, 0.44), bound=bound, bound.sigma2=c(0.5, 1.5),
#'                         bootsCI=FALSE, B=1000, seed=231111, verbose=TRUE)
#' binary.re0$result
#' print(binary.re0)
#' autoplot(object=binary.re0, point.size=2.75, width.SI=1, width.CI=0.6,
#'          axis.title.x.size=15, axis.text.size=16, legend.text.size=15,
#'          myxlim=NULL)
#'
#' \dontrun{
#' ## Sensitivity analysis with bootstrap percentile confidence interval
#' binary.re1 <- GSAMU.alt(data=data_r4, outcome="triglyceride_b", outcome.type="binary",
#'                         link="logit", hazard.model=NULL,
#'                         confounder=c("sex", "race", "age"),
#'                         exposure=c("Retinyl palmitate", "Retinol", "trans-beta-carotene",
#'                                    "alpha-Tocopherol", "gamma-Tocopherol"),
#'                         delta=c(0.11, 0.22, 0.33, 0.44), bound=bound, bound.sigma2=c(0.5, 1.5),
#'                         bootsCI=TRUE, B=1000, seed=231111, verbose=TRUE)
#' binary.re1$result
#' print(binary.re1)
#' autoplot(object=binary.re1, point.size=2.75, width.SI=1.55, width.CI=0.6,
#'          axis.title.x.size=15, axis.text.size=16, legend.text.size=15,
#'          myxlim=c(-2, 4))
#' }
#'
#' @seealso
#'  \code{\link[GSAMU]{print.GSAMU}}, \code{\link[GSAMU]{autoplot.GSAMU}}
#'
#' @references
#' Lee S, Jeong B, Lee D, Lee W (2024+):
#' Sensitivity analysis for effects of multiple exposures in the presence of unmeasured confounding: non-Gaussian and time-to-event outcomes
#' \emph{xxx}. DOI: xxx.
#'
#' @keywords Methods
#'
#' @export
GSAMU.alt <- function(data, outcome, outcome.type, link=NULL, hazard.model=NULL,
                      confounder, exposure,
                      delta, bound, bound.sigma2,
                      bootsCI=TRUE, B=1000, seed=231111, verbose=TRUE) {
  ##############################################################################
  ##
  if (outcome.type == "binary") {
    if (!(link %in% c("logit", "probit"))) {
      stop("For binary outcome, only logit or probit link function is allowed.")
    }
    if (!is.null(hazard.model)) {
      warning(paste0("For binary outcome, \"hazard.model\" does not used (i.e., ", hazard.model," is not used)."))
    }
  } else if (outcome.type == "count") {
    if (link != "log") {
      stop("For count outcome, only log link function is allowed.")
    }
    if (!is.null(hazard.model)) {
      warning(paste0("For count outcome, \"hazard.model\" does not used (i.e., ", hazard.model," is not used)."))
    }
  } else if (outcome.type == "continuous") {
    if (link != "identity") {
      stop("For continuous outcome, only identity link function is allowed.")
    }
    if (!is.null(hazard.model)) {
      warning(paste0("For continuous outcome, \"hazard.model\" does not used (i.e., ", hazard.model," is not used)."))
    }
  } else if (outcome.type == "timetoevent") {
    if (!(hazard.model %in% c("coxph", "ah"))) {
      stop("For time-to-event outcome, only coxph (Cox PH model) or ah (additive hazard model) is allowed.")
    }
    if (!is.null(link)) {
      warning(paste0("For time-to-event outcome, \"link\" function does not used (i.e., ", link," is not used)."))
    }
  } else {
    stop("\"outcome.type\" is incorrect. This must be one of \"continuous\", \"count\", \"binary\", or \"timetoevent\". Set again. See the Arguments explanation.")
  }
  ##
  if (!(outcome.type %in% c("continuous", "count", "binary", "timetoevent"))) {
    stop("For other types of outcome, GSAMU() is currently undergoing development.")
  }
  if (outcome.type %in% c("continuous", "count", "binary")) {
    if (!(link %in% c("identity", "log", "logit", "probit"))) {
      stop("For other types of link function, GSAMU() is currently undergoing development.")
    }
  }
  if (outcome.type == "timetoevent") {
    if (length(outcome) != 2) {
      stop("For time-to-event outcome, \"outcome\" must be the variable name of time and status.")
    }
  }
  ##
  if (sum(is.element(el=exposure, set=colnames(data))) != length(exposure)) {
    stop("The specified exposure contains variable(s) that are not in \"data\".")
  }
  if (sum(is.element(el=confounder, set=colnames(data))) != length(confounder)) {
    stop("The specified confounder contains avariable(s) that are not in \"data\".")
  }
  ##
  k <- length(confounder)
  p <- length(exposure)
  if (k == 0) {
    stop("There must be at least one confounder.")
  }
  if (p == 0) {
    stop("There must be at least one exposure.")
  }
  if (length(bound) != (p)*2) {
    stop("\"bound\" is incorrect. Set again. See the Arguments explanation.")
  } else if (sum(bound[1:(p)] <= bound[(p+1):(2*p)]) != p) {
    stop("There is a case in \"bound\" where the lower bound is greater than the upper bound.
         Reset the order of \"bound\". See the Arguments explanation.")
  }
  if (!(length(bound.sigma2) %in% c(1,2))) {
    stop("\"bound.sigma2\" is incorrect. Set again. See the Arguments explanation.")
  }
  ##
  if (length(delta) < 1) {
    stop("The length of \"delta\" should be larger than 0. Set again.")
  } else if (length(delta) > 5) {
    stop("The length of \"delta\" should be smaller than 6. Set again.")
  }

  ##############################################################################
  ## Rearrange data
  if (outcome.type == "timetoevent") {
    # Time-to-event outcome (time, status)
    dataX <- data[, c(confounder, exposure)]
    scaled.dataX <- scale(dataX)
    dataY <- data[, colnames(data) %in% outcome]

    if (!is.numeric(dataY[,1])) {
      stop("Time variable is not numeric.")
    } else if (sum(dataY[,1] >= 0) != length(dataY[,1])) {
      stop("Time variable must be positive.")
    } else if (length(dataY[,1]) != length(dataY[,2])) {
      stop("Time and status are different lengths.")
    } else if (!all.equal(sort(unique(dataY[,2])), c(0,1))) {
      stop("Status indicator must be 0 (for censored) and 1 (for event).")
    }

    data <- data.frame("time"=dataY[,1], "status"=dataY[,2], scaled.dataX)

    if (hazard.model == "coxph") {
      ## Fitting the Cox PH model
      fitmodel <- coxph(formula=Surv(time, status) ~ ., data=data)
    } else {
      ## Fitting the additive hazard model
      # fitmodel <- ah(Surv(time, status) ~ ., data=data) # using addhazard package
      formula_str <- as.formula(paste("Surv(time, status) ~ ",
                                      paste(c(paste0("const(`", confounder, "`)"),
                                              paste0("const(`", exposure, "`)")), collapse=" + ")))
      fitmodel <- aalen(formula=formula_str, data=data) # using timereg package
    }

  } else {
    # Continuous, count, or binary outcome
    dataX <- data[, c(confounder, exposure)]
    scaled.dataX <- scale(dataX)
    dataY <- data[, colnames(data) == outcome]

    if (outcome.type == "binary") {
      if (!all.equal(sort(unique(dataY)), c(0,1))) {
        stop("For binary outcome, \"outcome\" must be 0, 1.")
      }
    } else if (outcome.type == "count") {
      if (sum(dataY >= 0) != length(dataY)) {
        stop("For count outcome, \"outcome\" must be positive.")
      }
    } else if (outcome.type == "continuous") {
      if (!is.numeric(dataY)) {
        stop("For continuous outcome, \"outcome\" must be numeric.")
      }
    }

    data <- data.frame("Y"=dataY, scaled.dataX)
    if (link == "log") {
      ## Fitting the regression model
      fitmodel <- glm(formula=Y ~ ., family=poisson(link=link), data=data)
    } else if (link == "identity") {
      ## Fitting the regression model (identity)
      fitmodel <- lm(formula=Y ~ ., data=data)
    } else {
      ## Fitting the regression model (probit or logit)
      fitmodel <- glm(formula=Y ~ ., family=binomial(link=link), data=data)
    }
  }

  ##############################################################################
  ## Sensitivity analysis
  sens.results <- GSAMU_alt(data=data, fitmodel=fitmodel, exposure=exposure,
                            delta=delta, k=k, p=p, bound=bound, bound.sigma2=bound.sigma2,
                            link=link, hazard.model=hazard.model)

  ## Bootstrap confidence interval
  if (bootsCI == TRUE) {
    ## Check whether data are scaled
    if (sum(round(apply(dataX, 2, mean), 5) == 0) == ncol(dataX) | sum(round(apply(dataX, 2, sd), 5) == 1) == ncol(dataX)) {
      stop("Scaled data cannot be used when calculating bootstrap confidence interval. Put unscaled data.")
    }

    ## Bootstrap results
    set.seed(seed)
    sens.results.boots <- list()
    for(ii in 1:B){
      ##
      data_boot <- data[sample.int(nrow(data), nrow(data), replace=TRUE), ]

      if (outcome.type == "binary") {
        data_boot[, c(2:ncol(data_boot))] <- scale(data_boot[, c(2:ncol(data_boot))])
        if (fitmodel$family$link == "probit") {
          fitmodel.boots <- glm(formula=formula(fitmodel), family=binomial(link="probit"), data=data_boot)
        } else if (fitmodel$family$link == "logit") {
          fitmodel.boots <- glm(formula=formula(fitmodel), family=binomial(link="logit"), data=data_boot)
        }
      } else if (outcome.type %in% c("count", "timetoevent")) {
        if (outcome.type == "count") {
          data_boot[, c(2:ncol(data_boot))] <- scale(data_boot[, c(2:ncol(data_boot))])
          fitmodel.boots <- glm(formula=formula(fitmodel), family=poisson(link="log"), data=data_boot)
        } else if (outcome.type == "continuous") {
          data_boot[, c(2:ncol(data_boot))] <- scale(data_boot[, c(2:ncol(data_boot))])
          fitmodel.boots <- lm(formula=formula(fitmodel), data=data_boot)
        } else {
          data_boot[, c(3:ncol(data_boot))] <- scale(data_boot[, c(3:ncol(data_boot))])
          if (hazard.model == "coxph") {
            fitmodel.boots <- coxph(formula=formula(fitmodel), data=data_boot)
          } else {
            # fitmodel.boots <- ah(formula(fitmodel), data=data_boot)
            formula_str <- as.formula(paste("Surv(time, status) ~ ",
                                            paste(c(paste0("const(`", confounder, "`)"),
                                                    paste0("const(`", exposure, "`)")), collapse=" + ")))
            fitmodel.boots <- aalen(formula=formula_str, data=data_boot) # using timereg package
          }
        }
      }
      re.boots <- GSAMU_alt(data=data_boot, fitmodel=fitmodel.boots, exposure=exposure,
                            delta=delta, k=k, p=p, bound=bound, bound.sigma2=bound.sigma2,
                            link=link, hazard.model=hazard.model)

      sens.results.boots[[ii]] <- re.boots

      rm(data_boot, fitmodel.boots, re.boots)

      if(verbose == TRUE & ii %% 100 == 0){
        cat(ii, "th end! \n")
      }
    }

    ## Summary
    for(i in 1:B){
      if(i == 1){
        result.temp <- report.f(dd=sens.results.boots, index=i)
        results.mat <- result.temp
      } else {
        result.temp <- report.f(dd=sens.results.boots, index=i)
        results.mat <- rbind(results.mat, result.temp)
      }
    }

    # NOTE: nullify variables used for non-standard evaluation for tidyverse/ggplot2 below
    label <- delta  <- Lower.SI <- Upper.SI <- NULL
    results.mat <- as.data.frame(results.mat)
    rownames(results.mat) <- NULL
    results.mat2 <- results.mat %>%
      group_by(label, delta) %>%
      summarise(Lower.CI=quantile(Lower.SI, 0.025),
                Upper.CI=quantile(Upper.SI, 0.975),
                .groups="keep")

    ## Merge the bootstrap results with sensitivity analysis results
    sens.results <- merge(sens.results, results.mat2, by=c("label", "delta"), sort=FALSE)
  }

  ## Rearrange results
  # levels.exposure <- c(gsub(c(" "), ".", exposure), "Joint effect")
  # levels.exposure <- gsub(c("-"), ".", levels.exposure)
  levels.exposure <- c(exposure, "Joint effect")
  sens.results$label <- factor(sens.results$label, levels=levels.exposure)
  sens.results <- sens.results[order(sens.results$label), ]
  rownames(sens.results) <- NULL

  ##
  if (bootsCI == FALSE) {
    returns <- list(result=sens.results, result.boots=NULL, outcome.type=outcome.type,
                    link=link, hazard.model=hazard.model,
                    n=nrow(data), confounder=confounder, exposure=exposure,
                    bound=bound, bound.sigma2=bound.sigma2)
  } else {
    returns <- list(result=sens.results, result.boots=sens.results.boots, outcome.type=outcome.type,
                    link=link, hazard.model=hazard.model,
                    n=nrow(data), confounder=confounder, exposure=exposure,
                    bound=bound, bound.sigma2=bound.sigma2)
  }
  structure(returns, class=c("GSAMU", class(returns)))
  # invisible(x=returns)
}
