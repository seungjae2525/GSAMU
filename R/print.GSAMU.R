#' @title Print for \code{GSAMU} objects
#'
#' @description Print the sensitivity analysis results for object of class \code{GSAMU}.
#'
#' @param x An object for class \code{GSAMU}.
#' @param digits Print digits. Default: max(1L, getOption("digits") - 3L).
#' @param ... Further arguments (currently not used).
#'
#' @details Print the sensitivity analysis results for object (\code{GSAMU}) of class \code{GSAMU}.
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
#'            # sex   race   age
#'            0.0, 0.0, 0.0,
#'            # Retinyl palmitate    Retinol  trans-beta-carotene
#'            rep(-0.19, 3),
#'            # alpha-Tocopherol  gamma-Tocopherol
#'            rep(-0.19, 2),
#'
#'            ## upper bound
#'            # sex   race   age
#'            0.32, 0.5, 0.29,
#'            # Retinyl palmitate    Retinol  trans-beta-carotene
#'            rep(0.15, 3),
#'            # alpha-Tocopherol  gamma-Tocopherol
#'            rep(0.10, 2))
#'
#' ## Sensitivity analysis
#' binary.re0 <- GSAMU(data=data_r4, outcome="triglyceride_b", outcome.type="binary",
#'                     link="logit", hazard.model=NULL,
#'                     confounder=c("sex", "race", "age"),
#'                     exposure=c("Retinyl palmitate", "Retinol", "trans-beta-carotene",
#'                                "alpha-Tocopherol", "gamma-Tocopherol"),
#'                     delta=c(0.11, 0.22, 0.33, 0.44), bound=bound,
#'                     bootsCI=FALSE, B=1000, seed=231111, verbose=TRUE)
#' print(binary.re0)
#'
#' \dontrun{
#' ## Sensitivity analysis with bootstrap percentile confidence interval
#' binary.re1 <- GSAMU(data=data_r4, outcome="triglyceride_b", outcome.type="binary",
#'                     link="logit", hazard.model=NULL,
#'                     confounder=c("sex", "race", "age"),
#'                     exposure=c("Retinyl palmitate", "Retinol", "trans-beta-carotene",
#'                                "alpha-Tocopherol", "gamma-Tocopherol"),
#'                     delta=c(0.11, 0.22, 0.33, 0.44), bound=bound,
#'                     bootsCI=TRUE, B=1000, seed=231111, verbose=TRUE)
#' print(binary.re1)
#' }
#'
#' @keywords Print
#'
#' @seealso
#'  \code{\link[GSAMU]{GSAMU}}, \code{\link[GSAMU]{GSAMU.alt}}
#'
#' @export

print.GSAMU <- function (x, digits=max(1L, getOption("digits") - 3L), ...){

  if (!inherits(x, "GSAMU")){
    stop("Argument 'x' must be an object of class \"GSAMU\".")
  }

  xx <- x

  cat("\n")
  cat("Data information: \n")
  cat("Sample size =", xx$n, "\n")
  cat("Number of confounders =", length(xx$confounder), "\n")
  cat("Number of exposures =", length(xx$exposure), "\n")
  cat("Outcome type =", ifelse(xx$outcome.type == "binary", "Binary",
                               ifelse(xx$outcome.type == "count", "Count",
                                      ifelse(xx$outcome.type == "continuous", "Continuous", "Time to event"))), "\n")

  if (xx$outcome.type %in% c("count", "binary")) {
    cat("Link function =", xx$link, "\n")
  } else if (xx$outcome.type == "continuous") {
    cat("Link function =", "identity", "\n")
  } else {
    cat("Hazard model =", xx$hazard.model, "\n")
  }

  ##
  if (is.null(xx$bound.sigma2)) {
    AA <- matrix(xx$bound, nrow=length(c(xx$confounder, xx$exposure)), ncol=2, byrow=FALSE)
    rownames(AA) <- c(paste0("phi_", "(", xx$confounder, ")"),
                      paste0("rho_", "(", xx$exposure, ")"))
    colnames(AA) <- c("Lower bound", "Upper bound")
    cat("\n")
    cat("Range of sensitivity paramters: \n")
    cat("For confounders \n")
    print(AA[c(1:length(xx$confounder)), ])
    cat("\n")
    cat("For exposures \n")
    print(AA[c((1+length(xx$confounder)):nrow(AA)), ])
  } else {
    AA <- matrix(xx$bound, nrow=length(c(xx$exposure)), ncol=2, byrow=FALSE)
    rownames(AA) <- c(paste0("rho_", "(", xx$exposure, ")"))
    colnames(AA) <- c("Lower bound", "Upper bound")
    cat("\n")
    cat("Range of sensitivity paramters: \n")
    cat("For exposures \n")
    print(AA[c(1:nrow(AA)), ])

    cat("\n")
    cat("Range of variance of distribution of U | (L,X): \n")
    cat(paste0("[", xx$bound.sigma2[1], ", ", xx$bound.sigma2[2], "]"))
  }

  ##
  delta <- unique(xx$result$delta)
  which.0.delta <- which(delta == 0)

  ##
  if (xx$outcome.type == "count") {
    effect.name <- "Conditional log-relative ratio"
  } else if (xx$outcome.type == "timetoevent") {
    if (xx$hazard.model == "coxph") {
      effect.name <- "Conditional log-hazard ratio"
    } else {
      effect.name <- "Conditional hazard difference"
    }
  } else {
    if (xx$link == "logit") {
      effect.name <- "Conditional log-odds ratio"
    } else {
      effect.name <- "Conditional exposure effects"
    }
  }

  ##
  cat("\n")
  cat("Sensitivity analysis results: \n")
  unique.label <- unique(xx$result$label)
  if (is.null(xx$result.boots)) {
    for (i in 1:length(unique.label)) {
      cat(paste0(i, ". ", c(xx$exposure, "Joint effect")[i], " \n"))
      xx.temp <- xx$result[xx$result$label == unique.label[i], ]
      xx.temp <- xx.temp[-which.0.delta, ] # remove row with delta=0
      tmp1 <- cbind(sprintf(paste0("%.", digits, "f"), xx.temp$coef),
                    paste0("[", sprintf(paste0("%.", digits, "f"), xx.temp$Lower.SI), ", ",
                           sprintf(paste0("%.", digits, "f"), xx.temp$Upper.SI), "]"))
      tmp1 <-data.frame(tmp1)
      colnames(tmp1) <- c(effect.name, 'Sensitivity interval')
      rownames(tmp1) <- c(paste0(rep(expression(delta), times=nrow(xx.temp)), '=', delta[-which.0.delta]))
      #
      if (max(delta[-which.0.delta]) < 0) {
        tmp1 <- tmp1[rev(1:nrow(tmp1)), , drop=FALSE]
      }
      print(tmp1)
      cat("\n")
    }
  } else {
    for (i in 1:length(unique.label)) {
      cat(paste0(i, ". ", c(xx$exposure, "Joint effect")[i], " \n"))
      xx.temp <- xx$result[xx$result$label == unique.label[i], ]
      xx.temp <- xx.temp[-which.0.delta, ] # remove row with delta=0
      tmp1 <- cbind(sprintf(paste0("%.", digits, "f"), xx.temp$coef),
                    paste0("[", sprintf(paste0("%.", digits, "f"), xx.temp$Lower.SI), ", ",
                           sprintf(paste0("%.", digits, "f"), xx.temp$Upper.SI), "]"),
                    paste0("[", sprintf(paste0("%.", digits, "f"), xx.temp$Lower.CI), ", ",
                           sprintf(paste0("%.", digits, "f"), xx.temp$Upper.CI), "]"))
      tmp1 <- data.frame(tmp1)
      colnames(tmp1) <- c(effect.name, 'Sensitivity interval', 'Confidence interval')
      rownames(tmp1) <- c(paste0(rep(expression(delta), times=nrow(xx.temp)), '=', delta[-which.0.delta]))

      print(tmp1)
      cat("\n")
    }
  }

  # invisible(x)
  invisible(NULL)
}
