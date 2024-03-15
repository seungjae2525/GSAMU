#' @importFrom ggplot2 autoplot
#' @export
ggplot2::autoplot

#' @title Plot for sensitivity analysis results
#'
#' @param object An object for class \code{GSAMU}.
#' @param point.size Point size of the the lower and upper bound of the sensitivity range. Default: 1.4.
#' @param width.SI Width of horizon line which represents the sensitivity interval. Default: 2.75.
#' @param width.CI Width of horizon line which represents the confidence interval of the population sensitivity interval. Default: 1.55.
#' @param axis.title.x.size Size of x axis title. Default: 15.
#' @param axis.text.size Size of x and y axis text. Default: 12.
#' @param legend.text.size label text size. Default: 4.
#' @param myxlim The minimum and maximum values of x-axis.
#' @param ... Further arguments (currently not used).
#'
#' @examples
#' ############################################
#' ## Real data analysis for binary outcome ##
#' ############################################
#' ### Log transformed of exposures. continuous variables are standardized
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
#' # ## Scaled dataset
#' # data_r4 <- usedata
#' # data_r4[,c(2:9)] <- scale(data_r4[,c(2:9)])
#'
#' # change delta up to 0.44
#' uprimeupperrange <- 0.44
#'
#' ## glm fitting (working regression model)
#' fit.glm <- glm(triglyceride_b~., family=binomial(link="logit"), data=data_r4)
#' summary(fit.glm)
#' round(cbind("Est"=exp(coef(fit.glm)), exp(confint(fit.glm)),
#'       "P value"=summary(fit.glm)$coefficients[,4]), 3)
#'
#' ##
#' beta <- glm(triglyceride_b~., family=binomial(link="logit"), data=data_r4)$coefficient
#'
#' ##
#' k <- 3
#' p <- 5
#'
#' ##
#' bound <- c(## upper bound
#'   # male   race   age
#'   0.32, 0.5, 0.29,
#'   # Retinyl palmitate    Retinol  trans-beta-carotene
#'   rep(0.15, 3),
#'   # alpha-Tocopherol  gamma-Tocopherol
#'   rep(0.10, 2),
#'
#'   ## lower bound
#'   # male   race   age
#'   0.0, 0.0, 0.0,
#'   # Retinyl palmitate    Retinol  trans-beta-carotene
#'   rep(-0.19, 3),
#'   # alpha-Tocopherol  gamma-Tocopherol
#'   rep(-0.19, 2))
#'
#' ## Sensitivity analysis
#' binary.re0 <- GSAMU(data=data_r4, outcome="triglyceride_b", outcome.type="binary", link="logit",
#'                     confounder=c("sex", "race", "age"),
#'                     exposure=c("Retinyl palmitate", "Retinol", "trans-beta-carotene",
#'                                "alpha-Tocopherol", "gamma-Tocopherol"),
#'                     delta.range=c(0.11, 0.44), delta.diff=0.11, k=k, p=p, bound=bound,
#'                     bootsCI=FALSE, B=1000, seed=231111, verbose=TRUE,
#'                     report.result=TRUE, decimal.p=3)
#' autoplot(object=binary.re0, point.size=2.75, width.SI=1, width.CI=0.6,
#'          axis.title.x.size=15, axis.text.size=16, legend.text.size=15,
#'          myxlim=NULL)
#' \dontrun{
#' ## Sensitivity analysis with bootstrap percentile confidence interval
#' binary.re1 <- GSAMU(data=data_r4, outcome="triglyceride_b", outcome.type="binary", link="logit",
#'                     confounder=c("sex", "race", "age"),
#'                     exposure=c("Retinyl palmitate", "Retinol", "trans-beta-carotene",
#'                                "alpha-Tocopherol", "gamma-Tocopherol"),
#'                     delta.range=c(0.11, 0.44), delta.diff=0.11, k=k, p=p, bound=bound,
#'                     bootsCI=TRUE, B=1000, seed=231111, verbose=TRUE,
#'                     report.result=TRUE, decimal.p=3)
#' autoplot(object=binary.re1, point.size=2.75, width.SI=1.55, width.CI=0.6,
#'          axis.title.x.size=15, axis.text.size=16, legend.text.size=15,
#'          myxlim=NULL)
#' }
#'
#' @seealso
#'  \code{\link[GSAMU]{GSAMU}}
#'
#' @references
#' Lee S, Jeong B, Lee D, Lee W (2024+):
#' Sensitivity analysis for effects of multiple exposures in the presence of unmeasured confounding: non-Gaussian and time-to-event outcomes
#' \emph{xxx}. DOI: xxx.
#'
#' @importFrom ggplot2 autoplot
#'
#' @rdname autoplot.GSAMU
#'
#' @export
#'
autoplot.GSAMU <- function(object, point.size=2.75, width.SI=1.55, width.CI=0.6,
                           axis.title.x.size=15, axis.text.size=12, legend.text.size=12,
                           myxlim=NULL, ...){
  autoplot_GSAMU(sens.result=object, point.size=point.size, width.SI=width.SI, width.CI=width.CI,
                 axis.title.x.size=axis.title.x.size, axis.text.size=axis.text.size, legend.text.size=legend.text.size,
                 myxlim=myxlim, ...=...)
}

#########---------------------------------------------------------------------------------------------
autoplot_GSAMU <- function(sens.result, point.size=2.75, width.SI=1.55, width.CI=0.6,
                           axis.title.x.size=15, axis.text.size=12, legend.text.size=12,
                           myxlim=NULL, ...){
  outcome.type <- sens.result$outcome.type
  if (outcome.type == "count") {
    xaxis.name <- "Conditional log-relative ratio"
  } else if (outcome.type == "timetoevent") {
    xaxis.name <- "Conditional log-hazard ratio"
  } else {
    if (sens.result$link == "logit") {
      xaxis.name <- "Conditional log-odds ratio"
    } else {
      xaxis.name <- "Conditional single- and joint-exposure effects"
    }
  }

  ##
  sens.result <- sens.result$result

  ##
  if (ncol(sens.result) == 5) {
    ### Without bootstrap results
    ## Remove delta=0 results
    sens.result <- sens.result[sens.result$delta != 0, ]
  }
  ##
  alldelta <- unique(sens.result$delta)
  if (alldelta[1] == 0) {
    c0upper <- alldelta[length(alldelta)]
  } else {
    c0upper <- alldelta[1]
  }
  ##
  if (length(alldelta) == 1) {
    myColors <- c("#01579B")
  } else if (length(alldelta) == 2) {
    if (c0upper < 0) {
      myColors <- c("#01579B","#81D4FA")
    } else {
      myColors <- c("#81D4FA","#01579B")
    }
  } else if (length(alldelta) == 3) {
    if (c0upper < 0) {
      myColors <- c("#01579B","#29B6F6","#81D4FA")
    } else {
      myColors <- c("#81D4FA","#29B6F6","#01579B")
    }
  } else if (length(alldelta) == 4) {
    if (c0upper < 0) {
      myColors <- c("#01579B","#0288D1","#29B6F6","#81D4FA")
    } else {
      myColors <- c("#81D4FA","#29B6F6","#0288D1","#01579B")
    }
  } else if (length(alldelta) == 5) {
    if (c0upper < 0){
      myColors <- c("#01309b", "#01579B", "#0288D1", "#29B6F6", "#81D4FA")
    } else {
      myColors <- c("#81D4FA", "#29B6F6", "#0288D1", "#01579B", "#01309b")
    }
  } else {
    stop("Warnings!")
  }

  ## Reset the levels
  if (c0upper < 0) {
    names(myColors) <- levels(factor(sens.result$delta))
    trev <- levels(factor(sens.result$delta))
    sens.result$delta <- factor(sens.result$delta, levels=trev)
    sens.result$label <- factor(sens.result$label, levels=rev(levels(sens.result$label)))
  } else {
    names(myColors) <- levels(factor(sens.result$delta))
    trev <- rev(levels(factor(sens.result$delta)))
    sens.result$delta <- factor(sens.result$delta, levels=trev)
    sens.result$label <- factor(sens.result$label, levels=rev(levels(sens.result$label)))
  }

  # NOTE: nullify variables used for non-standard evaluation for tidyverse/ggplot2 below
  label <- delta  <- Lower.SI <- Upper.SI  <- Lower.CI <- Upper.CI <- NULL

  if (ncol(sens.result) == 5) {
    ### Without bootstrap results
    fp <- ggplot(sens.result) +
      geom_linerange(aes(x=coef, y=label, xmin=Lower.SI, xmax=Upper.SI, color=factor(delta)),
                     position=position_dodgev(height=0.7),
                     linewidth=width.SI) + # line width
      geom_point(aes(x=coef, y=label, color=factor(delta)), size=point.size, # point size
                 position=position_dodgev(height=0.7)) +
      scale_colour_manual(name=expression(bold(italic("\u03B4")~"(delta)")), values=myColors,
                          breaks=rev(levels(factor(sens.result$delta))),
                          guide=guide_legend(override.aes=list(shape=NA))) +
      # coord_cartesian(xlim=c(-0.25, 2)) +  # xlim
      geom_vline(xintercept=0, linetype="solid", linewidth=0.85, colour="black")+
      ggplot2::theme_minimal() +
      ggplot2::theme(rect=element_blank(),
                     text=element_text(colour="black", face="bold", family="serif")) %+replace%
      ggplot2::theme(axis.text.y=element_text(colour="black", size=axis.text.size, face="bold", family="serif"),
                     axis.text.y.right=element_text(hjust=1, size=axis.text.size, face="bold", family="serif"),
                     axis.text.x=element_text(colour="black", size=axis.text.size, face="bold", family="serif"),
                     axis.title.x=element_text(colour="black", size=axis.title.x.size, face="bold", family="serif"),
                     panel.border=element_blank(),
                     strip.text=element_text(hjust=0, face="bold", family="serif"),
                     strip.background=element_rect(colour="white", fill="aliceblue"),
                     panel.grid.major.x=element_line(colour="gray70", linewidth=0.8, linetype=2),
                     panel.grid.minor.x=element_line(colour="gray90", linewidth=0.7, linetype=2),
                     panel.grid.major.y=element_blank(),
                     panel.grid.minor.y=element_blank(),
                     legend.title=element_text(size=legend.text.size+1, face="bold", family="serif"),
                     legend.text=element_text(size=legend.text.size, face="bold", family="serif"),
                     legend.key.size=unit(0.75,"cm"),
                     legend.key.height=unit(0.55,"cm")) + # height between legends
      xlab(xaxis.name)+
      ylab("")
  } else {
    ## With bootstrap results
    fp <- ggplot(sens.result) +
      geom_linerange(aes(x=coef, y=label, xmin=Lower.CI, xmax=Upper.CI, color=factor(delta)),
                     position=position_dodgev(height=0.7),
                     linewidth=width.CI) + # line width
      guides(colour=guide_legend(override.aes=list(size=2, # legend point size
                                                   linewidth=1))) + # legend line width
      geom_linerange(aes(x=coef, y=label, xmin=Lower.SI, xmax=Upper.SI, color=factor(delta)),
                     position=position_dodgev(height=0.7),
                     linewidth=width.SI) + # line width
      geom_point(aes(x=coef, y=label, color=factor(delta)), size=point.size, # point size
                 position=position_dodgev(height=0.7)) +
      scale_colour_manual(name=expression(bold(italic("\u03B4")~"(delta)")), values=myColors,
                          breaks=rev(levels(factor(sens.result$delta))),
                          guide=guide_legend(override.aes=list(shape=NA))) +
      geom_vline(xintercept=0, linetype="solid", linewidth=0.85, colour="black")+
      ggplot2::theme_minimal() +
      ggplot2::theme(rect=element_blank(),
                     text=element_text(colour="black", face="bold", family="serif")) %+replace%
      ggplot2::theme(axis.text.y=element_text(colour="black", size=axis.text.size, face="bold", family="serif"),
                     axis.text.y.right=element_text(hjust=1, size=axis.text.size, face="bold", family="serif"),
                     axis.text.x=element_text(colour="black", size=axis.text.size, face="bold", family="serif"),
                     axis.title.x=element_text(colour="black", size=axis.title.x.size, face="bold", family="serif"),
                     panel.border=element_blank(),
                     strip.text=element_text(hjust=0, face="bold", family="serif"),
                     strip.background=element_rect(colour="white", fill="aliceblue"),
                     panel.grid.major.x=element_line(colour="gray70", linewidth=0.8, linetype=2),
                     panel.grid.minor.x=element_line(colour="gray90", linewidth=0.7, linetype=2),
                     panel.grid.major.y=element_blank(),
                     panel.grid.minor.y=element_blank(),
                     legend.title=element_text(size=legend.text.size+1, face="bold", family="serif"),
                     legend.text=element_text(size=legend.text.size, face="bold", family="serif"),
                     legend.key.size=unit(0.75,"cm"),
                     legend.key.height=unit(0.55,"cm")) + # height between legends
      xlab(xaxis.name)+
      ylab("")
  }

  if (!is.null(myxlim)) {
    fp <- fp + coord_cartesian(xlim=myxlim)
  }

  suppressWarnings(print(fp))
  invisible(x=fp)
}


