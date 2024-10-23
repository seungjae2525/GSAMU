##
GSAMU.binary <- function(data, fitmodel, exposure,
                         delta, k, p, bound){
  ##############################################################################
  ## calculate correlation matrix
  corrmatrix <- cor(data[,-1])
  invcor <- solve(corrmatrix)

  ## Coefficients
  beta <- fitmodel$coefficients[-1] # k+p vector

  ## Set the candidate of delta
  alldelta <- sort(c(0, delta))

  ## Probit or Logit
  if (fitmodel$family$link == "logit") {
    qq <- pi/8
  } else if (fitmodel$family$link == "probit") {
    qq <- 1
  } else {
    stop("Link function must be \"logit\" or \"probit\"")
  }

  ## Rearrange the bound (c(lower, upper) -> c(upper, lower))
  bound <- c(bound[(k+p+1):(k+p+k+p)], bound[1:(k+p)])

  ##############################################################################
  ##
  df_all <- c()
  for(ddd in 1:length(alldelta)){
    true_range <-c()
    delta.temp <- alldelta[ddd]

    for(s in 1:p){
      # linear constraints
      val <- c(1, bound)
      dir <- c(">=", rep("<=",(k+p)), rep(">=",(k+p)))
      Amat <- matrix(0, nrow=(2*(k+p)+1), ncol=((k+p)+1))
      Amat[-1,-1] <- as.matrix(rbind(diag(1,k+p), diag(1,k+p)))
      Amat[1,1] <- 1

      myQ3 <- invcor
      smallA <- matrix(0, nrow=(k+p+1), ncol=(k+p+1))
      smallA[1,1] <- 1
      smallA[-1,-1] <- (delta.temp^2*qq)*myQ3

      # linear objective function
      coef <- beta[(s+k)]

      a <- c(coef, -delta.temp*invcor[,(s+k)])

      mycop_max <- cop(f=linfun(a=a, name="Obj.function"),
                       lc=lincon(A=Amat, dir=dir, val=val, name=seq(1,dim(Amat)[1])),
                       qc=quadcon(smallA, d=0, dir="<=", val=(1+delta.temp^2*qq), use=TRUE),
                       max=TRUE)
      mycop_max$qc[[1]] <- quadcon(smallA, d=0, dir="<=", val=(1+delta.temp^2*qq), id=paste0("a",1:(k+p+1)))
      mycop_max$qc[[2]] <- quadcon(-smallA, d=0, dir="<=", val=-(1+delta.temp^2*qq), id=paste0("b",1:(k+p+1)))

      mycop_min <- cop(f=linfun(a=a, name="Obj.function"),
                       lc=lincon(A=Amat, dir=dir, val=val, name=seq(1,dim(Amat)[1])),
                       qc=quadcon(smallA, d=0, dir="<=", val=(1+delta.temp^2*qq), use=TRUE),
                       max=FALSE)
      mycop_min$qc[[1]] <- quadcon(smallA, d=0, dir="<=", val=(1+delta.temp^2*qq), id=paste0("a",1:(k+p+1)))
      mycop_min$qc[[2]] <- quadcon(-smallA, d=0, dir="<=", val=-(1+delta.temp^2*qq), id=paste0("b",1:(k+p+1)))

      res_max <- solvecop(mycop_max, quiet=TRUE)
      res_min <- solvecop(mycop_min, quiet=TRUE)
      maxbias <- validate(mycop_max, res_max, quiet=TRUE)$obj.fun
      minbias <- validate(mycop_min, res_min, quiet=TRUE)$obj.fun

      true_range <- rbind(true_range, c(minbias, maxbias))
    }

    ## for joint effect
    if (p == 1) {
      # For single exposure
      a <- c(sum(beta[-c(1:k)]), -delta.temp*invcor[-c(1:k),])
    } else {
      # For multiple exposures
      a <- c(sum(beta[-c(1:k)]), -delta.temp*colSums(invcor[-c(1:k),]))
    }

    mycop_max <- cop(f=linfun(a=a, name="Obj.function"),
                     lc=lincon(A=Amat, dir=dir, val=val, name=seq(1,dim(Amat)[1])),
                     qc=quadcon(smallA, d=0, dir="<=", val=(1+delta.temp^2*qq), use=TRUE),
                     max=TRUE)
    mycop_max$qc[[1]] <- quadcon(smallA, d=0, dir="<=", val=(1+delta.temp^2*qq), id=paste0("a",1:(k+p+1)))
    mycop_max$qc[[2]] <- quadcon(-smallA, d=0, dir="<=", val=-(1+delta.temp^2*qq), id=paste0("b",1:(k+p+1)))

    mycop_min <- cop(f=linfun(a=a, name="Obj.function"),
                     lc=lincon(A=Amat, dir=dir, val=val, name=seq(1,dim(Amat)[1])),
                     qc=quadcon(smallA, d=0, dir="<=", val=(1+delta.temp^2*qq), use=TRUE),
                     max=FALSE)
    mycop_min$qc[[1]] <- quadcon(smallA, d=0, dir="<=", val=(1+delta.temp^2*qq), id=paste0("a",1:(k+p+1)))
    mycop_min$qc[[2]] <- quadcon(-smallA, d=0, dir="<=", val=-(1+delta.temp^2*qq), id=paste0("b",1:(k+p+1)))

    res_max <- solvecop(mycop_max, quiet=TRUE)
    res_min <- solvecop(mycop_min, quiet=TRUE)
    maxbias <- validate(mycop_max, res_max, quiet=TRUE)$obj.fun
    minbias <- validate(mycop_min, res_min, quiet=TRUE)$obj.fun

    true_range <- rbind(true_range, c(minbias, maxbias))

    df_temp <- as.data.frame(cbind(true_range))
    df_temp$label <- c(colnames(data)[(1+k+1):(1+k+p)], "Joint effect")
    df_temp$coef <- c(beta[-c(1:(k))], sum(beta[-c(1:(k))]))
    colnames(df_temp) <- c("Lower.SI", "Upper.SI", "label", "coef")
    df_temp$delta <- rep(delta.temp, dim(df_temp)[1])
    df_temp <- df_temp[,c("label", "delta", "coef", "Lower.SI", "Upper.SI")]

    df_all <- rbind(df_all, df_temp)
  }


  df_all2 <- c()
  exposure2 <- c(names(beta[-c(1:k)]), "Joint effect")
  for (i in 1:length(exposure2)) {
    df_temp <- df_all[df_all$label == exposure2[i], ]
    df_all2 <- rbind(df_all2, df_temp)
  }

  df_all <- df_all2
  rownames(df_all) <- NULL

  # df_all$label <- factor(df_all$label, levels=unique(df_all$label))
  return(df_all)
}

##
GSAMU.count.hazard <- function(data, fitmodel, exposure,
                               delta, k, p, bound){
  ##############################################################################
  ## Coefficients of cox or glm model
  if (is(fitmodel)[1] == "coxph") {
    # Time to event outcome (time, status)
    # h*(L,Xs+1,X-s) - h*(L,Xs,X-s)
    beta <- fitmodel$coefficients # k+p vector

    ## calculate correlation matrix
    corrmatrix <- cor(data[,-c(1,2)])
    invcor <- solve(corrmatrix)
  } else if (is(fitmodel)[1] == "aalen") {
    # Time to event outcome (time, status)
    # h*(L,Xs+1,X-s) - h*(L,Xs,X-s)
    # beta <- fitmodel$coef # k+p vector; using af
    beta <- coef(fitmodel, digits=Inf)[,1] # k+p vector; using aalen

    ## calculate correlation matrix
    corrmatrix <- cor(data[,-c(1,2)])
    invcor <- solve(corrmatrix)
  } else {
    # Count outcome
    # h*(L,Xs+1,X-s) - h*(L,Xs,X-s)
    beta <- fitmodel$coefficients[-1] # k+p vector

    ## calculate correlation matrix
    corrmatrix <- cor(data[,-1])
    invcor <- solve(corrmatrix)
  }

  ##############################################################################
  ## Set the candidate of delta
  alldelta <- sort(c(0, delta))

  ##############################################################################
  ## Rearrange the bound (c(lower, upper) -> c(upper, lower))
  bound <- c(bound[(k+p+1):(k+p+k+p)], bound[1:(k+p)])

  ##
  val <- bound
  dir <- c(rep("<=", (k+p)), rep(">=", (k+p)))
  Amat <- as.matrix(rbind(diag(1,(k+p)), diag(1,(k+p))))

  myQ3 <- invcor

  ##
  singlebias <- apply(invcor, 1, function(x) {
    biascal_new(invcor=x, Amat=Amat, dir=dir, val=val, k=k, p=p, Q=myQ3)
  })

  ## Joint bias
  forjointmyq3 <- myQ3[,-c(1:k)]
  if (p == 1) {
    # Max
    mycop_max <- cop(f=linfun(a=forjointmyq3, name="Obj.function"),
                     lc=lincon(A=Amat, dir=dir, val=val, name=seq(1, (k+p)*2)),
                     qc=quadcon(Q=myQ3, d=0, dir="<=", val=1, use=TRUE),
                     max=TRUE)
    # Min
    mycop_min <- cop(f=linfun(a=forjointmyq3, name="Obj.function"),
                     lc=lincon(A=Amat, dir=dir, val=val, name=seq(1, (k+p)*2)),
                     qc=quadcon(Q=myQ3, d=0, dir="<=", val=1, use=TRUE),
                     max=FALSE)
  } else {
    # Max
    mycop_max <- cop(f=linfun(a=rowSums(forjointmyq3), name="Obj.function"),
                     lc=lincon(A=Amat, dir=dir, val=val, name=seq(1, (k+p)*2)),
                     qc=quadcon(Q=myQ3, d=0, dir="<=", val=1, use=TRUE),
                     max=TRUE)
    # Min
    mycop_min <- cop(f=linfun(a=rowSums(forjointmyq3), name="Obj.function"),
                     lc=lincon(A=Amat, dir=dir, val=val, name=seq(1, (k+p)*2)),
                     qc=quadcon(Q=myQ3, d=0, dir="<=", val=1, use=TRUE),
                     max=FALSE)
  }

  res_max <- solvecop(op=mycop_max, solver="alabama", quiet=TRUE)
  res_min <- solvecop(op=mycop_min, solver="alabama", quiet=TRUE)
  maxbias <- validate(op=mycop_max, sol=res_max, quiet=TRUE)$obj.fun
  minbias <- validate(op=mycop_min, sol=res_min, quiet=TRUE)$obj.fun

  Joint <- c(minbias, maxbias)

  resultdf <- as.data.frame(t(cbind(singlebias, Joint)))

  coeff <- beta[-c(1:k)]

  # label <- c(names(coeff), "Joint effect")
  label <- c(exposure, "Joint effect")

  jointeffcoef <- sum(coeff)

  if (k==0) {
    resultdf <- resultdf
  } else {
    resultdf <- resultdf[-c(1:k),]
  }

  df.init <- data.frame(cbind(label, resultdf, c(coeff, jointeffcoef)))
  colnames(df.init) <- c("label", "bias_low", "bias_upper", "coef")

  for (i in 1:length(alldelta)) {
    tempc0 <- alldelta[i]

    df_temp <- df.init

    c_1 <- df_temp$coef - tempc0*df_temp$bias_upper
    c_2 <- df_temp$coef - tempc0*df_temp$bias_low

    df_temp$cond_min <- pmin(c_1, c_2)
    df_temp$cond_max <- pmax(c_1, c_2)

    df_temp$c0 <- rep(tempc0, dim(df_temp)[1])

    if (i == 1) {
      df <- df_temp
    } else {
      df <- rbind(df, df_temp)
    }
  }

  colnames(df) <- c("label", "bias_low", "bias_upper", "coef", "Lower.SI", "Upper.SI", "delta")
  df_all <- data.frame(df[,c("label", "delta", "coef", "Lower.SI", "Upper.SI")])
  rownames(df_all) <- NULL

  df_all2 <- c()
  exposure2 <- c(exposure, "Joint effect")
  for (i in 1:length(exposure2)) {
    df_temp <- df_all[df_all$label == exposure2[i], ]
    df_all2 <- rbind(df_all2, df_temp)
  }

  df_all <- df_all2
  rownames(df_all) <- NULL

  return(df_all)
}

##
GSAMU_alt <- function(data, fitmodel, exposure, delta, k, p,
                      bound, bound.sigma2, link, hazard.model){
  ##############################################################################
  ## Coefficients of cox or glm model
  if (is(fitmodel)[1] == "coxph") {
    # Time to event outcome (time, status)
    # h*(L,Xs+1,X-s) - h*(L,Xs,X-s)
    beta <- fitmodel$coefficients # k+p vector

    # ## calculate correlation matrix
    # corrmatrix <- cor(data[,-c(1,2)])
    # invcor <- solve(corrmatrix)
  } else if (is(fitmodel)[1] == "aalen") {
    # Time to event outcome (time, status)
    # h*(L,Xs+1,X-s) - h*(L,Xs,X-s)
    # beta <- fitmodel$coef # k+p vector; using af
    beta <- coef(fitmodel, digits=Inf)[,1] # k+p vector; using aalen

    # ## calculate correlation matrix
    # corrmatrix <- cor(data[,-c(1,2)])
    # invcor <- solve(corrmatrix)
  } else {
    # Count or binary outcome
    # h*(L,Xs+1,X-s) - h*(L,Xs,X-s)
    beta <- fitmodel$coefficients[-1] # k+p vector

    # ## calculate correlation matrix
    # corrmatrix <- cor(data[,-1])
    # invcor <- solve(corrmatrix)
  }

  ##############################################################################
  ## Set the candidate of delta
  alldelta <- sort(c(0, delta))

  ## Upper and lower bounds
  lower_bound <- bound[1:p]
  upper_bound <- bound[(p+1):length(bound)]
  resultdf <- matrix(c(lower_bound, sum(lower_bound), upper_bound, sum(upper_bound)),
                     nrow=p+1, ncol=2, byrow=FALSE)
  rownames(resultdf) <- c(exposure, "Joint effect")

  ##
  coef <- beta[c((k+1):(k+p))]
  jointeffcoef <- sum(coef)
  label <- c(exposure, "Joint effect")
  df.init <- data.frame(label, resultdf, c(coef, jointeffcoef))
  colnames(df.init) <- c("label", "bias_low", "bias_upper", "coef")

  ##
  if (!is.null(link) & is.null(hazard.model)) {
    ## binary outcome
    if (link == "logit") {
      qq.min <- sqrt(1 + pi/8 * alldelta^{2} * min(bound.sigma2))
      qq.max <- sqrt(1 + pi/8 * alldelta^{2} * max(bound.sigma2))
    } else if (link == "probit") {
      qq.min <- sqrt(1 + 1 * alldelta^{2} * min(bound.sigma2))
      qq.max <- sqrt(1 + 1 * alldelta^{2} * max(bound.sigma2))
    }
  } else {
    ## other outcomes
    qq.min <- rep(1, length(alldelta))
    qq.max <- rep(1, length(alldelta))
  }

  ##############################################################################
  ##
  df <- c()
  for(ddd in 1:length(alldelta)){
    delta.temp <- alldelta[ddd]

    df_temp <- df.init

    c_1 <- df_temp$coef * qq.min[ddd] - delta.temp*df_temp$bias_low
    c_2 <- df_temp$coef * qq.max[ddd] - delta.temp*df_temp$bias_low
    c_3 <- df_temp$coef * qq.min[ddd] - delta.temp*df_temp$bias_upper
    c_4 <- df_temp$coef * qq.max[ddd] - delta.temp*df_temp$bias_upper

    df_temp$cond_min <- pmin(c_1, c_2, c_3, c_4)
    df_temp$cond_max <- pmax(c_1, c_2, c_3, c_4)

    df_temp$c0 <- rep(delta.temp, dim(df_temp)[1])

    if (ddd == 1) {
      df <- df_temp
    } else {
      df <- rbind(df, df_temp)
    }
  }

  colnames(df) <- c("label", "bias_low", "bias_upper", "model_output", "cond_min", "cond_max", "delta")
  sensresult <- data.frame(df[,c("label", "model_output", "cond_min", "cond_max", "delta")])
  rownames(sensresult) <- NULL

  ##
  df_all <- c()
  exposure2 <- c(exposure, "Joint effect")
  for (i in 1:length(exposure2)) {
    df_temp <- sensresult[sensresult$label == exposure2[i], ]
    df_all <- rbind(df_all, df_temp)
  }

  df_all <- as.data.frame(df_all)[, c(1, 5, 2, 3, 4)]
  colnames(df_all) <- c("label", "delta", "coef", "Lower.SI", "Upper.SI")
  rownames(df_all) <- NULL

  # df_all$label <- factor(df_all$label, levels=unique(df_all$label))
  return(df_all)
}
