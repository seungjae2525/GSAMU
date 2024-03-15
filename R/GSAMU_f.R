##
GSAMU.binary <- function(data, fitmodel,
                         delta.range, delta.diff, k, p, bound){

  ##############################################################################
  ## calculate correlation matrix
  corrmatrix <- cor(data[,-1])
  invcor <- solve(corrmatrix)

  ## Coefficients
  beta <- fitmodel$coefficients[-1] # k+p vector

  ## Set the candidate of delta
  if (delta.range[2] >= 0) {
    delta.range <- sort(delta.range, decreasing=FALSE)
  } else {
    delta.range <- sort(delta.range, decreasing=TRUE)
  }
  alldelta <- c(0, seq(delta.range[1], delta.range[2], by=delta.diff))
  alldelta <- sort(alldelta)

  ## Probit or Logit
  if (fitmodel$family$link == "logit") {
    qq <- pi/8
  } else if (fitmodel$family$link == "probit") {
    qq <- 1
  } else {
    stop("Link function must be \"logit\" or \"probit\"")
  }

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
    a <- c(sum(beta[-c(1:k)]), -delta.temp*colSums(invcor[-c(1:k),]))
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

  # df_all$label <- factor(df_all$label, levels=unique(df_all$label))
  return(df_all)
}

##
GSAMU.count.cox <- function(data, fitmodel,
                            delta.range, delta.diff, k, p, bound){
  ##############################################################################
  ## Coefficients of cox or glm model
  if (is(fitmodel)[1] == "coxph") {
    # Time to event outcome (time, status)
    # h*(L,Xs+1,X-s) - h*(L,Xs,X-s)
    beta <- fitmodel$coefficients# k+p vector

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
  if (delta.range[2] >= 0) {
    delta.range <- sort(delta.range, decreasing=FALSE)
  } else {
    delta.range <- sort(delta.range, decreasing=TRUE)
    delta.diff <- -delta.diff
  }
  alldelta <- c(0, seq(delta.range[1], delta.range[2], by=delta.diff))
  alldelta <- sort(alldelta)

  ##############################################################################
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

  res_max <- solvecop(op=mycop_max, solver="alabama", quiet=TRUE)
  res_min <- solvecop(op=mycop_min, solver="alabama", quiet=TRUE)
  maxbias <- validate(op=mycop_max, sol=res_max, quiet=TRUE)$obj.fun
  minbias <- validate(op=mycop_min, sol=res_min, quiet=TRUE)$obj.fun

  Joint <- c(minbias, maxbias)

  resultdf <- as.data.frame(t(cbind(singlebias, Joint)))

  coeff <- beta[-c(1:k)]
  label <- c(names(coeff), "Joint effect")

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

  return(df_all)
}
