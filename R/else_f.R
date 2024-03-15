##
biascal_new <- function(invcor, Amat, dir, val, k, p, Q){
  aprime <- invcor

  ## Max
  mycop_max <- cop(f=linfun(a=aprime, name="Obj.function"),
                   lc=lincon(A=Amat, dir=dir, val=val, name=seq(1,(k+p)*2)),
                   qc=quadcon(Q=Q, d=0, dir="<=", val=1, use=TRUE),
                   max=TRUE)
  ## Min
  mycop_min <- cop(f=linfun(a=aprime, name="Obj.function"),
                   lc=lincon(A=Amat, dir=dir, val=val, name=seq(1,(k+p)*2)),
                   qc=quadcon(Q=Q, d=0, dir="<=", val=1, use=TRUE),
                   max=FALSE)

  res_max <- solvecop(op=mycop_max, solver="alabama", quiet=TRUE)
  res_min <- solvecop(op=mycop_min, solver="alabama", quiet=TRUE)
  maxbias <- validate(op=mycop_max, sol=res_max, quiet=TRUE)$obj.fun
  minbias <- validate(op=mycop_min, sol=res_min, quiet=TRUE)$obj.fun

  return(c(minbias, maxbias))
}

# ##
# makedat_plot <- function(sens.result) {
#   object <- sens.result$result
#   object$outcome.type <- sens.result$outcome.type
#   object$link <- sens.result$link
#   object
# }

##
report.f <- function(dd, index){
  results.mat <- dd[[index]]
  results.mat$No <- index
  return(results.mat)
}
