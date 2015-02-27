phylocurve <- function(formula,tree,data,ymin=.01,ymax=.99,ylength=100,species_identifier="species",verbose=FALSE)
{
  if(!(species_identifier %in% colnames(data))) stop("Add species names column to data.")
  dat <- model.frame(formula,data=data)
  dat$species <- data$species
  data <- dat
  nspecies <- length(tree$tip.label)
  y_name <- colnames(data)[1]
  x_name <- colnames(data)[2]
  taxa <- tree$tip.label
  tip_coefficients <- matrix(NA,nrow = nspecies,ncol = 2,dimnames = list(taxa,c("(Intercept)",x_name)))
  if(verbose) cat("Phase 1: Fitting species curves: ")
  counter <- 0
  for(i in 1:nspecies)
  {
    if(verbose)
    {
      if((ceiling(i/nspecies*100)>counter))
      {
        cat(paste(counter,"%   ",sep=""))
        counter <- counter + 10
      }
      if(i==nspecies) cat(paste(100,"%   ",sep=""))

    }
    tip_coefficients[taxa[i],] <- coef(glm(formula,family = quasibinomial("logit"),data = data[data$species==taxa[i],]))
  }
  if(verbose) cat("\n")
  pgls_curve(tree = tree,tip_coefficients = tip_coefficients,ymin = ymin,ymax = ymax,ylength = ylength,verbose = verbose)
}

pgls_curve <- function(tree,tip_coefficients,varAY,vals_only=FALSE,ymin,ymax,ylength,verbose)
{
  if(is.null(rownames(tip_coefficients)))
  {
    warning("No row names for tip_coefficients. Assuming rows match up with tree tips.")
    rownames(tip_coefficients) <- tree$tip.label
  }
  if(missing(varAY))
  {
    D <- dist.nodes(tree)
    nspecies <- length(tree$tip.label)
    Nnode <- tree$Nnode
    MRCA <- mrca(tree, full = TRUE)
    M <- D[as.character(nspecies + 1), MRCA]
    dim(M) <- rep(sqrt(length(M)), 2)
    varAY <- M[(nspecies+1):(nspecies+Nnode), 1:nspecies]
    varA <- M[(nspecies+1):(nspecies+Nnode), (nspecies+1):(nspecies+Nnode)]
    colnames(varAY) <- tree$tip.label
  }
  nspecies <- length(tree$tip.label)
  Nnode <- tree$Nnode
  X <- matrix(1,nspecies,1)
  y <- seq(ymin,ymax,length=ylength)
  rownames(X) <- tree$tip.label
  tip_coefficients <- tip_coefficients[tree$tip.label,]
  Yinv <- t(apply(tip_coefficients,1,function(X) logit_inv(X,y)))
  #root <- as.double(solve(t(X)%*%solveC%*%X)%*%(t(X)%*%solveC%*%Yinv))
  if(verbose) cat("Phase 2: Reconstructing root curve","\n")
  root_three_point <- three.point.compute(phy = tree,P = X,Q = Yinv)
  root <- solve(root_three_point$PP) %*% root_three_point$QP
  root <- (matrix(rep(root,nspecies),length(root),nspecies))
  #anc_vals <- t(varAY%*%solveC%*%(Yinv-t(root)))+root[,1:Nnode]
  if(verbose) cat("Phase 3: Reconstructing internal node curves: ")
  #return(list(tree=tree,P=t(varAY),Q=(Yinv-t(root))))
  P <- t(varAY)
  Q <- (Yinv-t(root))
  var_BM <- apply(Q,2,function(X) t(three.point.compute(tree,X)$PP)/(nspecies-1))
  anc_three_point <- matrix(NA,ncol(Q),ncol(P))
  DD <- D1 <- matrix(NA,ncol(P),1)
  for(i in 1:ceiling(ncol(P)/50))
  {
    if(verbose) cat(paste(round(i/ceiling(ncol(P)/50)*100),"%   ",sep=""))
    temp <- three.point.compute(phy = tree,P = P[,(1+((i-1)*50)):min(50+(i-1)*50,ncol(P))],Q = Q)
    DD[(1+((i-1)*50)):min(50+(i-1)*50,ncol(P))] <- diag(temp$PP)
    D1[(1+((i-1)*50)):min(50+(i-1)*50,ncol(P))] <- temp$P1
    anc_three_point[,(1+((i-1)*50)):min(50+(i-1)*50,ncol(P))] <- temp$QP
  }
  vec11 <- temp$vec11
  var_nodes <- t((diag(varA) - DD + diag((1 - D1) %*% solve(vec11) %*% t(1 - D1))))
  var_nodes <- matrix(rep(var_nodes,ylength),ncol=ylength)
  var_BM <- matrix(rep(var_BM,tree$Nnode),ncol=tree$Nnode)
  CI <- t(var_nodes)*(var_BM)*1.96
  #anc_three_point <- three.point.compute(phy = tree,P = P,Q = Q)$QP
  anc_vals <- (anc_three_point)+root[,1:Nnode]
  
  if(vals_only) return(t(anc_vals))
  if(verbose) cat("\nPhase 4: Fitting ancesral curves")
  ret <- t(apply(anc_vals,2,function(X) coef(glm(y~X,family=quasibinomial("logit")))))
  rownames(ret) <- (nspecies+1):(nspecies+Nnode)
  rownames(CI) <- rownames(anc_vals) <- seq(ymin,ymax,length=ylength)
  colnames(CI) <- colnames(anc_vals) <- (nspecies+1):(nspecies+tree$Nnode)
  lower_CI <- anc_vals - CI
  upper_CI <- anc_vals + CI
  return(list(node_coefficients=ret,fitted_x=anc_vals,lower_CI_x=lower_CI,upper_CI_x=upper_CI,y_vals=seq(ymin,ymax,length=ylength)))
}

# logit function
# first parameter is the glm intercept
# second parameter the glm slope
logit_fx <- function(theta,x)
{
  b <- theta[1]
  m <- theta[2]
  exp(m*x+b)/(1+exp(m*x+b))
}

# inverse logit function
logit_inv <- function(theta,y)
{
  b <- theta[1]
  m <- theta[2]
  ((log(-y/(y-1))-b)/m)
}

# transform logit glm parameters to uncorrelated slope parameter
logit_slope_func <- function(theta)
{
  b <- theta[1]
  m <- theta[2]
  m/4
}

# transform logit glm parameters to uncorrelated EC50 parameter
logit_ec50_func <- function(theta)
{
  b <- theta[1]
  m <- theta[2]
  -b/m
}

# back-calculates glm intercept parameter
logit_b_func <- function(slope,ec50)
{
  m <- slope*4
  b <- -ec50 * m
  return(b)
}

# back-calculates glm slope parameter
logit_m_func <- function(slope,ec50)
{
  m <- slope*4
  return(m)
}

simcurves <- function(nspecies = 30,x_length=20,nreps=20,startree=FALSE,lambda=1,seed)
{
  if(!missing(seed)) set.seed(seed)
  x <- seq(0,15,length=x_length) # environmental gradient
  tree <- pbtree(n=nspecies)
  if(!missing(lambda)) simtree <- rescale(tree,"lambda",lambda=lambda) else simtree <- tree
  if(startree) simtree <- starTree(tree$tip.label,rep(1,length(tree$tip.label)))
  
  logit_ec50_true <- fastBM(simtree,5,bounds=c(2,20),internal = TRUE)
  logit_slope_true <- fastBM(simtree,.5,bounds=c(.2,1),internal = TRUE)
  
  logit_ec50 <- logit_ec50_true[1:nspecies]
  logit_slope <- logit_slope_true[1:nspecies]
  
  logit_b <- logit_b_func(slope=logit_slope,ec50=logit_ec50)
  logit_m <- logit_m_func(slope=logit_slope,ec50=logit_ec50)
  logit_coefs <- cbind(logit_b,logit_m)
  logit_sim_coefs <- logit_coefs
  
  logit_y <- data.frame(species =   as.character(t(matrix(rep(simtree$tip.label,length(x)),nrow=nspecies))),
                        x =  rep(x,nspecies),y=rep(0,nspecies*length(x)))
  for(i in 1:nspecies)
  {
    logit_y[1:length(x)+length(x)*(i-1),3] <- logit_fx(c(logit_b[i],logit_m[i]),x)
  }
  ret <- list(data = logit_y,tree = simtree,
              true_coefs = data.frame(Intercept=logit_b_func(logit_slope_true,logit_ec50_true),slope=logit_m_func(logit_slope_true,logit_ec50_true),row.names = names(logit_ec50_true))
  )
  ret
}
