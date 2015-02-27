phylocurve <- function(formula,tree,data,ymin=.01,ymax=.99,ylength=30,tip_coefficients,species_identifier="species",verbose=FALSE)
{
  if(missing(tip_coefficients))
  {
    if(!(species_identifier %in% colnames(data))) stop("Add species names column to data.")
    tip_coefficients <- get_tip_coefficients(formula = formula,tree = tree,data = data,ymin = ymin,ymax = ymax,ylength = ylength,species_identifier = species_identifier,verbose = verbose)
  } else if(verbose) cat("Phase 1: Tip estimates already provided\n")
  pgls_curve(tree = tree,tip_coefficients = tip_coefficients,ymin = ymin,ymax = ymax,ylength = ylength,verbose = verbose)
}

get_tip_coefficients <- function(formula,tree,data,ymin=.01,ymax=.99,ylength=30,species_identifier="species",verbose=FALSE)
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
  tip_coefficients
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
  return(list(node_coefficients=ret,fitted_x=anc_vals,lower_CI_x=lower_CI,upper_CI_x=upper_CI,y_vals=seq(ymin,ymax,length=ylength),tip_coefficients=tip_coefficients))
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

simcurves <- function(nspecies = 30,x_length=20,startree=FALSE,lambda=1,seed)
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

# wrapper for physignal
phylocurve.signal <- function(tip_coefficients,tree,ymin=.01,ymax=.99,ylength=30,iter=1000)
{
  x <- tip_coefficients
  x <- t(apply(x,1,function(X) logit_inv(X,seq(ymin,ymax,length=ylength))))
  physignal_no_plot(tree,x,iter)
}

# faster version of physignal from geomorph package by avoiding multiple calls to solve()
# also avoids plotting
physignal_no_plot <- function (phy, A, iter = 249, method = c("Kmult", "SSC")) 
{
  method <- match.arg(method)
  if (any(is.na(A)) == T) {
    stop("Data matrix contains missing values. Estimate these first (see 'estimate.missing').")
  }
  if (length(dim(A)) == 3) {
    if (is.null(dimnames(A)[[3]])) {
      stop("Data matrix does not include taxa names as dimnames for 3rd dimension.")
    }
    x <- two.d.array(A)
  }
  if (length(dim(A)) == 2) {
    if (is.null(rownames(A))) {
      stop("Data matrix does not include taxa names as dimnames for rows.")
    }
    x <- A
  }
  if (class(phy) != "phylo") 
    stop("tree must be of class 'phylo.'")
  if (!is.binary.tree(phy)) 
    stop("tree is not fully bifurcating.")
  N <- length(phy$tip.label)
  if (N != dim(x)[1]) {
    stop("Number of taxa in data matrix and tree are not not equal.")
  }
  if (length(match(rownames(x), phy$tip.label)) != N) 
    stop("Data matrix missing some taxa present on the tree.")
  if (length(match(phy$tip.label, rownames(x))) != N) 
    stop("Tree missing some taxa in the data matrix.")
  if (any(is.na(match(sort(phy$tip.label), sort(rownames(x))))) == 
        T) {
    stop("Names do not match between tree and data matrix.")
  }
  x <- x[phy$tip.label, ]
  if (is.null(dim(x)) == TRUE) {
    x <- matrix(x, dimnames = list(names(x)))
  }
  if (method == "Kmult") {
    Kmult <- function(x, phy,C,solveC,D.mat,K.denom) {
      x <- as.matrix(x)
      N <- length(phy$tip.label)
      ones <- array(1, N)
      C <- C[row.names(x), row.names(x)]
      solveC <- solveC[row.names(x), row.names(x)]
      D.mat <- D.mat[row.names(x), row.names(x)]
      a.obs <- colSums(solveC) %*% x/sum(solveC)
      distmat <- as.matrix(dist(rbind(as.matrix(x), a.obs)))
      MSEobs.d <- sum(distmat[(1:N), (N + 1)]^2)
      dist.adj <- as.matrix(dist(rbind((D.mat %*% (x - (ones %*% a.obs))), 0)))
      MSE.d <- sum(dist.adj[(1:N), (N + 1)]^2)
      K.stat <- (MSEobs.d/MSE.d)/K.denom
      return(K.stat)
    }
    C <- vcv.phylo(phy)
    solveC <- solve(C)
    eigC <- eigen(C)
    D.mat <- solve(eigC$vectors %*% diag(sqrt(eigC$values)) %*% 
                     t(eigC$vectors))
    colnames(D.mat) <- rownames(D.mat) <- colnames(C)
    ones <- array(1, N)
    K.denom <- (sum(diag(C)) - N * solve(t(ones) %*% 
                                           solveC %*% ones))/(N - 1)
    K.obs <- Kmult(x, phy,C,solveC,D.mat,K.denom)
    P.val <- 1
    K.val <- rep(0, iter)
    for (i in 1:iter) {
      x.r <- as.matrix(x[sample(nrow(x)), ])
      rownames(x.r) <- rownames(x)
      K.rand <- Kmult(x.r, phy,C,solveC,D.mat,K.denom)
      P.val <- ifelse(K.rand >= K.obs, P.val + 1, P.val)
      K.val[i] <- K.rand
    }
    P.val <- P.val/(iter + 1)
    K.val[iter + 1] = K.obs
    if (dim(x)[2] > 1) {
      #plotGMPhyloMorphoSpace(phy, A, ancStates = FALSE)
    }
    return(list(phy.signal = K.obs, pvalue = P.val))
  }
  if (method == "SSC") {
    anc.states <- NULL
    options(warn = -1)
    for (i in 1:ncol(x)) {
      tmp <- as.vector(fastAnc(phy, x[, i]))
      anc.states <- cbind(anc.states, tmp)
    }
    colnames(anc.states) <- NULL
    dist.mat <- as.matrix(dist(rbind(as.matrix(x), as.matrix(anc.states)))^2)
    SSC.o <- 0
    for (i in 1:nrow(phy$edge)) {
      SSC.o <- SSC.o + dist.mat[phy$edge[i, 1], phy$edge[i, 
                                                         2]]
    }
    P.val <- 1
    SSC.val <- rep(0, iter)
    for (ii in 1:iter) {
      x.r <- x[sample(nrow(x)), ]
      if (is.null(dim(x.r)) == TRUE) {
        x.r <- matrix(x.r)
      }
      row.names(x.r) <- row.names(x)
      anc.states.r <- NULL
      options(warn = -1)
      for (i in 1:ncol(x.r)) {
        tmp <- as.vector(fastAnc(phy, x.r[, i]))
        anc.states.r <- cbind(anc.states.r, tmp)
      }
      colnames(anc.states.r) <- NULL
      dist.mat.r <- as.matrix(dist(rbind(as.matrix(x.r), 
                                         as.matrix(anc.states.r)))^2)
      SSC.r <- 0
      for (i in 1:nrow(phy$edge)) {
        SSC.r <- SSC.r + dist.mat.r[phy$edge[i, 1], phy$edge[i, 
                                                             2]]
      }
      P.val <- ifelse(SSC.r <= SSC.o, P.val + 1, P.val)
      SSC.val[ii] <- SSC.r
    }
    P.val <- P.val/(iter + 1)
    SSC.val[iter + 1] = SSC.o
    if (dim(x)[2] > 1) {
      #plotGMPhyloMorphoSpace(phy, A, ancStates = FALSE)
    }
    return(list(phy.signal = SSC.o, pvalue = P.val))
  }
}

# wrapper for procD.pgls
phylocurve.pgls <- function(tip_coefficients,univariate_trait,tree,ymin=.01,ymax=.99,ylength=30,iter=1000)
{
  x <- tip_coefficients
  y <- t(apply(x,1,function(X) logit_inv(X,seq(ymin,ymax,length=ylength))))
  procd <- procD.pgls2(f1 = y~univariate_trait,phy = tree,iter = iter, y=y)
  #list(procd=procd,univariate_trait=univariate_trait)
  procd
}

procD.pgls2 <- function (f1, phy, iter = 999, int.first = FALSE, RRPP = FALSE, 
          verbose = FALSE,y) 
{
  data = NULL
  form.in <- formula(f1)
  if (int.first == TRUE) 
    ko = TRUE
  else ko = FALSE
  Terms <- terms(form.in, keep.order = ko)
  k <- length(attr(Terms, "term.labels"))
  mf <- model.frame(form.in)
  Y <- as.matrix(mf[1])
  if (length(dim(Y)) != 2) {
    stop("Response matrix (shape) not a 2D array. Use 'two.d.array' first.")
  }
  if (any(is.na(Y)) == T) {
    stop("Response data matrix (shape) contains missing values. Estimate these first (see 'estimate.missing').")
  }
  if (is.null(dimnames(Y)[[1]])) {
    stop("No species names with Y-data")
  }
  N <- length(phy$tip.label)
  if (length(match(rownames(Y), phy$tip.label)) != N) 
    stop("Data matrix missing some taxa present on the tree.")
  if (length(match(phy$tip.label, rownames(Y))) != N) 
    stop("Tree missing some taxa in the data matrix.")
  C <- vcv.phylo(phy)
  C <- C[rownames(Y), rownames(Y)]
  eigC <- eigen(C)
  lambda <- zapsmall(eigC$values)
  if (any(lambda == 0)) {
    warning("Singular phylogenetic covariance matrix. Proceed with caution")
    lambda = lambda[lambda > 0]
  }
  eigC.vect = eigC$vectors[, 1:(length(lambda))]
  Pcor <- solve(eigC.vect %*% diag(sqrt(lambda)) %*% t(eigC.vect))
  PY <- Pcor %*% Y
  Xs = geomorph:::mod.mats(mf)
  anova.parts.obs <- geomorph:::anova.pgls.parts(form.in, X = NULL, Pcor, 
                                      Yalt = "observed", keep.order = ko)
  anova.tab <- anova.parts.obs$table
  df <- anova.parts.obs$df[1:k]
  dfE <- anova.parts.obs$df[k + 1]
  P <- array(0, c(k, 1, iter + 1))
  P[, , 1] <- SS.obs <- anova.parts.obs$F[1:k]
  for (i in 1:iter) {
    if (RRPP == TRUE) {
      SS.ran <- geomorph:::SS.pgls.random(Y, Xs, SS = SS.obs, Pcor, 
                               Yalt = "RRPP")
    }
    else SS.ran <- geomorph:::SS.pgls.random(Y, Xs, Pcor, SS = SS.obs, 
                                  Yalt = "resample")
    SS.r <- SS.ran$SS
    Yr <- SS.ran$Y
    SSE.r <- SS.ran$SSE
    Fs.r <- (SS.r/df)/(SSE.r/dfE)
    P[, , i + 1] <- Fs.r
  }
  P.val <- geomorph:::Pval.matrix(P)
  Z <- geomorph:::Effect.size.matrix(P)
  anova.tab <- data.frame(anova.tab, Z = c(Z, NA, NA), P.value = c(P.val, 
                                                                   NA, NA))
  anova.title = "\nRandomization of Raw Values used\n"
  attr(anova.tab, "heading") <- paste("\nType I (Sequential) Sums of Squares and Cross-products\n", 
                                      anova.title)
  class(anova.tab) <- c("anova", class(anova.tab))
  if (verbose == TRUE) {
    list(anova.table = anova.tab, call = match.call(), SS.rand = P)
  }
  else anova.tab
}