.onAttach <- function(libname, pkgname) {
  #version <- read.dcf(file=system.file("DESCRIPTION", package=pkgname),fields="Version")
  packageStartupMessage(pkgname, " is beta software. Please report any issues to the package author.")
}

phylocurve <- function(formula,tree,data,ymin=.01,ymax=.99,ylength=30,tip_coefficients,species_identifier="species",verbose=FALSE)
{
  if(missing(tip_coefficients))
  {
    if(!(species_identifier %in% colnames(data))) stop("Add species names column to data.")
    tip_coefficients <- get_tip_coefficients(formula = formula,tree = tree,data = data,ymin = ymin,ymax = ymax,ylength = ylength,species_identifier = species_identifier,verbose = verbose)
  } else if(verbose) cat("Phase 1: Tip estimates already provided\n")
  pgls_curve(tree = tree,tip_coefficients = tip_coefficients,ymin = ymin,ymax = ymax,ylength = ylength,verbose = verbose)
}

phylocurve.generalized <- function(tree,X,Y)
{
  nspecies <- length(tree$tip.label)
  p <- fast_anc_hand(x = X,Y = Y,tree = tree)
  tree1 <- multi2di(tree)
  ancx <- apply(t(matrix(p$aligned_coordinates$x,ncol = nspecies,dimnames = list(NULL,tree1$tip.label))),2,function(X) ace(X,tree1,method="pic")$ace[1])
  ancy <- apply(t(matrix(p$aligned_coordinates$y,ncol = nspecies,dimnames = list(NULL,tree1$tip.label))),2,function(X) ace(X,tree1,method="pic")$ace[1])
  p$nr <- nrow(p$aligned_coordinates)/nspecies
  
  ret <- c(p,list(anc_X=ancx,anc_Y=ancy,tree=tree))
  ret
}

polynomial.fit <- function(data, x_variable, y_variable, method = "BIC",
                           nterms = 2,min_x = -Inf,max_x = Inf,min_y = 0, max_y = Inf,eval_length=30)
{
  if(is.null(colnames(data))) stop("data variable must have column names")
  speciesInd <- match("species", colnames(data))
  xInd <- match(x_variable, colnames(data))
  yInd <- match(y_variable, colnames(data))
  if(is.na(speciesInd)) stop("data variable must have column named 'species'")
  if(is.na(xInd)) stop("data variable must have a column name matching x_variable")
  if(is.na(yInd)) stop("data variable must have a column name matching y_variable")
  dat <- data[,c(speciesInd,xInd,yInd)]
  colnames(dat) <- c("species","x","y")
  dat <- dat[dat$x>=min_x & dat$x<=max_x,]
  new_x <- unique(dat[,2])
  new_x <- seq(min(new_x),max(new_x),length=max(length(new_x),eval_length))
  species <- unique(as.character(dat[,1]))
  nspecies <- length(species)
  
  mod_list <- vector("list",nspecies)
  names(mod_list) <- unique(as.character(dat$species))
  
  ret <- matrix(NA,nrow = nspecies,ncol = length(new_x))
  colnames(ret) <- new_x
  rownames(ret) <- species
  mod <- "y ~ x"
  for(i in 1:nterms)
  {
    mod <- paste(mod," + I(x^",i,") ",sep="",collapse="")
  }
  mod <- formula(mod)
  
  for(i in 1:nspecies)
  {
    temp_species <- species[i]
    x <- dat[dat[,1]==temp_species,2]
    y <- dat[dat[,1]==temp_species,3]
    complete_cases_y <- which(complete.cases(y))
    x <- x[complete_cases_y]
    y <- y[complete_cases_y]
    n <- length(complete_cases_y)
    temp_lm <- lm(mod)
    k <- if(method=="AIC") 2 else log(n)
    mod_list[[i]] <- step(temp_lm,direction = "both",k = k,trace = 0)
    if(length(coef(mod_list[[i]]))==1) mod_list[[i]] <- lm(y~x)
    pred_y <- predict.lm(mod_list[[i]],newdata = data.frame(x=new_x))
    pred_y[pred_y<min_y] <- min_y
    pred_y[pred_y>max_y] <- max_y
    ret[i,] <- pred_y
  }
  ret <- ret[complete.cases(ret),]
  list(X=new_x,Y=ret)
}

nonlinear.fit <- function(data, x_variable, y_variable, fct = LL2.3(),
                          min_x = -Inf,max_x = Inf,min_y = 0, max_y = Inf,eval_length = 30,...)
{
  if(is.null(colnames(data))) stop("data variable must have column names")
  speciesInd <- match("species", colnames(data))
  xInd <- match(x_variable, colnames(data))
  yInd <- match(y_variable, colnames(data))
  if(is.na(speciesInd)) stop("data variable must have column named 'species'")
  if(is.na(xInd)) stop("data variable must have a column name matching x_variable")
  if(is.na(yInd)) stop("data variable must have a column name matching y_variable")
  dat <- data[,c(speciesInd,xInd,yInd)]
  colnames(dat) <- c("species","x","y")
  dat <- dat[dat$x>=min_x & dat$x<=max_x,]
  new_x <- unique(dat[,2])
  new_x <- seq(min(new_x),max(new_x),length=max(length(new_x),eval_length))
  species <- unique(as.character(dat[,1]))
  nspecies <- length(species)
  
  mod_list <- vector("list",nspecies)
  names(mod_list) <- unique(as.character(dat$species))
  
  ret <- matrix(NA,nrow = nspecies,ncol = length(new_x))
  colnames(ret) <- new_x
  rownames(ret) <- species
  
  for(i in 1:nspecies)
  {
    temp_species <- species[i]
    x <- dat[dat[,1]==temp_species,2]
    y <- dat[dat[,1]==temp_species,3]
    complete_cases_y <- which(complete.cases(y))
    x <- x[complete_cases_y]
    y <- y[complete_cases_y]
    temp_lm <- try(drm(y~x,fct = fct,...),silent=TRUE)
    if(class(temp_lm)=="try-error")
    {
      warning("Convergence failed for ",temp_species,". Using simple linear regression.",immediate. = TRUE)
      pred_y <- predict(lm(y~x),newdata = data.frame(x=new_x))
    } else
    {
      mod_list[[i]] <- temp_lm
      pred_y <- predict(mod_list[[i]],newdata = data.frame(x=new_x))
    }
    pred_y[pred_y<min_y] <- min_y
    pred_y[pred_y>max_y] <- max_y
    ret[i,] <- pred_y
  }
  
  ret <- ret[complete.cases(ret),]
  list(X=new_x,Y=ret)
}

GP.fit <- function(data, x_variable, y_variable,
                   min_x = -Inf,max_x = Inf,min_y = 0, max_y = Inf,eval_length = 30,...)
{
  if(is.null(colnames(data))) stop("data variable must have column names")
  speciesInd <- match("species", colnames(data))
  xInd <- match(x_variable, colnames(data))
  yInd <- match(y_variable, colnames(data))
  if(is.na(speciesInd)) stop("data variable must have column named 'species'")
  if(is.na(xInd)) stop("data variable must have a column name matching x_variable")
  if(is.na(yInd)) stop("data variable must have a column name matching y_variable")
  dat <- data[,c(speciesInd,xInd,yInd)]
  colnames(dat) <- c("species","x","y")
  dat <- dat[dat$x>=min_x & dat$x<=max_x,]
  new_x <- unique(dat[,2])
  new_x <- seq(min(new_x),max(new_x),length=max(length(new_x),eval_length))
  species <- unique(as.character(dat[,1]))
  nspecies <- length(species)
  
  mod_list <- vector("list",nspecies)
  names(mod_list) <- unique(as.character(dat$species))
  
  ret <- matrix(NA,nrow = nspecies,ncol = length(new_x))
  colnames(ret) <- new_x
  rownames(ret) <- species
  
  for(i in 1:nspecies)
  {
    temp_species <- species[i]
    x <- dat[dat[,1]==temp_species,2]
    y <- dat[dat[,1]==temp_species,3]
    complete_cases_y <- which(complete.cases(y))
    x <- x[complete_cases_y]
    y <- y[complete_cases_y]
    temp_lm <- try(GP_fit(X = normalize_to_01(c(range(new_x),x))[-(1:2)],Y = y,...),silent=TRUE)
    if(class(temp_lm)=="try-error")
    {
      warning("Convergence failed for ",temp_species,". Using simple linear regression.",immediate. = TRUE)
      pred_y <- predict(lm(y~x),newdata = data.frame(x=new_x))
    } else
    {
      mod_list[[i]] <- temp_lm
      pred_y <- predict.GP(mod_list[[i]],xnew = normalize_to_01(new_x))[[1]]
    }
    pred_y[pred_y<min_y] <- min_y
    pred_y[pred_y>max_y] <- max_y
    ret[i,] <- pred_y
  }
  ret <- ret[complete.cases(ret),]
  list(X=new_x,Y=ret)
}

phylocurve.trim <- function(phylocurve.generalized,min_Y = -Inf,max_Y = Inf,min_X = -Inf,max_X = Inf)
{
  p <- phylocurve.generalized
  tree <- p$tree
  ancx <- p$anc_X
  ancy <- p$anc_Y
  range_x <- which(ancx>min_X & ancx<max_X)
  range_y <- which(ancy>min_Y & ancy<max_Y)
  range <- intersect(range_x,range_y)
  nspecies <- length(tree$tip.label)
  nr <- nrow(p$aligned_coordinates)/nspecies
  p$aligned_coordinates <- p$aligned_coordinates[as.double(sapply((1:nspecies-1)*max(nr),function(X) X + range)),]
  p$aligned_data <- p$aligned_data[,c(1,range+1,range+nr+1)]
  p$aligned_X <- p$aligned_X[,range]
  p$aligned_Y <- p$aligned_Y[,range]
  p$anc_X <- p$anc_X[range]
  p$anc_Y <- p$anc_Y[range]
  p$nr <- nr
  p
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

get_aligned_function_data <- function(tip_coefficients,ylength=30,ymin=.01,ymax=.99)
{
  Y <- t(apply(tip_coefficients,1,logit_inv,y=seq(ymin,ymax,length=ylength)))
  data.frame(species=rownames(Y),Y)
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
  nvar <- ncol(P)
  var_BM <- apply(Q,2,function(X) t(three.point.compute(tree,X)$PP)/(nspecies-1))
  anc_three_point <- matrix(NA,ncol(Q),ncol(P))
  DD <- D1 <- matrix(NA,nvar,1)
  for(i in 1:ceiling(nvar/50))
  {
    if(verbose) cat(paste(round(i/ceiling(nvar/50)*100),"%   ",sep=""))
    temp <- three.point.compute(phy = tree,P = P[,(1+((i-1)*50)):min(50+(i-1)*50,nvar)],Q = Q)
    DD[(1+((i-1)*50)):min(50+(i-1)*50,nvar)] <- diag(temp$PP)
    D1[(1+((i-1)*50)):min(50+(i-1)*50,nvar)] <- temp$P1
    anc_three_point[,(1+((i-1)*50)):min(50+(i-1)*50,nvar)] <- temp$QP
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
  return(list(node_coefficients=ret,fitted_x=anc_vals,lower_CI_x=lower_CI,upper_CI_x=upper_CI,y_vals=seq(ymin,ymax,length=ylength),tip_coefficients=tip_coefficients,tip_X=Yinv))
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

sim.curves <- function(nspecies = 30,x_length=20,startree=FALSE,lambda=1,seed)
{
  if(!missing(seed)) set.seed(seed)
  x <- seq(0,10,length=x_length) # environmental gradient
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

threepoint_prepare <- function (phy, P, Q = NULL) 
{
  if (!inherits(phy, "phylo")) 
    stop("object \"phy\" is not of class \"phylo\".")
  phy <- reorder(phy, "pruningwise")
  n <- length(phy$tip.label)
  N <- dim(phy$edge)[1]
  ROOT <- n + 1L
  root.edge <- if (is.null(phy$root.edge)) 0 else phy$root.edge
  anc <- phy$edge[, 1]
  des <- phy$edge[, 2]
  flag <- 0
  if (is.null(Q)) {
    flag <- 1
    Q <- rep(1, n)
    names(Q) <- phy$tip.label
  }
  P <- as.matrix(P)
  Q <- as.matrix(Q)
  if (is.null(rownames(P))) stop("P needs to have row names.")
  ordr <- match(phy$tip.label, rownames(P))
  if (sum(is.na(ordr)) > 0) stop("row names of P do not match the tree tip labels.")
  P <- P[ordr, , drop = F]
  if (is.null(rownames(Q))) stop("Q needs to have row names.")
  ordr <- match(phy$tip.label, rownames(Q))
  if (sum(is.na(ordr)) > 0) stop("row names of Q do not match the tree tip labels.")
  Q <- Q[ordr, , drop = F]
  if (nrow(P) != n) 
    stop("the number of rows in P needs to be the same as the number of tips in the tree.")
  if (nrow(Q) != n) 
    stop("the number of rows in Q needs to be the same as the number of tips in the tree.")
  P <- cbind(rep(1, n), P)
  Q <- cbind(rep(1, n), Q)
  colP <- ncol(P)
  colQ <- ncol(Q)
  nout <- 2 + colP + colP^2 + colQ + colQ^2 + colP * colQ
  return(list(arg1=as.integer(N), arg2=as.integer(n), arg3=as.integer(phy$Nnode), 
              arg4=as.integer(colP), arg5=as.integer(colQ), arg6=as.integer(ROOT), 
              arg7=as.double(root.edge), arg8=as.double(phy$edge.length), arg9=as.integer(des), 
              arg10=as.integer(anc), arg11=as.double(as.vector(P)), arg12=as.double(as.vector(Q)), 
              result = double(nout),colP=colP,colQ=colQ))
}

threepoint_direct <- function(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10, arg11, arg12, result,colP,colQ)
{
  tmp <- .C("threepoint", arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10, arg11, arg12, result=result)$result
  PP <- matrix(tmp[2 + colP + 1:(colP^2)], colP, colP)
  QQ <- matrix(tmp[2 + colP + colP^2 + colQ + 1:(colQ^2)], colQ, 
               colQ)
  QP <- matrix(tmp[2 + colP + colP^2 + colQ + colQ^2 + 1:(colP * 
                                                            colQ)], colQ, colP)
  return(list(vec11 = PP[1, 1], P1 = PP[1, -1], PP = PP[-1, 
                                                        -1], Q1 = QQ[1, -1], QQ = QQ[-1, -1], QP = QP[-1, -1], logd = tmp[1]))
}

environment(threepoint_direct) <- asNamespace('phylolm')

sim.mult <- function(nspecies,R,error,nreps=1,nmissing=0,model,parameters,anc,tree,seed,nsims=1)
{
  if(nreps>1 & missing(error)) error <- rep(.5,ncol(R))
  if(length(R)==1) R <- matrix(R)
  if(!missing(seed)) set.seed(seed)
  if(missing(tree)) tree <- pbtree(n=nspecies)
  if(missing(anc)) anc <- rep(0,ncol(R)) else if(length(anc)==1) anc <- rep(anc,ncol(R))
  if(!missing(model))
  {
    atree <- tree
    args <- list(x=atree,model=model[1],lambda=1)
    for(i in 1:length(model))
    {
      args[[2]] <- model[i]
      names(args)[3] <- names(parameters)[i]
      args[[3]] <- parameters[[i]]
      atree <- do.call(rescale,args)
      #atree <- transf.branch.lengths(phy = atree,model = model[i],parameters = parameters[i])$tree
    }
  } else  atree <- tree
  ret <- vector("list",nsims)
  for(j in 1:nsims)
  {
    if(length(R)>1) Y <- sim.corrs(tree = atree,vcv = R,anc = anc,internal = FALSE) else Y <- as.matrix(cbind(fastBM(atree,a = anc,sig2 = R[1,1])))
    
    if(nreps>1)
    {
      Y_raw <- matrix(0,nrow = nspecies*nreps,ncol = ncol(R))
      for(i in 1:nreps)
      {
        Y_raw[1:nspecies+(nspecies*(i-1)),] <- matrix(rnorm(n = length(Y),mean = Y,sd = sqrt(as.double(t(matrix(rep(error,nspecies),nrow = length(error),ncol = nspecies))))),nrow = nspecies,ncol = ncol(R))
      }
    } else Y_raw <- Y
    if(nmissing>0)
    {
      Y_raw[sample(x = length(Y_raw),size = nmissing)] <- NA
    }
    ### reorder species by tree tips
    Y_raw <- data.frame(species=as.character(rep(tree$tip.label,nreps)),Y_raw)
    colnames(Y_raw)[1:ncol(R)+1] <- colnames(Y) <- paste("V",1:ncol(R),sep="")
    Y_means <- apply(Y_raw[,1:ncol(R)+1,drop=F],2,function(X) ave(X,Y_raw$species,FUN = function(Y) mean(Y,na.rm = TRUE)))
    Y_means <- Y_means[tree$tip.label,,drop=F]
    Y_array <- array(t(Y_means),dim = c(ncol(R),1,nspecies),dimnames = list(colnames(Y),NULL,rownames(Y_means))) 
    Y_means <- data.frame(species=as.character(rownames(Y_means)),Y_means)
    ret[[j]] <- list(Y_means=Y_means,Y_raw=Y_raw,Y_array=Y_array,Y_original=Y,tree=tree,tree_sim=atree)
  }
  if(nsims==1) ret[[1]] else ret
}

rate.mult <- function(tree=tree,Y=Y,type=c("mult","diag","all"),
                      method=c("REML","ML"),error=c("none","estimate","supply"),
                      error_n=20,error_supply,model="BM",fixed_sigma2,fixed_model_pars)
{
  rate.mult.args <- as.list(match.call())
  rate.mult.args <- rate.mult.args[2:length(rate.mult.args)]
  if(is.null(rate.mult.args$error)) rate.mult.args$error <- "none"
  if(is.null(rate.mult.args$model)) rate.mult.args$model <- "BM"
  if(missing(type)) rate.mult.args$type <- type[1]
  error <- error[1]
  if(missing(fixed_sigma2)) fixed_sigma2 <- NA
  if(missing(fixed_model_pars)) fixed_model_pars <- NA
  nspecies <- length(tree$tip.label)
  if(error=="estimate")
  {
    mean_data <- apply(Y[,2:ncol(Y),drop=FALSE],2,function(X) ave(X,Y[,1],FUN = function(Y) mean(Y,na.rm=TRUE)))
    mean_data <- data.frame(species=rownames(mean_data),mean_data)
    vars <- apply(Y[,2:ncol(Y),drop=FALSE],2,function(X) ave(X,Y[,1],FUN = function(Y) var(Y,na.rm=TRUE)))
    vars[is.na(vars)] <- 0
    ns <- apply(Y[,2:ncol(Y),drop=FALSE],2,function(X) ave(X,Y[,1],FUN = function(Y) length(which(!is.na(Y)))))
    num <- colSums(vars * (ns-1))
    temp_ns <- ns
    temp_ns[temp_ns<1] <- NA
    denom <- as.double(colSums(ns,na.rm = TRUE)) - as.double(apply(temp_ns,2,function(X) length(which(!is.na(X)))))
    pooled <- num / denom
    pooled[is.na(pooled)] <- 0
    pooled <- matrix(1,nspecies) %*% pooled
    #pooled <- matrix(1,nspecies) %*% colSums((ns-1)*vars/apply(ns,2,function(X) sum(X[X>1 & !is.na(X)])-length(X[X>1 & !is.na(X)])),na.rm=TRUE)
    vars[ns<error_n] <- pooled[ns<error_n]
    temp_ns[is.na(temp_ns)] <- 1
    errs <- vars[tree$tip.label,,drop=FALSE] / temp_ns
    Y <- mean_data[tree$tip.label,,drop=FALSE]
  } else if(error=="supply")
  {
    if(length(unlist(error_supply))==(ncol(Y)-1))
    {
      error_supply <- matrix(error_supply,1)
      vars <- matrix(1,nspecies) %*% error_supply
      rownames(vars) <- tree$tip.label
    } else vars <- error_supply
    errs <- vars[tree$tip.label,,drop=FALSE]
  }
  
  if(is.array(Y))
  {
    var_names <- names(Y[,,1])
    Y <- data.frame(species=names(Y[1,,]),t(matrix(Y,ncol = nspecies)),row.names=names(Y[1,,]))
    colnames(Y)[2:ncol(Y)] <- var_names
  }
  if(nrow(Y)>length(tree$tip.label))
  {
    mean_data <- apply(Y[,2:ncol(Y),drop=FALSE],2,function(X) ave(X,Y[,1],FUN = function(Y) mean(Y,na.rm=TRUE)))
    mean_data <- data.frame(species=rownames(mean_data),mean_data)
    Y <- mean_data[tree$tip.label,,drop=FALSE]
  }
  type <- type[1]
  method <- method[1]
  nvar <- ncol(Y)-1
  ones <- rbind(setNames(rep(1,nspecies),nm = tree$tip.label))
  anc <- sigma2 <- SS <- rbind(setNames(numeric(nvar),colnames(Y)[1:nvar+1]))
  P <- Y[,1:nvar+1,drop=F]
  
  ########## Below code adapted from phylolm package
  OU <- c("OUrandomRoot","OUfixedRoot")
  tol <- 1e-10  
  tree <- reorder(tree,"pruningwise")
  N <- dim(tree$edge)[1]
  ROOT <- nspecies + 1L
  anc <- tree$edge[, 1]
  des <- tree$edge[, 2]
  externalEdge <- des<=nspecies
  
  dis <- pruningwise.distFromRoot(tree)[1:nspecies]
  Tmax <- mean(dis)
  ## Default bounds
  bounds.default <- matrix(c(1e-7/Tmax,50/Tmax,1e-7/Tmax,50/Tmax,1e-7,1,1e-6,1,1e-5,3,-3/Tmax,0), ncol=2, byrow=TRUE)
  #bounds.default <- matrix(c(1/Tmax,50/Tmax,1/Tmax,50/Tmax,1e-7,1,1e-6,1,1e-5,3,-3/Tmax,0), ncol=2, byrow=TRUE)
  rownames(bounds.default) <- c("alpha","alpha","lambda","kappa","delta","rate")
  colnames(bounds.default) <- c("min","max")
  ## Default starting values
  starting.values.default <- c(0.5/Tmax,0.5/Tmax,0.5,0.5,0.5,-1/Tmax)
  starting.values.default[1:5] <- log(starting.values.default[1:5])
  bounds.default[1:5,] <- log(bounds.default[1:5,])
  names(starting.values.default) <- c("OUrandomRoot","OUfixedRoot","lambda","kappa","delta","EB")
  model_i <- match(model,names(starting.values.default))
  bounds.default <- bounds.default[model_i,,drop=FALSE]
  starting.values.default <- starting.values.default[model_i]
  names(starting.values.default) <- rownames(bounds.default)
  if(model[1]=="BM") starting.values.default <- bounds.default <- NULL
  
  flag <- 0 # flag and D are used for OU model if tree is not ultrametric:
  D <- NULL # for the generalized 3-point structure
  ## preparing for OU model
  if (any(model %in% OU)) {
    D <- numeric(nspecies)
    if (!is.ultrametric(tree)) flag <- 1
    dis <- pruningwise.distFromRoot(tree)
    Tmax <- max(dis[1:nspecies])
    D = Tmax - dis[1:nspecies]
    D = D - mean(D)
    tree$edge.length[externalEdge] <- tree$edge.length[externalEdge] + D[des[externalEdge]]
    ## phy is now ultrametric, with height Tmax:
    Tmax <- Tmax + min(D)
  }
  ########## Above code adapted from phylolm package
  
  temp_tree <- tree
  if(any(is.na(Y))) missing_data <- TRUE else missing_data <- FALSE
  if(!missing_data & error=="none" & model[1]=="BM") # basic mult model: no missing data, no within-species variation, BM model
  {
    if(!is.na(fixed_sigma2[1]))
    {
      ret_val <- est_BM(Y = Y,tree = tree,type = type,method = method,nspecies = nspecies,nvar = nvar,ones = ones,anc = anc,sigma2 = sigma2,SS = SS,fixed_sigma2 = fixed_sigma2)
    } else
    {
      ret_val <- est_BM(Y = Y,tree = tree,type = type,method = method,nspecies = nspecies,nvar = nvar,ones = ones,anc = anc,sigma2 = sigma2,SS = SS)
    }
    logL <- ret_val$logL
    sigma2 <- ret_val$sigma2
  } else
  {
    f <- function(theta,temp_Y,temp_tree,model,temp_error,temp_nspecies,temp_nvar,temp_ones,temp_anc,temp_sigma2,temp_SS)
    {
      theta[names(theta)!="rate"] <- exp(theta[names(theta)!="rate"])
      if(flag)
      {
        temp_Y <- exp(-theta["alpha"]*D) * temp_Y
      }
      for(j in 1:length(model))
      {
        if(model[1]!="BM")
        {
          if(is.na(fixed_model_pars[1]))
          {
            temp_tree <- transf.branch.lengths(phy = temp_tree,model = model[j],parameters = as.list(theta[is.na(fixed_sigma2[1])+j]),check.ultrametric = FALSE,D = D,check.names = FALSE)$tree
          } else
          {
            temp_tree <- transf.branch.lengths(phy = temp_tree,model = model[j],parameters = as.list(fixed_model_pars[j]),check.ultrametric = FALSE,D = D,check.names = FALSE)$tree
          }
        }
      }
      if(any(OU%in%model) & !is.na(fixed_sigma2[1])) temp_tree$edge.length <- temp_tree$edge.length / 2/if(is.na(fixed_model_pars[1])) theta[which(OU%in%model)] else fixed_model_pars[which(OU%in%model)]
      if(is.na(fixed_sigma2[1])) temp_tree$edge.length <- temp_tree$edge.length * theta[1] else temp_tree$edge.length <- temp_tree$edge.length * fixed_sigma2
      if(error!="none") temp_tree$edge.length[temp_tree$edge[,2]<=temp_nspecies] <- temp_tree$edge.length[temp_tree$edge[,2]<=temp_nspecies] + temp_error[temp_tree$tip.label[temp_tree$edge[temp_tree$edge[,2]<=temp_nspecies,2]],drop=F]
      logL <- try(est_BM(Y = temp_Y,tree = temp_tree,type = type,method = method,nspecies = temp_nspecies,
                         nvar = temp_nvar,ones = temp_ones,anc = temp_anc,sigma2 = temp_sigma2,SS = temp_SS,calc_sigma = FALSE)$logL,silent=TRUE)
      if(class(logL)=="try-error") return(-Inf)
      if(is.na(logL)) -Inf else logL
    }
    if(type=="mult" & !missing_data & error=="none" & model[1]!="BM") # basic mult model but not BM model
    {
      init <- c(sigma2=log(est_BM(Y = Y,tree = tree,type = type,method = method,nspecies = nspecies,
                                  nvar = nvar,ones = ones,anc = anc,sigma2 = sigma2,SS = SS)$sigma2),starting.values.default)
      if(!is.na(fixed_sigma2[1]))
      {
        if(is.na(fixed_model_pars[1]))
        {
          ret <- optim(init[2:length(init)],f,method = "L-BFGS-B",lower = bounds.default[,1],upper = bounds.default[,2],control=list(fnscale=-1),
                       temp_Y=Y,temp_tree=temp_tree,model=model,temp_error=error,temp_nspecies=nspecies,temp_nvar=nvar,temp_ones=ones,temp_anc=anc,temp_sigma2=sigma2,temp_SS=SS)
          ret <- list(par=c(log(fixed_sigma2),ret$par),val=ret$val)
        } else
        {
          logL <- f(fixed_model_pars,temp_Y=Y,temp_tree=temp_tree,model=model,temp_error=error,temp_nspecies=nspecies,temp_nvar=nvar,temp_ones=ones,temp_anc=anc,temp_sigma2=sigma2,temp_SS=SS)
          fixed_model_pars[names(fixed_model_pars)!="rate"] <- log(fixed_model_pars[names(fixed_model_pars)!="rate"])
          ret <- list(par=c(log(fixed_sigma2),fixed_model_pars),val=logL)
        }
      } else if(!is.na(fixed_model_pars[1]))
      {
        ret <- optim(init[1],f,method = "BFGS",control=list(fnscale=-1),
                     temp_Y=Y,temp_tree=temp_tree,model=model,temp_error=error,temp_nspecies=nspecies,temp_nvar=nvar,temp_ones=ones,temp_anc=anc,temp_sigma2=sigma2,temp_SS=SS)
        fixed_model_pars[names(fixed_model_pars)!="rate"] <- log(fixed_model_pars[names(fixed_model_pars)!="rate"])
        ret <- list(par=c(ret$par,fixed_model_pars),val=ret$val)
      } else
      {
        ret <- optim(init,f,method = "L-BFGS-B",lower = c(-Inf,bounds.default[,1]),upper = c(Inf,bounds.default[,2]),control=list(fnscale=-1),
                     temp_Y=Y,temp_tree=temp_tree,model=model,temp_error=error,temp_nspecies=nspecies,temp_nvar=nvar,temp_ones=ones,temp_anc=anc,temp_sigma2=sigma2,temp_SS=SS)
      }
      logL <- ret$val
      sigma2 <- exp(ret$par[1])
      pars <- if(length(ret$par)>1) ret$par[2:length(ret$par)] else NA
      pars[names(pars)!="rate"] <- exp(pars[names(pars)!="rate"])
    } else if(type=="mult")# & (missing_data | error!="none")) # mult model with missing data and/or measurement error (with or without BM evolutionary model)
    {
      init <- 0
      if(is.na(fixed_sigma2[1]))
      {
        for(i in 1:nvar)
        {
          temp_tree <- tree
          temp_nspecies <- nspecies
          temp_ones <- ones
          temp_Y <- Y[,c(1,i+1),drop=F]
          if(missing_data)
          {
            if(any(is.na(Y[,i+1])))
            {
              keep <- which(!is.na(Y[,i+1]))
              temp_tree <- drop.tip(tree,tip = as.character(Y[-keep,1]))
              temp_nspecies <- length(keep)
              temp_Y <- Y[keep,c(1,i+1),drop=F]
              temp_ones <- temp_ones[,keep,drop=F]
            }
          }
          init <- init + (temp_nspecies-1)*est_BM(Y = temp_Y,tree = temp_tree,type = type,method = method,nspecies = nspecies,
                                                  nvar = 1,ones = temp_ones,anc = anc[i],sigma2 = sigma2[i],SS = SS[i])$sigma2
        }
        init <- init / (sum(apply(Y[,1:nvar+1,drop=F],2,function(X) ave(X,Y$species,FUN = function(Y) length(which(!is.na(Y))))),na.rm = TRUE) - nvar)
      } else init <- fixed_sigma2
      #if(!is.na(fixed_model_pars[1])) fixed_model_pars[names(fixed_model_pars)!="rate"] <- log(fixed_model_pars[names(fixed_model_pars)!="rate"])
      init <- c(log(init),if(is.na(fixed_model_pars[1])) starting.values.default else fixed_model_pars)
      names(init)[1] <- "sigma2"
      f_mult <- function(theta)
      {
        theta[names(theta)!="rate"] <- exp(theta[names(theta)!="rate"])
        logL <- 0
        temp_tree <- tree
        for(j in 1:length(model))
        {
          if(model[1]!="BM")
          {
            if(is.na(fixed_model_pars[1]))
            {
              temp_tree <- transf.branch.lengths(phy = temp_tree,model = model[j],parameters = as.list(theta[is.na(fixed_sigma2[1])+j]),check.ultrametric = FALSE,D = D,check.names = FALSE)$tree
            } else
            {
              temp_tree <- transf.branch.lengths(phy = temp_tree,model = model[j],parameters = as.list(fixed_model_pars[j]),check.ultrametric = FALSE,D = D,check.names = FALSE)$tree
            }
          }
        }
        
        if(any(OU%in%model) & !is.na(fixed_sigma2[1])) temp_tree$edge.length <- temp_tree$edge.length / 2/if(is.na(fixed_model_pars[1])) theta[which(OU%in%model)] else fixed_model_pars[which(OU%in%model)]
        if(is.na(fixed_sigma2[1])) temp_tree$edge.length <- temp_tree$edge.length * theta[1] else temp_tree$edge.length <- temp_tree$edge.length * fixed_sigma2
        if(flag)
        {
          Y <- exp(-theta["alpha"]*D) * Y
        }
        temp_nspecies <- length(temp_tree$tip.label)
        perm_tree <- temp_tree
        for(i in 1:nvar)
        {
          temp_tree <- perm_tree
          if(error!="none") temp_tree$edge.length[temp_tree$edge[,2]<=temp_nspecies] <- temp_tree$edge.length[temp_tree$edge[,2]<=temp_nspecies] + errs[temp_tree$edge[temp_tree$edge[,2]<=temp_nspecies,2],i,drop=F]
          
          temp_nspecies <- nspecies
          temp_ones <- ones
          temp_Y <- Y[,c(1,i+1),drop=F]
          temp_D <- D
          if(missing_data)
          {
            if(any(is.na(Y[,i+1])))
            {
              keep <- which(!is.na(Y[,i+1]))
              temp_tree <- drop.tip(temp_tree,tip = as.character(Y[-keep,1]))
              temp_nspecies <- length(keep)
              temp_Y <- Y[keep,c(1,i+1),drop=F]
              temp_ones <- temp_ones[,keep,drop=F]
              temp_D <- temp_D[keep]
            }
          }
          logL <- try(logL + 
                        est_BM(Y = temp_Y,tree = temp_tree,type = type,method = method,nspecies = temp_nspecies,
                               nvar = 1,ones = temp_ones,anc = 0,sigma2 = 0,SS = 0,calc_sigma = FALSE)$logL,silent=TRUE)
          
          if(class(logL)=="try-error") return(-Inf)
        }
        if(is.na(logL)) -Inf else logL
      }
      if(!is.na(fixed_sigma2[1]) & model[1]!="BM")
      {
        if(is.na(fixed_model_pars[1]))
        {
          ret <- optim(init[2:length(init)],f_mult,method = "L-BFGS-B",lower = bounds.default[,1],upper = bounds.default[,2],control=list(fnscale=-1))
          ret <- list(par=c(log(fixed_sigma2),ret$par),val=ret$val)
        } else
        {
          logL <- f_mult(fixed_model_pars)
          fixed_model_pars[names(fixed_model_pars)!="rate"] <- log(fixed_model_pars[names(fixed_model_pars)!="rate"])
          ret <- list(par=c(log(fixed_sigma2),fixed_model_pars),val=logL)
        }
      } else if(!is.na(fixed_model_pars[1]))
      {
        ret <- optim(init[1],f_mult,method = "BFGS",control=list(fnscale=-1))
        fixed_model_pars[names(fixed_model_pars)!="rate"] <- log(fixed_model_pars[names(fixed_model_pars)!="rate"])
        ret <- list(par=c(ret$par,fixed_model_pars),val=ret$val)
      } else
      {
        ret <- optim(init,f_mult,method = "L-BFGS-B",lower = c(-Inf,bounds.default[,1]),upper = c(Inf,bounds.default[,2]),control=list(fnscale=-1))
      }
      logL <- ret$val
      sigma2 <- exp(ret$par[1])
      if(model[1]=="BM")
      {
        pars <- NA 
      } else
      {
        if(is.na(fixed_model_pars))
        {
          pars <- ret$par[(is.na(fixed_sigma2[1]) + 1):length(ret$par)]
        } else pars <- fixed_model_pars
      }
      pars[names(pars)!="rate"] <- exp(pars[names(pars)!="rate"])
    } else if(type=="diag") # diag model with missing data and/or measurement error (with or without BM evolutionary model)
    {
      init <- rep(0,nvar)
      if(is.na(fixed_sigma2[1]))
      {
        for(i in 1:nvar)
        {
          temp_tree <- tree
          temp_nspecies <- nspecies
          temp_ones <- ones
          temp_Y <- Y[,c(1,i+1),drop=F]
          if(missing_data)
          {
            if(any(is.na(Y[,i+1])))
            {
              keep <- which(!is.na(Y[,i+1]))
              temp_tree <- drop.tip(tree,tip = as.character(Y[-keep,1]))
              temp_nspecies <- length(keep)
              temp_Y <- Y[keep,c(1,i+1),drop=F]
              temp_ones <- temp_ones[,keep,drop=F]
            }
          }
          init[i] <- est_BM(Y = temp_Y,tree = temp_tree,type = type,method = method,nspecies = nspecies,
                            nvar = 1,ones = temp_ones,anc = anc[i],sigma2 = sigma2[i],SS = SS[i])$sigma2
        }
      } else init <- fixed_sigma2
      init <- c(log(init),if(is.na(fixed_model_pars[1])) starting.values.default else fixed_model_pars)
      names(init)[1:nvar] <- "sigma2"
      f_diag <- function(theta,i,temp_tree)
      {
        theta[names(theta)!="rate"] <- exp(theta[names(theta)!="rate"])
        logL <- 0
        if(is.na(fixed_sigma2[1])) temp_tree$edge.length <- temp_tree$edge.length * theta[1] else temp_tree$edge.length <- temp_tree$edge.length * fixed_sigma2[i]
        
        if(flag)
        {
          Y <- exp(-theta["alpha"]*D) * Y
        }  
        temp_tree <- temp_tree
        temp_nspecies <- length(temp_tree$tip.label)
        
        if(error!="none") temp_tree$edge.length[temp_tree$edge[,2]<=temp_nspecies] <- temp_tree$edge.length[temp_tree$edge[,2]<=temp_nspecies] + errs[temp_tree$edge[temp_tree$edge[,2]<=temp_nspecies,2],i,drop=F]
        
        temp_nspecies <- nspecies
        temp_ones <- ones
        temp_Y <- Y[,c(1,i+1),drop=F]
        temp_D <- D
        if(missing_data)
        {
          if(any(is.na(Y[,i+1])))
          {
            keep <- which(!is.na(Y[,i+1]))
            temp_tree <- drop.tip(temp_tree,tip = as.character(Y[-keep,1]))
            temp_nspecies <- length(keep)
            temp_Y <- Y[keep,c(1,i+1),drop=F]
            temp_ones <- temp_ones[,keep,drop=F]
            temp_D <- temp_D[keep]
          }
        }
        logL <- try(est_BM(Y = temp_Y,tree = temp_tree,type = type,method = method,nspecies = temp_nspecies,
                           nvar = 1,ones = temp_ones,anc = 0,sigma2 = 0,SS = 0,calc_sigma = FALSE)$logL,silent=TRUE)
        if(class(logL)=="try-error") return(-Inf)
        if(is.na(logL)) -Inf else logL
      }
      
      if(model[1]=="BM")
      {
        #sigma2 <- log(init)
        sigma2 <- init
        logL <- 0
        for(i in 1:nvar)
        {
          if(is.na(fixed_sigma2[1]))
          {
            temp_ret <- optim(sigma2[i],f_diag,method = "BFGS",control=list(fnscale=-1),i=i,temp_tree=tree)
            logL <- logL + temp_ret$val
            sigma2[i] <- log(temp_ret$par)
          } else
          {
            logL <- logL + f_diag(0,i=i,temp_tree=tree)
            sigma2[i] <- fixed_sigma2[i]
          }
        }
        pars <- NA
      } else
      {
        f_diag_mod <- function(theta,ret_sigma2 = FALSE)
        {
          #sigma2 <- log(init[1:nvar])
          sigma2 <- init[1:nvar]
          theta[names(theta)!="rate"] <- exp(theta[names(theta)!="rate"])
          temp_tree <- tree
          #for(j in 1:length(model))
          #{
          #  if(model[1]!="BM")
          #  {
          #    temp_tree <- transf.branch.lengths(phy = temp_tree,model = model[j],parameters = as.list(theta[j]),check.ultrametric = FALSE,D = D,check.names = FALSE)$tree          
          #  }
          #}
          
          for(j in 1:length(model))
          {
            if(model[1]!="BM")
            {
              if(is.na(fixed_model_pars[1]))
              {
                temp_tree <- transf.branch.lengths(phy = temp_tree,model = model[j],parameters = as.list(theta[j]),check.ultrametric = FALSE,D = D,check.names = FALSE)$tree
              } else
              {
                temp_tree <- transf.branch.lengths(phy = temp_tree,model = model[j],parameters = as.list(fixed_model_pars[j]),check.ultrametric = FALSE,D = D,check.names = FALSE)$tree
              }
            }
          }
          
          if(any(OU%in%model) & !is.na(fixed_sigma2[1])) temp_tree$edge.length <- temp_tree$edge.length / 2/if(is.na(fixed_model_pars[1])) theta[which(OU%in%model)] else fixed_model_pars[which(OU%in%model)]
          
          model <- "BM"          
          logL <- 0
          for(i in 1:nvar)
          {  
            #temp_ret <- optim(c(sigma2[i]),f_diag,method = "BFGS",control=list(fnscale=-1),i=i,temp_tree=temp_tree)            
            #logL <- logL + temp_ret$val
            #sigma2[i] <- temp_ret$par
            
            if(is.na(fixed_sigma2[1]))
            {
              temp_ret <- optim(sigma2[i],f_diag,method = "BFGS",control=list(fnscale=-1),i=i,temp_tree=temp_tree)
              logL <- logL + temp_ret$val
              sigma2[i] <- temp_ret$par
            } else
            {
              logL <- logL + f_diag(fixed_sigma2[i],i=i,temp_tree=temp_tree)
              sigma2[i] <- fixed_sigma2[i]
            }
            
          }
          if(ret_sigma2 & is.na(fixed_sigma2[1])) return(sigma2) else if(ret_sigma2) return(log(fixed_sigma2))
          if(is.na(logL)) -Inf else logL
        }
        if(is.na(fixed_model_pars[1]))
        {
          ret <- optim(starting.values.default,fn = f_diag_mod,method = "L-BFGS-B",lower = bounds.default[,1],upper = bounds.default[,2],control=list(fnscale=-1))          
          pars <- ret$par
          logL <- ret$val  
        } else
        {
          pars <- fixed_model_pars
          logL <- f_diag_mod(fixed_model_pars)
          #fixed_model_pars[names(fixed_model_pars)!="rate"] <- log(fixed_model_pars[names(fixed_model_pars)!="rate"])
          #pars <- fixed_model_pars
        }
        if(is.na(fixed_sigma2[1])) sigma2 <- exp(f_diag_mod(pars,TRUE)) else sigma2 <- fixed_sigma2
        if(is.na(fixed_model_pars[1])) pars[names(pars)!="rate"] <- exp(pars[names(pars)!="rate"])
      }
    }
  }
  if(model[1]!="BM")
  {
    bounds.default[rownames(bounds.default)!="rate",] <- exp(bounds.default[rownames(bounds.default)!="rate",])
    if(any(abs(pars-bounds.default[,1,drop=F])<1e-10))
      warning(paste("the estimation of", names(pars)[abs(pars-bounds.default[,1])<1e-10], 'matches the lower bound for this parameter.'))
    if(any(abs(pars-bounds.default[,2,drop=F])<1e-10))
      warning(paste("the estimation of", names(pars)[abs(pars-bounds.default[,2])<1e-10], 'matches the upper bound for this parameter.'))
    if(flag) logL <- logL + pars["alpha"] * 2*sum(D)
    if(any(OU%in%model) & is.na(fixed_sigma2[1])) sigma2 <- 2*as.double(pars["alpha"])*sigma2
    ret_val <- list(sigma2=sigma2,pars=pars,logL=logL,method=method)
    logL <- ret_val$logL
    sigma2 <- ret_val$sigma2
  } else
  {
    ret_val <- list(sigma2=sigma2,logL=logL,method=method)
  }
  ret_val$model <- model
  ret_val$D <- D
  for(i in 1:length(rate.mult.args)) rate.mult.args[[i]] <- get(names(rate.mult.args)[i])
  if(error!="none")
  {
    rate.mult.args$error <- "supply"
    rate.mult.args$error_supply <- errs
  }
  ret_val$rate.mult.args <- rate.mult.args
  class(ret_val) <- "rate.mult"
  ret_val
}

est_BM <- function(Y,tree,type,method,nspecies,nvar,ones,anc,sigma2,SS,calc_sigma=TRUE,fixed_sigma2=NA)
{
  if(!is.na(fixed_sigma2[1])) calc_sigma <- TRUE
  temp_tree <- tree
  P <- Y[,1:nvar+1,drop=F]
  for(i in 1:ceiling(nvar/50))
  {
    rng <- (1+((i-1)*50)):min(50+(i-1)*50,nvar)
    temp <- three.point.compute(phy = temp_tree,P = P[,rng,drop=F])
    anc[rng] <- solve(temp$vec11)%*%temp$P1
  }
  P <- Y[,1:nvar+1,drop=F] - (t(ones)%*%anc)
  if(calc_sigma)
  {
    if(is.na(fixed_sigma2[1]))
    {
      for(i in 1:ceiling(nvar/50))
      {
        rng <- (1+((i-1)*50)):min(50+(i-1)*50,nvar)
        temp <- three.point.compute(phy = temp_tree,P = P[,rng,drop=F])
        sigma2[rng] <- if(length(temp$PP)==1) temp$PP else diag(temp$PP)
      }
      if(type=="mult")
      {
        sigma2 <- sum(sigma2)
      }
      if(method=="ML")
      {
        sigma2 <- sigma2 / nspecies
      } else if(method=="REML")
      {
        sigma2 <- sigma2 / (nspecies-1)
      }
    } else sigma2 <- fixed_sigma2
  }
  if(type=="mult")
  {
    if(calc_sigma) 
    { 
      if(is.na(fixed_sigma2[1])) sigma2 <- sigma2 / nvar
      temp_tree$edge.length <- tree$edge.length*sigma2
    }
    for(i in 1:ceiling(nvar/50))
    {
      rng <- (1+((i-1)*50)):min(50+(i-1)*50,nvar)
      temp <- three.point.compute(phy = temp_tree,P = P[,rng,drop=F])
      SS[rng] <- if(length(temp$PP)==1) temp$PP else diag(temp$PP)
    }
    logL <- -.5*(sum(SS) + nvar*(nspecies - if(method=="REML") 1 else 0)*log(2*pi) + nvar*(temp$logd + if(method=="REML") log(temp$vec11) else 0))
  } else if(type=="diag")
  {
    logL <- numeric(nvar)
    for(i in 1:nvar)
    {
      if(calc_sigma) temp_tree$edge.length <- tree$edge.length*sigma2[i]
      temp <- three.point.compute(phy = temp_tree,P = P[,i,drop=F])
      SS <- temp$PP
      nvar <- 1
      logL[i] <- -.5*(sum(SS) + nvar*(nspecies - if(method=="REML") 1 else 0)*log(2*pi) + nvar*(temp$logd + if(method=="REML") log(temp$vec11) else 0))
    }
    logL <- sum(logL)
  } else if(type=="all") return(est_all(tree,Y,method))
  list(sigma2=sigma2,logL=logL,method=method)
}

est_all <- function(tree,Y,method,estimate_logl=TRUE)
{
  nvar <- ncol(Y) - 1
  tr <- multi2di(tree)
  anc <- (matrix(1,nrow(Y),1)%*%apply(Y[,1:nvar+1,drop=F],2,function(X) ace(x = setNames(X,Y$species),phy=tr,method="pic")$ace[1]))[1,]
  
  #temp <- three.point.compute(phy = tree,P = data.frame(Y[,1:nvar+1,drop=F],row.names = Y[,1]))
  #anc <- solve(temp$vec11)%*%temp$P1
  nspecies <- length(tree$tip.label)
  ones <- cbind(setNames(rep(1,nspecies),tree$tip.label))
  P <- Y[,1:nvar+1,drop=F] - (ones%*%anc)
  SS <- three.point.compute(phy = tree,P = P)$PP
  sigma2 <- SS / (nspecies - (as.integer(method=="REML")))
  if(estimate_logl)
  {
    SS <- t(unlist(P)) %*% solve(sigma2 %x% vcv(tree)) %*% unlist(P) # INEFFICIENT -- implement tree transversal function
    logL <- -.5*(SS + (as.double(determinant(sigma2,logarithm = TRUE)$modulus)*nspecies + nvar*three.point.compute(tree,ones)$logd) + (nspecies - (as.integer(method=="REML")))*log(2*pi))
    list(sigma2=sigma2,logL=logL,method=method)
  } else
    list(sigma2=sigma2,logL=NA,method=method)
}

compare.rate.mult <- function(rate.mult.fitted,groups,fit_individual=FALSE)
{
  full_tree <- rate.mult.fitted$rate.mult.args$tree
  ngroups <- nlevels(groups)
  group_levels <- levels(groups)
  rmf <- rate.mult.fitted$rate.mult.args
  rate.mult.args <- rep(list(rmf),ngroups)
  
  alt_model <- null_model <- vector("list",length = ngroups)
  names(rate.mult.args) <- names(alt_model) <- names(null_model) <- group_levels
  alt.logL <- null.logL <- 0
  
  if(fit_individual | (rate.mult.fitted$model[1]=="BM"))
  {
    for(i in 1:ngroups)
    {
      keep_taxa <- names(groups)[which(groups==group_levels[i])]
      drop_taxa <- name.check(phy = full_tree,data.names = keep_taxa)$tree_not_data
      rate.mult.args[[i]]$tree <- drop.tip(full_tree,tip = drop_taxa)
      rate.mult.args[[i]]$Y <- rate.mult.args[[i]]$Y[keep_taxa,]
      if(rate.mult.args[[i]]$error=="supply") rate.mult.args[[i]]$error_supply <- rate.mult.args[[i]]$error_supply[keep_taxa,]
      alt_model[[i]] <- do.call(rate.mult,rate.mult.args[[i]])
      rate.mult.args[[i]]$fixed_sigma2 <- rate.mult.fitted$sigma2
      if(rate.mult.fitted$model[1]!="BM") rate.mult.args[[i]]$fixed_model_pars <- rate.mult.fitted$pars
      null_model[[i]] <- do.call(rate.mult,rate.mult.args[[i]])
      null.logL <- null.logL + null_model[[i]]$logL
      alt.logL <- alt.logL + alt_model[[i]]$logL
    }
  } else
  {
    for(i in 1:ngroups)
    {
      keep_taxa <- names(groups)[which(groups==group_levels[i])]
      drop_taxa <- name.check(phy = full_tree,data.names = keep_taxa)$tree_not_data
      rate.mult.args[[i]]$tree <- drop.tip(full_tree,tip = drop_taxa)
      rate.mult.args[[i]]$Y <- rate.mult.args[[i]]$Y[keep_taxa,]
      rate.mult.args[[i]]$fixed_model_pars <- rate.mult.fitted$pars
      if(rate.mult.args[[i]]$error=="supply") rate.mult.args[[i]]$error_supply <- rate.mult.args[[i]]$error_supply[keep_taxa,]
      alt_model[[i]] <- do.call(rate.mult,rate.mult.args[[i]])
      rate.mult.args[[i]]$fixed_sigma2 <- rate.mult.fitted$sigma2
      null_model[[i]] <- do.call(rate.mult,rate.mult.args[[i]])
      null.logL <- null.logL + null_model[[i]]$logL
      alt.logL <- alt.logL + alt_model[[i]]$logL
    }
  }
  null.pars <- length(rate.mult.fitted$sigma2)
  alt.pars <- null.pars*ngroups
  if(rate.mult.fitted$model[1]!="BM")
  {
    null.pars <- null.pars + length(rate.mult.fitted$pars)
    if(fit_individual)
    {
      alt.pars <- alt.pars + (length(rate.mult.fitted$pars)*ngroups)
    } else alt.pars <- alt.pars + length(rate.mult.fitted$pars)
  }
  df <- as.integer(alt.pars-null.pars)
  chi_sq=-2*null.logL+2*alt.logL
  p <- pchisq(chi_sq,df = as.integer(alt.pars - null.pars),lower.tail = FALSE)
  ret <- list(null.logL=null.logL,null.pars=as.integer(null.pars),alt.logL=alt.logL,alt.pars=as.integer(alt.pars),df=df,chi_sq=chi_sq,p=p,method=rate.mult.fitted$method,null_model_list=null_model,alt_model_list=alt_model)
  class(ret) <- "lr.test"
  ret
}

compare.multivar.rate.mult <- function(null_model,alt_model_list)
{
  ntraits <- length(alt_model_list)
  alt.logL <- 0
  alt.pars <- 0
  for(i in 1:ntraits)
  {
    alt.logL <- alt.logL + alt_model_list[[i]]$logL
    alt.pars <- alt.pars + length(alt_model_list[[i]]$sigma2)
    if(alt_model_list[[i]]$rate.mult.args$model[1]!="BM") alt.pars <- alt.pars + length(alt_model_list[[i]]$pars)
  }
  null.logL <- null_model$logL
  null.pars <- length(null_model$sigma2)
  if(null_model$rate.mult.args$model[1]!="BM") null.pars <- null.pars + length(null_model$pars)
  df <- as.integer(alt.pars-null.pars)
  chi_sq=-2*null.logL+2*alt.logL
  p <- pchisq(chi_sq,df = as.integer(alt.pars - null.pars),lower.tail = FALSE)
  ret <- list(null.logL=null.logL,null.pars=as.integer(null.pars),alt.logL=alt.logL,alt.pars=as.integer(alt.pars),df=df,chi_sq=chi_sq,p=p,method=null_model$method)
  class(ret) <- "lr.test"
  ret
}

K.mult <- function(rate.mult.fitted,iter=1000)
{
  tree <- reorder(rate.mult.fitted$rate.mult.args$tree,"pruningwise")
  nspecies <- length(tree$tip.label)
  nvar <- ncol(rate.mult.fitted$rate.mult.args$Y) - 1
  if(any(is.na(rate.mult.fitted$rate.mult.args$Y))) missing_data <- TRUE else missing_data <- FALSE
  if(rate.mult.fitted$rate.mult.args$error!="none")
  {
    error <- TRUE
    errs <- rate.mult.fitted$rate.mult.args$error_supply
  } else error <- FALSE
  D <- numeric(nspecies)
  if(rate.mult.fitted$rate.mult.args$model[1]=="BM") BM <- TRUE else
  {
    BM <- FALSE
    model <- rate.mult.fitted$rate.mult.args$model
    pars <- rate.mult.fitted$pars
    if("alpha" %in% names(pars) & !is.ultrametric(tree))
    {
      D <- rate.mult.fitted$rate.mult.args$D
      rate.mult.fitted$rate.mult.args$Y[,1:nvar+1] <- exp(-pars["alpha"]*rate.mult.fitted$rate.mult.args$D) * rate.mult.fitted$rate.mult.args$Y[,1:nvar+1]
    }
    for(j in 1:length(model))
    {
      temp_tree <- transf.branch.lengths(phy = tree,model = model[j],parameters = as.list(pars[j]),check.ultrametric = FALSE,D = D,check.names = FALSE)$tree
    }
  }
  Y <- rate.mult.fitted$rate.mult.args$Y
  Y <- Y[tree$tip.label,]
  type <- rate.mult.fitted$rate.mult.args$type
  sigma2 <- rate.mult.fitted$sigma2
  
  K.mult_f <- function(Y)
  {
    K.num.sum <- 0
    K.denom.sum <- 0
    for(i in 1:nvar)
    {
      temp_sigma2 <- if(type=="mult") sigma2 else sigma2[i]
      temp_tree <- tree
      temp_dist_from_root <- dist_from_root
      if(!BM)
      { 
        if(any("alpha"==names(pars)))
        {
          temp_tree$edge.length <- temp_tree$edge.length / 2 / pars[which(names(pars)=="alpha")]
          temp_dist_from_root <- temp_dist_from_root / 2 / pars[which(names(pars)=="alpha")]
          tp_prep[[7]] <- temp_tree$root.edge * temp_sigma2
        }
      }
      temp_tree$edge.length <- temp_tree$edge.length * temp_sigma2
      temp_dist_from_root <- temp_dist_from_root * temp_sigma2
      temp_nspecies <- nspecies
      if(error)
      {
        temp_tree$edge.length[temp_tree$edge[,2]<=temp_nspecies] <- temp_tree$edge.length[temp_tree$edge[,2]<=temp_nspecies] + errs[temp_tree$edge[temp_tree$edge[,2]<=temp_nspecies,2],i,drop=F]
        temp_dist_from_root <- temp_dist_from_root + errs[temp_tree$edge[temp_tree$edge[,2]<=temp_nspecies,2],i,drop=F]
      }
      temp_Y <- Y[,i+1,drop=F]
      temp_ones <- ones
      if(missing_data)
      {
        keep <- which(!is.na(temp_Y))
        temp_tree <- reorder(drop.tip(tree,tip = as.character(Y[-keep,1])),"pruningwise")
        temp_nspecies <- length(keep)
        temp_Y <- temp_Y[keep,,drop=F]
        temp_ones <- temp_ones[keep]
        temp_dist_from_root <- temp_dist_from_root[keep]
      }
      if(nspecies!=temp_nspecies)
      {
        tp_prep <- threepoint_prepare(temp_tree,temp_Y)
        x_span <- (temp_nspecies+1):length(tp_prep[[11]])
      }
      temp_dist_from_root <- sum(temp_dist_from_root)
      temp_Y <- temp_Y[,1]
      tp_prep[[11]][x_span] <- temp_Y
      tp_prep[[8]] <- as.double(temp_tree$edge.length)
      tp <- do.call(threepoint_direct,tp_prep)
      K.denom <- (temp_dist_from_root - (temp_nspecies / tp$vec11))/(temp_nspecies - 1)
      a.obs <- t(tp$P1/tp$vec11)
      x_a <- temp_Y - temp_ones %*% a.obs
      
      tp_prep[[11]][x_span] <- as.double(x_a)
      tp <- do.call(threepoint_direct,tp_prep)
      K.denom <- tp$PP * K.denom
      K.num <- sum(x_a^2)
      K.num.sum <- K.num.sum + K.num
      K.denom.sum <- K.denom.sum + K.denom
    }
    return(K.num.sum/K.denom.sum)
  }
  
  dist_from_root <- pruningwise.distFromRoot(reorder(tree,"pruningwise"))[1:nspecies]
  tp_prep <- threepoint_prepare(tree,Y[,2,drop=F])
  ones <- rep(1,nspecies)
  
  tp_prep <- threepoint_prepare(tree,Y[,2,drop=F])
  x_span <- (nspecies+1):length(tp_prep[[11]])
  K.obs <- K.mult_f(Y)
  P.val <- 1
  K.val <- rep(0, iter)
  for (i in 1:iter) {
    rand_Y <- Y[sample(nrow(Y)), ]
    rownames(rand_Y) <- rand_Y[,1] <- rownames(Y)
    K.rand <- K.mult_f(rand_Y)
    P.val <- ifelse(K.rand >= K.obs, P.val + 1, P.val)
    K.val[i] <- K.rand
  }
  P.val <- P.val/(iter + 1)
  K.val[iter + 1] = K.obs
  return(list(phy.signal = K.obs, pvalue = P.val))
}

pgls.mult <- function(rate.mult.fitted,X)
{
  tree <- reorder(rate.mult.fitted$rate.mult.args$tree,"pruningwise")
  nspecies <- length(tree$tip.label)
  type <- rate.mult.fitted$rate.mult.args$type
  sigma2 <- rate.mult.fitted$sigma2
  
  if(class(X)=="data.frame") X <- as.matrix(X)
  if(length(X)==nspecies | !is.matrix(X))
  {
    X <- X[tree$tip.label,]
    X <- matrix(cbind(1,as.matrix(X)),ncol = 2,dimnames=list(tree$tip.label,c("Intercept","X")))
  } else
  {
    X <- X[tree$tip.label,]
    X <- cbind(Intercept=1,X)
  }
  
  nvar <- ncol(rate.mult.fitted$rate.mult.args$Y) - 1
  if(any(is.na(rate.mult.fitted$rate.mult.args$Y))) missing_data <- TRUE else missing_data <- FALSE
  if(rate.mult.fitted$rate.mult.args$error!="none")
  {
    error <- TRUE
    errs <- rate.mult.fitted$rate.mult.args$error_supply
  } else error <- FALSE
  D <- numeric(nspecies)
  if(rate.mult.fitted$rate.mult.args$model[1]=="BM") BM <- TRUE else
  {
    BM <- FALSE
    model <- rate.mult.fitted$rate.mult.args$model
    pars <- rate.mult.fitted$pars
    if("alpha" %in% names(pars) & !is.ultrametric(tree))
    {
      D <- rate.mult.fitted$rate.mult.args$D
      rate.mult.fitted$rate.mult.args$Y[,1:nvar+1] <- exp(-pars["alpha"]*rate.mult.fitted$rate.mult.args$D) * rate.mult.fitted$rate.mult.args$Y[,1:nvar+1]
    }
    for(j in 1:length(model))
    {
      temp_tree <- transf.branch.lengths(phy = tree,model = model[j],parameters = as.list(pars[j]),check.ultrametric = FALSE,D = D,check.names = FALSE)$tree
    }
  }
  
  Y <- rate.mult.fitted$rate.mult.args$Y
  tp_prep <- threepoint_prepare(tree,X,Y[,2,drop=F])
  x_span <- (nspecies+1):length(tp_prep[[11]])
  y_span <- (nspecies+1):length(tp_prep[[12]])
  ones <- rep(1,nspecies)
  k <- ncol(X)
  rdf <- nspecies - k
  SS_total <- SS_reg <- SS_res <- SS_res1 <- 0
  SSEs <- numeric(k)
  for(i in 1:nvar)
  {
    temp_nspecies <- nspecies
    temp_sigma2 <- if(type=="mult") sigma2 else sigma2[i]
    temp_tree <- tree
    if(!BM)
    { 
      if(any("alpha"==names(pars)))
      {
        temp_tree$edge.length <- temp_tree$edge.length / 2 / pars[which(names(pars)=="alpha")]
        tp_prep[[7]] <- temp_tree$root.edge * temp_sigma2
      }
    }
    temp_tree$edge.length <- temp_tree$edge.length * temp_sigma2 / temp_sigma2
    temp_nspecies <- nspecies
    if(error)
    {
      temp_tree$edge.length[temp_tree$edge[,2]<=temp_nspecies] <- temp_tree$edge.length[temp_tree$edge[,2]<=temp_nspecies] + errs[temp_tree$edge[temp_tree$edge[,2]<=temp_nspecies,2],i,drop=F]
    }
    temp_Y <- Y[,i+1,drop=F]
    temp_X <- X
    temp_ones <- ones
    if(missing_data)
    {
      keep <- which(!is.na(temp_Y))
      temp_tree <- reorder(drop.tip(tree,tip = as.character(Y[-keep,1])),"pruningwise")
      temp_nspecies <- length(keep)
      temp_Y <- temp_Y[keep,,drop=F]
      temp_X <- X[keep,,drop=F]
      temp_ones <- temp_ones[keep]
    }
    if(nspecies!=temp_nspecies)
    {
      tp_prep <- threepoint_prepare(temp_tree,temp_X,temp_Y)
      x_span <- (temp_nspecies+1):length(tp_prep[[11]])
      y_span <- (temp_nspecies+1):length(tp_prep[[12]])
    }
    temp_Y <- temp_Y[,1]
    tp_prep[[11]][x_span] <- temp_X
    tp_prep[[12]][y_span] <- temp_Y
    tp_prep[[8]] <- as.double(temp_tree$edge.length)
    tp <- do.call(threepoint_direct,tp_prep)
    
    a.obs <- t(tp$Q1/tp$vec11)
    temp_beta <- solve(tp$PP)%*%tp$QP
    
    y_ones_diff <- temp_Y - temp_ones %*% a.obs
    tp_prep[[12]][y_span] <- as.double(y_ones_diff)
    tp <- do.call(threepoint_direct,tp_prep)
    SS_total <- SS_total + tp$QQ
    
    y_hat_diff <- temp_Y - temp_X%*% temp_beta
    tp_prep[[12]][y_span] <- as.double(y_hat_diff)
    tp <- do.call(threepoint_direct,tp_prep)
    SS_res <- SS_res + tp$QQ
    
    for(j in 2:k)
    {
      y_hat_diff <- temp_Y - temp_X[,c(1:j)]%*%temp_beta[c(1:j)]
      tp_prep[[12]][y_span] <- as.double(y_hat_diff)
      tp <- do.call(threepoint_direct,tp_prep)
      SSEs[j] <- SSEs[j] + tp$QQ
    }
    
    SS_reg <- SS_total - SS_res
    var <- SS_reg / (temp_nspecies - 1)
  }
  SSEs[1] <- SS_total
  SS.tmp <- c(SSEs[-1],SSEs[k])
  df <- 1:k
  df.tmp <- c(df[-1],df[k])
  df <- (df.tmp - df)[1:(k-1)]
  SS <- (SSEs - SS.tmp)[1:(k-1)]
  
  MS <- SS/df
  R2 <- SS/SS_total
  SSE.model <- SS_total - sum(SS)
  dfE <- nrow(Y)-(sum(df)+1)
  MSE <- SSE.model/dfE
  Fs <- MS/MSE
  df <- c(df,dfE,nrow(Y)-1)
  SS <- c(SS,SSE.model, SS_total)
  MS <- c(MS,MSE,NA)
  R2 <- c(R2,NA,NA)
  Fs <- c(Fs,NA,NA)
  ps <- as.double(pf(Fs, k - 1, rdf, lower.tail = FALSE))
  
  R2_all <- SS_reg / SS_total
  R2adj <- 1 - (1 - R2_all) * (nspecies - 1)/(rdf)
  Fstat <- rdf / (k-1)*R2_all/(1-R2_all)
  pval <- as.double(pf(Fstat, k - 1, rdf, lower.tail = FALSE))
  
  df <- c(df,NA,sum(df[1:(k-1)]))
  SS <- c(SS,NA,SS_reg)
  MS <- c(MS,NA,SS_reg/rdf)
  R2 <- c(R2,NA,R2_all)
  Fs <- c(Fs,NA,Fstat)
  ps <- c(ps,NA,pval)
  
  a.tab <- data.frame(df=df,SS=SS,MS=MS,Rsq=R2,F=Fs,P.value = ps)
  rownames(a.tab)[1:(k-1)] <- colnames(X)[2:k]
  rownames(a.tab)[k:(k+1)] <- c("Residuals","Total")
  rownames(a.tab)[(k+2):(k+3)] <- c("-----","Model")
  colnames(a.tab) <- c("df","SS","MS","Rsq","F","P.value")
  class(a.tab) <- c("anova", class(a.tab))
  attr(a.tab, "heading") <- paste("\nType I (Sequential) Sums of Squares and Cross-products\n","\nAnalysis of Variance Table\n")
  
  
  a.tab
}

print.rate.mult <- function(x,...)
{
  cat("\nEvolutionary rate(s):\n")
  cat(x$sigma2)
  cat("\n\nLog-likelihood: ",x$logL,"\n")
  cat("Method: ",x$method,"\n")
  cat("\nEvolutionary model: ")
  if(x$model[1]!="BM")
  {
    cat("\n")
    print(data.frame(Parameter=names(x$pars),Value=x$pars,row.names = x$model))
  } else cat("BM")
  cat("\n")
}

print.lr.test <- function(x,...)
{
  cat("Null model:\n")
  cat("\t",if(x$method=="REML") "Restricted log-" else "Log-") 
  cat("likelihood = ",x$null.logL)
  cat("\t","Number of parameters = ",x$null.pars)
  
  cat("\nAlternative model:\n")
  cat("\t",if(x$method=="REML") "Restricted log-" else "Log-") 
  cat("likelihood = ",x$alt.logL)
  cat("\t","Number of parameters = ",x$alt.pars)
  
  cat("\n\nDegrees of freedom = ",x$df) 
  cat("\nChi-square value = ",x$chi_sq)
  cat("\np-value = ",x$p)
}

# calculates Felsenstein's ancestral state reconstruction values
ace_hand <- function(x,Y,tree,gpr_fit=TRUE,gp_list)
{
  gp_missing <- missing(gp_list)
  nspecies <- length(tree$tip.label)
  if(gp_missing)
  {
    gp_list <- vector("list",nspecies)
    names(gp_list) <- tree$tip.label
  }
  if(gpr_fit) x <- normalize_to_01((x-min(x))/(max(x)-min(x)))
  
  Y <- Y[tree$tip.label,]
  Y <- rbind(Y,matrix(0,nrow=tree$Nnode,ncol=ncol(Y)))
  v1_lookup <- v2_lookup <- n1 <- n2 <- d1 <- d2 <- double(1)
  corrected_d <- tree$edge.length
  for(i in (tree$Nnode+nspecies):(nspecies+1))
  {
    v1_lookup <- which(tree$edge[,1]==i)[1]
    v2_lookup <- which(tree$edge[,1]==i)[2]
    n1 <- tree$edge[v1_lookup,2] # node 1 number
    n2 <- tree$edge[v2_lookup,2] # node 2 number
    d1 <- (1-(corrected_d[v1_lookup]/(corrected_d[v1_lookup]+corrected_d[v2_lookup])))
    d2 <- (1-(corrected_d[v2_lookup]/(corrected_d[v1_lookup]+corrected_d[v2_lookup])))
    if(gpr_fit)
    {
      if(!(gp_missing)) gp1 <- gp_list[[n1]] else gp1 <- GP_fit(x,Y[n1,])
      if(!(gp_missing)) gp2 <- gp_list[[n2]] else gp2 <- GP_fit(x,Y[n2,])
      if(n1<=nspecies) gp_list[[n1]] <- gp1
      if(n2<=nspecies) gp_list[[n2]] <- gp2
      time_warp <- dtw(predict(gp1,x)$Y_hat,predict(gp2,x)$Y_hat)
    } else
    {
      gp1 <- coef(glm(Y[n1,]~x,family=quasibinomial("logit")))
      gp2 <- coef(glm(Y[n2,]~x,family=quasibinomial("logit")))
      time_warp <- dtw(logit_fx(gp1,x),logit_fx(gp2,x))
    }
    
    new_x1 <- (((max(x)-min(x))*(time_warp$index1-min(time_warp$index1))) / (max(time_warp$index1) - min(time_warp$index1))) + min(x)
    new_x2 <- (((max(x)-min(x))*(time_warp$index2-min(time_warp$index2))) / (max(time_warp$index2) - min(time_warp$index2))) + min(x)
    if(gpr_fit)
    {
      new_y1 <- predict(gp1,new_x1)$Y_hat
      new_y2 <- predict(gp2,new_x2)$Y_hat
    } else
    {
      new_y1 <- logit_fx(gp1,new_x1)
      new_y2 <- logit_fx(gp2,new_x2)
    }
    anc_x <- (d1*new_x1) + (d2*new_x2)
    anc_y <- (d1*new_y1) + (d2*new_y2)
    if(gpr_fit)
    {
      anc_gp <- GP_fit(anc_x,anc_y)
      Y[i,] <- predict(anc_gp,x)$Y_hat
    } else
    {
      anc_gp <- coef(glm(anc_y~anc_x,family=quasibinomial("logit")))
      Y[i,] <- logit_fx(anc_gp,x)
    }
    corrected_d[which(tree$edge[,2]==i)] <- corrected_d[which(tree$edge[,2]==i)]+((corrected_d[v1_lookup]*corrected_d[v2_lookup])/(corrected_d[v1_lookup]+corrected_d[v2_lookup]))
  }
  anc_y <- Y[(length(tree$tip.label)+1):(tree$Nnode+length(tree$tip.label)),]
  return(list(anc_y=anc_y,gp_list=gp_list,anc_gp=anc_gp))
}

# calculates ancestral state reconstruction with pairwise dynamic time warping and 
# the fast PIC method (adapted from fastAnc in phytools)
fast_anc_hand <- function(x,Y,tree,root_only=TRUE,gpr_fit=TRUE)
{
  if(is.null(rownames(Y)))
  {
    rownames(Y) <- tree$tip.label
    warning("No taxa labels attached to coefficients. Assuming order matches tips on tree.")
  } else
  {
    temp <- match(tree$tip.label,rownames(Y))
    Y <- Y[temp,]
  }
  
  tip_gp <- vector("list",length = nrow(Y))
  
  M <- tree$Nnode
  N <- length(tree$tip.label)
  a <- tree
  node_vals <- matrix(0,M,ncol(Y))
  # calculate root state
  for(i in 1:M + N)
  {
    a <- multi2di(root(tree,node=(i)))
    if(i==(1+N))
    {
      res <- ace_hand(x,Y,tree=a,gpr_fit = gpr_fit)
      gp_list <- res$gp_list
      anc_Y <- res[[1]][1,]
      anc_gp <- res$anc_gp
    } else
    {
      res <- ace_hand(x,Y,tree=a,gpr_fit = gpr_fit,gp_list = gp_list)
    }
    node_vals[i-N,] <- res[[1]][1,]
    if(root_only)
    {
      #return(node_vals[1,])
      break
    }
  }
  rownames(node_vals) <- 1:M + N
  
  alignments <- vector("list",nrow(Y))
  names(alignments) <- names(res$gp_list)
  max_dtw <- rep(1,length(x))
  for(i in 1:nrow(Y))
  {
    alignments[[i]] <- dtw(Y[i,],anc_Y)
    temp_max <- ave(alignments[[i]]$index2,alignments[[i]]$index2,FUN = function(X) length(X))
    max_dtw[temp_max > max_dtw] <- temp_max[temp_max > max_dtw]
  }
  anc_align <- numeric()
  for(i in 1:length(x))
  {
    anc_align <- c(anc_align,rep(i,max_dtw[i]))
  }
  
  x_align <- matrix(0,nrow = nrow(Y),ncol = length(anc_align))
  
  for(j in 1:nrow(Y))
  {
    low <- 1
    for(i in 1:length(x))
    {
      temp_align <- alignments[[j]]$index1
      temp_ref <- alignments[[j]]$index2
      temp_align <- temp_align[temp_ref==i]
      temp_range <- (round(seq(min(temp_align),max(temp_align),length=max_dtw[i])))
      high <- low + max_dtw[i] - 1
      x_align[j,low:high] <-  temp_range
      low <- high + 1
    }
  }
  convert_x <- (x-min(x))/(max(x)-min(x))
  new_x <- apply(x_align,1,function(X) x[X])
  new_convert_x <- (new_x-min(x))/(max(x)-min(x))
  new_Y <- matrix(0,nrow = nrow(Y),ncol = nrow(new_x))
  for(i in 1:nrow(Y))
  {
    new_Y[i,] <- predict.GP(gp_list[[i]],xnew = new_convert_x[,i])$Y_hat
  }
  new_x <- t(new_x)
  rownames(new_Y) <- rownames(new_x) <- tree$tip.label
  aligned_data <- data.frame(species=tree$tip.label,new_x,new_Y,row.names=tree$tip.label)
  colnames(aligned_data) <- c("species",paste("x",1:ncol(new_x),sep=""),paste("y",1:ncol(new_Y),sep=""))
  aligned_coordinates <- data.frame(species = as.character(t(matrix(rep(tree$tip.label,ncol(new_x)),ncol=ncol(new_x)))),
                                    x = as.double(t(new_x)),
                                    y = as.double(t(new_Y)))
  return(list(aligned_data=aligned_data,aligned_coordinates=aligned_coordinates,aligned_X=new_x,aligned_Y=new_Y))
}

normalize_to_01 <- function(x)
{
  max_x <- max(x)
  min_x <- min(x)
  ret <- ((x-min_x) / (max_x - min_x))
  ret[ret<0] <- 0
  ret[ret>1] <- 1
  ret
}