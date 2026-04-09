Tobit_model <- function(metabolites, meta_data, formula){
  X <- model.matrix(as.formula(formula), data = meta_data)
  p <- nrow(metabolites)
  k <- ncol(X)
  n <- nrow(X)
  
  # Initial values for beta
  beta_init <- matrix(NA, p, k)
  for(i in 1:p){
    y <- as.numeric(metabolites[i,])
    y[y==0] <- min(y[y != 0])/2
    y <- log(y)
    l <- lm(y~0+X)
    beta_init[i,] <- l$coefficients
  }
  colnames(beta_init) <- colnames(X)
  
  # Setting detection limits an initial values for sigma
  d <- apply(metabolites, 1, function(x) min(log(x)[x!=0])-0.1)
  sigma_init <- apply(log(metabolites), 1, function(x) sd(x[x!=-Inf]))
  sigma_init_log <- log(sigma_init)
  sigma_init_log[sigma_init_log<log(0.2)] <- log(0.2)
  sigma_init_log[is.na(sigma_init_log)] <- log(0.5)
  
  # Estimating beta and sigma
  sigma_est <- rep(NA, nrow(metabolites))
  likelihoods <- rep(NA, nrow(metabolites))
  beta_est <- matrix(NA, nrow(metabolites), k)
  for(i in 1:nrow(metabolites)){
    cat(i, "\r")
    y <- log(metabolites[i,])
    likelihood <- function(param){
      sigma <- exp(param[1])
      beta <- param[-1]
      zero_idx <- y == -Inf
      l_nonzero <- sum( log(1/sigma*dnorm((y-X%*%beta)[!zero_idx]/sigma)) )
      l_zero <- sum( log(pnorm((d[i]-X%*%beta)[zero_idx]/sigma)) )
      l <- l_nonzero + l_zero
      return(-l)
    }
    
    tryCatch(
      optim(par = c(sigma_init_log[i], beta_init[i,]), fn = likelihood),
      error = function(e) {
        message("An error occurred. Skipping function call.")
        return(NULL)
      }
    ) -> opt
    
    if(!is.null(opt)){
      sigma_est[i] <- exp(opt$par[1])
      beta_est[i,] <- opt$par[-1]
      likelihoods[i] <- opt$value
    }
  }
  
  colnames(beta_est) <- colnames(X)
  rownames(beta_est) <- rownames(metabolites)
  res <- data.frame(beta_est, sigma = sigma_est,
                    likelihood = likelihoods, d = d)
  
  out <- list(res = res, X = X, metabolites = metabolites, beta_init = beta_init)
  
  return(out)
}

# LRT ----
Tobit_LRT <- function(TT, var, p_adjust = "fdr"){
  X <- TT$X
  cols_to_remove <- which(colnames(X) %in% var)
  X <- X[,-cols_to_remove, drop = F]
  k <- ncol(X)
  n <- nrow(X)
  
  d <- TT$res$d
  beta_init <- TT$res[,1:ncol(TT$X)][,-cols_to_remove, drop=F]
  sigma_init_log <- log(TT$res$sigma)
  
  sigma_est <- rep(NA, nrow(TT$metabolites))
  likelihoods <- rep(NA, nrow(TT$metabolites))
  beta_est <- matrix(NA, nrow(TT$metabolites), k)
  for(i in 1:nrow(TT$metabolites)){
    cat(i, "\r")
    if(!is.na(TT$res$likelihood[i])){
      y <- log(TT$metabolites[i,])
      # y <- log(RA[,i])
      likelihood <- function(param){
        sigma <- exp(param[1])
        beta <- param[-1]
        zero_idx <- y == -Inf
        l_nonzero <- sum( log(1/sigma*dnorm((y-X%*%beta)[!zero_idx]/sigma)) )
        l_zero <- sum( log(pnorm((d[i]-X%*%beta)[zero_idx]/sigma)) )
        l <- l_nonzero + l_zero
        return(-l)
      }
      # opt <- optim(par = c(sigma_init_log[i], TT$beta_init[-cols_to_remove,i]),
      #              fn = likelihood)
      tryCatch(
        optim(par = c(sigma_init_log[i], beta_init[i,]), fn = likelihood),
        error = function(e) {
          message("An error occurred. Skipping function call.")
          return(NULL)
        }
      ) -> opt
      # opt <- optim(par = c(sigma_init_log[i], c(mean(d), 0, 0.72, 1.35,-0.81,-11)),
      #              fn = likelihood)
      if(!is.null(opt)){
        sigma_est[i] <- exp(opt$par[1])
        beta_est[i,] <- opt$par[-1]
        likelihoods[i] <- opt$value
      }
    }
  }
  LR_stat <- -2*(TT$res$likelihood - likelihoods)
  p_vals <- 1-pchisq(LR_stat, df = ncol(TT$X)-ncol(X))
  p_vals[is.na(p_vals)] <- 1
  
  res <- cbind(
    TT$res, p_vals = p_vals, p_adj = p.adjust(p_vals, method = p_adjust)
  )
  
  return(res)
}
