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
    
    # meta_data$y <- as.numeric(y)
    # v <- vglm(y~Genotype, family = tobit(Lower = d[i]), data = meta_data)

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
  res <- data.frame(beta_est, sigma = sigma_est, likelihood = likelihoods, d = d)
  
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

Tobit_model2 <- function(
  metabolites,
  meta_data,
  formula,
  ref = "WT",
  test = "two-sided",
  presence_absence = TRUE,
  # p_adjust = "fdr",
  seed
) {
  X <- model.matrix(as.formula(formula), data = meta_data)
  p <- nrow(metabolites)
  k <- ncol(X)
  n <- nrow(X)
  meta_data <- data.frame(sample_id = rownames(meta_data), meta_data)

  # Setting detection limits an initial values for sigma
  d <- apply(metabolites, 1, function(x) min(log(x)[x != 0]) - 0.1)

  empty_matrix <- function(metabolites, X) {
    k <- ncol(X)
    empty_mat <- matrix(NA, nrow(metabolites), k)
    colnames(empty_mat) <- colnames(X)
    rownames(empty_mat) <- rownames(metabolites)
    return(empty_mat)
  }

  beta_est_pa <- NULL
  p_vals_pa <- NULL

  if (presence_absence) {
    beta_est_pa <- empty_matrix(metabolites, X)
    p_vals_pa <- empty_matrix(metabolites, X)
  }

  # Estimating beta and sigma
  sigma_est <- rep(NA, nrow(metabolites))
  beta_est <- empty_matrix(metabolites, X)
  p_vals <- empty_matrix(metabolites, X)
  # resid_mat <- matrix(NA, nrow(metabolites), ncol(metabolites))
  resid_list <- list()

  if (!missing(seed)) {
    set.seed(seed)
  }

  for (i in 1:p) {
    cat(i, "\r")
    y <- t(log(metabolites[i, ]))
    y <- data.frame(sample_id = rownames(y), y)
    df <- base::merge(meta_data, y, by = "sample_id")
    colnames(df)[ncol(df)] <- "intensity"

    dt <- data.table(df)
    dt <- dt[order(Genotype)]
    dt[, pres := ifelse(intensity == -Inf, 0, 1)]
    check <- dt[, .(all_zero = mean(intensity == -Inf) == 1), Genotype]

    if (presence_absence) {
      pres_abs <- logistf::logistf(pres ~ Genotype, data = dt)
      beta_est_pa[i, ] <- pres_abs$coefficients
      p_vals_pa[i, ] <- pres_abs$prob
    }

    if (check$all_zero[1] == FALSE & all(check$all_zero[-1]) == FALSE) {
      pres_gt <- check$Genotype[!check$all_zero]
      dt <- dt[Genotype %in% pres_gt]

      v <- vglm(
        intensity ~ Genotype,
        data = dt,
        family = tobit(Lower = d[i], imethod = 2)
      )
      res <- coef(summary(v))
      rownames(res) <- gsub(":1", "", rownames(res))
      idx <- rownames(res)[-2]

      if (test == "one-sided") {
        p_vals_tobit <- pt(-res[idx, 3], df = nrow(dt) - ncol(X))
      } else if(test == "two-sided"){
        p_vals_tobit <- 2*pt(-abs(res[idx, 3]), df = nrow(dt) - ncol(X))
      }

      sigma_est[i] <- exp(res[2, 1])
      beta_est[i, idx] <- res[idx, 1]
      p_vals[i, idx] <- p_vals_tobit

      eta <- predict(v)
      mu <- eta[, 1]
      sigma <- exp(eta[, 2])
      y_obs <- y[, 2]

      resid_rquant <- pnorm(y_obs, mean = mu, sd = sigma)
      cens_idx <- y_obs == -Inf
      resid_rquant[y_obs == -Inf] <- runif(
        n = sum(y_obs == -Inf),
        0,
        pnorm(d[i], mean = mu, sd = sigma)
      )
      resid_rquant <- qnorm(resid_rquant, mean = 0, sd = 1)

      resid_list[[i]] <- resid_rquant
    }
  }

  result <- data.frame(beta_est, sigma = sigma_est, d = d)

  out <- list(
    res = result,
    p_vals = p_vals,
    beta_est_pa = beta_est_pa,
    p_vals_pa = p_vals_pa,
    resid = resid_list
  )

  return(out)
}

# Experiment with brtobit
Tobit_model3 <- function(
  metabolites,
  meta_data,
  formula,
  ref = "WT",
  test = "two-sided",
  # p_adjust = "fdr",
  seed
) {
  X <- model.matrix(as.formula(formula), data = meta_data)
  p <- nrow(metabolites)
  k <- ncol(X)
  n <- nrow(X)
  meta_data <- data.frame(sample_id = rownames(meta_data), meta_data)

  # Setting detection limits an initial values for sigma
  d <- apply(metabolites, 1, function(x) min(log(x)[x != 0]) - 0.1)

  empty_matrix <- function(metabolites, X) {
    k <- ncol(X)
    empty_mat <- matrix(NA, nrow(metabolites), k)
    colnames(empty_mat) <- colnames(X)
    rownames(empty_mat) <- rownames(metabolites)
    return(empty_mat)
  }

  # Estimating beta and sigma
  sigma_est <- rep(NA, nrow(metabolites))
  beta_est <- empty_matrix(metabolites, X)
  p_vals <- empty_matrix(metabolites, X)
  # resid_mat <- matrix(NA, nrow(metabolites), ncol(metabolites))
  resid_list <- list()

  if (!missing(seed)) {
    set.seed(seed)
  }

  for (i in 1:p) {
    cat(i, "\r")
    y <- t(log(metabolites[i, ]))
    y <- data.frame(sample_id = rownames(y), y)
    df <- base::merge(meta_data, y, by = "sample_id")
    colnames(df)[ncol(df)] <- "intensity"

    dt <- data.table(df)
    dt <- dt[order(Genotype)]
    check <- dt[, .(all_zero = mean(intensity == -Inf) == 1), Genotype]

    if (!any(check$all_zero)) {
      v <- vglm(
        intensity ~ Genotype,
        data = dt,
        family = tobit(Lower = d[i], imethod = 2)
      )
      res <- coef(summary(v))
      rownames(res) <- gsub(":1", "", rownames(res))

      sigma <- exp(res[2, 1])
      mu <- predict(v)[, 1]
      res <- res[-2, ]
    } else {
      dt[, intensity := intensity - d[i]]
      dt[intensity == -Inf, intensity := 0]

      b <- tryCatch(
        brtobit(
          intensity ~ Genotype,
          data = dt
        ),
        error = function(e) {
          return(NULL)
        }
      )

      if (is.null(b)) {
        m <- lm.fit(model.matrix(~Genotype, dt), dt$intensity)
        b <- brtobit(
          intensity ~ Genotype,
          data = dt,
          start = c(m$coefficients, "(Variance)" = var(dt$intensity))
        )
      }

      res <- coef(summary(b))

      sigma <- sqrt(res["(Variance)", 1])
      mu <- predict(b) + d[i]
      res <- res[-nrow(res), ]

      res[1, 1] <- res[1, 1] + d[i]
      res[1, 3] <- res[1, 1] / res[1, 2]
      res[1, 4] <- 2 * pnorm(-abs(res[1, 3]), sd = res[1, 2])
    }

    idx <- rownames(res)

    if (test == "one-sided") {
      p_vals_tobit <- pnorm(-res[idx, 3])
    } else if (test == "two-sided") {
      p_vals_tobit <- res[idx, 4]
    }

    sigma_est[i] <- sigma
    beta_est[i, idx] <- res[idx, 1]
    p_vals[i, idx] <- p_vals_tobit[idx]

    y_obs <- y[, 2]
    resid_rquant <- pnorm(y_obs, mean = mu, sd = sigma)
    cens_idx <- y_obs == -Inf
    resid_rquant[y_obs == -Inf] <- runif(
      n = sum(y_obs == -Inf),
      0,
      pnorm(d[i], mean = mu, sd = sigma)
    )
    resid_rquant <- qnorm(resid_rquant, mean = 0, sd = 1)

    resid_list[[i]] <- resid_rquant
  }

  result <- data.frame(beta_est, sigma = sigma_est, d = d)

  out <- list(
    res = result,
    p_vals = p_vals,
    resid = resid_list
  )

  return(out)
}
