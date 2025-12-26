# Function to compute additional sample size needed to adjust for covariates, eq (26)
size_covadjust <- function(ntot, pibar, vpi, p, q, bb) {
  b <- bb / sum(bb)
  d <- 1 + q / (ntot * pibar - 2)
  e <- numeric(p)
  for (j in 1:p) {
    e[j] <- sum(d[j] - (vpi[1:j] * d[1:j]) / vpi[j])
  }
  size <- 0
  for (j in 1:p) {
    size <- size + b[j] * (d[j] + e[j] / (ntot * pibar[j] - j + 1))
  }
  return(ntot * (size - 1))
}

# Function to compute estimated KR adjusted variance and degrees of freedom (f), eq. (21) 
KR_var_df <- function(ntot, pibar, vpi, l2sig, qs, p, q) {
  m <- ntot * pibar
  varx <- (vpi / ntot) * (1 + q / (ntot * pibar - q - 3))
  c <- numeric(p)
  for (j in 1:p) {
    fac <- 1 - (j - 1) / (m[j] - qs)
    sum_term <- 0
    if (j < p) {
      for (k in (j+1):p) {
        sum_term <- sum_term + l2sig[k] / (m[k] - qs - k)
      }
    }
    c[j] <- fac * (l2sig[j] + sum_term)
  }
  num <- sum(c * varx) ^ 2
  denom1 <- 0
  for (j in 2:p) {
    inner_sum <- 0
    for (t in 1:(j-1)) {
      inner_sum <- inner_sum + c[t] * varx[t] ^ 2
    }
    denom1 <- denom1 + c[j] * inner_sum / (m[j] - qs - j)
  }
  denom1 <- 2 * denom1
  denom2 <- sum((c ^ 2) * (varx ^ 2) / (m - qs))
  f <- num / (denom1 + denom2)
  vartau <- sum(c * varx)
  for (j in 2:p) {
    inner_sum <- 0
    for (t in 1:(j-1)) {
      inner_sum <- inner_sum + (varx[j] - varx[t])
    }
    vartau <- vartau + 2 * c[j] * inner_sum / (m[j] - qs)
  }
  return(c(vartau, f))
}

# Function to compute sample size / power for MMRM
power <- function(covrentdata, alpha, power, ntot, gamma0, 
                  p, q, taup, type, M0, Ml, Mu) {
  
  # Set Ml to -Mu if missing (for equivalence margin)
  if (is.na(Ml)) {Ml <- -Mu}
  
  # Check input parameters
  if (!(toupper(type) %in% c("SUP", "NI", "EQUI"))) {
    stop("ERROR: Type should be SUP, NI, or EQUI.")
  }
  else if ((!is.na(power) & (power <= 0 | power >= 1)) |
      (!is.na(ntot) & ntot <= 0) |
      (is.na(power) & is.na(ntot)) |
      (!is.na(power) & !is.na(ntot))) {
    stop("ERROR: Either 0 < power < 1 with ntot = NA, OR ntot > 0 and power = NA")
  }
  else if (p <= 0 | as.integer(p) != p) {
    stop("ERROR: Number of visits (p) must be a positive integer.")
  }
  else if (q < 0 | as.integer(q) != q) {
    stop("ERROR: Number of covariates (q) must be a non-negative integer.")
  }
  else if (toupper(type) == "SUP" & taup == 0) {
    stop("ERROR: Effect at last visit (taup) should not be 0 in a superiority trial.")
  }
  else if (toupper(type) == "EQUI" &
      (is.na(Mu) | 
       is.na(Ml) | 
       taup >= Mu |
       taup <= Ml |
       Ml >= Mu)) {
    stop("ERROR: Equivalence margin must satisfy Ml < taup < Mu with Mu not missing.")
  }
  else if (toupper(type) == "NI" & (is.na(M0) | taup == M0)) {
    stop("ERROR: NI margin must satisfy M0 not missing and M0 != taup.")
  }
  
  # Define variables
  covpi <- covrentdata
  nvis <- nrow(covpi)
  sigma <- covpi[,1:nvis]
  rent0 <- covpi[,nvis+1]
  rent1 <- covpi[,nvis+2]
  rootb <- t(chol(sigma))
  dia <- diag(rootb)
  L <- rootb %*% diag(1 / dia)
  sig <- diag(dia ^ 2)
  qs <- q + 2
  gamma1 <- 1 - gamma0
  pibar <- gamma0 * rent0 + gamma1 * rent1
  vpi <- 1 / (gamma0 * rent0) + 1 / (gamma1 * rent1)
  l2sig <- L[p,]^2 %*% sig
  bb <- l2sig * vpi
  vtau <- sum(bb)
  type <- tolower(type)
  ntot_given <- ntot
  power_given <- power
  
  # Option 1: compute sample size if ntot is missing
  if (is.na(ntot_given)) {
    if (type == "sup") {
      effect <- taup
    }
    else if (type == "ni") {
      effect <- taup - M0
    }
    else if (type == "equi") {
      effect1 <- taup - Ml
      effect2 <- Mu - taup
      effect <- min(effect1, effect2)
    }
    # Asymptotic normal sample size approximation, eq (24)
    norm_asy <- (qnorm(power_given) + qnorm(1 - alpha/2)) ^ 2 * vtau / (effect ^ 2)
    # Normal approximation adjusted for covariates, eq (25)
    norm_adj <- norm_asy + size_covadjust(ntot=norm_asy, pibar, vpi, p, q, bb)
    # Covariate adjusted sample size based on T distribution, eq (5)
    varparm <- KR_var_df(ntot=norm_adj, pibar, vpi, l2sig, qs, p, q)
    rho <- varparm[2] / (norm_adj * pibar[1] - qs)
    t_adj <- norm_adj + qnorm(1 - alpha/2) ^ 2 / (2 * rho)
    # Store results - currently using normal approximation
    norm_adj <- ceiling(norm_adj)
    t_adj <- ceiling(t_adj)
    final_ntot <- norm_adj # (if t-dist. approx. preferred, replace 'norm_adj' with 't_adj')
    # If 1:1 allocation, return total sample size as even number
    if (gamma0 == 0.5) {
      if (final_ntot %% 2 != 0) {final_ntot <- final_ntot + 1}
    }
    final_power <- power_given
  }
  
  # Option 2: compute power if power is missing
  else if (is.na(power_given)) {
    varparm <- KR_var_df(ntot=ntot_given, pibar, vpi, l2sig, qs, p, q)
    vartau <- varparm[1]
    df <- varparm[2]
    if (type == "sup") {
      effect <- taup
      noncen <- abs(effect) / sqrt(vartau)
      tcrit <- qt(1 - alpha/2, df)
      power <- 1 - pt(tcrit, df, ncp = noncen)
    }
    else if (type == "ni") {
      effect <- taup - M0
      noncen <- abs(effect) / sqrt(vartau)
      tcrit <- qt(1 - alpha/2, df)
      power <- 1 - pt(tcrit, df, ncp = noncen)
    }
    else if (type == "equi") {
      effect1 <- Mu - taup
      effect2 <- taup - Ml
      noncen1 <- abs(effect1) / sqrt(vartau)
      noncen2 <- abs(effect2) / sqrt(vartau)
      tcrit <- qt(1 - alpha/2, df)
      power <- 1 - pt(tcrit, df, ncp = noncen1) - pt(tcrit, df, ncp = noncen2)
    }
    # Store results
    final_ntot <- ntot_given
    final_power <- round(power, 3) # round to 3 decimal places
  }
  
  # Create dataframe for results
  if (is.na(ntot_given)) {
    power <- final_power
    ntot_z <- norm_adj
    ntot_t <- t_adj
    ntot <- final_ntot
    note <- "Sample size computed"
    res <- data.frame(power, ntot, note)
  }
  else if (is.na(power_given)) {
    power <- final_power
    ntot <- final_ntot
    note <- "Power computed"
    res <- data.frame(power, ntot, note)
  }
  
  # Output dataframe
  return(res)
  
}