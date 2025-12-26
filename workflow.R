#########################################################################################
## Workflow: Power / Sample Size Calculations for MMRM in Unstratified Clinical Trials. 
## Reference: Tang (2018), Statistics in Medicine.
#########################################################################################

# Read in defined functions
source("cov_matrix.R")
source("power_calc.R")
source("plot_power.R")

# Wrapper function for workflow
run_workflow <- function(

# Enter covariance parameters
times = c(8, 18, 28, 40, 52),           # Post-baseline visit times
stds = c(175, 200, 225, 250, 275),      # Standard deviations at each visit
cov_type = "AR",                        # Covariance structure (AR or CS)
rho = 0.8,                              # Correlation parameter (0 < rho < 1)

# Enter scenarios to be considered for trial (can list multiple)
dropout_list = c(0, 0.25, 0.3),         # Dropout rates  
margin_list = c(0),                     # Non-inferiority margins (<0) or equivalence margins (>0) or set to 0 for superiority
diff_list = c(50, 55, 60),              # Treatment effects at last visit
power_list = c(0.8, 0.9),               # Study powers or set to NA if want computed
ntot_list = c(NA),                      # Total sample sizes or set to NA if want computed

# Enter additional trial parameters
type = "sup",                            # Type of trial (sup, ni, or equi)
q = 3,                                  # Number of model covariates (excluding intercept and treatment)
alpha = 0.05,                           # Type I error rate
gamma0 = 0.5,                           # Proportion of total participants allocated to control group

# Enter plotting parameters
plot_yn = "Y",                          # Plot results? (Y/N)
title = NA,                             # Enter custom title for plots or NA for default
xlabel = NA,                            # Enter custom x-axis label or NA for default
ylabel = NA                             # Enter custom y-axis label or NA for default

) {

#########################################################################################

  # Covariance matrix
  cov <- cov_matrix(times=times, stds=stds, type=cov_type, rho=rho)
  res_cov <- cov
  
  # Create useful variables
  p <- length(times)
  time_col <- p + 1
  last_time <- times[p]
  
  # Construct grid of scenarios
  param <- expand.grid(dropout = dropout_list, 
                       margin = margin_list, 
                       diff = diff_list, 
                       pow = power_list, 
                       ntot = ntot_list)
  nparm <- nrow(param)
  
  # Initialize data frame for results
  powall <- as.data.frame(matrix(NA, nrow = nparm, ncol = 6))
  colnames(powall) <- c("dropout", "margin", "diff", 
                        "power", "ntot", "note")
  
  # Perform sample size / power calculation for each scenario
  for (i in 1:nparm) {
    
    # Calculate retention data based on exponential dropout rate
    lambda <- -log(1 - param$dropout[i]) / last_time
    rent0 <- exp(-lambda * times)
    rent1 <- rent0
    covtemp <- cbind(cov[,1:p], rent0, rent1)
    
    # Call function to compute sample size or power
    powtemp <- power(covrentdata=covtemp, alpha=alpha, power=param$pow[i],
                     gamma0=gamma0, p=p, q=q, taup=param$diff[i],
                     type=type, ntot=param$ntot[i], M0=param$margin[i],
                     Ml=NA, Mu=param$margin[i])
    powtemp <- data.frame(dropout=param$dropout[i], margin=param$margin[i], 
                          diff=param$diff[i], powtemp)
    powall[i,] <- powtemp
    
  }
  
  # Results table
  res_table <- powall
  
  # Plot results
  if (toupper(plot_yn) == "Y") {
    
    # Y-axis variable and label
    if (all(is.na(ntot_list))) {
      yvar <- "ntot"
      if (is.na(ylabel)) {
        ylabel <- "Total sample size for MMRM"
      }
    }
    else if (all(is.na(power_list))) {
      yvar <- "power"
      if (is.na(ylabel)) {
        ylabel <- "Power for MMRM"
      }
    }
    
    # Initialize plot list
    p_list <- list()
    
    # Superiority trials
    if (toupper(type) == "SUP") {
      
      # X-axis variable and label
      xvar <- "diff"
      if (is.na(xlabel)) {
        xlabel <- paste("Difference at week", last_time, sep = " ")
      }
      
      # Plot sample size results
      if (all(is.na(ntot_list))) {
        n_pow <- length(power_list)
        for (k in 1:n_pow) {
          power_tmp <- power_list[k]
          if (is.na(title)) {
            title_tmp <- paste("Target power ",
                               power_tmp, 
                               " for superiority trial, ",
                               toupper(cov_type), "(", rho, ") covariance",
                               sep = "")
          }
          else {title_tmp <- title}
          p <- plot_power(data=powall, xvar=xvar, yvar=yvar,
                     pow=power_tmp, n=NA,
                     xlabel=xlabel, ylabel=ylabel, title=title_tmp)
          p_list[[length(p_list)+1]] <- p
        }
      }
      
      # Plot power results
      else if (all(is.na(power_list))) {
        n_ntot <- length(ntot_list)
        for (k in 1:n_ntot) {
          ntot_tmp <- ntot_list[k]
          if (is.na(title)) {
            title_tmp <- paste("Total sample size ", ntot_tmp,
                               " for superiority trial, ",
                               toupper(cov_type), "(",
                               rho, ") covariance", sep = "")
          }
          else {title_tmp <- title}
          p <- plot_power(data=powall, xvar=xvar, yvar=yvar,
                     pow=NA, n=ntot_tmp,
                     xlabel=xlabel, ylabel=ylabel, title=title_tmp)
          
          p_list[[length(p_list)+1]] <- p
        }
      }
      
  
      
    }
    
    # Non-inferiority trials
    else if (toupper(type) == "NI") {
      
      # X-axis variable and label
      xvar <- "margin"
      if (is.na(xlabel)) {
        xlabel <- "Non-inferority margin"
      }
    
      # Plot sample size results
      if (all(is.na(ntot_list))) {
        n_pow <- length(power_list)
        n_diff <- length(diff_list)
        for (k in 1:n_pow) {
          power_tmp <- power_list[k]
          for (j in 1:n_diff) {
            diff_tmp <- diff_list[j]
            if (is.na(title)) {
              title_tmp <- paste("Target power ", power_tmp, 
                                 " for non-inferiority trial, difference of ",
                                 diff_tmp, ", ", toupper(cov_type), "(",
                                 rho, ")", sep = "")
            }
            else {title_tmp <- title}
            powall_tmp <- powall |> filter(diff == diff_tmp)
            p <- plot_power(data=powall_tmp, xvar=xvar, yvar=yvar,
                       pow=power_tmp, n=NA,
                       xlabel=xlabel, ylabel=ylabel, title=title_tmp)
            p_list[[length(p_list)+1]] <- p
          }
        }
      }
    
      # Plot power results
      else if (all(is.na(power_list))) {
        n_ntot <- length(ntot_list)
        n_diff <- length(diff_list)
        for (k in 1:n_ntot) {
          ntot_tmp <- ntot_list[k]
          for (j in 1:n_diff) {
            diff_tmp <- diff_list[j]
            if (is.na(title)) {
              title_tmp <- paste("Total sample size ", ntot_tmp,
                                 " for non-inferiority trial, difference of ",
                                 diff_tmp, ", ", toupper(cov_type), "(", rho, ") covariance",
                                 sep = "")
            }
            else {title_tmp <- title}
            powall_tmp <- powall |> filter(diff == diff_tmp)
            p <- plot_power(data=powall_tmp, xvar=xvar, yvar=yvar, 
                       pow=NA, n=ntot_tmp, 
                       xlabel=xlabel, ylabel=ylabel, title=title_tmp)
            p_list[[length(p_list)+1]] <- p
          }
        }
      }
      
    }
    
    # Equivalence trials
    else if (toupper(type) == "EQUI") {
    
      # X-axis variable and label
      xvar <- "margin"
      if (is.na(xlabel)) {
        xlabel <- "Equivalence margin"
      }
    
      # Plot sample size results
      if (all(is.na(ntot_list))) {
        n_pow <- length(power_list)
        n_diff <- length(diff_list)
        for (k in 1:n_pow) {
          power_tmp <- power_list[k]
          for (j in 1:n_diff) {
            diff_tmp <- diff_list[j]
            if (is.na(title)) {
              title_tmp <- paste("Target power ", power_tmp, " for equivalence trial, difference of ", 
                                 diff_tmp, ", ", toupper(cov_type), "(", rho, ") covariance", sep = "")
            }
            else {title_tmp <- title}
            powall_tmp <- powall |> filter(diff == diff_tmp)
            p <- plot_power(data=powall_tmp, xvar=xvar, yvar=yvar,
                       pow=power_tmp, n=NA,
                       xlabel=xlabel, ylabel=ylabel, title=title_tmp)
            p_list[[length(p_list)+1]] <- p
          }
        }
      }
    
      # Plot power results
      else if (all(is.na(power_list))) {
        n_ntot <- length(ntot_list)
        n_diff <- length(diff_list)
        for (k in 1:n_ntot) {
          ntot_tmp <- ntot_list[k]
          for (j in 1:n_diff) {
            diff_tmp <- diff_list[j]
            if (is.na(title)) {
              title_tmp <- paste("Total sample size ", ntot_tmp, " for equivalence trial, difference of ",
                                 diff_tmp, ", ", toupper(cov_type), "(", rho, ") covariance", sep = "")
            }
            else {title_tmp <- title}
            powall_tmp <- powall |> filter(diff == diff_tmp)
            p <- plot_power(data=powall_tmp, xvar=xvar, yvar=yvar,
                       pow=NA, n=ntot_tmp,
                       xlabel=xlabel, ylabel=ylabel, title=title_tmp)
            p_list[[length(p_list)+1]] <- p
          }
        }
      }
      
    }

  }
  
  if (toupper(plot_yn) == "Y") {
    return(list(
      covmat = res_cov,
      table = res_table,
      plots = p_list
    ))
  } else {
    return(list(
      covmat = res_cov,
      table = res_table
    ))
  }
  
}

## Print results
# x <- run_workflow()
# print(x$covmat)
# print(x$table)
# print(x$plots)
