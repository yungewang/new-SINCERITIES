function(DATA, distance = 1, method = 1, noDIAG = 0, SIGN = 1, normalize_distance = TRUE, standardize_distance = TRUE) {
  
  if (!distance %in% c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11)) {
    stop("Choose distance metric with 1-7: 
         1 for Kolmogorov-Smirnov (KM), 
         2 for Cramer-von Mises (CM), 
         3 for Anderson-Darling (AD), 
         4 for Mean Expression Difference, 
         5 for Earth Mover's Distance (EMD), 
         6 for Bhattacharyya Distance (BD), 
         7 for Kullback-Leibler Divergence (KL), 
         8 for Jensen-Shannon Divergence
         9 for Total Variation
         10 for Hellinger Distance
         11 for Pearson Divergence")
  }
  
  if (!method %in% c(1, 2, 3, 4, 5)) {
    stop("Choose regularization regression strategy with 1-4: 
         1 for RIDGE, 
         2 for ELASTIC NET with automatic detection of optimal alpha parameter, 
         3 for LASSO, 
         4 for ELASTIC NET with manual selection of alpha parameter
         5 for MCP")
  }
  
  if (!noDIAG %in% c(0, 1)) {
    stop("noDIAG should be either 0 or 1")
  }
  
  if (!SIGN %in% c(0, 1)) {
    stop("SIGN should be either 0 or 1")
  }
  
  # Initialization
  single_cell_data <- DATA$singleCELLdata
  time <- DATA$time
  numGENES <- DATA$numGENES
  num_time_points <- length(time)
  
  if (num_time_points < 5) {
    stop("** DATA with a number of time points < 5. Please run SINCERITIES_CROSS_VALIDATION function **")
  }
  
  # Distribution Distance
  DISTANCE_matrix <- matrix(data = 0, nrow = num_time_points - 1, ncol = numGENES)
  
  cmtest2 <- dget("/Users/yungewang/Downloads/SINCERITIES-R v2.0 2/SINCERITIES functions/cmtest2.R")
  
  for (ti in 1:(num_time_points - 1)) {
    data_ti <- t(single_cell_data[[ti]])
    data_ti_plus1 <- t(single_cell_data[[ti + 1]])
    
    for (gi in 1:numGENES) {
      p1 <- data_ti[gi, ]
      p2 <- data_ti_plus1[gi, ]
      
      # Compute distances based on the selected metric
      if (distance == 1) {
        test.stat <- ks.test(p1, p2)
        DISTANCE_matrix[ti, gi] <- test.stat$statistic
      } else if (distance == 2) {
        DISTANCE_matrix[ti, gi] <- cmtest2(p1, p2)$CM_limiting_stat
      } else if (distance == 3) {
        test.stat <- ad.test(p1, p2)
        DISTANCE_matrix[ti, gi] <- test.stat$ad[2, 1]
      } else if (distance == 4) {
        DISTANCE_matrix[ti, gi] <- abs(mean(p1) - mean(p2))
      } else if (distance == 5) {
        DISTANCE_matrix[ti, gi] <- emd::emd2d(p1, p2)
      } else if (distance == 6) {
        DISTANCE_matrix[ti, gi] <- -log(sum(sqrt(p1 * p2)))
      } else if (distance == 7) {
        p1_safe <- (p1 / sum(p1)) + 1e-10
        p2_safe <- (p2 / sum(p2)) + 1e-10
        DISTANCE_matrix[ti, gi] <- sum(p1_safe * log(p1_safe / p2_safe))
      } else if (distance == 8) {
        p1_safe <- (p1 / sum(p1)) + 1e-10
        p2_safe <- (p2 / sum(p2)) + 1e-10
        M <- 0.5 * (p1_safe + p2_safe)
        DISTANCE_matrix[ti, gi] <- 0.5 * sum(p1_safe * log(p1_safe / M) + p2_safe * log(p2_safe / M))
      } else if (distance == 9) {  # Total Variation Distance
        p1_safe <- p1 / sum(p1)
        p2_safe <- p2 / sum(p2)
        DISTANCE_matrix[ti, gi] <- 0.5 * sum(abs(p1_safe - p2_safe))
      } else if (distance == 10) {  # Hellinger Distance
        p1_safe <- (p1 / sum(p1)) + 1e-10
        p2_safe <- (p2 / sum(p2)) + 1e-10
        DISTANCE_matrix[ti, gi] <- sqrt(1 - sum(sqrt(p1_safe * p2_safe)))
      } else if (distance == 11) {  # Pearson Divergence
        p1_safe <- (p1 / sum(p1)) + 1e-10
        p2_safe <- (p2 / sum(p2)) + 1e-10
        DISTANCE_matrix[ti, gi] <- sum(p2_safe * ((p1_safe / p2_safe - 1)^2))
      }
    }
  }
  
  # Optional Normalization
  if (normalize_distance) {
    deltaT <- replicate(dim(DISTANCE_matrix)[2], time[2:length(time)] - time[1:(length(time) - 1)])
    DISTANCE_matrix <- DISTANCE_matrix / deltaT
  }
  
  # Optional Standardization
  if (standardize_distance) {
    mean_val <- mean(DISTANCE_matrix, na.rm = TRUE)
    sd_val <- sd(DISTANCE_matrix, na.rm = TRUE)
    DISTANCE_matrix <- (DISTANCE_matrix - mean_val) / sd_val
  }
  
  # Generate Y and X_matrix for glmnet
  if(method == 1) {
    alphas <- 0
  } else if(method == 2) {
    alphas <- seq(0, 1, 0.1)
  } else if(method == 3) {
    alphas <- 1
  } else if(method == 4) {
    input <- readline(" *** Please input manually the alpha values (between 0 and 1) separated by comma: ")
    alphas <- as.numeric(unlist(strsplit(input, ",")))
  }
  
  X_matrix <- DISTANCE_matrix[1:(num_time_points-2),]
  
  nfold <- dim(X_matrix)[1]
  foldid <- 1:nfold
  keep <- TRUE
  pred_lambda_min <- matrix(0, nrow=numGENES, ncol=numGENES)
  
  lambda_res <- vector()
  alpha_res <- vector()
  
  for (gi in 1:numGENES) {
    lambda <- vector()
    cvERROR <- vector()
    beta <- matrix(data=0, nrow=dim(X_matrix)[2], ncol=length(alphas))
    
    for (test in 1:length(alphas)) {
      Y_vector <- DISTANCE_matrix[2:(num_time_points-1), gi]
      if(noDIAG == 1) {
        CV_results <- cv.glmnet(X_matrix, Y_vector, alpha=alphas[test], exclude=gi, nfolds=nfold, foldid=foldid, 
                                keep=keep, lower.limits=0, upper.limits=Inf, grouped=FALSE)
      } else {
        CV_results <- cv.glmnet(X_matrix, Y_vector, alpha=alphas[test], nfolds=nfold, foldid=foldid, 
                                keep=keep, lower.limits=0, upper.limits=Inf, grouped=FALSE)
      }
      lambda[test] <- CV_results$lambda.min
      cvERROR[test] <- CV_results$cvm[CV_results$lambda == CV_results$lambda.min]
      coef.CV_results <- coef(CV_results, s='lambda.min')
      beta[coef.CV_results@i[-1], test] <- coef.CV_results@x[-1]
      
      # Print lambda for each alpha value
      cat(paste("Gene", gi, "- Alpha:", alphas[test], "- Lambda:", lambda[test], "\n"))
    }
    
    minIdx <- max(which(cvERROR == min(cvERROR)))
    lambda_res[gi] <- lambda[minIdx]
    alpha_res[gi] <- alphas[minIdx]
    pred_lambda_min[, gi] <- beta[, minIdx]
    
    cat("Summary of Lambda Values:\n")
    print(lambda_res)
  }
  
  if(SIGN == 1) {
    parcorr_matrix <- pcor(DATA$totDATA, method='pearson')$estimate
    pred_lambda_min <- pred_lambda_min * sign(parcorr_matrix)
  }
  
  result <- list(DISTANCE_matrix=DISTANCE_matrix, adj_matrix=pred_lambda_min)
  return(result)
  
  
}
