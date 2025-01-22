library(glmnet)
library(ppcor)
library(igraph)

process_window <- function(window_indices, DATA, distance_matrix, noDIAG, window_label) {
  num_genes <- DATA$numGENES

  window_data <- do.call(rbind, DATA$singleCELLdata[window_indices])
  partial_corr <- pcor(window_data, method = "pearson")$estimate

  alpha_matrix <- matrix(0, nrow = num_genes, ncol = num_genes)
  lambda_values <- numeric(num_genes)
  
  for (gi in 1:num_genes) {

    X_matrix <- distance_matrix[1:(nrow(distance_matrix) - 1), ]
    Y_vector <- distance_matrix[2:nrow(distance_matrix), gi]
    
    CV_results <- cv.glmnet(as.matrix(X_matrix), Y_vector, alpha = 1,
                            nfolds = min(3, nrow(X_matrix)), lower.limits = 0, upper.limits = Inf, grouped = FALSE)

    best_lambda <- CV_results$lambda.min
    lambda_values[gi] <- best_lambda
    best_coefficients <- as.numeric(coef(CV_results, s = "lambda.min"))[-1]

    if (noDIAG == 1) {
      best_coefficients[gi] <- 0
    }
    best_coefficients <- best_coefficients * sign(partial_corr[, gi])
    alpha_matrix[, gi] <- best_coefficients
  }

  cat(paste("Lambda Values for", window_label, ":\n"))
  print(lambda_values)
  
  cat(paste("Alpha Matrix for", window_label, ":\n"))
  print(alpha_matrix)
  
  max_value <- max(alpha_matrix, na.rm = TRUE)
  adj_matrix_normalized <- alpha_matrix / max_value  # Normalize adjacency matrix
  
  genes <- sprintf("G%d", 1:nrow(adj_matrix_normalized))  

  table <- final_ranked_predictions(adj_matrix_normalized, genes, SIGN = 1, 
                                    fileNAME = paste0("results_", window_label), saveFile = FALSE)
  
  graph <- graph_from_adjacency_matrix(t(adj_matrix_normalized), mode = "directed", weighted = TRUE)
  
  E(graph)$color <- ifelse(E(graph)$weight > 0, "black", "red")  # Activation: Black, Repression: Red
  E(graph)$lty <- ifelse(E(graph)$color == "red", 2, 1)  # Red edges dashed, black edges solid

  layout <- layout_in_circle(graph)

  plot(
    graph,
    layout = layout,
    vertex.size = 20,               # Size of nodes
    vertex.color = "skyblue",       # Color of nodes
    vertex.label = V(graph)$name,   # Label nodes with gene names (G1, G2, ...)
    vertex.label.cex = 0.8,         # Label size
    edge.arrow.size = 0.3,          # Arrow size
    edge.color = E(graph)$color,    # Edge colors
    edge.lty = E(graph)$lty,        # Edge line type
    main = window_label
  )
  
  return(list(alpha_matrix = alpha_matrix, lambda_values = lambda_values))
}
