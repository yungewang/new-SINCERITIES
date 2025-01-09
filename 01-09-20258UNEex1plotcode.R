## 1. 01/08/2025 example code for the example 1 in 8UNE dataset, since this datset has 20 examples(networks),8 time points:
#the golden standard plot for this ex1 using shape package:
library(igraph)
library(shape)


# Load the golden standard network for Example 1
numEXAMPLES <- 1  # Example index
mat <- readMat('In silico single cell data/20_nets_10genes_8UNEVENtime_sigma01B_no_initial_points2.mat')
true_network <- mat$networks[,,numEXAMPLES]

# Identify activation and repression edges
gold_activation_edges <- which(true_network > 0, arr.ind = TRUE)
gold_repression_edges <- which(true_network < 0, arr.ind = TRUE)

# Handle self-loops for activation and repression
self_activation <- gold_activation_edges[gold_activation_edges[, 1] == gold_activation_edges[, 2], , drop = FALSE]
self_repression <- gold_repression_edges[gold_repression_edges[, 1] == gold_repression_edges[, 2], , drop = FALSE]

# Create edge lists for activation and repression
activation_edges_list <- cbind(gold_activation_edges[, 1], gold_activation_edges[, 2])
repression_edges_list <- cbind(gold_repression_edges[, 1], gold_repression_edges[, 2])

# Combine all edges into a single list
all_edges_list <- rbind(activation_edges_list, repression_edges_list)

# Create the graph object
g <- graph_from_edgelist(all_edges_list, directed = TRUE)

# Use a circle layout for visualization
layout_coords <- layout_in_circle(g)

# Plot the graph
plot(
  g,
  vertex.label = paste0("Gene ", 1:vcount(g)),  # Label nodes as Gene 1, Gene 2, etc.
  edge.arrow.size = 0.5,                         # Suppress default arrows
  edge.color = NA,                             # Hide default edges
  vertex.color = "skyblue",                    # Node color
  vertex.size = 20,                            # Node size
  main = paste("Gene Regulatory Network - Example", numEXAMPLES),
  layout = layout_coords
)

# Add activation edges (green arrows)
if (nrow(activation_edges_list) > 0) {
  for (i in 1:nrow(activation_edges_list)) {
    from <- activation_edges_list[i, 1]
    to <- activation_edges_list[i, 2]
    Arrows(
      layout_coords[from, 1], layout_coords[from, 2],  
      layout_coords[to, 1], layout_coords[to, 2],     
      arr.type = "triangle", col = "green",           
      arr.width = 0.1, arr.length = 0.2             
    )
  }
}

# Add repression edges 
if (nrow(repression_edges_list) > 0) {
  for (i in 1:nrow(repression_edges_list)) {
    from <- repression_edges_list[i, 1]
    to <- repression_edges_list[i, 2]
    x1 <- layout_coords[from, 1]
    y1 <- layout_coords[from, 2]
    x2 <- layout_coords[to, 1]
    y2 <- layout_coords[to, 2]
    segments(x1, y1, x2, y2, col = "red", lwd = 1.5, lty = 2)  # Dashed red line
    filledcircle(
      mid = c(x2, y2), 
      r1 = 0.03,         # Circle radius
      col = "red"        # Circle color
    )
  }
}

# Add self-activation edges 
if (nrow(self_activation) > 0) {
  for (i in 1:nrow(self_activation)) {
    node <- self_activation[i, 1]
    Arrows(
      layout_coords[node, 1], layout_coords[node, 2], 
      layout_coords[node, 1] + 0.1, layout_coords[node, 2] + 0.1,  
      arr.type = "triangle", col = "green", arr.width = 0.1, arr.length = 0.2
    )
  }
}

# Add self-repression edges 
if (nrow(self_repression) > 0) {
  for (i in 1:nrow(self_repression)) {
    node <- self_repression[i, 1]
    filledcircle(
      mid = c(layout_coords[node, 1] + 0.1, layout_coords[node, 2] + 0.1), 
      r1 = 0.03,
      col = "red"
    )
  }
}
#########using igraph only to plot the golden standard(This is what I used):
library(igraph)

# Adjust edge coordinates to position arrows outside nodes
adjust_arrow_position <- function(layout, from, to, node_radius) {
  dx <- layout[to, 1] - layout[from, 1]
  dy <- layout[to, 2] - layout[from, 2]
  distance <- sqrt(dx^2 + dy^2)
  adjust_x <- node_radius * dx / distance
  adjust_y <- node_radius * dy / distance
  list(
    from = c(layout[from, 1] + adjust_x, layout[from, 2] + adjust_y),
    to = c(layout[to, 1] - adjust_x, layout[to, 2] - adjust_y)
  )
}

# Create the graph object
g <- graph_from_edgelist(rbind(activation_edges_list, repression_edges_list), directed = TRUE)
V(g)$name <- paste0("G", 1:vcount(g))  # Label nodes as G1, G2, etc.

# Get layout for the graph
layout_coords <- layout_in_circle(g)  # Circular layout for clarity

# Plot the graph without edges first
plot(
  g,
  layout = layout_coords,
  vertex.size = 20,
  vertex.color = "skyblue",
  vertex.label = V(g)$name,
  vertex.label.cex = 0.8,
  edge.arrow.size = 0,  # Suppress default arrows
  edge.color = NA,      # Hide default edges
  main = "Gene Regulatory Network - Example 1"
)

# Node radius for adjusting arrows
node_radius <- 0.15  # Adjust based on node size and layout scaling
(green solid arrows)
for (i in 1:nrow(activation_edges_list)) {
  from <- activation_edges_list[i, 1]
  to <- activation_edges_list[i, 2]
  adjusted <- adjust_arrow_position(layout_coords, from, to, node_radius)
  arrows(
    x0 = adjusted$from[1], y0 = adjusted$from[2],
    x1 = adjusted$to[1], y1 = adjusted$to[2],
    col = "green", lwd = 1.5, length = 0.1  # Reduced arrow size
  )
}

# Add repression edges
for (i in 1:nrow(repression_edges_list)) {
  from <- repression_edges_list[i, 1]
  to <- repression_edges_list[i, 2]
  adjusted <- adjust_arrow_position(layout_coords, from, to, node_radius)
  arrows(
    x0 = adjusted$from[1], y0 = adjusted$from[2],
    x1 = adjusted$to[1], y1 = adjusted$to[2],
    col = "red", lty = 2, lwd = 1.5, length = 0.1  # Reduced arrow size and dashed lines
  )
}


########2. time 1-4 use ridge and KS distance, plot the network:
library(kSamples)
library(glmnet)
library(ppcor)
library(R.matlab)
library(igraph)

# Data loading
set.seed(12345)
SIGN <- 1
mat <- readMat('In silico single cell data/20_nets_10genes_8UNEVENtime_sigma01B_no_initial_points2.mat')
time <- as.vector(mat$time.points)
numGENES <- as.vector(mat$n)

# Load SINCERITIES_PLUS and supporting functions
SINCERITIES_PLUS <- dget("SINCERITIES functions/SINCERITIES_PLUS.R")
final_ranked_predictions <- dget("SINCERITIES functions/final_ranked_predictions.R")

# Process only the first 4 time points
data_time_series <- mat$data.tot.array[[1]][[1]]
singleCELLdata <- list()
for (i in 1:4) {
  singleCELLdata[[i]] <- data_time_series[, i, ]
}

totDATA <- matrix(nrow = 0, ncol = dim(data_time_series)[3])
for (i in 1:4) {
  totDATA <- rbind(totDATA, data_time_series[, i, ])
}

DATA <- list(
  time = time[1:4],        # Use only time points 1-4
  numGENES = numGENES,
  singleCELLdata = singleCELLdata,
  genes = sprintf("G%d", 1:numGENES),  # Node names: G1, G2, ...
  totDATA = totDATA
)

# Run SINCERITIES_PLUS
result <- SINCERITIES_PLUS(DATA, noDIAG = 1, SIGN = 1, CV_nfolds = 5)
adj_matrix <- result$adj_matrix / max(result$adj_matrix)

# Generate ranked predictions
table <- final_ranked_predictions(adj_matrix, DATA$genes, SIGN = SIGN, 
                                  fileNAME = 'results_time_1_4', saveFile = FALSE)

# # Filter out "no regulation" edges and apply the threshold
# threshold <- 0.85
# filtered_data <- table[table$Edges != "no regulation" & table$Interaction >= threshold, ]
# 
# # Ensure all nodes are included
# all_nodes <- sprintf("G%d", 1:numGENES)
# graph <- graph_from_data_frame(filtered_data, directed = TRUE, vertices = data.frame(name = all_nodes))
# 
# # Assign colors and line types to edges based on interaction type
# E(graph)$color <- ifelse(filtered_data$Edges == "activation", "green", "red")
# E(graph)$lty <- ifelse(E(graph)$color == "red", 2, 1)  # Red edges dashed, green edges solid
# 
# # Generate layout
# layout <- layout_in_circle(graph)
# 
# # Plot the graph
# plot(
#   graph,
#   layout = layout,
#   vertex.size = 20,               # Size of nodes
#   vertex.color = "skyblue",       # Color of nodes
#   vertex.label = V(graph)$name,   # Label nodes with gene names (G1, G2, ...)
#   vertex.label.cex = 0.8,         # Label size
#   edge.arrow.size = 0.3,          # Arrow size
#   edge.color = E(graph)$color,    # Edge colors
#   edge.lty = E(graph)$lty,        # Edge line type
#   main = "Gene Regulatory Network (Time 1-4, Threshold 0.85)"
# )
# Calculate the top 12% threshold for positive interactions
top_12_percent_threshold <- quantile(table$Interaction, 0.88)

# Filter interactions based on the threshold
filtered_data <- table[table$Edges != "no regulation" & table$Interaction >= top_12_percent_threshold, ]

# Ensure all nodes are included
all_nodes <- sprintf("G%d", 1:numGENES)
graph <- graph_from_data_frame(filtered_data, directed = TRUE, vertices = data.frame(name = all_nodes))

# Assign colors and line types based on interaction type
E(graph)$color <- ifelse(filtered_data$Edges == "activation", "green", "red")
E(graph)$lty <- ifelse(E(graph)$color == "red", 2, 1)  # Red edges dashed, green edges solid

# Generate layout
layout <- layout_in_circle(graph)

# Plot the graph
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
  main = "Gene Regulatory Network (Time 1-4, Top 12%)"
)

# Print the threshold used for verification
cat("Top 12% threshold (Interaction values):", top_12_percent_threshold, "\n")

#calculate the AUROC and AUPR:# Load the true network (gold standard)
true_network <- mat$networks[,,1]  # Example 1 from the dataset
diag(true_network) <- 0  # Remove self-loops

# Prepare predicted adjacency matrix from `adj_matrix`
predicted_1_4 <- adj_matrix  # Adjacency matrix from SINCERITIES_PLUS
diag(predicted_1_4) <- 0  # Ensure no self-loops

# Function to compute AUROC and AUPR
library(pROC)
library(PRROC)

compute_metrics <- function(predicted_interactions, true_adj_matrix) {
  # Flatten matrices
  predicted_scores <- as.vector(predicted_interactions)
  true_labels <- as.vector(true_adj_matrix)
  
  # Ensure binary labels: 1 for edges, 0 for no edges
  true_labels <- ifelse(true_labels != 0, 1, 0)
  
  # Remove NA or irrelevant values
  valid_indices <- !is.na(predicted_scores) & !is.na(true_labels)
  predicted_scores <- predicted_scores[valid_indices]
  true_labels <- true_labels[valid_indices]
  
  # Compute AUROC
  roc_obj <- roc(true_labels, predicted_scores, levels = c(0, 1))
  auroc <- auc(roc_obj)
  
  # Compute AUPR
  pr_curve <- pr.curve(scores.class0 = predicted_scores, weights.class0 = true_labels, curve = TRUE)
  aupr <- pr_curve$auc.integral
  
  return(list(AUROC = auroc, AUPR = aupr))
}

# Compute AUROC and AUPR for Time 1-4
metrics_1_4 <- compute_metrics(predicted_1_4, true_network)
cat("Time 1-4: AUROC =", metrics_1_4$AUROC, "AUPR =", metrics_1_4$AUPR, "\n")


######## time 3-6:
library(kSamples)
library(glmnet)
library(ppcor)
library(R.matlab)
library(igraph)

# Data loading
set.seed(12345)
SIGN <- 1
mat <- readMat('In silico single cell data/20_nets_10genes_8UNEVENtime_sigma01B_no_initial_points2.mat')
time <- as.vector(mat$time.points)
numGENES <- as.vector(mat$n)

# Load SINCERITIES_PLUS and supporting functions
SINCERITIES_PLUS <- dget("SINCERITIES functions/SINCERITIES_PLUS.R")
final_ranked_predictions <- dget("SINCERITIES functions/final_ranked_predictions.R")

# Process time points 3-6
data_time_series <- mat$data.tot.array[[1]][[1]]
singleCELLdata <- list()
for (i in 3:6) {
  singleCELLdata[[i - 2]] <- data_time_series[, i, ]
}

totDATA <- matrix(nrow = 0, ncol = dim(data_time_series)[3])
for (i in 3:6) {
  totDATA <- rbind(totDATA, data_time_series[, i, ])
}

DATA <- list(
  time = time[3:6],        # Use only time points 3-6
  numGENES = numGENES,
  singleCELLdata = singleCELLdata,
  genes = sprintf("G%d", 1:numGENES),  # Node names: G1, G2, ...
  totDATA = totDATA
)

# Run SINCERITIES_PLUS
result <- SINCERITIES_PLUS(DATA, noDIAG = 1, SIGN = 1, CV_nfolds = 5)
adj_matrix <- result$adj_matrix / max(result$adj_matrix)

# Generate ranked predictions
table <- final_ranked_predictions(adj_matrix, DATA$genes, SIGN = SIGN, 
                                  fileNAME = 'results_time_3_6', saveFile = FALSE)

# # Filter out "no regulation" edges and apply the threshold
# threshold <- 0.85
# filtered_data <- table[table$Edges != "no regulation" & table$Interaction >= threshold, ]
# 
# # Ensure all nodes are included
# all_nodes <- sprintf("G%d", 1:numGENES)
# graph <- graph_from_data_frame(filtered_data, directed = TRUE, vertices = data.frame(name = all_nodes))
# 
# # Assign colors and line types to edges based on interaction type
# E(graph)$color <- ifelse(filtered_data$Edges == "activation", "green", "red")
# E(graph)$lty <- ifelse(E(graph)$color == "red", 2, 1)  # Red edges dashed, green edges solid
# 
# # Generate layout
# layout <- layout_in_circle(graph)
# 
# # Plot the graph
# plot(
#   graph,
#   layout = layout,
#   vertex.size = 20,               # Size of nodes
#   vertex.color = "skyblue",       # Color of nodes
#   vertex.label = V(graph)$name,   # Label nodes with gene names (G1, G2, ...)
#   vertex.label.cex = 0.8,         # Label size
#   edge.arrow.size = 0.3,          # Arrow size
#   edge.color = E(graph)$color,    # Edge colors
#   edge.lty = E(graph)$lty,        # Edge line type
#   main = "Gene Regulatory Network (Time 3-6, Threshold 0.85)"
# )
#we can use 12% threshold here:
# Calculate the top 12% threshold for positive interactions
top_12_percent_threshold <- quantile(table$Interaction, 0.88)

# Filter interactions based on the threshold
filtered_data <- table[table$Interaction >= top_12_percent_threshold, ]

# Ensure all nodes are included
all_nodes <- sprintf("G%d", 1:numGENES)
graph <- graph_from_data_frame(filtered_data, directed = TRUE, vertices = data.frame(name = all_nodes))

# Assign colors and line types based on the 'Edges' column
E(graph)$color <- ifelse(filtered_data$Edges == "activation", "green", "red")
E(graph)$lty <- ifelse(E(graph)$color == "red", 2, 1)  # Red edges dashed, green edges solid

# Generate layout
layout <- layout_in_circle(graph)

# Plot the graph
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
  main = "Gene Regulatory Network (Time 3-6, Top 12%)"
)
#
# Map the filtered data to a predicted adjacency matrix
predicted_3_6 <- matrix(0, nrow = numGENES, ncol = numGENES)
for (i in 1:nrow(filtered_data)) {
  source <- as.numeric(sub("G", "", filtered_data$Source[i]))
  target <- as.numeric(sub("G", "", filtered_data$Target[i]))
  predicted_3_6[source, target] <- filtered_data$Interaction[i]
}

# Compute AUROC and AUPR for Time 3-6
metrics_3_6 <- compute_metrics(predicted_3_6, true_network)
cat("Time 3-6: AUROC =", metrics_3_6$AUROC, "AUPR =", metrics_3_6$AUPR, "\n")

#################### time 5-8：
# Process only the time points 5-8
data_time_series <- mat$data.tot.array[[1]][[1]]
singleCELLdata <- list()
for (i in 5:8) {
  singleCELLdata[[i - 4]] <- data_time_series[, i, ]
}

totDATA <- matrix(nrow = 0, ncol = dim(data_time_series)[3])
for (i in 5:8) {
  totDATA <- rbind(totDATA, data_time_series[, i, ])
}

DATA <- list(
  time = time[5:8],        # Use only time points 5-8
  numGENES = numGENES,
  singleCELLdata = singleCELLdata,
  genes = sprintf("G%d", 1:numGENES),  # Node names: G1, G2, ...
  totDATA = totDATA
)

# Run SINCERITIES_PLUS
result <- SINCERITIES_PLUS(DATA, noDIAG = 1, SIGN = 1, CV_nfolds = 5)
adj_matrix <- result$adj_matrix / max(result$adj_matrix)

# Generate ranked predictions
table <- final_ranked_predictions(adj_matrix, DATA$genes, SIGN = SIGN, 
                                  fileNAME = 'results_time_5_8', saveFile = FALSE)

# Calculate the top 12% threshold for positive interactions
top_12_percent_threshold <- quantile(table$Interaction, 0.88)

# Filter interactions based on the threshold
filtered_data <- table[table$Edges != "no regulation" & table$Interaction >= top_12_percent_threshold, ]

# Ensure all nodes are included
all_nodes <- sprintf("G%d", 1:numGENES)
graph <- graph_from_data_frame(filtered_data, directed = TRUE, vertices = data.frame(name = all_nodes))

# Assign colors and line types based on interaction type
E(graph)$color <- ifelse(filtered_data$Edges == "activation", "green", "red")
E(graph)$lty <- ifelse(E(graph)$color == "red", 2, 1)  # Red edges dashed, green edges solid

# Generate layout
layout <- layout_in_circle(graph)

# Plot the graph
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
  main = "Gene Regulatory Network (Time 5-8, Top 12%)"
)

# Print the threshold used for verification
cat("Top 12% threshold (Interaction values):", top_12_percent_threshold, "\n")

# Map the filtered data to a predicted adjacency matrix
predicted_5_8 <- matrix(0, nrow = numGENES, ncol = numGENES)
for (i in 1:nrow(filtered_data)) {
  source <- as.numeric(sub("G", "", filtered_data$Source[i]))
  target <- as.numeric(sub("G", "", filtered_data$Target[i]))
  predicted_5_8[source, target] <- filtered_data$Interaction[i]
}

# Compute AUROC and AUPR for Time 5-8
metrics_5_8 <- compute_metrics(predicted_5_8, true_network)
cat("Time 5-8: AUROC =", metrics_5_8$AUROC, "AUPR =", metrics_5_8$AUPR, "\n")




