#' Network Visualization Functions
#'
#' Comprehensive visualization tools for protein-protein interaction networks

#' Visualize PPI Network
#'
#' Create various visualizations of protein-protein interaction networks
#'
#' @param network An igraph object representing the PPI network
#' @param input_proteins Vector of protein names to extract sub-network. If NULL, visualizes the whole network.
#' @param layout_type Layout algorithm for network visualization ("fr", "kk", "circle", "spring", "tree")
#' @param node_color Color scheme for nodes ("degree", "community", "source", "uniform")
#' @param edge_color Color scheme for edges ("weight", "source", "uniform")
#' @param node_size Size of nodes ("degree", "uniform", "betweenness")
#' @param show_labels Whether to display node labels
#' @param highlight_nodes Vector of node names to highlight with red color and star shape
#' @param title Plot title
#' @param save_path Optional path to save the plot
#' @param width Plot width in inches
#' @param height Plot height in inches
#' @param interactive Create interactive plot using visNetwork
#'
#' @return A ggplot object or interactive visNetwork object
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Basic network visualization
#' networks <- build_ppi_network(data_sources = c("STRING", "BIOGRID"))
#' merged <- merge_ppi_networks(networks, merge_method = "union")
#' viz <- visualize_ppi_network(merged)
#'
#' # Sub-network visualization
#' sub_viz <- visualize_ppi_network(merged, input_proteins = c("TP53", "BRCA1"))
#'
#' # Highlight specific proteins
#' highlight_viz <- visualize_ppi_network(merged, highlight_nodes = c("TP53", "BRCA1"))
#'
#' # Interactive visualization
#' interactive_viz <- visualize_ppi_network(merged, interactive = TRUE)
#'
#' # Save to file
#' visualize_ppi_network(merged, save_path = "ppi_network.png", width = 12, height = 8)
#' }
visualize_ppi_network <- function(network,
                                  input_proteins = NULL,
                                  layout_type = "fr",
                                  node_color = "degree",
                                  edge_color = "weight",
                                  node_size = "degree",
                                  show_labels = TRUE,
                                  highlight_nodes = NULL,
                                  title = "PPI Network",
                                  save_path = NULL,
                                  width = 10,
                                  height = 8,
                                  interactive = FALSE) {
  if (!igraph::is.igraph(network)) {
    stop("Input must be an igraph object")
  }

  if (igraph::vcount(network) == 0) {
    warning("Empty network provided")
    return(NULL)
  }

  # Handle sub-network extraction based on input_proteins
  if (!is.null(input_proteins)) {
    input_proteins <- toupper(input_proteins)

    # Find nodes that exist in the network
    network_nodes <- toupper(igraph::V(network)$name)
    valid_proteins <- input_proteins[input_proteins %in% network_nodes]

    if (length(valid_proteins) == 0) {
      warning("None of the input proteins found in the network")
      return(NULL)
    }

    # Extract sub-network including neighbors
    sub_network <- igraph::induced_subgraph(
      network,
      vids = igraph::V(network)[toupper(igraph::V(network)$name) %in% valid_proteins]
    )

    # Add first-order neighbors if available
    neighbors <- unique(unlist(igraph::neighborhood(sub_network, order = 1)))
    if (length(neighbors) > 0) {
      sub_network <- igraph::induced_subgraph(network, vids = neighbors)
    }

    network <- sub_network

    if (igraph::vcount(network) == 0) {
      warning("No sub-network could be extracted")
      return(NULL)
    }

    # Update title to indicate sub-network
    if (is.null(title) || title == "PPI Network") {
      title <- paste("Sub-network with", length(valid_proteins), "seed proteins")
    }
  }

  if (interactive) {
    # Create interactive visualization using visNetwork
    if (!requireNamespace("visNetwork", quietly = TRUE)) {
      warning("visNetwork package not available, falling back to static plot")
      interactive <- FALSE
    } else {
      return(.create_interactive_network(network, node_color, node_size, title, highlight_nodes))
    }
  }

  # Static visualization
  return(.create_static_network(
    network, layout_type, node_color, edge_color,
    node_size, show_labels, highlight_nodes, title,
    save_path, width, height
  ))
}

#' Create Static Network Visualization
#' @noRd
.create_static_network <- function(network, layout_type, node_color, edge_color,
                                   node_size, show_labels, highlight_nodes, title,
                                   save_path, width, height) {
  # Prepare layout
  layout <- switch(layout_type,
    "fr" = igraph::layout_with_fr(network),
    "kk" = igraph::layout_with_kk(network),
    "circle" = igraph::layout_in_circle(network),
    "spring" = igraph::layout_with_spring(network),
    "tree" = igraph::layout_as_tree(network),
    igraph::layout_with_fr(network)
  )

  # Prepare node data
  igraph::V(network)$degree <- igraph::degree(network)
  igraph::V(network)$betweenness <- igraph::betweenness(network, normalized = TRUE)
  igraph::V(network)$community <- as.factor(igraph::membership(igraph::cluster_fast_greedy(network)))

  # Get edge attributes
  edge_weights <- igraph::E(network)$weight
  if (is.null(edge_weights)) {
    edge_weights <- rep(1, igraph::ecount(network))
  }

  # Prepare node colors
  node_colors <- switch(node_color,
    "degree" = colorRampPalette(c("lightblue", "darkblue"))(length(unique(igraph::V(network)$degree)))[as.numeric(cut(igraph::V(network)$degree, breaks = length(unique(igraph::V(network)$degree))))],
    "community" = RColorBrewer::brewer.pal(max(3, min(8, length(unique(igraph::V(network)$community)))), "Set3")[as.numeric(igraph::V(network)$community)],
    "uniform" = rep("skyblue", igraph::vcount(network)),
    rep("skyblue", igraph::vcount(network))
  )

  # Handle highlighted nodes
  if (!is.null(highlight_nodes)) {
    highlight_nodes <- toupper(highlight_nodes)
    highlight_indices <- igraph::V(network)$name %in% highlight_nodes
    node_colors[highlight_indices] <- "red"
  }

  # Prepare node sizes
  node_sizes <- switch(node_size,
    "degree" = sqrt(igraph::V(network)$degree) * 2 + 3,
    "betweenness" = sqrt(igraph::V(network)$betweenness) * 10 + 3,
    "uniform" = rep(5, igraph::vcount(network)),
    rep(5, igraph::vcount(network))
  )



  edge_widths <- sqrt(edge_weights) * 0.5 + 0.5

  # Create plot
  plot <- igraph::plot.igraph(network,
    layout = layout,
    vertex.color = node_colors,
    vertex.size = node_sizes,
    vertex.label = ifelse(show_labels, igraph::V(network)$name, NA),
    vertex.label.cex = 0.7,
    vertex.label.color = "black",
    edge.width = edge_widths,
    edge.arrow.size = 0.3,
    main = title
  )

  # Add legend if needed
  if (node_color == "degree") {
    legend("topright",
      legend = c("Low", "Medium", "High"),
      col = colorRampPalette(c("lightblue", "darkblue"))(3),
      pch = 19, title = "Node Degree"
    )
  }

  # Save if requested
  if (!is.null(save_path)) {
    if (!requireNamespace("Cairo", quietly = TRUE)) {
      png(save_path, width = width * 300, height = height * 300, res = 300)
    } else {
      Cairo::CairoPNG(save_path, width = width * 300, height = height * 300, res = 300)
    }
    plot
    dev.off()
    message("Network visualization saved to: ", save_path)
  }

  return(plot)
}

#' Create Interactive Network Visualization
#' @noRd
.create_interactive_network <- function(network, node_color, node_size, title, highlight_nodes = NULL) {
  # Prepare node data
  node_df <- data.frame(
    id = igraph::V(network)$name,
    label = igraph::V(network)$name,
    title = igraph::V(network)$name,
    degree = igraph::degree(network),
    betweenness = igraph::betweenness(network, normalized = TRUE),
    community = as.factor(igraph::membership(igraph::cluster_fast_greedy(network)))
  )

  # Prepare edge data
  edge_df <- igraph::as_data_frame(network, what = "edges")

  # Determine node colors
  if (node_color == "degree") {
    node_df$color <- colorRampPalette(c("lightblue", "darkblue"))(length(unique(node_df$degree)))[as.numeric(cut(node_df$degree, breaks = length(unique(node_df$degree))))]
  } else if (node_color == "community") {
    colors <- RColorBrewer::brewer.pal(max(3, min(8, length(unique(node_df$community)))), "Set3")
    node_df$color <- colors[as.numeric(node_df$community)]
  } else {
    node_df$color <- "skyblue"
  }

  # Determine node sizes
  if (node_size == "degree") {
    node_df$size <- sqrt(node_df$degree) * 5 + 10
  } else if (node_size == "betweenness") {
    node_df$size <- sqrt(node_df$betweenness) * 20 + 10
  } else {
    node_df$size <- 15
  }

  # Handle highlighted nodes
  if (!is.null(highlight_nodes)) {
    highlight_nodes <- toupper(highlight_nodes)
    node_df$borderWidth <- ifelse(toupper(node_df$id) %in% highlight_nodes, 3, 1)
    node_df$borderColor <- ifelse(toupper(node_df$id) %in% highlight_nodes, "red", "black")
    node_df$shape <- ifelse(toupper(node_df$id) %in% highlight_nodes, "star", "dot")
  }

  # Create interactive network
  visNetwork::visNetwork(
    nodes = node_df,
    edges = edge_df,
    width = "100%",
    height = "600px"
  ) %>%
    visNetwork::visNodes(
      shape = "dot",
      scaling = list(min = 10, max = 30),
      font = list(size = 12, color = "black")
    ) %>%
    visNetwork::visEdges(
      width = 1,
      smooth = list(type = "continuous")
    ) %>%
    visNetwork::visLayout(randomSeed = 42) %>%
    visNetwork::visIgraphLayout(layout = "layout_with_fr") %>%
    visNetwork::visOptions(
      highlightNearest = list(enabled = TRUE, degree = 1, hover = TRUE),
      nodesIdSelection = TRUE
    ) %>%
    visNetwork::visInteraction(
      tooltipDelay = 200,
      hideEdgesOnDrag = TRUE
    )
}

#' Plot Network Summary
#'
#' Create summary plots for network analysis
#'
#' @param network An igraph object representing the PPI network
#' @param save_path Optional path to save the plots
#' @param width Plot width in inches
#' @param height Plot height in inches
#'
#' @return List of ggplot objects
#'
#' @examples
#' \dontrun{
#' # Create a sample network
#' library(igraph)
#' g <- erdos.renyi.game(50, 0.1)
#' V(g)$name <- paste0("GENE", 1:50)
#'
#' # Generate summary plots
#' plots <- plot_network_summary(g)
#'
#' # Display degree distribution
#' print(plots$degree_distribution)
#' }
#'
#' @export
plot_network_summary <- function(network, save_path = NULL, width = 8, height = 6) {
  if (!igraph::is.igraph(network)) {
    stop("Input must be an igraph object")
  }

  if (igraph::vcount(network) == 0) {
    warning("Empty network provided")
    return(NULL)
  }

  plots <- list()

  # Degree distribution
  degrees <- igraph::degree(network)
  plots$degree_distribution <- ggplot2::ggplot(
    data.frame(degree = degrees),
    ggplot2::aes(x = degree)
  ) +
    ggplot2::geom_histogram(binwidth = 1, fill = "skyblue", color = "black") +
    ggplot2::labs(title = "Node Degree Distribution", x = "Degree", y = "Count") +
    ggplot2::theme_minimal()

  # Betweenness distribution
  betweenness <- igraph::betweenness(network, normalized = TRUE)
  plots$betweenness_distribution <- ggplot2::ggplot(
    data.frame(betweenness = betweenness),
    ggplot2::aes(x = betweenness)
  ) +
    ggplot2::geom_histogram(fill = "lightcoral", color = "black", bins = 20) +
    ggplot2::labs(title = "Node Betweenness Distribution", x = "Betweenness Centrality", y = "Count") +
    ggplot2::theme_minimal()

  # Edge weight distribution (if available)
  edge_weights <- igraph::E(network)$weight
  if (!is.null(edge_weights)) {
    plots$weight_distribution <- ggplot2::ggplot(
      data.frame(weight = edge_weights),
      ggplot2::aes(x = weight)
    ) +
      ggplot2::geom_histogram(fill = "lightgreen", color = "black", bins = 20) +
      ggplot2::labs(title = "Edge Weight Distribution", x = "Weight", y = "Count") +
      ggplot2::theme_minimal()
  }

  # Community size distribution
  communities <- igraph::cluster_fast_greedy(network)
  comm_sizes <- table(igraph::membership(communities))
  plots$community_sizes <- ggplot2::ggplot(
    data.frame(size = as.numeric(comm_sizes)),
    ggplot2::aes(x = size)
  ) +
    ggplot2::geom_bar(fill = "lightsteelblue", color = "black", stat = "count") +
    ggplot2::labs(title = "Community Size Distribution", x = "Community Size", y = "Count") +
    ggplot2::theme_minimal()

  # Save plots if requested
  if (!is.null(save_path)) {
    base_path <- tools::file_path_sans_ext(save_path)
    ext <- tools::file_ext(save_path)

    for (plot_name in names(plots)) {
      filename <- paste0(base_path, "_", plot_name, ".", ext)
      ggplot2::ggsave(filename, plots[[plot_name]], width = width, height = height)
    }
    message("Network summary plots saved with base name: ", base_path)
  }

  return(plots)
}

#' Plot Network Comparison
#'
#' Compare multiple networks side by side
#'
#' @param networks_list List of igraph objects to compare
#' @param network_names Names for the networks
#' @param save_path Optional path to save the comparison plot
#' @param width Plot width in inches
#' @param height Plot height in inches
#'
#' @return ggplot object
#'
#' @examples
#' \dontrun{
#' # Create sample networks
#' library(igraph)
#' net1 <- erdos.renyi.game(50, 0.1)
#' net2 <- erdos.renyi.game(60, 0.15)
#'
#' # Compare networks
#' plots <- plot_network_comparison(
#'   networks_list = list(net1, net2),
#'   network_names = c("Network 1", "Network 2")
#' )
#'
#' # Display comparison
#' print(plots$density)
#' }
#'
#' @export
plot_network_comparison <- function(networks_list, network_names = NULL,
                                    save_path = NULL, width = 12, height = 8) {
  if (length(networks_list) == 0) {
    stop("No networks provided for comparison")
  }

  if (is.null(network_names)) {
    network_names <- paste0("Network", seq_along(networks_list))
  }

  # Calculate network metrics
  metrics <- data.frame(
    Network = network_names,
    Nodes = vapply(networks_list, igraph::vcount, integer(1)),
    Edges = vapply(networks_list, igraph::ecount, integer(1)),
    Density = vapply(networks_list, igraph::edge_density, numeric(1)),
    Avg_Degree = vapply(networks_list, function(x) mean(igraph::degree(x)), numeric(1)),
    Diameter = vapply(networks_list, function(x) igraph::diameter(x, directed = FALSE), numeric(1)),
    Clustering = vapply(networks_list, igraph::transitivity, numeric(1), type = "global")
  )

  # Create comparison plots
  plots <- list()

  # Nodes vs Edges
  plots$nodes_edges <- ggplot2::ggplot(metrics, ggplot2::aes(x = Nodes, y = Edges, label = Network)) +
    ggplot2::geom_point(size = 3) +
    ggplot2::geom_text(vjust = -0.5) +
    ggplot2::labs(title = "Nodes vs Edges", x = "Number of Nodes", y = "Number of Edges") +
    ggplot2::theme_minimal()

  # Density comparison
  plots$density <- ggplot2::ggplot(metrics, ggplot2::aes(x = Network, y = Density)) +
    ggplot2::geom_bar(stat = "identity", fill = "steelblue") +
    ggplot2::labs(title = "Network Density Comparison", y = "Density") +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  # Average degree
  plots$avg_degree <- ggplot2::ggplot(metrics, ggplot2::aes(x = Network, y = Avg_Degree)) +
    ggplot2::geom_bar(stat = "identity", fill = "coral") +
    ggplot2::labs(title = "Average Degree Comparison", y = "Average Degree") +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  # Combined metrics heatmap
  metrics_scaled <- as.data.frame(scale(metrics[, -1]))
  metrics_scaled$Network <- metrics$Network

  metrics_long <- tidyr::pivot_longer(metrics_scaled, -Network, names_to = "Metric", values_to = "Value")

  plots$heatmap <- ggplot2::ggplot(metrics_long, ggplot2::aes(x = Network, y = Metric, fill = Value)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
    ggplot2::labs(title = "Network Metrics Heatmap") +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  # Save if requested
  if (!is.null(save_path)) {
    base_path <- tools::file_path_sans_ext(save_path)
    ext <- tools::file_ext(save_path)

    for (plot_name in names(plots)) {
      filename <- paste0(base_path, "_", plot_name, ".", ext)
      ggplot2::ggsave(filename, plots[[plot_name]], width = width, height = height)
    }
    message("Network comparison plots saved with base name: ", base_path)
  }

  return(plots)
}
