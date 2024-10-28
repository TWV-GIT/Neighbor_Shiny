library(shiny)
library(MASS)
library(mvtnorm)
library(ggplot2)
library(dplyr)
library(tidyr)
library(data.table)

# UI ------------------------
ui <- fluidPage(
  tags$head(
    tags$style(HTML("
      .loading-spinner {
        position: fixed;
        top: 50%;
        left: 50%;
        transform: translate(-50%, -50%);
        background: rgba(255, 255, 255, 0.9);
        padding: 20px;
        border-radius: 5px;
        box-shadow: 0 0 10px rgba(0,0,0,0.1);
        z-index: 1000;
        text-align: center;
      }
      .loading-spinner .spinner-border {
        width: 3rem;
        height: 3rem;
      }
    "))
  ),
  
  titlePanel("Neighbor Analysis Simulation"),
  
  sidebarLayout(
    sidebarPanel(
      numericInput("n_centers_dist", "Number of recombination clusters:", 
                   value = 19, min = 1, max = 20),
      numericInput("n_convergent_centers", "Number of convergent clusters:", 
                   value = 5, min = 1, max = 20),
      numericInput("n_samples_sampled", "Number of sampled CDR3s:", 
                   value = 1000, min = 100, max = 1000),
      numericInput("sample_ratio", "Ratio of sampled to background CDR3s:", 
                   value = 10, min = 1, max = 20),
      numericInput("radius", "Radius for neighborhood analysis:", 
                   value = 0.2, min = 0.1, max = 2, step = 0.1),
      sliderInput("intensity_range", "Relative intensity range:",
                  min = 0, max = 1, value = c(0.3, 0.5)),
      numericInput("seed", "Random seed:", value = 3),
      actionButton("run", "Run Analysis", class = "btn-primary",
                   style = "width: 100%")
    ),
    mainPanel(
      conditionalPanel(
        condition = "input.run != 0 && $('html').hasClass('shiny-busy')",
        tags$div(class = "loading-spinner",
                 tags$div(class = "spinner-border text-primary", role = "status"),
                 tags$br(),
                 tags$span("Calculating... Please wait"))
      ),
      
      # First row: Distribution, Neighborhood, and Enrichment Map
      fluidRow(
        column(6, h3("Recombination pGen Plot"), plotOutput("distribution_plot")),
        column(6, h3("Logratio Plot"), plotOutput("neighborhood_plot"))
      ),
      
      # Second row: Analysis Plot and p-value plot
      fluidRow(
        column(6, h3("Neighbour enrichment Plot"), plotOutput("analysis_plot")),
        column(6, h3("p-value Plot"), plotOutput("enrichment_plot"))
      ),
      
      # Third row: summary statistics (remove?)
      fluidRow(
        column(12, h3("Statistics"), verbatimTextOutput("stats"))
      )
    )
  )
)

# Calculations ------------------------

server <- function(input, output) {
  results <- reactiveVal(NULL)
  
  # Center generation --> add option to define centers?
  generate_centers <- function(n, min_x = -2, max_x = 2, min_y = -2, max_y = 2) {
    means <- matrix(runif(n * 2, c(min_x, min_y), c(max_x, max_y)), ncol = 2)
    scales <- matrix(runif(n * 2, 0.1, 0.5), ncol = 2)
    thetas <- runif(n, 0, 2*pi)
    
    centers <- vector("list", n)
    for(i in 1:n) {
      R <- matrix(c(cos(thetas[i]), -sin(thetas[i]), 
                    sin(thetas[i]), cos(thetas[i])), 2, 2)
      S <- diag(scales[i,]^2)
      centers[[i]] <- list(
        mean = means[i,],
        sigma = R %*% S %*% t(R)
      )
    }
    return(centers)
  }
  
  calculate_neighborhoods <- function(samples_all, samples_dist, radius) {
    dt_all <- as.data.table(samples_all)
    dt_dist <- as.data.table(samples_dist)
    setnames(dt_all, c("x", "y"))
    setnames(dt_dist, c("x", "y"))
    
    total_sampled <- nrow(dt_all)
    total_background <- nrow(dt_dist)
    
    calc_neighbors <- function(point, dt, radius) {
      sqrt((dt$x - point[1])^2 + (dt$y - point[2])^2) <= radius
    }
    
    results <- dt_all[, {
      neighbors_sampled <- sum(calc_neighbors(c(x, y), dt_all, radius)) - 1 # -1 to remove self
      neighbors_background <- sum(calc_neighbors(c(x, y), dt_dist, radius))
    
      ratio <- neighbors_sampled/((neighbors_background +1)/10) # Based on formula of Fabio A.
      
      # Calculate hypergeometric test p-value (models sampling without replacement)
      pvalue <- tryCatch({
        if (neighbors_sampled > 0) {  # Only calculate p-value if we found neighbors in set 2
          phyper(
            q = neighbors_sampled - 1,           # Subtract 1 to get strictly greater than
            m = total_sampled - 1,              # Total sampled CDR3s minus the query point (total_sampled - 1)
            n = total_background,                  # Total background CDR3s (total_background)
            k = neighbors_background + neighbors_sampled, # total number of neighbor CDR3s found (neighbors_background + neighbors_sampled - 1)
            lower.tail = FALSE               # We want P(X > x)
          )
        } else {
          1.0  # No enrichment if no sampled CDR3s neighbors
        }
      }, error = function(e) 1.0)
      
      # Ensure p-value is valid but keep original ratio
      pvalue <- pmin(pmax(pvalue, .Machine$double.xmin), 1)
      
      list(
        points_in_radius_sampled = neighbors_sampled,
        points_in_radius_background = neighbors_background,
        ratio = ratio,
        pvalue = pvalue,
        logratio = log10(ratio),
        neg_log10_pvalue = -log10(pvalue)
      )
    }, by = .(x, y)]
    
    return(results)
  }

    observeEvent(input$run, {
    withProgress(message = 'Performing calculations', value = 0, {
      set.seed(input$seed)
      
      # Generate centers
      incProgress(0.1, detail = "Generating centers...")
      centers_dist <- generate_centers(input$n_centers_dist)
      centers_convergent <- generate_centers(input$n_convergent_centers)
      all_centers <- c(centers_dist, centers_convergent)
      
      # Calculate weights using vectorization
      incProgress(0.2, detail = "Calculating weights...")
      weights_dist <- runif(input$n_centers_dist, 
                            input$intensity_range[1], 
                            input$intensity_range[2])
      weights_dist <- weights_dist / sum(weights_dist)
      
      weights_all <- runif(input$n_centers_dist + input$n_convergent_centers, 
                           input$intensity_range[1], 
                           input$intensity_range[2])
      weights_all <- weights_all / sum(weights_all)
      
      # Generate distribution using vectorization
      incProgress(0.3, detail = "Generating distribution...")
      x <- seq(-3, 3, length.out = 100)
      y <- seq(-3, 3, length.out = 100)
      grid <- as.matrix(expand.grid(x = x, y = y))
      
      # Vectorized distribution calculation
      z <- matrix(0, nrow = length(x), ncol = length(y))
      for (i in 1:input$n_centers_dist) {
        z_component <- matrix(
          weights_dist[i] * dmvnorm(grid, 
                                    mean = centers_dist[[i]]$mean, 
                                    sigma = centers_dist[[i]]$sigma),
          nrow = length(x)
        )
        z <- z + z_component
      }
      
      # Generate samples using vectorization
      incProgress(0.5, detail = "Generating samples...")
      n_samples_background <- input$n_samples_sampled * input$sample_ratio
      
      # Vectorized sample generation
      components_dist <- sample(1:input$n_centers_dist, n_samples_background, 
                                replace = TRUE, prob = weights_dist)
      samples_dist <- t(sapply(components_dist, function(i) {
        mvrnorm(1, mu = centers_dist[[i]]$mean, Sigma = centers_dist[[i]]$sigma)
      }))
      
      components_all <- sample(1:(input$n_centers_dist + input$n_convergent_centers), 
                               input$n_samples_sampled, replace = TRUE, prob = weights_all)
      samples_all <- t(sapply(components_all, function(i) {
        mvrnorm(1, mu = all_centers[[i]]$mean, Sigma = all_centers[[i]]$sigma)
      }))
      
      # Calculate neighborhood statistics using optimized function
      incProgress(0.7, detail = "Calculating neighborhood statistics...")
      neighborhood_stats <- calculate_neighborhoods(samples_all, samples_dist, input$radius)
      
      # Prepare centers data
      centers_df <- data.frame(
        x = sapply(all_centers, function(c) c$mean[1]),
        y = sapply(all_centers, function(c) c$mean[2]),
        type = c(rep("recombined", input$n_centers_dist), 
                 rep("convergent", input$n_convergent_centers)),
        weight = weights_all
      )
      
      # Store results
      results(list(
        plot_data = data.table(expand.grid(x = x, y = y), z = as.vector(z)),
        centers_df = centers_df,
        samples_dist = samples_dist,
        neighborhood_stats = neighborhood_stats
      ))
      
      incProgress(1, detail = "Done!")
    })
  })

# Plot ------------------------
  output$distribution_plot <- renderPlot({
    req(results())
    ggplot(results()$plot_data, aes(x = x, y = y)) +
      geom_raster(aes(fill = z)) +
      scale_fill_viridis_c() +
      geom_point(data = results()$centers_df[results()$centers_df$type == "recombined",],
                 aes(size = weight), shape = 4, color = "red") +
      scale_size_continuous(range = c(3, 8)) +
      labs(x = "Dim 1", y = "Dim 2") +
      theme_minimal() +
      coord_fixed(ratio = 1) +
      theme(legend.position = "none")
  })
  
  output$neighborhood_plot <- renderPlot({
    req(results())
    ggplot() +
      geom_point(data = as.data.frame(results()$samples_dist), 
                 aes(x = V1, y = V2), 
                 color = "blue", alpha = 0.2) +
      geom_point(data = results()$neighborhood_stats,
                 aes(x = x, y = y, color = ratio), 
                 size = 3) +
      scale_color_gradient(low = "yellow", high = "red",
                           name = "logratio") +
      geom_point(data = results()$centers_df,
                 aes(x = x, y = y, shape = type),
                 color = "red", size = 3) +
      scale_shape_manual(values = c("convergent" = 3, "recombined" = 4)) +
      labs(x = "Dim 1", y = "Dim 2") +
      theme_minimal() +
      coord_fixed(ratio = 1)
  })
  
  output$enrichment_plot <- renderPlot({
    req(results())
    ggplot() +
      geom_point(data = as.data.frame(results()$samples_dist), 
                 aes(x = V1, y = V2), 
                 color = "gray80", 
                 alpha = 0.2) +
      geom_point(data = results()$neighborhood_stats,
                 aes(x = x, y = y, color = neg_log10_pvalue), 
                 size = 3) +
      scale_color_gradientn(
        colors = c("black", "orange", "red", "purple", "blue"),
        name = "-log10(p-value)"
      ) +
      geom_point(data = dplyr::filter(results()$centers_df, type == 'convergent'),
                 aes(x = x, y = y),
                 color = "red", 
                 size = 4,
                 shape = 3,
                 stroke = 3) +
      labs(x = "Dim 1", 
           y = "Dim 2") +
      theme_minimal() +
      coord_fixed(ratio = 1) +
      theme(
        legend.position = "right",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "gray90")
      )
  })
  
  output$analysis_plot <- renderPlot({
    req(results())
    ggplot(results()$neighborhood_stats, 
           aes(x = logratio, 
               y = points_in_radius_sampled)) +
      geom_vline(xintercept = c(-1, 0, 1, 2), 
                 linetype = "dashed", 
                 color = "gray80") +
      geom_point(aes(color = neg_log10_pvalue), 
                 size = 2) +
      expand_limits(x = 2.5)+
      scale_color_gradientn(colors = c("black", "orange", "red", "purple", "blue")) +
      labs(x = 'logratio', 
           y = 'neighbors',
           subtitle = paste0("radius = ", input$radius,
                             ", number of clusters = ", input$n_centers_dist),
           color = '-log10(p-value)') +
      theme_bw() +
      theme(
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "gray90"),
        plot.title = element_text(size = 12),
        axis.title = element_text(size = 10),
        legend.position = "right"
      ) +
      geom_vline(xintercept = 0, color = "black", linetype = "dashed")
  })
  
  output$stats <- renderPrint({
    req(results())
    stats <- results()$neighborhood_stats[, .(
      mean_points_sampled = mean(points_in_radius_sampled),
      mean_points_background = mean(points_in_radius_background),
      mean_ratio = mean(ratio),
      median_ratio = median(ratio)
    )]
    
    cat("\nNeighborhood Analysis Summary:\n")
    cat("Mean points in radius (Sampled CDR3s):", stats$mean_points_sampled, "\n")
    cat("Mean points in radius (Background CDR3s):", stats$mean_points_background, "\n")
    cat("Mean ratio:", stats$mean_ratio, "\n")
    cat("Median ratio:", stats$median_ratio, "\n")
  })
}

shinyApp(ui = ui, server = server)  
