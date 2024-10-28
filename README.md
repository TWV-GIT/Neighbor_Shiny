# Neighbor Analysis Modeling - Shiny App

This repository contains a Shiny app written in R, designed to model neighbor analysis with customizable parameters. 

## Overview

The app generates random background Complementarity-Determining Region 3 (CDR3) clusters and 'convergent' CDR3 clusters, then performs neighbor analysis to produce common plots and statistical outputs. This allows for the visualization of how cluster properties impact the resulting p-values and plot characteristics.

The app is hosted at https://twverdonckt.shinyapps.io/neighbor_shiny/

## Key Features

- **Random Background and Convergent Cluster Generation**: The app creates random CDR3 clusters alongside 'convergent' clusters to model various neighborhood structures.
- **Neighbor Analysis**: Performs neighbor analysis to evaluate spatial relationships within generated clusters.
- **Plot Generation**: Produces commonly used plots to help interpret the effects of:
  - Cluster size and count
  - Topology
  - Other customizable parameters

## Important Notes

- **No TCRdist3 Metrics**: The app does not employ TCRdist3 or related distance metrics.
- **No Real CDR3 Sequences**: It operates purely on generated clusters and does not simulate actual CDR3 sequences.

## Getting Started

1. Clone this repository:
    ```bash
    git clone https://github.com/TWV-GIT/Neighbor_Shiny.git
    ```
2. Open the Shiny app in R or RStudio and run `shiny::runApp("path/to/your/app")` to start the application.

## Dependencies

Make sure the following R packages are installed:
- **Shiny**
- **dplyr**
- **ggplot2**
- **MASS**
- **mvtnorm**
- **tidyr**
- **data.table**

Install dependencies with:
```r
install.packages(c("shiny", "dplyr", "ggplot2", "MASS", "mvtnorm", "tidyr", "data.table"))
