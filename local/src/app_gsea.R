#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(shiny)
  library(DT)
  library(ggplot2)
  library(fgsea)
})

parse_args <- function(args) {
  values <- list(
    host = "127.0.0.1",
    port = 3838L,
    rnk_dir = ".",
    gmt_dir = ".",
    results = "multi_GSEA.gz"
  )
  names_map <- c(
    "--host" = "host", "--port" = "port", "--rnk-dir" = "rnk_dir",
    "--gmt-dir" = "gmt_dir", "--results" = "results"
  )
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (!key %in% names(names_map) || i == length(args)) {
      stop("Usage: app_gsea.R [--host HOST] [--port PORT] [--rnk-dir DIR] ",
           "[--gmt-dir DIR] [--results FILE]")
    }
    values[[names_map[[key]]]] <- args[[i + 1L]]
    i <- i + 2L
  }
  values$port <- suppressWarnings(as.integer(values$port))
  if (is.na(values$port) || values$port < 1L || values$port > 65535L) {
    stop("--port must be an integer between 1 and 65535")
  }
  values
}

read_ranking <- function(path) {
  ranking <- read.delim(path, header = FALSE, stringsAsFactors = FALSE)
  if (ncol(ranking) < 2L) stop("Ranking must have at least two columns: ", path)
  ranking <- ranking[, 1:2]
  names(ranking) <- c("gene", "score")
  ranking$gene <- trimws(ranking$gene)
  ranking$score <- suppressWarnings(as.numeric(ranking$score))
  ranking <- ranking[nzchar(ranking$gene) & is.finite(ranking$score), ]
  if (!nrow(ranking)) stop("No valid gene scores in: ", path)
  if (anyDuplicated(ranking$gene)) {
    stop("Gene identifiers must be unique in: ", path)
  }
  sort(setNames(ranking$score, ranking$gene), decreasing = TRUE)
}

read_gsea_results <- function(path) {
  if (!file.exists(path)) stop("GSEA results not found: ", path)
  has_header <- grepl("header_added", basename(path), fixed = TRUE)
  results <- read.delim(path, header = has_header, stringsAsFactors = FALSE,
                        check.names = FALSE)
  if (!has_header) {
    expected <- c("contrast", "msigdb_type", "pathway", "pval", "padj",
                  "ES", "NES", "nMoreExtreme", "size", "leadingEdge")
    if (ncol(results) != length(expected)) {
      stop("Expected 10 columns in headerless GSEA results, found ", ncol(results))
    }
    names(results) <- expected
  }
  required <- c("contrast", "msigdb_type", "pathway", "pval", "padj",
                "NES", "size")
  missing <- setdiff(required, names(results))
  if (length(missing)) stop("Missing GSEA columns: ", paste(missing, collapse = ", "))
  for (column in c("pval", "padj", "ES", "NES", "size")) {
    if (column %in% names(results)) {
      results[[column]] <- suppressWarnings(as.numeric(results[[column]]))
    }
  }
  results
}

normalise_id <- function(x) sub("\\.rnk$", "", basename(as.character(x)))

options <- parse_args(commandArgs(trailingOnly = TRUE))
rnk_files <- list.files(options$rnk_dir, pattern = "\\.rnk$", full.names = TRUE)
gmt_files <- list.files(options$gmt_dir, pattern = "\\.gmt$", full.names = TRUE)
if (!length(rnk_files)) stop("No .rnk files found in: ", options$rnk_dir)
if (!length(gmt_files)) stop("No .gmt files found in: ", options$gmt_dir)

rankings <- setNames(lapply(rnk_files, read_ranking), normalise_id(rnk_files))
pathways <- list()
for (gmt_file in gmt_files) pathways <- c(pathways, fgsea::gmtPathways(gmt_file))
gsea_table <- read_gsea_results(options$results)
gsea_table$contrast_id <- normalise_id(gsea_table$contrast)
available_rankings <- intersect(names(rankings), unique(gsea_table$contrast_id))
if (!length(available_rankings)) {
  stop("No ranking name matches the contrasts in ", options$results,
       ". Rankings: ", paste(names(rankings), collapse = ", "),
       "; contrasts: ", paste(unique(gsea_table$contrast_id), collapse = ", "))
}
rankings <- rankings[available_rankings]

ui <- fluidPage(
  titlePanel("SnakeGSEA local explorer"),
  sidebarLayout(
    sidebarPanel(
      selectInput("ranking", "Ranking", choices = names(rankings)),
      selectInput("collection", "MSigDB collection",
                  choices = c("ALL", sort(unique(gsea_table$msigdb_type)))),
      sliderInput("padj", "Maximum FDR (padj)", min = 0.001, max = 1,
                  value = 0.05, step = 0.01),
      radioButtons("direction", "Enrichment direction",
                   choices = c("All" = "all", "NES > 0" = "positive",
                               "NES < 0" = "negative")),
      helpText("Select a table row to display its enrichment curve and leading edge.")
    ),
    mainPanel(
      DTOutput("results_table"),
      h4(textOutput("selection_title")),
      plotOutput("enrichment_plot", height = "380px"),
      h4("Leading-edge genes"),
      DTOutput("leading_edge_table")
    )
  )
)

server <- function(input, output, session) {
  filtered_results <- reactive({
    data <- gsea_table[gsea_table$contrast_id == input$ranking, , drop = FALSE]
    if (input$collection != "ALL") {
      data <- data[data$msigdb_type == input$collection, , drop = FALSE]
    }
    data <- data[is.finite(data$padj) & data$padj <= input$padj, , drop = FALSE]
    if (input$direction == "positive") data <- data[data$NES > 0, , drop = FALSE]
    if (input$direction == "negative") data <- data[data$NES < 0, , drop = FALSE]
    data
  })

  output$results_table <- renderDT({
    data <- filtered_results()
    shown <- intersect(c("msigdb_type", "pathway", "NES", "pval", "padj", "size"),
                       names(data))
    datatable(data[, shown, drop = FALSE], selection = "single", rownames = FALSE,
              options = list(pageLength = 10, scrollX = TRUE))
  })

  selected_row <- reactive({
    data <- filtered_results()
    selected <- input$results_table_rows_selected
    if (!length(selected) && nrow(data)) selected <- 1L
    req(length(selected), nrow(data) >= selected)
    data[selected, , drop = FALSE]
  })

  output$selection_title <- renderText(selected_row()$pathway[[1]])

  output$enrichment_plot <- renderPlot({
    row <- selected_row()
    pathway <- row$pathway[[1]]
    genes <- pathways[[pathway]]
    req(length(genes))
    fgsea::plotEnrichment(genes, rankings[[input$ranking]]) +
      ggtitle(pathway) + theme_bw(base_size = 13)
  })

  output$leading_edge_table <- renderDT({
    row <- selected_row()
    genes <- character()
    if ("leadingEdge" %in% names(row) && !is.na(row$leadingEdge[[1]])) {
      genes <- trimws(strsplit(as.character(row$leadingEdge[[1]]), ",")[[1]])
    }
    scores <- rankings[[input$ranking]][genes]
    data <- data.frame(gene = genes, score = unname(scores), stringsAsFactors = FALSE)
    data <- data[is.finite(data$score), , drop = FALSE]
    datatable(data, rownames = FALSE, options = list(pageLength = 10))
  })
}

message("SnakeGSEA local explorer listening on http://", options$host, ":", options$port)
runApp(shinyApp(ui, server), host = options$host, port = options$port,
       launch.browser = FALSE)
