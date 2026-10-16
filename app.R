library(shiny)
library(GEOquery)
library(limma)
library(ggplot2)
library(pheatmap)
library(ggrepel)
library(pROC)
flca_tools <- new.env(parent = globalenv())
sys.source("flca_ms_sil_module.R", envir = flca_tools)
sys.source("renderSSplot.R", envir = flca_tools)
sys.source("kano.R", envir = flca_tools)
source("flca_process.R", local = TRUE)
`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L) y else x
options(shiny.maxRequestSize = 100 * 1024^2)

default_samples <- paste(c(
  "SampleID\tCondition",
  "GSM2109532\tCase",
  "GSM2109533\tControl",
  "GSM2109534\tCase",
  "GSM2109535\tControl",
  "GSM2109536\tCase",
  "GSM2109537\tControl",
  "GSM2109538\tCase",
  "GSM2109539\tControl",
  "GSM2109540\tCase",
  "GSM2109541\tControl",
  "GSM2109542\tCase",
  "GSM2109543\tControl",
  "GSM2109544\tCase",
  "GSM2109545\tControl",
  "GSM2109546\tCase",
  "GSM2109547\tControl",
  "GSM2109548\tCase",
  "GSM2109549\tControl",
  "GSM2109550\tCase",
  "GSM2109551\tControl"
), collapse = "\n")

ui <- fluidPage(
  tags$head(tags$style(HTML("
    .upload-input .btn-file {
      background: #18794e; border-color: #18794e; color: #fff;
      font-weight: 700;
    }
    .upload-input .btn-file:hover,
    .upload-input .btn-file:focus {
      background: #12623f; border-color: #12623f; color: #fff;
    }
    #clear_samples {
      background: #b45309; border-color: #b45309; color: #fff;
      font-weight: 700;
    }
    #clear_samples:hover, #clear_samples:focus {
      background: #92400e; border-color: #92400e; color: #fff;
    }
    #run { font-weight: 700; }
  "))),
  titlePanel("Differential expression: volcano, heatmap, ROC, and CCC clusters"),
  sidebarLayout(
    sidebarPanel(
      radioButtons("data_source", "Expression source",
                   choices = c("GEO series" = "geo", "Upload expression file" = "upload", "Online GitHub link" = "link"),
                   selected = "geo", inline = TRUE),
      textInput("accession", "GEO series accession", value = "GSE79973"),
      tags$div(class = "upload-input",
               fileInput("expression_file", "Input data file (CSV, TSV, TXT, XLSX, or gzipped text)",
                         accept = c(".csv", ".tsv", ".txt", ".xlsx", ".gz"),
                         buttonLabel = "Browse input data")),
      textInput("online_expression_url", "Online GitHub data URL", value = "",
                placeholder = "https://github.com/owner/repo/blob/main/expression.csv"),
      p("Choose Online GitHub link for a public CSV, TSV, TXT, XLSX, or gzipped file. A GitHub blob URL is converted to the matching raw download. A ?url= link selects this third source automatically."),
      p("For uploads: feature IDs in the first column, sample IDs as column headings; an optional GeneSymbol column is supported. Use normalized log2 expression values. Uploaded IDs and values are used as supplied, without probe annotation lookup."),
      textAreaInput(
        "sample_info", "Sample information (SampleID and Condition)",
        value = default_samples, rows = 22, width = "100%"
      ),
      actionButton("clear_samples", "Clear sample information"),
      p("With this box empty, upload a nurse-ratio-style CSV (node names in column 1, numeric observations in the remaining columns, optional Profile) and run FLCA Process."),
      numericInput("p_cutoff", "Significance cutoff", value = 0.05,
                   min = 0, max = 1, step = 0.01),
      selectInput("p_metric", "Significance measure",
                  choices = c("Raw P-value" = "raw", "Adjusted FDR" = "fdr"),
                  selected = "raw"),
      numericInput("fc_cutoff", "Absolute log2 fold-change cutoff",
                   value = 1, min = 0, step = 0.1),
      selectInput("heatmap_direction", "Heatmap genes",
                  choices = c("Top 10 up- and top 10 down-regulated genes" = "both",
                              "Top 10 up-regulated genes" = "up",
                              "Top 10 down-regulated genes" = "down"),
                  selected = "both"),
      selectInput("rank_by", "Rank selected genes by",
                  choices = c("Selected significance measure" = "p",
                              "Fold change" = "fc"),
                  selected = "p"),
      sliderInput("label_font_size", "Network and Kano name font size",
                  min = 2, max = 7, value = 3.5, step = 0.25),
      actionButton("run", "Run analysis / FLCA Process", class = "btn-primary"),
      br(), br(),
      downloadButton("download_results", "Download all results")
    ),
    mainPanel(
      conditionalPanel("input.sample_info && input.sample_info.trim().length > 0",
      p("Positive log2 fold change means higher expression in Case than Control."),
      textOutput("summary"),
      plotOutput("volcano", height = "600px"),
      h4("Heatmap of selected genes"),
      p("Each row shows expression standardized across the selected samples (blue = lower, red = higher). The most significant probe represents each gene."),
      imageOutput("heatmap", height = "700px"),
      h4("Expression distributions of selected genes"),
      plotOutput("density", height = "350px"),
      plotOutput("sample_boxplot", height = "400px"),
      h4("Gene ROC curve"),
      p("AUC and 96% DeLong confidence intervals use the same samples used to select the genes. Perfect in-sample AUC can produce a 1.00 to 1.00 interval; independent validation is needed before predictive use."),
      selectInput("roc_feature", "Gene for ROC curve", choices = character(0)),
      plotOutput("roc_plot", height = "450px"),
      h4("AUC for selected up- and down-regulated genes"),
      tableOutput("roc_table"),
      downloadButton("download_roc", "Download ROC summary"),
      h4("Lin's CCC clusters of selected genes"),
      p("The top 10 up- and top 10 down-regulated genes are nodes. Lin's CCC compares their expression across the selected samples. Positive CCC pairs form the full network; the strongest positive link for each gene forms the reduced network. Clusters are connected components of the reduced network."),
      p("Bubble size shows the sum of each gene's expression scores. Edge width shows CCC; bubble color shows cluster and shape shows regulation. Node positions use the unweighted network topology."),
      textOutput("ccc_summary"),
      plotOutput("ccc_full_plot", height = "650px"),
      plotOutput("ccc_reduced_plot", height = "650px"),
      imageOutput("ccc_heatmap", height = "650px"),
      h4("FLCA Process plots for selected genes"),
      p("value = sum of expression scores across sample columns; value2 = sum of incident positive CCC links. Displayed values use two decimals."),
      plotOutput("ccc_sankey", height = "850px"),
      plotOutput("ccc_ssplot", height = "1050px"),
      h4("SSplot quality by cluster"),
      p("SS is mean silhouette width. Qw and Qu are each cluster's contributions to weighted and unweighted modularity; their contributions sum to the OVERALL value."),
      tableOutput("ccc_cluster_quality"),
      tags$div(style = "margin-bottom: 48px; clear: both;",
               plotOutput("ccc_kano", height = "850px")),
      tags$div(style = "margin-top: 32px; clear: both;",
               plotOutput("ccc_chord", height = "650px")),
      h4("CCC cluster membership"),
      tableOutput("ccc_table"),
      downloadButton("download_ccc", "Download CCC clusters"),
      h4("Quality: comparison of clustering algorithms"),
      tableOutput("ccc_algorithm_quality"),
      downloadButton("download_ccc_quality", "Download algorithm comparison"),
      h4("Top 20 features by P-value"),
      tableOutput("results_table")
      ),
      conditionalPanel("!input.sample_info || input.sample_info.trim().length === 0",
        h3("FLCA Process: uploaded node time series"),
        p("value = sum of scores across observation columns; value2 = sum of incident positive CCC links. Displayed values use two decimals."),
        p("Bubble size shows value. Edge width shows CCC; bubble color shows cluster. Node positions use the unweighted network topology."),
        textOutput("flca_summary"),
        h4("Step 2. Positive pairwise Lin's CCC network"),
        plotOutput("flca_full_plot", height = "650px"),
        h4("Step 3. Strongest-link FLCA network"),
        plotOutput("flca_reduced_plot", height = "650px"),
        h4("Step 4. CCC heatmap"),
        imageOutput("flca_ccc_heatmap", height = "650px"),
        h4("Step 5. Sankey plot"),
        plotOutput("flca_sankey", height = "850px"),
        h4("Step 6. SSplot cluster summary"),
        plotOutput("flca_ssplot", height = "1050px"),
        h4("SSplot quality by cluster"),
        p("SS is mean silhouette width. Qw and Qu are each cluster's contributions to weighted and unweighted modularity; their contributions sum to the OVERALL value."),
        tableOutput("flca_cluster_quality"),
        h4("Step 6b. Kano plot"),
        tags$div(style = "margin-bottom: 48px; clear: both;",
                 plotOutput("flca_kano", height = "850px")),
        h4("Step 6c. Chord plot"),
        tags$div(style = "margin-top: 32px; clear: both;",
                 plotOutput("flca_chord", height = "650px")),
        h4("FLCA clusters"),
        tableOutput("flca_table"),
        downloadButton("download_flca", "Download FLCA clusters"),
        h4("Quality: comparison of clustering algorithms"),
        tableOutput("flca_algorithm_quality"),
        downloadButton("download_flca_quality", "Download algorithm comparison")
      )
    )
  )
)

make_volcano <- function(x) {
  labeled <- x$results[!is.na(x$results$GeneSymbol) &
                         nzchar(x$results$GeneSymbol) &
                         x$results$group != "Not significant", , drop = FALSE]
  labeled <- do.call(rbind, lapply(c("Up-regulated", "Down-regulated"),
                                   function(group) {
    rows <- labeled[labeled$group == group, , drop = FALSE]
    head(rows[!duplicated(rows$GeneSymbol), , drop = FALSE], 10)
  }))
  ggplot(x$results, aes(x = logFC, y = logP, color = group)) +
    geom_point(alpha = 0.4, size = 1.75, na.rm = TRUE) +
    scale_color_manual(values = c("Down-regulated" = "blue",
                                  "Not significant" = "grey",
                                  "Up-regulated" = "red"), drop = FALSE) +
    geom_vline(xintercept = c(-x$fc_cutoff, x$fc_cutoff),
               linetype = "dotted") +
    geom_hline(yintercept = -log10(x$p_cutoff), linetype = "dotted") +
    labs(title = paste("Volcano plot:", x$accession),
         x = "log2 fold change (Case vs Control)",
         y = if (x$p_metric == "fdr") "-log10 adjusted FDR"
             else "-log10 P-value",
         color = NULL) +
    ggrepel::geom_text_repel(
      data = labeled, aes(label = GeneSymbol), size = 3,
      max.overlaps = Inf, show.legend = FALSE, seed = 1
    ) +
    theme_bw(base_size = 15) +
    theme(legend.position = "bottom",
          axis.title = element_text(face = "bold"),
          axis.text = element_text(face = "bold"))
}

bold_heatmap_axes <- function(heatmap) {
  for (name in c("row_names", "col_names", "row_annotation_names",
                 "col_annotation_names")) {
    index <- which(heatmap$gtable$layout$name == name)
    if (length(index)) {
      gp <- heatmap$gtable$grobs[[index[1L]]]$gp
      gp$fontface <- NULL
      gp$font <- 2
      heatmap$gtable$grobs[[index[1L]]]$gp <- gp
    }
  }
  heatmap
}

get_gene_symbols <- function(eset, probe_ids) {
  symbols <- rep(NA_character_, length(probe_ids))
  names(symbols) <- probe_ids
  features <- Biobase::fData(eset)
  symbol_column <- grep("^(gene[._ ]?symbol|symbol)$", names(features),
                        ignore.case = TRUE, value = TRUE)
  if (length(symbol_column) > 0) {
    symbols <- as.character(features[probe_ids, symbol_column[1]])
    names(symbols) <- probe_ids
  }
  if (Biobase::annotation(eset) == "GPL570" &&
      requireNamespace("hgu133plus2.db", quietly = TRUE) &&
      requireNamespace("AnnotationDbi", quietly = TRUE)) {
    missing <- is.na(symbols) | !nzchar(symbols)
    if (any(missing)) {
      symbols[missing] <- unname(AnnotationDbi::mapIds(
        hgu133plus2.db::hgu133plus2.db,
        keys = probe_ids[missing], column = "SYMBOL",
        keytype = "PROBEID", multiVals = "first"
      ))
    }
  }
  trimws(symbols)
}

read_expression_upload <- function(file_info, sample_ids) {
  path <- file_info$datapath
  file_name <- file_info$name
  lower_name <- tolower(file_name)
  if (grepl("\\.xlsx$", lower_name)) {
    validate(need(requireNamespace("readxl", quietly = TRUE),
                  "Install the readxl package to upload XLSX files."))
    data <- as.data.frame(readxl::read_xlsx(path, .name_repair = "minimal"),
                          check.names = FALSE)
  } else {
    first_line <- local({
      header_connection <- if (grepl("\\.gz$", lower_name))
        gzfile(path) else file(path)
      on.exit(close(header_connection))
      readLines(header_connection, n = 1, warn = FALSE)
    })
    validate(need(length(first_line) > 0, "The uploaded file is empty."))
    connection <- if (grepl("\\.gz$", lower_name)) gzfile(path) else file(path)
    separator <- if (grepl("\\t", first_line)) "\\t" else if (
      grepl(",", first_line)) "," else ""
    data <- read.table(connection, header = TRUE, sep = separator,
                       check.names = FALSE, stringsAsFactors = FALSE,
                       comment.char = "", quote = "\"", fill = FALSE)
  }
  validate(need(ncol(data) >= 3 && nrow(data) > 0,
                "Upload a matrix with feature IDs in column 1 and samples in the other columns."))
  names(data) <- trimws(names(data))
  validate(need(!anyDuplicated(names(data)),
                "Uploaded column headings must be unique."))
  feature_ids <- trimws(as.character(data[[1]]))
  validate(need(!anyNA(feature_ids) && all(nzchar(feature_ids)),
                "Feature IDs in the first column must be filled in."))
  sample_columns <- match(sample_ids, names(data))
  validate(need(!anyNA(sample_columns) && all(sample_columns > 1),
                "Every SampleID in the textarea must match a sample column in the upload."))

  raw_values <- as.matrix(data[, sample_columns, drop = FALSE])
  values <- suppressWarnings(matrix(
    as.numeric(raw_values), nrow = nrow(raw_values),
    dimnames = list(feature_ids, sample_ids)
  ))
  invalid <- is.na(values) & !is.na(raw_values) &
    !(trimws(as.character(raw_values)) %in% c("", "NA", "NaN"))
  validate(need(!any(invalid), "Sample columns must contain numeric expression values."),
           need(any(is.finite(values)), "The uploaded matrix has no usable expression values."))

  symbol_column <- which(tolower(names(data)) %in%
                           c("genesymbol", "gene_symbol", "gene symbol", "symbol"))
  symbol_column <- setdiff(symbol_column, 1)
  symbols <- if (length(symbol_column)) {
    trimws(as.character(data[[symbol_column[1]]]))
  } else {
    rep(NA_character_, length(feature_ids))
  }
  symbols[is.na(symbols) | !nzchar(symbols)] <-
    feature_ids[is.na(symbols) | !nzchar(symbols)]

  if (anyDuplicated(feature_ids)) {
    counts <- rowsum(1L * (!is.na(values)), feature_ids, reorder = FALSE)
    sums <- rowsum(replace(values, is.na(values), 0),
                   feature_ids, reorder = FALSE)
    values <- sums / counts
    values[counts == 0] <- NA_real_
    symbols <- symbols[match(rownames(values), feature_ids)]
  }
  names(symbols) <- rownames(values)
  list(expression = values, symbols = symbols, label = file_name)
}

# Download a public GitHub-hosted expression file for the Online GitHub link
# mode.  Restricting hosts prevents a shared URL from accessing local networks.
github_expression_raw_url <- function(source_url) {
  source_url <- trimws(as.character(source_url)[1L])
  if (is.na(source_url) || !nzchar(source_url)) stop("Provide a GitHub data URL.", call. = FALSE)
  source_url <- sub("#.*$", "", source_url)
  blob <- regexec("^https://github\\.com/([^/]+)/([^/]+)/blob/([^/]+)/(.+\\.(?:csv|tsv|txt|xlsx|gz))(?:\\?.*)?$", source_url, perl = TRUE)
  parts <- regmatches(source_url, blob)[[1L]]
  if (length(parts)) return(paste0("https://raw.githubusercontent.com/", parts[2L], "/", parts[3L], "/", parts[4L], "/", parts[5L]))
  if (grepl("^https://raw\\.githubusercontent\\.com/[^/]+/[^/]+/[^/]+/.+\\.(?:csv|tsv|txt|xlsx|gz)(?:\\?.*)?$", source_url, perl = TRUE)) return(source_url)
  stop("Use an HTTPS GitHub blob or raw URL ending in .csv, .tsv, .txt, .xlsx, or .gz.", call. = FALSE)
}

download_online_expression <- function(source_url) {
  raw_url <- github_expression_raw_url(source_url)
  url_path <- sub("\\?.*$", "", raw_url)
  extension <- sub("^.*(\\.(?:csv|tsv|txt|xlsx|gz))$", "\\1", url_path, ignore.case = TRUE)
  if (!grepl("^\\.(csv|tsv|txt|xlsx|gz)$", extension, ignore.case = TRUE)) extension <- ".csv"
  dest <- tempfile(fileext = extension)
  status <- suppressWarnings(utils::download.file(raw_url, dest, mode = "wb", quiet = TRUE, method = "libcurl"))
  if (!identical(status, 0L) || !file.exists(dest) || file.size(dest) == 0L)
    stop("The online GitHub data file could not be downloaded.", call. = FALSE)
  list(datapath = dest, name = basename(url_path), source_url = raw_url)
}
select_regulated <- function(x, direction = "both", rank_by = "p") {
  if (is.null(rank_by) || length(rank_by) == 0) rank_by <- "p"
  top_rows <- function(group) {
    rows <- x$results[x$results$group == group, , drop = FALSE]
    if (rank_by == "fc") {
      rows <- rows[order(if (group == "Up-regulated") -rows$logFC
                         else rows$logFC, rows$P.Value), , drop = FALSE]
    } else {
      rank_p <- if (x$p_metric == "fdr") rows$adj.P.Val else rows$P.Value
      rows <- rows[order(rank_p, rows$P.Value), , drop = FALSE]
    }
    key <- rows$GeneSymbol
    key[is.na(key) | !nzchar(key)] <- rows$ProbeID[is.na(key) | !nzchar(key)]
    head(rows[!duplicated(key), , drop = FALSE], 10)
  }
  up <- top_rows("Up-regulated")
  down <- top_rows("Down-regulated")
  switch(direction, up = up, down = down,
         both = rbind(up, down), rbind(up, down))
}

make_heatmap <- function(x, direction = "both", rank_by = "p") {
  selected <- select_regulated(x, direction, rank_by)
  probes <- selected$ProbeID
  if (length(probes) == 0) return(NULL)

  values <- x$expression[probes, , drop = FALSE]
  row_sd <- apply(values, 1, sd, na.rm = TRUE)
  values <- values[is.finite(row_sd) & row_sd > 0, , drop = FALSE]
  if (nrow(values) == 0) return(NULL)
  row_mean <- rowMeans(values, na.rm = TRUE)
  for (i in seq_len(nrow(values))) {
    values[i, is.na(values[i, ])] <- row_mean[i]
  }
  z <- t(scale(t(values)))
  z[z < -2.5] <- -2.5
  z[z > 2.5] <- 2.5

  group <- x$condition[colnames(z)]
  annotation_col <- data.frame(Group = group, row.names = colnames(z))
  regulation <- as.character(selected$group[match(rownames(z), selected$ProbeID)])
  annotation_row <- data.frame(Regulation = regulation,
                               row.names = rownames(z))
  annotation_colors <- list(
    Group = c(Control = "#35C5C4", Case = "#EE9B99"),
    Regulation = c("Up-regulated" = "#D73027",
                   "Down-regulated" = "#4575B4")
  )
  row_labels <- x$results$GeneSymbol[match(rownames(z), x$results$ProbeID)]
  row_labels[is.na(row_labels) | !nzchar(row_labels)] <-
    rownames(z)[is.na(row_labels) | !nzchar(row_labels)]
  row_labels <- make.unique(row_labels, sep = " #")
  bold_heatmap_axes(pheatmap(
    z, color = colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(101),
    breaks = seq(-2.5, 2.5, length.out = 102),
    annotation_col = annotation_col,
    annotation_row = annotation_row,
    annotation_colors = annotation_colors,
    cluster_rows = nrow(z) > 1, cluster_cols = ncol(z) > 1,
    labels_row = row_labels,
    fontsize_row = 9, fontsize_col = 8, angle_col = 45,
    border_color = NA, main = paste("Top regulated genes:", x$accession),
    silent = TRUE
  ))
}

selected_expression_long <- function(x, rank_by = "p") {
  selected <- select_regulated(x, "both", rank_by)
  if (nrow(selected) == 0) return(NULL)
  values <- x$expression[selected$ProbeID, , drop = FALSE]
  data.frame(
    Expression = as.vector(values),
    SampleID = rep(colnames(values), each = nrow(values)),
    Condition = factor(rep(x$condition[colnames(values)],
                           each = nrow(values)), levels = c("Control", "Case"))
  )
}

calculate_roc <- function(x, rank_by = "p") {
  selected <- select_regulated(x, "both", rank_by)
  if (nrow(selected) == 0) return(list(table = data.frame(), curves = list()))
  curves <- vector("list", nrow(selected))
  names(curves) <- selected$ProbeID
  summary_rows <- lapply(seq_len(nrow(selected)), function(i) {
    feature_id <- selected$ProbeID[i]
    score <- as.numeric(x$expression[feature_id, ])
    group <- factor(x$condition[colnames(x$expression)],
                    levels = c("Control", "Case"))
    keep <- is.finite(score) & !is.na(group)
    curve <- NULL
    confidence <- c(NA_real_, NA_real_, NA_real_)
    if (all(table(group[keep]) >= 2) && length(unique(score[keep])) >= 2) {
      curve <- tryCatch(
        pROC::roc(group[keep], score[keep],
                  levels = c("Control", "Case"),
                  direction = if (selected$logFC[i] >= 0) "<" else ">",
                  quiet = TRUE), error = function(e) NULL
      )
      if (!is.null(curve)) {
        confidence <- tryCatch(
          as.numeric(pROC::ci.auc(curve, conf.level = 0.96,
                                  method = "delong")),
          error = function(e) confidence
        )
      }
    }
    curves[[i]] <<- curve
    data.frame(
      GeneSymbol = selected$GeneSymbol[i],
      ProbeID = feature_id,
      Regulation = as.character(selected$group[i]),
      logFC = selected$logFC[i],
      AUC = if (is.null(curve)) NA_real_ else as.numeric(pROC::auc(curve)),
      CI96_Lower = confidence[1], CI96_Upper = confidence[3],
      Samples = sum(keep)
    )
  })
  list(table = do.call(rbind, summary_rows), curves = curves)
}

lin_ccc <- function(a, b) {
  keep <- is.finite(a) & is.finite(b)
  if (sum(keep) < 3L) return(NA_real_)
  a <- a[keep]; b <- b[keep]
  ma <- mean(a); mb <- mean(b)
  den <- mean((a - ma)^2) + mean((b - mb)^2) + (ma - mb)^2
  if (!is.finite(den) || den <= 0) return(NA_real_)
  2 * mean((a - ma) * (b - mb)) / den
}

network_components <- function(labels, edges) {
  membership <- setNames(rep(NA_integer_, length(labels)), labels)
  cluster <- 0L
  for (node in labels) {
    if (!is.na(membership[node])) next
    cluster <- cluster + 1L
    queue <- node
    membership[node] <- cluster
    while (length(queue)) {
      current <- queue[1L]
      queue <- queue[-1L]
      neighbors <- c(edges$To[edges$From == current],
                     edges$From[edges$To == current])
      neighbors <- neighbors[is.na(membership[neighbors])]
      if (length(neighbors)) {
        membership[neighbors] <- cluster
        queue <- c(queue, neighbors)
      }
    }
  }
  membership
}

calculate_ccc_clusters <- function(x, rank_by = "p") {
  selected <- select_regulated(x, "both", rank_by)
  if (nrow(selected) < 2L) return(NULL)
  ids <- selected$ProbeID
  values <- x$expression[ids, , drop = FALSE]
  labels <- selected$GeneSymbol
  labels[is.na(labels) | !nzchar(labels)] <- ids[is.na(labels) | !nzchar(labels)]
  labels <- make.unique(labels, sep = " #")
  n <- length(ids)
  ccc <- matrix(NA_real_, n, n, dimnames = list(labels, labels))
  diag(ccc) <- 1
  pairs <- utils::combn(seq_len(n), 2L)
  for (j in seq_len(ncol(pairs))) {
    i1 <- pairs[1L, j]; i2 <- pairs[2L, j]
    ccc[i1, i2] <- ccc[i2, i1] <- lin_ccc(values[i1, ], values[i2, ])
  }
  positive <- is.finite(ccc[cbind(pairs[1L, ], pairs[2L, ])]) &
    ccc[cbind(pairs[1L, ], pairs[2L, ])] > 0
  full <- data.frame(From = labels[pairs[1L, positive]],
                     To = labels[pairs[2L, positive]],
                     CCC = ccc[cbind(pairs[1L, positive], pairs[2L, positive])])
  reduced <- full[FALSE, , drop = FALSE]
  if (nrow(full)) {
    strongest <- character(0)
    for (node in labels) {
      incident <- which(full$From == node | full$To == node)
      if (length(incident)) {
        best <- incident[order(-full$CCC[incident],
                               full$From[incident], full$To[incident])][1L]
        strongest <- c(strongest, best)
      }
    }
    reduced <- full[sort(unique(as.integer(strongest))), , drop = FALSE]
  }
  membership <- network_components(labels, reduced)
  edge_value <- ccc_incident_strength(labels, full)
  coordinates <- network_topology_layout(labels, reduced)
  nodes <- data.frame(
    Gene = labels, FeatureID = ids,
    Regulation = as.character(selected$group),
    Cluster = as.integer(membership[labels]),
    X = coordinates[, 1L], Y = coordinates[, 2L],
    value = rowSums(values, na.rm = TRUE),
    value2 = edge_value,
    stringsAsFactors = FALSE
  )
  nodes$Cluster <- factor(nodes$Cluster)
  best_ccc <- vapply(labels, function(node) {
    incident <- reduced$CCC[reduced$From == node | reduced$To == node]
    if (length(incident)) max(incident) else NA_real_
  }, numeric(1))
  summary <- nodes[c("Gene", "FeatureID", "Regulation", "Cluster")]
  summary$value <- sprintf("%.2f", nodes$value)
  summary$value2 <- sprintf("%.2f", nodes$value2)
  summary$Strongest_CCC <- ifelse(is.finite(best_ccc),
                                 sprintf("%.2f", unname(best_ccc)), "")
  list(nodes = nodes, full = full, reduced = reduced,
       matrix = ccc, table = summary)
}

make_ccc_network_plot <- function(result, reduced = FALSE, label_size = 3.5) {
  if (is.null(label_size) || length(label_size) != 1L ||
      !is.finite(label_size)) label_size <- 3.5
  edges <- if (reduced) result$reduced else result$full
  nodes <- result$nodes
  node_kind <- if (all(nodes$Regulation == "Node")) "nodes" else "genes"
  edges$X <- edges$Y <- edges$Xend <- edges$Yend <- numeric(nrow(edges))
  if (nrow(edges)) {
    edges$X <- nodes$X[match(edges$From, nodes$Gene)]
    edges$Y <- nodes$Y[match(edges$From, nodes$Gene)]
    edges$Xend <- nodes$X[match(edges$To, nodes$Gene)]
    edges$Yend <- nodes$Y[match(edges$To, nodes$Gene)]
  }
  ggplot() +
    geom_segment(data = edges,
                 aes(x = X, y = Y, xend = Xend, yend = Yend,
                     linewidth = CCC), color = "#67788B",
                 alpha = if (reduced) 0.7 else 0.28) +
    scale_linewidth(range = if (reduced) c(0.5, 2.5) else c(0.2, 1.2),
                    limits = c(0, 1)) +
    geom_point(data = nodes, aes(x = X, y = Y, fill = Cluster,
                                 shape = Regulation, size = value),
               color = "#25364A") +
    scale_size_continuous(name = "Sum of observation scores (value)",
                          range = c(3.5, 10),
                          labels = function(x) sprintf("%.2f", x)) +
    scale_shape_manual(values = c("Up-regulated" = 21,
                                  "Down-regulated" = 24,
                                  "Node" = 21)) +
    ggrepel::geom_text_repel(data = nodes, aes(x = X, y = Y, label = Gene),
                             size = label_size, fontface = "bold",
                             max.overlaps = Inf, seed = 1,
                             box.padding = 0.4) +
    coord_equal() +
    labs(title = if (reduced) "Strongest positive CCC link per gene (FLCA reduced)"
         else "Positive pairwise Lin's CCC network",
         subtitle = sprintf("%d %s, %d positive links, %d reduced-network clusters",
                            nrow(nodes), node_kind, nrow(edges),
                            length(unique(nodes$Cluster))),
         x = NULL, y = NULL, linewidth = "Lin's CCC") +
    theme_void(base_size = 13) +
    theme(legend.position = "bottom", legend.box = "vertical",
          plot.title = element_text(face = "bold"),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          plot.background = element_rect(fill = "white", color = NA),
          panel.background = element_rect(fill = "white", color = NA))
}

make_ccc_heatmap <- function(result) {
  z <- result$matrix
  z[!is.finite(z)] <- 0
  annotation <- data.frame(
    Cluster = result$nodes$Cluster,
    Regulation = result$nodes$Regulation,
    row.names = result$nodes$Gene
  )
  bold_heatmap_axes(pheatmap::pheatmap(
    z, color = colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(101),
    breaks = seq(-1, 1, length.out = 102),
    cluster_rows = nrow(z) > 2, cluster_cols = ncol(z) > 2,
    annotation_row = annotation, annotation_col = annotation,
    fontsize_row = 9, fontsize_col = 9, border_color = NA,
    main = if (all(result$nodes$Regulation == "Node"))
      "Lin's CCC across uploaded nodes" else
      "Lin's CCC across selected genes", silent = TRUE
  ))
}

server <- function(input, output, session) {
  auto_link_run <- reactiveVal(0L)
  pending_link_url <- reactiveVal(NULL)
  clear_link_samples_once <- reactiveVal(FALSE)
  force_link_mode_once <- reactiveVal(FALSE)

  run_link_as_flca <- function(url = NULL) {
    if (!is.null(url) && nzchar(trimws(url))) pending_link_url(trimws(url))
    clear_link_samples_once(TRUE)
    force_link_mode_once(TRUE)
    updateTextAreaInput(session, "sample_info", value = "")
    if (!is.null(pending_link_url()) || nzchar(trimws(input$online_expression_url %||% "")))
      auto_link_run(auto_link_run() + 1L)
    else
      shiny::showNotification("Enter an Online GitHub data URL to run the FLCA workflow.", type = "message", duration = 7)
  }

  observeEvent(input$clear_samples, {
    updateTextAreaInput(session, "sample_info", value = "")
  })
  observeEvent(input$expression_file, {
    updateRadioButtons(session, "data_source", selected = "upload")
  })
  observeEvent(input$data_source, {
    if (identical(input$data_source, "link")) run_link_as_flca()
  }, ignoreInit = TRUE)
  observeEvent(session$clientData$url_search, {
    query <- shiny::parseQueryString(session$clientData$url_search)
    if (is.null(query$url) || !nzchar(trimws(query$url))) return(invisible(NULL))
    pending_link_url(query$url)
    updateTextInput(session, "online_expression_url", value = query$url)
    updateRadioButtons(session, "data_source", selected = "link")
    run_link_as_flca(query$url)
  }, once = TRUE, ignoreNULL = TRUE)

  analysis <- eventReactive(list(input$run, auto_link_run()), {
    source_mode <- if (isTRUE(force_link_mode_once())) "link" else input$data_source %||% "geo"
    accession <- toupper(trimws(input$accession))
    link_url <- trimws(input$online_expression_url %||% "")
    if (!nzchar(link_url)) link_url <- pending_link_url() %||% ""
    use_empty_samples <- identical(source_mode, "link") && isTRUE(clear_link_samples_once())
    sample_text <- if (use_empty_samples) "" else input$sample_info
    clear_link_samples_once(FALSE)
    force_link_mode_once(FALSE)
    online_file <- NULL
    if (identical(source_mode, "link")) {
      validate(need(nzchar(link_url), "Enter an online GitHub data URL."))
      online_file <- tryCatch(download_online_expression(link_url), error = function(e) e)
      validate(need(!inherits(online_file, "error"),
                    if (inherits(online_file, "error")) conditionMessage(online_file) else "Could not download the online data file."))
      on.exit(unlink(online_file$datapath), add = TRUE)
    }
    input_file <- if (identical(source_mode, "link")) online_file else input$expression_file

    if (source_mode %in% c("upload", "link") && !nzchar(trimws(sample_text))) {
      validate(need(!is.null(input_file), "Choose or link a nurse-ratio-style CSV for FLCA Process."))
      values <- read_flca_upload(input_file)
      return(list(mode = "flca", flca = calculate_flca_upload(values, input_file$name)))
    }
    if (identical(source_mode, "geo")) {
      validate(need(grepl("^GSE[0-9]+$", accession),
                    "Enter a GEO series accession such as GSE79973."))
    }
    samples <- tryCatch(
      read.table(text = sample_text, header = TRUE, sep = "",
                 stringsAsFactors = FALSE, check.names = FALSE,
                 comment.char = "", quote = ""),
      error = function(e) NULL
    )
    validate(need(!is.null(samples) &&
                    identical(names(samples), c("SampleID", "Condition")),
                  "Paste a two-column table headed SampleID and Condition."))
    samples$SampleID <- trimws(samples$SampleID)
    samples$Condition <- trimws(samples$Condition)
    validate(
      need(nrow(samples) >= 4 && !anyNA(samples) &&
             all(nzchar(samples$SampleID)) &&
             all(samples$Condition %in% c("Case", "Control")),
           "Provide at least four samples, each labeled Case or Control."),
      need(!anyDuplicated(samples$SampleID), "Sample IDs must be unique."),
      need(all(table(factor(samples$Condition,
                            levels = c("Control", "Case"))) >= 2),
           "Provide at least two samples in each condition."),
      need(is.finite(input$p_cutoff) && input$p_cutoff > 0 &&
             input$p_cutoff < 1, "P-value cutoff must be between 0 and 1."),
      need(is.finite(input$fc_cutoff) && input$fc_cutoff >= 0,
           "Fold-change cutoff must be nonnegative.")
    )

    if (source_mode %in% c("upload", "link")) {
      validate(need(!is.null(input_file), "Choose or link an expression matrix file."))
      uploaded <- read_expression_upload(input_file, samples$SampleID)
      expression <- uploaded$expression
      symbols <- uploaded$symbols
      source_label <- uploaded$label
    } else {
      local_matrix <- file.path("data", paste0(accession, "_series_matrix.txt.gz"))
      geo_data <- withProgress(message = paste("Loading", accession), value = 0, {
        if (file.exists(local_matrix)) {
          list(getGEO(filename = local_matrix, getGPL = FALSE))
        } else {
          getGEO(accession, GSEMatrix = TRUE, getGPL = FALSE)
        }
      })
      platform_matches <- which(vapply(
        geo_data,
        function(x) all(samples$SampleID %in% colnames(exprs(x))),
        logical(1)
      ))
      validate(need(length(platform_matches) > 0,
                    "The sample IDs were not all found on one GEO platform."))
      eset <- geo_data[[platform_matches[1]]]
      expression <- exprs(eset)[, samples$SampleID, drop = FALSE]
      symbols <- get_gene_symbols(eset, rownames(expression))
      source_label <- accession
    }
    validate(need(is.numeric(expression) && nrow(expression) > 0 &&
                    any(is.finite(expression)),
                  "The selected GEO series has no usable numeric expression matrix."))
    if (identical(source_mode, "geo") &&
        quantile(expression, 0.99, na.rm = TRUE) > 100) {
      validate(need(all(expression >= 0, na.rm = TRUE),
                    "Expression values include negatives; automatic log2 transformation is unsafe."))
      expression <- log2(expression + 1)
    }
    expression[!is.finite(expression)] <- NA_real_

    condition <- factor(samples$Condition, levels = c("Control", "Case"))
    design <- model.matrix(~ condition)
    rownames(design) <- samples$SampleID
    fit <- eBayes(lmFit(expression, design))
    results <- topTable(fit, coef = "conditionCase", number = Inf,
                        adjust.method = "BH", sort.by = "P")
    results <- data.frame(ProbeID = rownames(results), results,
                          row.names = NULL, check.names = FALSE)
    results$GeneSymbol <- unname(symbols[results$ProbeID])
    p_metric <- if (identical(input$p_metric, "fdr")) "fdr" else "raw"
    selected_p <- if (p_metric == "fdr") results$adj.P.Val else results$P.Value
    results$logP <- -log10(pmax(selected_p, .Machine$double.xmin))
    results$group <- "Not significant"
    significant <- !is.na(selected_p) & selected_p < input$p_cutoff
    results$group[significant & results$logFC > input$fc_cutoff] <-
      "Up-regulated"
    results$group[significant & results$logFC < -input$fc_cutoff] <-
      "Down-regulated"
    results$group <- factor(results$group,
                            levels = c("Down-regulated", "Not significant",
                                       "Up-regulated"))
    list(mode = "expression", results = results, expression = expression,
         condition = setNames(as.character(condition), samples$SampleID),
         accession = source_label, p_metric = p_metric,
         p_cutoff = input$p_cutoff, fc_cutoff = input$fc_cutoff)
  }, ignoreNULL = FALSE)

  expression_analysis <- reactive({
    x <- analysis()
    req(identical(x$mode, "expression"))
    x
  })
  flca_analysis <- reactive({
    x <- analysis()
    req(identical(x$mode, "flca"))
    x$flca
  })

  output$summary <- renderText({
    x <- expression_analysis()
    counts <- table(x$results$group)
    sprintf("%s: %d features; %d up-regulated and %d down-regulated at the selected cutoffs.",
            x$accession, nrow(x$results), counts["Up-regulated"],
            counts["Down-regulated"])
  })

  output$volcano <- renderPlot({
    make_volcano(expression_analysis())
  })

  output$heatmap <- renderImage({
    heatmap <- make_heatmap(expression_analysis(), input$heatmap_direction, input$rank_by)
    validate(need(!is.null(heatmap),
                  "No eligible probes for this heatmap at the selected cutoffs."))
    image_file <- tempfile(fileext = ".png")
    png(image_file, width = 1800, height = 1300, res = 160)
    on.exit(dev.off(), add = TRUE)
    grid::grid.newpage()
    grid::grid.draw(heatmap$gtable)
    list(src = image_file, contentType = "image/png", width = "100%",
         alt = "Heatmap of the selected up- and down-regulated genes")
  }, deleteFile = TRUE)

  output$density <- renderPlot({
    long <- selected_expression_long(expression_analysis(), input$rank_by)
    validate(need(!is.null(long), "No regulated genes meet the selected cutoffs."))
    long <- long[is.finite(long$Expression), , drop = FALSE]
    ggplot(long, aes(x = Expression, fill = Condition)) +
      geom_density(alpha = 0.4) +
      scale_fill_manual(values = c(Control = "#35C5C4", Case = "#EE9B99")) +
      labs(title = "Expression density of selected genes") +
      theme_bw(base_size = 13) +
      theme(axis.title = element_text(face = "bold"),
            axis.text = element_text(face = "bold"))
  })

  output$sample_boxplot <- renderPlot({
    long <- selected_expression_long(expression_analysis(), input$rank_by)
    validate(need(!is.null(long), "No regulated genes meet the selected cutoffs."))
    long <- long[is.finite(long$Expression), , drop = FALSE]
    ggplot(long, aes(x = SampleID, y = Expression, fill = Condition)) +
      geom_boxplot(outlier.size = 0.7) +
      scale_fill_manual(values = c(Control = "#35C5C4", Case = "#EE9B99")) +
      labs(title = "Selected-gene expression by sample", x = "Sample ID") +
      theme_bw(base_size = 12) +
      theme(axis.title = element_text(face = "bold"),
            axis.text = element_text(face = "bold"),
            axis.text.x = element_text(angle = 60, hjust = 1,
                                       face = "bold"))
  })

  roc_summary <- reactive(calculate_roc(expression_analysis(), input$rank_by))
  observeEvent(roc_summary(), {
    rows <- roc_summary()$table
    if (nrow(rows) > 0) {
      labels <- paste0(rows$GeneSymbol, " (", rows$ProbeID, ")")
      updateSelectInput(session, "roc_feature",
                        choices = setNames(rows$ProbeID, labels),
                        selected = rows$ProbeID[1])
    }
  })

  output$roc_plot <- renderPlot({
    roc_result <- roc_summary()
    validate(need(nrow(roc_result$table) > 0,
                  "No regulated genes meet the selected cutoffs."))
    feature_id <- input$roc_feature
    if (is.null(feature_id) ||
        !feature_id %in% roc_result$table$ProbeID) {
      feature_id <- roc_result$table$ProbeID[1]
    }
    curve <- roc_result$curves[[feature_id]]
    validate(need(!is.null(curve),
                  "This gene has insufficient expression values for an ROC curve."))
    row <- roc_result$table[match(feature_id, roc_result$table$ProbeID), ]
    points <- data.frame(FPR = 1 - as.numeric(curve$specificities),
                         Sensitivity = as.numeric(curve$sensitivities))
    points <- points[order(points$FPR, points$Sensitivity), ]
    ggplot(points, aes(x = FPR, y = Sensitivity)) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                  color = "grey60") +
      geom_path(linewidth = 1.1, color = "#2166AC") +
      coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
      labs(title = paste("ROC:", row$GeneSymbol),
           subtitle = sprintf("%s | AUC %.3f\n96%% CI %.3f to %.3f",
                              row$ProbeID, row$AUC,
                              row$CI96_Lower, row$CI96_Upper),
           x = "1 - specificity", y = "Sensitivity") +
      theme_bw(base_size = 14) +
      theme(axis.title = element_text(face = "bold"),
            axis.text = element_text(face = "bold"))
  })

  output$roc_table <- renderTable({
    roc_summary()$table
  }, digits = 3)

  output$download_roc <- downloadHandler(
    filename = function() paste0(gsub("[^A-Za-z0-9_-]", "_",
                                    expression_analysis()$accession), "_ROC_96CI.csv"),
    content = function(file) write.csv(roc_summary()$table, file,
                                       row.names = FALSE)
  )

  ccc_clusters <- reactive(calculate_ccc_clusters(expression_analysis(), input$rank_by))
  ccc_cluster_metrics <- reactive({
    ccc <- ccc_clusters()
    validate(need(!is.null(ccc), "No regulated nodes are available."))
    flca_cluster_metrics(ccc)
  })
  ccc_algorithm_quality <- reactive({
    ccc <- ccc_clusters()
    validate(need(!is.null(ccc), "No regulated nodes are available."))
    flca_quality_comparison(ccc)
  })
  output$ccc_summary <- renderText({
    ccc <- ccc_clusters()
    validate(need(!is.null(ccc),
                  "At least two regulated genes are needed for CCC clustering."))
    sprintf("%d genes; %d positive CCC links; %d strongest links; %d FLCA clusters.",
            nrow(ccc$nodes), nrow(ccc$full), nrow(ccc$reduced),
            length(unique(ccc$nodes$Cluster)))
  })
  output$ccc_full_plot <- renderPlot({
    ccc <- ccc_clusters()
    validate(need(!is.null(ccc),
                  "At least two regulated genes are needed for CCC clustering."))
    make_ccc_network_plot(ccc, label_size = input$label_font_size)
  }, width = 1100, height = 650, res = 110)
  output$ccc_reduced_plot <- renderPlot({
    ccc <- ccc_clusters()
    validate(need(!is.null(ccc),
                  "At least two regulated genes are needed for CCC clustering."))
    make_ccc_network_plot(ccc, reduced = TRUE,
                          label_size = input$label_font_size)
  }, width = 1100, height = 650, res = 110)
  output$ccc_heatmap <- renderImage({
    ccc <- ccc_clusters()
    validate(need(!is.null(ccc),
                  "At least two regulated genes are needed for a CCC heatmap."))
    heatmap <- make_ccc_heatmap(ccc)
    image_file <- tempfile(fileext = ".png")
    png(image_file, width = 1600, height = 1400, res = 150)
    on.exit(dev.off(), add = TRUE)
    grid::grid.newpage()
    grid::grid.draw(heatmap$gtable)
    list(src = image_file, contentType = "image/png", width = "100%",
         alt = "Lin's CCC matrix of selected genes")
  }, deleteFile = TRUE)
  output$ccc_sankey <- renderPlot({
    ccc <- ccc_clusters()
    validate(need(!is.null(ccc), "No regulated nodes are available."))
    plot_flca_sankey(ccc)
  }, width = 1400, height = 850, res = 110)
  output$ccc_ssplot <- renderPlot({
    ccc <- ccc_clusters()
    validate(need(!is.null(ccc), "No regulated nodes are available."))
    plot_flca_ssplot(ccc)
  }, width = 1400, height = 1050, res = 110)
  output$ccc_cluster_quality <- renderTable({
    format_quality_table(ccc_cluster_metrics()[c("Cluster", "n", "SS",
                                               "Qw", "Qu")])
  })
  output$ccc_kano <- renderPlot({
    ccc <- ccc_clusters()
    validate(need(!is.null(ccc), "No regulated nodes are available."))
    plot_flca_kano(ccc, label_size = input$label_font_size)
  }, width = 1100, height = 850, res = 110)
  output$ccc_chord <- renderPlot({
    ccc <- ccc_clusters()
    validate(need(!is.null(ccc), "No regulated nodes are available."))
    plot_flca_chord(ccc)
  }, width = 850, height = 650, res = 110)
  output$ccc_table <- renderTable({
    ccc <- ccc_clusters()
    validate(need(!is.null(ccc), "No CCC clusters are available."))
    ccc$table
  }, digits = 3)
  output$download_ccc <- downloadHandler(
    filename = function() paste0(gsub("[^A-Za-z0-9_-]", "_",
                                    expression_analysis()$accession), "_CCC_clusters.csv"),
    content = function(file) write.csv(ccc_clusters()$table, file,
                                       row.names = FALSE)
  )
  output$ccc_algorithm_quality <- renderTable({
    format_quality_table(ccc_algorithm_quality())
  })
  output$download_ccc_quality <- downloadHandler(
    filename = function() paste0(gsub("[^A-Za-z0-9_-]", "_",
                                    expression_analysis()$accession),
                                  "_algorithm_quality.csv"),
    content = function(file) write.csv(format_quality_table(ccc_algorithm_quality()),
                                       file, row.names = FALSE)
  )

  output$results_table <- renderTable({
    head(expression_analysis()$results[c("ProbeID", "GeneSymbol", "logFC", "P.Value",
                              "adj.P.Val", "group")], 20)
  }, digits = 4)

  output$download_results <- downloadHandler(
    filename = function() paste0(gsub("[^A-Za-z0-9_-]", "_",
                                    expression_analysis()$accession), "_limma_results.csv"),
    content = function(file) write.csv(expression_analysis()$results, file, row.names = FALSE)
  )

  output$flca_summary <- renderText({
    x <- flca_analysis()
    sprintf("%s: %d nodes across %d observations; %d positive CCC links; %d strongest links; %d FLCA clusters.",
            x$label, nrow(x$nodes), length(x$observations),
            nrow(x$full), nrow(x$reduced),
            length(unique(x$nodes$Cluster)))
  })
  output$flca_full_plot <- renderPlot({
    make_ccc_network_plot(flca_analysis(), label_size = input$label_font_size) +
      labs(title = "Step 2: positive pairwise Lin's CCC network")
  }, width = 1100, height = 650, res = 110)
  output$flca_reduced_plot <- renderPlot({
    make_ccc_network_plot(flca_analysis(), reduced = TRUE,
                          label_size = input$label_font_size) +
      labs(title = "Step 3: strongest-link FLCA network")
  }, width = 1100, height = 650, res = 110)
  output$flca_ccc_heatmap <- renderImage({
    heatmap <- make_ccc_heatmap(flca_analysis())
    image_file <- tempfile(fileext = ".png")
    png(image_file, width = 1600, height = 1400, res = 150)
    on.exit(dev.off(), add = TRUE)
    grid::grid.newpage()
    grid::grid.draw(heatmap$gtable)
    list(src = image_file, contentType = "image/png", width = "100%",
         alt = "Lin's CCC heatmap of uploaded nodes")
  }, deleteFile = TRUE)
  output$flca_sankey <- renderPlot({ plot_flca_sankey(flca_analysis()) },
                                   width = 1400, height = 850, res = 110)
  output$flca_ssplot <- renderPlot({ plot_flca_ssplot(flca_analysis()) },
                                   width = 1400, height = 1050, res = 110)
  output$flca_cluster_quality <- renderTable({
    format_quality_table(flca_cluster_metrics(flca_analysis())[
      c("Cluster", "n", "SS", "Qw", "Qu")])
  })
  output$flca_kano <- renderPlot({
    plot_flca_kano(flca_analysis(), label_size = input$label_font_size)
  },
                                 width = 1100, height = 850, res = 110)
  output$flca_chord <- renderPlot({ plot_flca_chord(flca_analysis()) },
                                  width = 850, height = 650, res = 110)
  output$flca_table <- renderTable({ flca_analysis()$table }, digits = 3)
  output$download_flca <- downloadHandler(
    filename = function() paste0(gsub("[^A-Za-z0-9_-]", "_",
                                    flca_analysis()$label), "_FLCA_clusters.csv"),
    content = function(file) write.csv(flca_analysis()$table, file,
                                       row.names = FALSE)
  )
  flca_algorithm_quality <- reactive(flca_quality_comparison(flca_analysis()))
  output$flca_algorithm_quality <- renderTable({
    format_quality_table(flca_algorithm_quality())
  })
  output$download_flca_quality <- downloadHandler(
    filename = function() paste0(gsub("[^A-Za-z0-9_-]", "_",
                                    flca_analysis()$label),
                                  "_algorithm_quality.csv"),
    content = function(file) write.csv(format_quality_table(flca_algorithm_quality()),
                                       file, row.names = FALSE)
  )
}

shinyApp(ui, server, options = list(launch.browser = TRUE))
