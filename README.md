# jmirGEA
gene expression analysis
# jmirGEA
An R Shiny app for differential expression analysis and FLCA clustering of node observations. The app starts in GEO mode with `GSE79973` and a built-in Case/Control sample table. It can also analyze an uploaded expression matrix or run FLCA on a node-by-observation CSV.
## Run locally
Install R and the packages used by the app. `GEOquery`, `limma`, `Biobase`, `AnnotationDbi`, and `hgu133plus2.db` are Bioconductor packages; the others are available from CRAN.
```r
install.packages(c(
  "shiny", "ggplot2", "pheatmap", "ggrepel", "pROC",
  "dplyr", "tibble", "igraph", "circlize", "readxl"
))
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c(
  "GEOquery", "limma", "Biobase", "AnnotationDbi", "hgu133plus2.db"
))
```
From the project directory, run:
```r
shiny::runApp(".")
```
The four local R modules (`flca_process.R`, `flca_ms_sil_module.R`, `renderSSplot.R`, and `kano.R`) must stay beside `app.R`. An internet connection is needed to download GEO series or files supplied by GitHub URL. `readxl` is needed for XLSX input; `AnnotationDbi` and `hgu133plus2.db` provide a fallback symbol lookup for GPL570 probes. The chord plot uses `circlize`.
## Differential expression workflow
1. Choose **GEO series** (the default) and enter a `GSE` accession, or choose **Upload expression file** and select a CSV, TSV, TXT, XLSX, or gzipped text matrix.
2. Provide the sample table in the app's **Sample information** box. Its headings must be exactly `SampleID` and `Condition`; conditions must be `Case` or `Control`, with at least two samples in each group. For an upload, every `SampleID` must match a matrix column.
3. Set the significance measure (raw P-value or adjusted FDR), P-value cutoff, absolute log2 fold-change cutoff, gene direction, and ranking. Click **Run analysis / FLCA Process**.
Expression upload format:
```csv
FeatureID,GeneSymbol,Sample_A,Sample_B,Sample_C,Sample_D
F1,GENE1,8.2,8.5,6.1,6.3
F2,GENE2,5.4,5.2,7.8,7.6
```
```text
SampleID	Condition
Sample_A	Case
Sample_B	Case
Sample_C	Control
Sample_D	Control
```
Put feature IDs in the first column and sample IDs in the column headings. `GeneSymbol` is optional. Supply **normalized log2 expression values** for uploaded matrices; the app uses those values as given. GEO data may be automatically log2 transformed when its values indicate an unlogged scale. Positive log2 fold change means higher expression in Case than Control.
The app produces a limma results table, volcano plot, selected-gene heatmap and expression distributions, ROC curves with 96% DeLong confidence intervals, and Lin's concordance correlation coefficient (CCC) networks and clusters. It also displays FLCA Sankey, SSplot, Kano, and chord plots, plus cluster quality and clustering algorithm comparisons. The **Download all results** button saves the limma table as CSV; separate buttons save ROC and cluster summaries.
ROC estimates are calculated on the same samples used to select genes. Validate predictive performance on independent data before using the AUC for prediction.
## FLCA workflow for node observations
Clear **Sample information**, choose **Upload expression file**, and upload a CSV such as [`data/nurseratio.csv`](data/nurseratio.csv). Then click **Run analysis / FLCA Process**. The CSV needs unique node names in its first column and at least two columns of finite numeric observations. An optional `Profile` column is ignored. The app analyzes up to 20 nodes, selected by their mean observation value.
```csv
Node,Time1,Time2,Time3,Profile
Node_A,6.4,7.0,7.1,group1
Node_B,6.9,7.2,7.4,group2
```
FLCA builds a network of positive pairwise Lin's CCC values, retains the strongest positive link for each node in a reduced network, and uses connected components as clusters. It displays the networks, CCC heatmap, Sankey, SSplot, Kano, and chord plots, along with cluster and algorithm quality tables. Buttons download the cluster and comparison tables as CSV.
You can also select **Online GitHub link** and supply a public GitHub `blob` or `raw.githubusercontent.com` URL. A URL passed as `?url=` selects that source automatically. Selecting the link source clears the sample table and starts FLCA, which requires a **CSV** file. To analyze a linked expression matrix instead, enter its Case/Control sample table and run the analysis again; that path also accepts TSV, TXT, XLSX, and gzipped text files.
## Repository files
- [`app.R`](app.R): Shiny interface, differential expression analysis, ROC, and CCC views.
- [`flca_process.R`](flca_process.R): FLCA input handling, networks, plots, and quality comparison.
- [`flca_ms_sil_module.R`](flca_ms_sil_module.R), [`renderSSplot.R`](renderSSplot.R), [`kano.R`](kano.R): local FLCA visualization and quality modules.
- [`data/`](data/): example node-observation CSV files.
- `test_*.R`: local test scripts for upload and FLCA behavior.
Run a test script from the project directory with `Rscript test_upload.R` or `Rscript test_flca.R` after installing the packages above.
