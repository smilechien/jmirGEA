# jmirGEA
gene expression analysis
For SampleID & Condition in genecluster7.csv
https://www.raschonline.com/raschonline/cbpabc.asp?cbp=ROCtop10eachcluster
# jmirGEA
# jmirGEA

An R Shiny app for differential expression analysis and FLCA clustering of node observations. The app starts in GEO mode with `GSE79973` and a built-in Case/Control sample table. It can also analyze an uploaded expression matrix or run FLCA on a node-by-observation CSV.

## Run locally

Install R and the packages used by the app. `GEOquery`, `limma`, `Biobase`, `AnnotationDbi`, and `hgu133plus2.db` are Bioconductor packages; the others are available from CRAN.

```r
install.packages(c(
  "shiny", "ggplot2", "pheatmap", "ggrepel", "pROC",
  "dplyr", "tibble", "igraph", "circlize", "readxl", "psych"
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

The sample information box accepts the two columns as tab- or space-separated text or as a Markdown table with `|` borders. Text uploads encoded in UTF-8 or CP950 are supported. The factor analysis panel requires more complete, variable feature rows than sample columns; it cannot estimate sample loadings for a matrix with 31 rows and 202 samples. Nutrition measurements are not gene expression, and raw nutrient amounts are not normalized log2 values, so interpret any differential results in that context.

The app produces a limma results table, volcano plot, selected-gene heatmap and expression distributions, ROC curves with 96% DeLong confidence intervals, and Lin's concordance correlation coefficient (CCC) networks, clusters, and a vertical hierarchical gene dendrogram. A bar aligned beneath the dendrogram's gene labels marks up-regulated genes in red and down-regulated genes in blue. It also displays FLCA Sankey, SSplot, Kano, and chord plots, plus cluster quality and clustering algorithm comparisons. A two-factor analysis treats complete gene rows as observations and sample columns as variables, then compares each sample's factor loadings with its Case/Control label. The **Download all results** button saves the limma table as CSV; separate buttons save ROC and cluster summaries.

ROC estimates are calculated on the same samples used to select genes. Validate predictive performance on independent data before using the AUC for prediction.

### Two-factor sample comparison

The factor analysis uses all complete, variable gene rows, with **sample columns as variables**. It standardizes each gene across samples, then fits a two-factor minimum-residual model with varimax rotation using `psych`. The sample table shows Factor 1 and Factor 2 loadings beside each sample's Case/Control label. Each sample is assigned to the factor with the larger absolute loading; a condition-by-factor table, best in-sample match percentage, and adjusted Rand index summarize agreement. A parallel-analysis plot and the first three PCA eigenvalues are also shown. These comparisons are exploratory and do not validate a classifier.

## FLCA workflow for node observations

Clear **Sample information**, choose **Upload expression file**, and upload a CSV such as [`data/nurseratio.csv`](data/nurseratio.csv). Then click **Run analysis / FLCA Process**. The CSV needs unique node names in its first column and at least two columns of finite numeric observations. An optional `Profile` column is ignored. The app analyzes up to 20 nodes, selected by their mean observation value.

```csv
Node,Time1,Time2,Time3,Profile
Node_A,6.4,7.0,7.1,group1
Node_B,6.9,7.2,7.4,group2
```

FLCA builds a network of positive pairwise Lin's CCC values, retains the strongest positive link for each node in a reduced network, and uses connected components as clusters. It displays the networks, CCC heatmap, Sankey, SSplot, Kano, and chord plots, along with cluster and algorithm quality tables. Buttons download the cluster and comparison tables as CSV.

While FLCA runs, a percentage progress bar shows file reading, pairwise CCC calculation, network construction, cluster assignment, and summary preparation. The selected-gene CCC/FLCA calculation uses the same progress display.

You can also select **Online GitHub link** and supply a public GitHub `blob` or `raw.githubusercontent.com` URL. A URL passed as `?url=` selects that source automatically. Selecting the link source clears the sample table and starts FLCA, which requires a **CSV** file. To analyze a linked expression matrix instead, enter its Case/Control sample table and run the analysis again; that path also accepts TSV, TXT, XLSX, and gzipped text files.

## Repository files

- [`app.R`](app.R): Shiny interface, differential expression analysis, ROC, and CCC views.
- [`flca_process.R`](flca_process.R): FLCA input handling, networks, plots, and quality comparison.
- [`flca_ms_sil_module.R`](flca_ms_sil_module.R), [`renderSSplot.R`](renderSSplot.R), [`kano.R`](kano.R): local FLCA visualization and quality modules.
- [`data/`](data/): example node-observation CSV files.
- `test_*.R`: local test scripts for upload and FLCA behavior.

Run a test script from the project directory with `Rscript test_upload.R` or `Rscript test_flca.R` after installing the packages above.

 - For nutritionvba.csf
- SampleID	Condition
H001	Case
H002	Case
H003	Case
H004	Case
H005	Case
H006	Case
H007	Case
H008	Case
H009	Case
H010	Case
H011	Case
H012	Case
H013	Case
H014	Case
H015	Case
H016	Case
H017	Case
H018	Case
H019	Case
H020	Case
H021	Case
H022	Case
H023	Case
H024	Case
H025	Case
H026	Case
H027	Case
H028	Case
H029	Case
H030	Case
H031	Case
H032	Case
H033	Case
H034	Case
H035	Case
H036	Case
H037	Case
H038	Case
H039	Case
H040	Case
H041	Case
H042	Case
H043	Case
H044	Case
H045	Case
H046	Case
H047	Case
H048	Case
H049	Case
H050	Case
H051	Case
H052	Case
H053	Case
H054	Case
H055	Case
H056	Case
H057	Case
H058	Case
H059	Case
H060	Case
H061	Case
H062	Case
H063	Case
H064	Case
H065	Case
H066	Case
H067	Case
H068	Case
H069	Case
H070	Case
H071	Case
H072	Case
H073	Case
H074	Case
H075	Case
H076	Case
H077	Case
H078	Case
H079	Case
H080	Case
H081	Case
H082	Case
H083	Case
H084	Case
H085	Case
H086	Case
H087	Case
H088	Case
H089	Case
H090	Case
H091	Case
H092	Case
H093	Case
H094	Case
H095	Case
H096	Case
H097	Case
H098	Case
H099	Case
H100	Case
H101	Case
H102	Case
N001	Control
N002	Control
N003	Control
N004	Control
N005	Control
N006	Control
N007	Control
N008	Control
N009	Control
N010	Control
N011	Control
N012	Control
N013	Control
N014	Control
N015	Control
N016	Control
N017	Control
N018	Control
N019	Control
N020	Control
N021	Control
N022	Control
N023	Control
N024	Control
N025	Control
N026	Control
N027	Control
N028	Control
N029	Control
N030	Control
N031	Control
N032	Control
N033	Control
N034	Control
N035	Control
N036	Control
N037	Control
N038	Control
N039	Control
N040	Control
N041	Control
N042	Control
N043	Control
N044	Control
N045	Control
N046	Control
N047	Control
N048	Control
N049	Control
N050	Control
N051	Control
N052	Control
N053	Control
N054	Control
N055	Control
N056	Control
N057	Control
N058	Control
N059	Control
N060	Control
N061	Control
N062	Control
N063	Control
N064	Control
N065	Control
N066	Control
N067	Control
N068	Control
N069	Control
N070	Control
N071	Control
N072	Control
N073	Control
N074	Control
N075	Control
N076	Control
N077	Control
N078	Control
N079	Control
N080	Control
N081	Control
N082	Control
N083	Control
N084	Control
N085	Control
N086	Control
N087	Control
N088	Control
N089	Control
N090	Control
N091	Control
N092	Control
N093	Control
N094	Control
N095	Control
N096	Control
N097	Control
N098	Control
N099	Control
N100	Control
<img width="179" height="4467" alt="image" src="https://github.com/user-attachments/assets/b141f636-3339-42d9-bff9-0a667313ad48" />

- [`flca_ms_sil_module.R`](flca_ms_sil_module.R), [`renderSSplot.R`](renderSSplot.R), [`kano.R`](kano.R): local FLCA visualization and quality modules.
- [`data/`](data/): example node-observation CSV files.
- `test_*.R`: local test scripts for upload and FLCA behavior.
Run a test script from the project directory with `Rscript test_upload.R` or `Rscript test_flca.R` after installing the packages above.
