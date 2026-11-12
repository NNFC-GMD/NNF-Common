# Reference run of scSLIDE in Seurat 6 on the Jost et al. 2020 Perturb-seq data.
#
# The Python notebook is written from the paper's Methods. This script runs
# Seurat itself on the same cells and with the same parameters, so the two
# results can be compared (compare.py). It calls documented functions, in the
# order of the scSLIDE vignette:
# https://satijalab.org/seurat/articles/scslide_vignette
#
# There is one exception, marked below: an internal function is called once, to
# write out two numbers per cell that FindmmNN computes and does not keep. We
# learned of them by reading Seurat's source (MIT licence) after the first
# comparison left one difference unexplained.
#
# Seurat 6 lives in an isolated library, .r-lib/ at the repository root, so the
# machine's own R library is left alone. Install it once with:
#
#   Rscript -e 'install.packages(c("Seurat", "pls", "RcppThread"), lib = ".r-lib", dependencies = FALSE)'
#
# Run from the repository root:
#
#   Rscript validation/scslide/seurat_reference.R --seed 123
#
# Outputs go to validation/scslide/out/, which git ignores.

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(flag, default) {
    position <- match(flag, args)
    if (is.na(position)) default else args[position + 1]
}
seed <- as.integer(arg_value("--seed", "123"))
data_dir <- arg_value("--data-dir", "data/GSE132080")
out_dir <- arg_value("--out", "validation/scslide/out")
r_lib <- arg_value("--r-lib", ".r-lib")

.libPaths(c(normalizePath(r_lib), .libPaths()))
suppressPackageStartupMessages({
    library(Matrix)
    library(Seurat)
})
stopifnot(packageVersion("Seurat") >= "6.0.0")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# The same parameters as notebooks/scslide/prepare.py.
n_variable <- 2000
n_pca <- 30
n_pls <- 20
n_landmarks <- 5000
k_landmarks <- 5

# --- Read the GEO files ------------------------------------------------------
geo <- function(name) file.path(data_dir, paste0("GSE132080_", name))
counts <- as(readMM(gzfile(geo("10X_matrix.mtx.gz"))), "CsparseMatrix")
barcodes <- readLines(gzfile(geo("10X_barcodes.tsv.gz")))
genes <- read.delim(gzfile(geo("10X_genes.tsv.gz")), header = FALSE, col.names = c("id", "symbol"))
rownames(counts) <- make.unique(genes$symbol)
colnames(counts) <- barcodes
identities <- read.csv(gzfile(geo("cell_identities.csv.gz")), stringsAsFactors = FALSE)

# --- Quality control ---------------------------------------------------------
assigned <- identities[
    identities$number_of_cells == 1 &
        identities$good_coverage %in% c("True", "TRUE", TRUE) &
        identities$guide_identity != "*",
]
stopifnot(nrow(assigned) == 19587, length(unique(assigned$guide_identity)) == 138)
rownames(assigned) <- assigned$cell_barcode

# Guides get short identifiers, because Seurat rewrites some characters in names.
guide_names <- sort(unique(assigned$guide_identity), method = "radix")
guides <- data.frame(
    id = sprintf("g%03d", seq_along(guide_names)),
    guide_identity = guide_names,
    stringsAsFactors = FALSE
)
assigned$guide <- guides$id[match(assigned$guide_identity, guides$guide_identity)]
assigned$gene <- ifelse(
    startsWith(assigned$guide_identity, "neg_ctrl_"),
    "non-targeting",
    sub("_.*$", "", assigned$guide_identity)
)

object <- CreateSeuratObject(
    counts = counts[, assigned$cell_barcode],
    meta.data = assigned[, c("guide", "gene", "guide_identity", "gemgroup")]
)
n_before <- ncol(object)
object <- subset(object, subset = nFeature_RNA >= 500 & nCount_RNA >= 1000)
cat(sprintf("%d cells with one assigned guide; the thresholds remove %d\n", n_before, n_before - ncol(object)))

# --- A standard Seurat object: normalise, variable genes, PCA ----------------
object <- NormalizeData(object, verbose = FALSE)
object <- FindVariableFeatures(object, nfeatures = n_variable, verbose = FALSE)
object <- ScaleData(object, verbose = FALSE)
object <- RunPCA(object, npcs = 50, verbose = FALSE)

# --- scSLIDE, step by step as in the vignette --------------------------------
# 1. Training cells: up to 1000 per guide, which here is every cell.
object <- SketchDataByGroup(
    object,
    assay = "RNA",
    ncells = 1000,
    sketched.assay = "TRAINING",
    method = "Uniform",
    group.by = "guide",
    verbose = FALSE
)
# 2. Landmark cells, by leverage score. The seed is the only thing that varies between runs.
object <- SketchData(
    object,
    assay = "TRAINING",
    ncells = n_landmarks,
    sketched.assay = "LANDMARK",
    method = "LeverageScore",
    seed = seed,
    verbose = FALSE
)
# 3. Genes correlated with the label, added to the variable genes.
label_genes <- QuickCorTest(
    object,
    assay = "TRAINING",
    layer = "data",
    group.by = NULL,
    Y = "gene",
    verbose = FALSE
)
pls_features <- union(VariableFeatures(object, assay = "RNA"), label_genes)
# 4. Canonical PLS against the target gene.
DefaultAssay(object) <- "TRAINING"
object <- ScaleData(object, features = pls_features, verbose = FALSE)
object <- RunPLS(
    object = object,
    assay = "TRAINING",
    features = pls_features,
    ncomp = n_pls,
    Y = c("gene"),
    pls.function = "cppls",
    verbose = FALSE
)
# 5. Project every cell onto the PLS components.
projected <- ProjectCellEmbeddings(
    query = object,
    query.assay = "RNA",
    reference = object,
    reference.assay = "TRAINING",
    reduction = "pls",
    dims = 1:n_pls,
    scale = TRUE,
    normalization.method = "LogNormalize",
    verbose = FALSE
)
object[["proj.pls"]] <- CreateDimReducObject(embeddings = projected, key = "pPLS_", assay = "RNA")
# 6. The joint graph from every cell to the landmarks, weights fixed at 0.5 and 0.5.
DefaultAssay(object) <- "RNA"
object <- FindmmNN(
    object,
    sketch.assay = "LANDMARK",
    reduction.list = list("pca", "proj.pls"),
    k.nn = k_landmarks,
    nn.name = "weighted.nn",
    dims.list = list(1:n_pca, 1:n_pls),
    fix.wnn.weights = c(0.5, 0.5)
)
# 7. The landmark-by-guide matrix, chi-square normalised.
sample_object <- GenerateSampleObject(
    object = object,
    sketch.assay = "LANDMARK",
    nn.name = "weighted.nn",
    group.by = "guide",
    k.nn = k_landmarks,
    normalization.method = "ChiSquared",
    add.meta.data = TRUE,
    return.seurat = TRUE
)

# --- Write what compare.py needs ----------------------------------------------
residuals <- as.matrix(LayerData(sample_object, assay = "LMC", layer = "data"))
landmark_counts <- as.matrix(LayerData(sample_object, assay = "LMC", layer = "counts"))
# Seurat may have rewritten the guide identifiers; map the columns back.
column_id <- guides$id[match(gsub("[^A-Za-z0-9]", "", colnames(residuals)), guides$id)]
stopifnot(!anyNA(column_id), !anyDuplicated(column_id))
colnames(residuals) <- column_id
colnames(landmark_counts) <- column_id
correlation <- cor(residuals[, guides$id])

tag <- sprintf("seed%d", seed)
write.csv(guides, file.path(out_dir, "guides.csv"), row.names = FALSE)
write.csv(correlation, file.path(out_dir, sprintf("correlation_%s.csv", tag)))
write.csv(
    landmark_counts[, guides$id],
    gzfile(file.path(out_dir, sprintf("landmark_counts_%s.csv.gz", tag)))
)

cells <- data.frame(
    cell_barcode = colnames(object),
    guide = object$guide,
    landmark = colnames(object) %in% colnames(object[["LANDMARK"]]),
    Embeddings(object, "pca")[, 1:n_pca],
    Embeddings(object, "proj.pls")[, 1:n_pls]
)
write.csv(cells, gzfile(file.path(out_dir, sprintf("cells_%s.csv.gz", tag))), row.names = FALSE)

# Every cell's landmarks as FindmmNN stored them, best first, so that compare.py
# can check them cell by cell. FindmmNN keeps 20 per cell and GenerateSampleObject
# reads the first k_landmarks of them. A neighbour is stored as a row number of
# the landmark-only object that FindmmNN builds, so that object's cell order is
# written next to it. The stored value is sqrt((1 - joint affinity) / 2).
neighbours <- object[["weighted.nn"]]
stopifnot(identical(Cells(neighbours), colnames(object)))
landmark_order <- colnames(subset(object, cells = Cells(object[["LANDMARK"]])))
neighbour_table <- data.frame(
    cell_barcode = Cells(neighbours),
    index = Indices(neighbours),
    value = Distances(neighbours)
)
write.csv(
    neighbour_table,
    gzfile(file.path(out_dir, sprintf("neighbours_%s.csv.gz", tag))),
    row.names = FALSE
)
writeLines(landmark_order, file.path(out_dir, sprintf("landmark_order_%s.txt", tag)))

# The one call to an internal function. To judge how near a landmark is,
# FindmmNN first measures every cell's surroundings among all cells: the
# distance to its nearest other cell (the offset) and how much farther its
# 20th nearest row lies (the width). It finds them with an approximate search
# and does not store them. The same internal call, with the arguments FindmmNN
# passes, gives them again, so compare.py can tell the rule from the search.
surroundings <- slot(
    Seurat:::FindModalityWeights(
        object = object,
        reduction.list = list("pca", "proj.pls"),
        dims.list = list(1:n_pca, 1:n_pls),
        fix.wnn.weights = c(0.5, 0.5),
        snn.far.nn = FALSE,
        l2.norm = TRUE,
        k.nn = 20,
        verbose = FALSE
    ),
    "params"
)
scale_table <- data.frame(
    cell_barcode = colnames(object),
    offset_pca = surroundings$nearest.dist[[1]],
    offset_pls = surroundings$nearest.dist[[2]],
    width_pca = surroundings$sigma.list[[1]],
    width_pls = surroundings$sigma.list[[2]]
)
write.csv(
    scale_table,
    gzfile(file.path(out_dir, sprintf("scales_%s.csv.gz", tag))),
    row.names = FALSE
)

summary_lines <- c(
    sprintf("seurat %s", packageVersion("Seurat")),
    sprintf("seed %d", seed),
    sprintf("cells %d", ncol(object)),
    sprintf("guides %d", ncol(residuals)),
    sprintf("landmarks %d", nrow(residuals)),
    sprintf("pls_features %d", length(pls_features)),
    sprintf("label_genes %d", length(label_genes)),
    sprintf("landmark_count_total %d", sum(landmark_counts))
)
writeLines(summary_lines, file.path(out_dir, sprintf("summary_%s.txt", tag)))
cat(summary_lines, sep = "\n")
