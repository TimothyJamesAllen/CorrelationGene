library(Seurat)
library(glmGamPoi)
library(Matrix)
library(tictoc)
library(Hmisc)
library(moments)
library(biomaRt)
obj = readRDS("/Users/artemis/Downloads/obj_new.rds")
randomgenelist <- c("FOXJ1",
                    "PPIE",
                    "C19orf52",
                    "PMP22",
                    "MGAM",
                    "TOM1",
                    "LINC01520",
                    "HIST2H2AB",
                    "SH2B3",
                    "DCAF6",
                    "ECA1",
                    "FCN1",
                    "ASPN",
                    "WDR20",
                    "SLC5A7",
                    "SNORD118",
                    "CYBA",
                    "UBXN2B",
                    "NAGA",
                    "FARP2",
                    "TRA-AGC2-1",
                    "C9orf69",
                    "OR13C2",
                    "FGF11",
                    "RBX1",
                    "DPRX",
                    "OR8J3",
                    "SHQ1",
                    "C10orf25",
                    "SLC12A5",
                    "SLC25A47",
                    "C4orf17",
                    "COLEC10",
                    "AMT",
                    "C5orf64")
correlationtest = function(df, goi, gene_list, PorS) {

  cat("Computing correlation matrix...\n")
  if (PorS == "P") {
    df <- rcorr(as.matrix(df), type = "pearson")
  } else if (PorS == "S") {
    df <- rcorr(as.matrix(df), type = "spearman")
  } else {
    cat("Error: PorS must be either P or S")
  }
  df_r = as.data.frame(df$r)
  df_p = as.data.frame(df$P)
  cat("extracting data...\n")
  goi_only_p = as.data.frame(df_p[goi], row.names = rownames(df_p))
  goi_only_r = as.data.frame(df_r[goi], row.names = rownames(df_r))
  cat("adding gene names...\n")
  goi_only_p$gene = rownames(goi_only_p)
  goi_only_r$gene = rownames(goi_only_r)
  merged = merge(goi_only_r, goi_only_p, by = "gene")
  colnames(merged)[2:3] <- c("R","P")  
  merged = merged[order(merged$R, decreasing = TRUE), ]
  merged$SFARI.Gene = ""
  cat("Merging data frames...\n")
  common_genes = intersect(merged$gene, gene_list)
  
  is_sfari <- ifelse(SFARI_genes$gene.symbol %in% common_genes, TRUE, FALSE)
  is_sfari <- is_sfari[1:nrow(merged)]

cat("Adding gene annotations...\n")
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
# Loop through each gene ID in the merged data frame
for (i in 1:nrow(merged)) {
  gene_id <- merged$gene[i]
  gene_info <- getBM(attributes = c("external_gene_name", "description"), filters = "external_gene_name", values = gene_id, mart = mart)
  if (nrow(gene_info) > 0) {
    merged$description[i] <- gene_info$description[1]
  } else {
    merged$description[i] <- "N/A"
  }
for (i in 1:nrow(merged)) {
  if (!is.na(is_sfari[i]) && is_sfari[i]) {
    merged[i, "SFARI.Gene"] <- "Y"
  }
  }

}

cat("Removing GOI...\n")

merged <- merged[-which(merged$gene == goi), ]


cat("Done.\n")
  return(merged)
}



Idents(obj) <- "celltype"
test = celltype_expression_RNA(obj, "Astrocytes","RFX4", xbox_rare_DNV_asd)
testc = correlation(test, "RFX4", xbox_rare_DNV_asd, "S")




print(result)

test_corr = rcorr(as.matrix(as), type = "spearman")
test_corr_r = as.data.frame(test_corr$r)
test_corr_p = as.data.frame(test_corr$P)

goi_only_p = as.data.frame(test_corr_p["RFX4"], row.names = rownames(test_corr_p))
goi_only_r = as.data.frame(test_corr_r["RFX4"], row.names = rownames(test_corr_r))

goi_only_p$gene = rownames(goi_only_p)
goi_only_r$gene = rownames(goi_only_r)

merged = merge(goi_only_r, goi_only_p, by = "gene")
colnames(merged)[2:3] <- c("R","P")
merged = merged[order(merged$R, decreasing = TRUE), ]
merged$SFARI.Gene = ""

common_genes = intersect(merged$gene, xbox_rare_DNV_asd)

is_sfari <- rep(FALSE, nrow(merged))
is_sfari[merged$gene %in% SFARI_genes$gene.symbol] <- TRUE

cat("Adding gene annotations...\n")
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
for (i in 1:nrow(merged)) {
  gene_id <- merged$gene[i]
  gene_info <- getBM(attributes = c("external_gene_name", "description"), filters = "external_gene_name", values = gene_id, mart = mart)
  if (nrow(gene_info) > 0) {
    merged$description[i] <- gene_info$description[1]
  } else {
    merged$description[i] <- "N/A"
  }
  if (is_sfari[i]) {
    merged[i, "SFARI.Gene"] <- "Y"
  }
}



library(tidyverse)
all_rfx3 = read.csv("/Users/artemis/Downloads/all_genes_corr_RFX3.csv")
all_rfx3 = as_tibble(all_rfx3)


VSMC_all_rfx3 <- all_rfx3 %>%
  filter(R > 0.25, celltype == "VSMC")


Oligodendrocytes_all_rfx3 = all_rfx3 %>%
  filter(R > 0.20 | R<0.20 , celltype == "Oligodendrocytes")

vsmc_gene_list = VSMC_all_rfx3$gene
cat(vsmc_gene_list, sep = " ,")
writeLines(vsmc_gene_list$V1, "vsmc_gene_list.txt")
vsmc_gene_list = readLines("vsmc_gene_list.txt")



writeLines(Oligodendrocytes_all_rfx3$gene, "oligodendrocytes_gene_list.txt")
oligodendrocytes_gene_list = readLines("oligodendrocytes_gene_list.txt")





library(WebGestaltR)




geneFile = "/Users/artemis/vsmc_gene_list.txt"
geneFile  = "/Users/artemis/oligodendrocytes_gene_list.txt"

res <- WebGestaltR(enrichMethod = "NTA",
                   organism = "hsapiens",
                   enrichDatabase = "network_PPI_BIOGRID",
                   enrichDatabaseType ="genesymbol",
                   interestGeneType = "genesymbol",
                   interestGeneFile = geneFile,
                   networkConstructionMethod="Network_Retrieval_Prioritization",
                   isOutput = TRUE)


library(BiocManager)
BiocManager::valid()


  BiocManager::install(c(
    "GenomeInfoDb", "GOfuncR", "MungeSumstats"
  ), update = TRUE, ask = FALSE, force = TRUE)

BiocManager::install("NeuCA")
library(NeuCA)


library(Seurat)
obj = readRDS("obj_new.RDS")
obj.se = as.SingleCellExperiment(obj)
saveRDS(obj.se, "obj_sce.RDS")
obj.se=readRDS("obj_sce.RDS")


obj_dev = readRDS("obj_preprocessed_umap.rds")
obj_dev.se = as.SingleCellExperiment(obj_dev)
saveRDS(obj_dev.se, "obj_dev_sce.RDS")
obj_dev.se=readRDS("obj_dev_sce.RDS")

library(SingleCellExperiment)
library(NeuCA)

true_cell_label = obj.se$celltype
ref_anno = data.frame(true_cell_label, row.names = colnames(obj.se@assays@data$counts))

ref = SingleCellExperiment(
  assays = list(counts = as.matrix(obj.se@assays@data$counts)),
  colData = ref_anno
)
rowData(ref)$feature_symbol =  rownames(ref)
ref =  ref[!duplicated(rownames(ref)),]
saveRDS(ref, "ref_sce.RDS")

predicted.label = NeuCA(train =  ref , test = obj_dev.se, model.size = "small", verbose = TRUE)



library(NeuCA)
library(SingleCellExperiment)
ref=readRDS("/Users/artemis/Downloads/ref_sce.RDS")
obj_dev.se=readRDS("/Users/artemis/Downloads/obj_dev_sce.RDS")
lapply(c("dplyr","Seurat","HGNChelper"), library, character.only = T)

library(Seurat)
obj_new = readRDS("/Users/artemis/Downloads/obj_new.RDS")

# load gene set preparation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
# load cell type annotation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
# DB file
db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx";
tissue = "Brain" 
gs_list = gene_sets_prepare(db_, tissue)
es.max = sctype_score(obj_new@assays$SCT@scale.data, scaled = TRUE, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)

BiocManager::install("EasyCellType")
library(EasyCellType)
library(org.Hs.eg.db)
library(AnnotationDbi)

markers = FindAllMarkers(obj_new, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

markers$entrezid <- mapIds(org.Hs.eg.db,
                           keys=markers$gene, #Column containing Ensembl gene ids
                           column="ENTREZID",
                           keytype="SYMBOL",
                           multiVals="first")
markers <- na.omit(markers)

library(dplyr)
markers_sort <- data.frame(gene=markers$gene, cluster=markers$cluster, 
                      score=markers$avg_log2FC) %>% 
  group_by(cluster) %>% 
  mutate(rank = rank(score),  ties.method = "random") %>% 
  arrange(desc(rank)) 
input.d <- as.data.frame(markers_sort[, 1:3])

annot.GSEA <- easyct(input.d, db="cellmarker", species="Human", 
                    tissue=c("Brain"),
                    test="GSEA", scoretype = "pos")

plot_dot(test= "GSEA", annot.GSEA)

#######################################


library(Seurat)
library("scPred")
library("magrittr")


reference = readRDS("obj_new.RDS")
query = readRDS("dev_processed.rds")


reference = getFeatureSpace(reference, "celltype")
reference <- trainModel(reference)
get_probabilities(reference) %>% head()
pdf("probbilities_reference.pdf")
plot_probabilities(reference)
dev.off()
save(reference, file = "reference.RDS")

query = NormalizeData(query)

query = scPredict(query,reference)

.make_names <- function(x){
  x <- gsub("\\+", "_plus", x)
  x <- gsub("\\-", "_minus", x)
  x <- make.names(x)
}

scPredict_edited <- function (new, reference, threshold = 0.55, max.iter.harmony = 20, 
    recompute_alignment = TRUE, seed = 66) 
{
    if (!(is(reference, "Seurat") | is(reference, "scPred"))) 
        stop("'object' must be of class 'scPred' or 'Seurat'")
    if (is(reference, "Seurat")) {
        spmodel <- reference@misc$scPred
    }
    else {
        spmodel <- reference
    }
    if (is.null(spmodel)) 
        stop("No feature space has been determined!")
    if (!length(spmodel@train)) 
        stop("No models have been trained!")
    if (!is(new, "Seurat")) 
        stop("New data must be a Seurat object")
    new <- project_query_edited(new, reference = spmodel, max.iter.harmony = max.iter.harmony, 
        recompute_alignment = recompute_alignment, seed = seed)
    new_embeddings_aligned <- Embeddings(new[["scpred"]])
    colnames(new_embeddings_aligned) <- colnames(spmodel@cell_embeddings)
    cellTypeModelNames <- names(spmodel@features)
    .predictCellClass <- function(cellType, spmodel, testEmbeddings) {
        features <- as.character(spmodel@features[[cellType]]$feature)
        model <- spmodel@train[[cellType]]
        prediction <- predict(model, newdata = scPred:::subsetMatrix(testEmbeddings, 
            features), type = "prob")
        rownames(prediction) <- rownames(testEmbeddings)
        prediction[, 1, drop = FALSE]
    }
    cat(crayon::green(cli::symbol$record, " Classifying cells...\n"))
    res <- sapply(cellTypeModelNames, .predictCellClass, spmodel, 
        new_embeddings_aligned)
    res <- as.data.frame(res)
    colnames(res) <- cellTypeModelNames
    rownames(res) <- colnames(new)
    classes <- cellTypeModelNames
    if (length(cellTypeModelNames) == 1) {
        metadata <- get_metadata(spmodel)
        cellClasses <- levels(metadata$pvar)
        res_comp <- 1 - res[, 1]
        negClass <- cellClasses[cellClasses != names(res)]
        res[[negClass]] <- res_comp
    }
    max_props <- as.data.frame(t(apply(res, 1, function(x) c(index = which.max(x), 
        max = x[which.max(x)]))))
    names(max_props) <- c("index", "max")
    max_props$generic_class <- names(res)[max_props$index]
    res <- cbind(res, max_props)
    pred <- ifelse(res$max > threshold, res$generic_class, "unassigned")
    names(pred) <- colnames(new)
    res$prediction <- pred
    res$index <- NULL
    res$no_rejection <- res$generic_class
    res$generic_class <- NULL
    names(res) <- .make_names(paste0("scpred_", names(res)))
    new <- AddMetaData(new, res)
    cat(crayon::green("DONE!\n"))
    new
}

project_query_edited <- function (new, reference, max.iter.harmony = 20, recompute_alignment = TRUE, 
    seed = 66, ...) 
{
    if (!(is(reference, "Seurat") | is(reference, "scPred"))) 
        stop("'object' must be of class 'scPred' or 'Seurat'")
    if (is(reference, "Seurat")) {
        spmodel <- reference@misc$scPred
    }
    else {
        spmodel <- reference
    }
    if (is.null(spmodel)) 
        stop("No feature space has been determined!")
    if (!is(new, "Seurat")) 
        stop("New data must be a Seurat object")
    if ("scpred" %in% names(new@reductions)) {
        if (recompute_alignment) {
            alignment <- TRUE
            cat(crayon::yellow(cli::symbol$figure_dash, "Data has already been aligned to a reference.\n"), 
                sep = "")
            cat(crayon::yellow(cli::symbol$sup_plus, "Skip data alignment using `recompute.alignment = FALSE`.\n"), 
                sep = "")
        }
        else {
            alignment <- FALSE
        }
    }
    else {
        alignment <- TRUE
    }
    if (alignment) {
        cat(crayon::green(cli::symbol$record, " Matching reference with new dataset...\n"))
        ref_loadings <- spmodel@feature_loadings
        ref_embeddings <- spmodel@cell_embeddings
        new_features <- rownames(new)
        reference_features <- rownames(ref_loadings)
        shared_features <- intersect(reference_features, new_features)
        cat(crayon::cyan("\t", cli::symbol$line, paste(length(reference_features), 
            "features present in reference loadings\n")))
        cat(crayon::cyan("\t", cli::symbol$line, paste(length(shared_features), 
            "features shared between reference and new dataset\n")))
        cat(crayon::cyan("\t", cli::symbol$line, paste0(round(length(shared_features)/length(reference_features) * 
            100, 2), "% of features in the reference are present in new dataset\n")))
        ref_loadings <- ref_loadings[shared_features, ]
        new_data <- GetAssayData(new, layer = "data")[shared_features, 
            ]
        means <- spmodel@scaling$means
        stdevs <- spmodel@scaling$stdevs
        new_data <- Matrix::t(new_data)
        names(means) <- names(stdevs) <- rownames(spmodel@scaling)
        means <- means[shared_features]
        stdevs <- stdevs[shared_features]
        i <- stdevs == 0
        if (any(i)) {
            warning(paste0(sum(i), " features have zero variance but are present in the feature loadings. \nDid you subset or integrated this data before?"))
            cat(crayon::yellow("Removing zero-variance genes from projection\n"))
            new_data <- new_data[, !i]
            ref_loadings <- ref_loadings[!i, ]
            means <- means[!i]
            stdevs <- stdevs[!i]
        }
        scaled_data <- scale(new_data, means, stdevs)
        new_embeddings <- scaled_data %*% ref_loadings
        dataset <- factor(c(rep("reference", nrow(ref_embeddings)), 
            rep("new", nrow(new_embeddings))), levels = c("reference", 
            "new"))
        rownames(ref_embeddings) <- paste0("ref_", rownames(ref_embeddings))
        rownames(new_embeddings) <- paste0("new_", rownames(new_embeddings))
        eigenspace <- as.data.frame(rbind(ref_embeddings, new_embeddings))
        meta_data <- data.frame(rownames(eigenspace), dataset = dataset)
        cat(crayon::green(cli::symbol$record, " Aligning new data to reference...\n"))
        set.seed(seed)
        harmony_embeddings <- harmony::HarmonyMatrix(eigenspace, meta_data, 
            "dataset", do_pca = FALSE, reference_values = "reference", 
            max.iter.harmony = max.iter.harmony, ...)
        new_embeddings_aligned <- harmony_embeddings[dataset == 
            "new", , drop = FALSE]
    }
    else {
        new_embeddings_aligned <- Embeddings(new, reduction = "scpred")
        colnames(new_embeddings_aligned) <- gsub("scpred_", spmodel@reduction_key, 
            colnames(new_embeddings_aligned))
    }
    rownames(new_embeddings_aligned) <- gsub("^new_", "", rownames(new_embeddings_aligned))
    new@reductions[["scpred"]] <- CreateDimReducObject(embeddings = new_embeddings_aligned, 
        key = "scpred_", assay = DefaultAssay(object = new))
    if (recompute_alignment) {
        rownames(new_embeddings) <- gsub("^new_", "", rownames(new_embeddings))
        new@reductions[["scpred_projection"]] <- CreateDimReducObject(embeddings = new_embeddings, 
            key = "Projection_", assay = DefaultAssay(object = new))
    }
    new
}

new = project_query_edited(query, reference, threshold = 0.2)
pdf("pred.pdf")
new@active.assay = "scpred"
DimPlot(new, group.by = "scpred_prediction", reduction = "scpred")
dev.off()



###### switch reference and query

reference = readRDS("dev_processed.rds")
query = readRDS("obj_new.RDS")

reference = getFeatureSpace(reference, "lineage")
reference <- trainModel(reference)
get_probabilities(reference) %>% head()
pdf("probbilities_reference2.pdf")
plot_probabilities(reference)
dev.off()





########

