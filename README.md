# CorrelationGene 0.1.5

CorrelationGene is a small R package that takes any Seurat object, gene of interest and gene list and calculates the Pearson correlations between the identities specified by the user. The output has each SFARI gene annotated, (https://gene.sfari.org/database/human-gene/).

Please raise any issues through Github or email me at timothy.james.allen@ki.se for questions/issues/contributions! 

# Future functions

- Take multiple genes of interests

- Further statistical functions on the correlations i.e. T test, Shapiro-test

- Statistical tests for differences in correlation matrices

- Make a vignette!

- Fix SFARI gene annotation ✓ (NEW!! in CorrelationGene 0.1.1)

- Add annotation for other NDDs 

- Add dependencies to the actual package

# Dependencies 

- Seurat
- Matrix
- tictoc
- Hmisc
- moments 
- org.Hs.eg.db
- bestNormalize

# Installation

>install.packages(devtools)

>library(devtools)

>devtools::install_github(CorrelationGene)

>library(CorrelationGene)

>library(Seurat)

>library(Matrix)

>library(tictoc)

>library(Hmisc)

>library(moments)

>library(org.Hs.eg.db)

>library("bestNormalize")


# Changelog
0.1.5 - Added scoring to the results that is based on a sigmoidal model that incorporates the R value, padj value, the maximum amount of cells in each celltype and orderNorm for normalization.

0.1.4 - Pipeline now suppports ensembl genes in the Seurat object, where you can specify if ensembl in the parameters. 
        The p.adj metric is now added in the result which uses the BH method.
        
        Removed biomaRt in favour for org.Hs.eg.db which adds annotations faster and more reliably.
        
0.1.3 - Gene descriptions are added using biomaRt automatically in the correlation() function. 

0.1.2 - You can now use spearman test in correlation() and express_cell_SCT/express_cell_RNA(). 

Example for 0.1.2:

>Spearmantest = express_cell_SCT(obj, "RFX3", xbox_genes$Gene.Name, "S")
