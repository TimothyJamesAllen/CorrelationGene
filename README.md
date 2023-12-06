# CorrelationGene 0.1.0

CorrelationGene is a small R package that takes any Seurat object, gene of interest and gene list and calculates the Pearson correlations between the identities specified by the user. The output has each SFARI gene annotated, (https://gene.sfari.org/database/human-gene/).

Please raise any issues through Github or email me at timothy.james.allen@ki.se for questions/issues/contributions! 

# Future functions

- Take multiple genes of interests

- Further statistical functions on the correlations i.e. T test, Shapiro-test, FDR

- Statistical tests for differences in correlation matrices

- Make a vignette!

- Fix SFARI gene annotation

- Add annotation for other NDDs 

- Add dependencies to the actual package

# Dependencies 

- Seurat
- Matrix
- tictoc
- Hmisc
- moments 
- org.Hs.eg.db
