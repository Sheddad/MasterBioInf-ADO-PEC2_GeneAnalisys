# Carreguem la llibreria
require(GEOquery)
# Amb l'identificador de l'estudi descarreguem les dades de GEO
gse <- getGEO("GSE2401", GSEMatrix=TRUE, AnnotGPL=TRUE)
esetFromGEO <- gse[[1]]
head(exprs(esetFromGEO))

gds <- getGEO("GDS1251", GSEMatrix=TRUE, AnnotGPL=TRUE)

# Boxplot
boxplot(filteredEset, which="all",las=2, main="Intensity distribution of RAW data", 
        cex.axis=0.6, col=sampleColor, names=sampleNames)

# HIERARQUICAL CLUSTERING
clust.euclid.average <- hclust(dist(t(exprs(filteredEset))),method="average")
plot(clust.euclid.average, labels=sampleNames, main="Hierarchical clustering of RawData", cex=0.7,  hang=-1)

clust.euclid.average <- hclust(dist(t(exprs(esetFromGEO))),method="average")
plot(clust.euclid.average, labels=sampleNames, main="Hierarchical clustering of Eset from GEO", cex=0.7,  hang=-1)


library(limma)
treat <- pData(esetFromGEO)$description
lev <- factor(treat, levels = unique(treat))
design <-model.matrix(~0+lev)
colnames(design) <- levels(lev)
rownames(design) <- pData(esetFromGEO)$geo_accession
print(design)

#COMPARISON
cont.matrix1 <- makeContrasts( 
  Hypotension.vs.Normal = Hypotesion-Normal,
  levels = design)
comparisonName <- "Efecte de la hipotensió induida"
print(cont.matrix1)


require("rae230a.db")  # Asegúrate de cambiar esto al paquete de anotación que estás usando
keys <- keys(rae230a.db, keytype = "ENTREZID")
head(keys)
tail(keys)

# Verificar si el identificador del gen está en la lista
identificador_a_buscar <- "1399148_s_at"
identificador_a_buscar %in% keys

require(org.Rn.eg.db)
keytypes(org.Rn.eg.db)
anotaciones<- AnnotationDbi::select (org.Rn.eg.db, keys=rownames(filteredData), columns=c( "SYMBOL"))


eset_metadata <- metadata(eset)
annotation(eset)

AnnotationDbi::select (rae230a.db, "1399148_s_at", "SYMBOL", "ENTREZID")





