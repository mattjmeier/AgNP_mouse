source("https://bioconductor.org/biocLite.R")
biocLite("biomformat")
biocLite("phyloseq")
biocLite("edgeR")
biocLite("DESeq2")
library("biomformat")
library("phyloseq")
library("scales")
library("edgeR")
library("DESeq2")
library('ggplot2')
library('dplyr')

### ANALYSIS FOR SHOTGUN DATA ANNOTATED BY MG-RAST

## Go to working directory
setwd("./")


# Sample data...
map <- read.table("./map.txt", header=T, sep="\t", row.names=1)

map$Group <- factor(map$Group, levels = c("Control","Antibiotics","AgNP 1", "AgNP 5", "Antbiotics plus AgNP 1", "Antbiotics plus AgNP 5"))
myLevels <- levels(map$Group)
# Limit based on R1
# map <- map[grep(".*R1", row.names(map)),]
dim(map)


# BIOM format from MG-RAST needs to have "type" changed to a compatible term such as "OTU table"
function.subsystems.biom <- read_biom("./All samples shotgun Subsystems function.biom")
function.subsystems.matrix <- as.matrix(biom_data(function.subsystems.biom))
dim(function.subsystems.matrix)
function.subsystems.matrix <- t(function.subsystems.matrix)
dim(function.subsystems.matrix)
#keep <- rowSums(function.subsystems.matrix) > 100000
#function.subsystems.matrix <- function.subsystems.matrix[keep,]
dim(function.subsystems.matrix)
colnames(function.subsystems.matrix) <- paste(observation_metadata(function.subsystems.biom)[[1]],
                                              ":",
                                              observation_metadata(function.subsystems.biom)[[2]],
                                              ":",
                                              observation_metadata(function.subsystems.biom)[[3]],
                                              ":",
                                              observation_metadata(function.subsystems.biom)[[4]])
# Limit to samples in map
function.subsystems.matrix <- function.subsystems.matrix[row.names(function.subsystems.matrix) %in% row.names(map),]
dim(function.subsystems.matrix)
my.matrix <- function.subsystems.matrix

function.KO.biom <- read_biom("./All samples shotgun KO function.biom")
function.KO.matrix <- as.matrix(biom_data(function.KO.biom))
dim(function.KO.matrix)
function.KO.matrix <- t(function.KO.matrix)
dim(function.KO.matrix)
keep <- rowSums(function.KO.matrix) > 100000
function.KO.matrix <- function.KO.matrix[keep,]
dim(function.KO.matrix)
colnames(function.KO.matrix) <- paste(observation_metadata(function.KO.biom)[[1]],
                                              ":",
                                              observation_metadata(function.KO.biom)[[2]],
                                              ":",
                                              observation_metadata(function.KO.biom)[[3]],
                                              ":",
                                              observation_metadata(function.KO.biom)[[4]])
# Limit to samples in map
function.KO.matrix <- function.KO.matrix[row.names(function.KO.matrix) %in% row.names(map),]
dim(function.KO.matrix)
my.matrix <- function.KO.matrix



# Limit based on matrix design
map <- map[row.names(map) %in% row.names(my.matrix), ]
dim(map)

namesOfAllSubsystems <- paste(observation_metadata(function.subsystems.biom)[[1]],
                              ":",
                              observation_metadata(function.subsystems.biom)[[2]],
                              ":",
                              observation_metadata(function.subsystems.biom)[[3]],
                              ":",
                              observation_metadata(function.subsystems.biom)[[4]])





#  EdgeR
# 
# for (coef in 2:length(levels(map$Group))) {
#   
# QL.results.function <- glmQLF.edgeR(x=map$Group, Y=my.matrix, coef=coef)
# # topTags(QL.results.function)
# write.table( topTags(QL.results.function, n=Inf)$table,
#             file=paste0(levels(map$Group)[coef],".edgeR_results.functions_subsystems.txt"),
#             sep="\t",
#             quote=F,
#             col.names=NA)
# }


# Get into Phyloseq
otumatrix = as(biom_data(function.subsystems.biom), "matrix")
OTU = otu_table(otumatrix, taxa_are_rows=TRUE)
taxmat = as.matrix(observation_metadata(function.subsystems.biom), rownames.force=TRUE)
TAX = tax_table(taxmat)
physeq = phyloseq(OTU, TAX, sample_data(map))
taxdf <- as.data.frame(tax_table(physeq))
taxdf$taxid <- row.names(taxdf)
ps.transformed <- transform_sample_counts(physeq, function(OTU) OTU/sum(OTU))

## Ordination
ps.ord <- ordinate(ps.transformed, "NMDS", "bray")
p1 = plot_ordination(ps.transformed, ps.ord, type="sample", shape="Antibiotics", color="Group", title="AgNP Mouse project: NMDS, Bray, samples plot")
p1 + geom_point(size=5)

p2 = plot_ordination(ps.transformed, ps.ord, type="taxa", color="functionalHierarchy1", title="AgNP Mouse project: NMDS, Bray, genes plot")
p2 + geom_point(size=0.5)

ps.ord <- ordinate(ps.transformed, "PCoA", "bray")
p3 = plot_ordination(ps.transformed, ps.ord, type="sample", shape="Antibiotics", color="Group", title="AgNP Mouse project: PCoA, Bray, samples plot")
p3 + geom_point(size=5)

metalresistance <-subset_taxa(ps.transformed, functionalHierarchy3=="Cobalt-zinc-cadmium_resistance")
ps.ord <- ordinate(metalresistance, "PCoA", "bray")
p4 = plot_ordination(metalresistance, ps.ord, type="sample", shape="Antibiotics", color="Group", title="AgNP Mouse project: PCoA, Bray, samples plot")
p4 + geom_point(size=5)


oxidativestress <-subset_taxa(ps.transformed, functionalHierarchy2=="Oxidative stress")
ps.ord <- ordinate(oxidativestress, "PCoA", "bray")
p5 = plot_ordination(oxidativestress, ps.ord, type="sample", shape="Antibiotics", color="Group", title="AgNP Mouse project: PCoA, Bray, samples plot")
p5 + geom_point(size=5)  +  geom_text(mapping = aes(label = Sample), size = 4)


resistance_anti_and_tox <- subset_taxa(ps.transformed, functionalHierarchy2=="Resistance to antibiotics and toxic compounds")
ABCtransporters <- subset_taxa(ps.transformed, functionalHierarchy2=="ABC transporters")
Transport_of_Nickel_and_Cobalt  <- subset_taxa(ps.transformed, functionalHierarchy3=="Transport_of_Nickel_and_Cobalt")
Multidrug_Resistance_Efflux_Pumps <- subset_taxa(ps.transformed, functionalHierarchy3=="Multidrug_Resistance_Efflux_Pumps")

plot_heatmap(metalresistance,
             sample.order=row.names(map),
             taxa.label="functionalHierarchy4",
             sample.label="Group",
             trans=log_trans(3))

plot_heatmap(resistance_anti_and_tox,
             sample.order=row.names(map),
             taxa.label="functionalHierarchy4",
             sample.label="Group",
             trans=log_trans(2))

plot_heatmap(resistance_anti_and_tox,
             method="PCoA",
             distance="bray",
             first.sample="Abio-1_R2",
             taxa.label="functionalHierarchy4",
             sample.label="Sample",
             trans=log_trans(2))

plot_heatmap(oxidativestress,
             sample.order=row.names(map),
             taxa.label="functionalHierarchy4",
             sample.label="Group",
             trans=log_trans(2))

plot_heatmap(ABCtransporters,
             sample.order=row.names(map),
             taxa.label="functionalHierarchy4",
             sample.label="Group",
             trans=log_trans(2))

plot_heatmap(Transport_of_Nickel_and_Cobalt,
             sample.order=row.names(map),
             taxa.label="functionalHierarchy4",
             sample.label="Group",
             trans=log_trans(2))

plot_heatmap(Multidrug_Resistance_Efflux_Pumps,
             sample.order=row.names(map),
             taxa.label="functionalHierarchy4",
             sample.label="Group",
             trans=log_trans(2))

plot_heatmap(ps.transformed,
             sample.order=row.names(map),
             taxa.label="functionalHierarchy4",
             sample.label="Group",
             trans=log_trans(2))


searchterm=c("nitro")
searched_names <- subset(taxdf, grepl(searchterm, taxdf$functionalHierarchy1, ignore.case=T) | grepl(searchterm, taxdf$functionalHierarchy2, ignore.case=T) | grepl(searchterm, taxdf$functionalHierarchy3, ignore.case=T)| grepl(searchterm, taxdf$functionalHierarchy4, ignore.case=T))
searched_tax_ids <- row.names(searched_names)
physeq_searched=prune_taxa(searched_tax_ids, ps_transformed)
physeq_searched_glom <- tax_glom(physeq_searched, taxrank="functionalHierarchy4")
plot_bar(physeq_searched_glom, x="Silver_concentration", fill="functionalHierarchy4")

siggenes <-subset_taxa(physeq, functionalHierarchy4=="Copper-translocating P-type ATPase (EC 3.6.3.4)" | functionalHierarchy4=="UDP-glucose 4-epimerase (EC 5.1.3.2)" | functionalHierarchy4=="Cation efflux system protein CusA" | functionalHierarchy4=="Cobalt-zinc-cadmium resistance protein CzcA" | functionalHierarchy4=="Cobalt-zinc-cadmium resistance protein" | functionalHierarchy4=="Coenzyme PQQ synthesis protein A" | functionalHierarchy4=="GTP cyclohydrolase I (EC 3.5.4.16) type 2" | functionalHierarchy4=="Dihydropteroate synthase (EC 2.5.1.15)" | functionalHierarchy4=="FIG004453: protein YceG like" | functionalHierarchy4=="Phosphocarrier protein, nitrogen regulation associated" | functionalHierarchy4=="Lipoprotein NlpD" | functionalHierarchy4=="Xylose isomerase (EC 5.3.1.5)" | functionalHierarchy4=="Chemotaxis regulator - transmits chemoreceptor signals to flagelllar motor components CheY" | functionalHierarchy4=="Dihydroorotate dehydrogenase (EC 1.3.3.1)" | functionalHierarchy4=="Serine acetyltransferase (EC 2.3.1.30)" | functionalHierarchy4=="Membrane-bound lytic murein transglycosylase B precursor (EC 3.2.1.-)" | functionalHierarchy4=="OpgC protein" | functionalHierarchy4=="Superoxide dismutase [Fe] (EC 1.15.1.1)" | functionalHierarchy4=="3-hydroxyisobutyrate dehydrogenase (EC 1.1.1.31)")

metal <-subset_taxa(ps_transformed, functionalHierarchy4=="Copper-translocating P-type ATPase (EC 3.6.3.4)" | functionalHierarchy4=="Cation efflux system protein CusA" | functionalHierarchy4=="Cobalt-zinc-cadmium resistance protein CzcA" | functionalHierarchy4=="Cobalt-zinc-cadmium resistance protein" | functionalHierarchy4=="Cobalt/zinc/cadmium efflux RND transporter, membrane fusion protein, CzcB family")


plot_bar(oxi, x="Silver_concentration", fill="functionalHierarchy4")

oxi<-subset_taxa(ps_transformed, functionalHierarchy4=="Coenzyme PQQ synthesis protein A" | functionalHierarchy4=="Peptide methionine sulfoxide reductase MsrB (EC 1.8.4.12)" | functionalHierarchy4=="4-hydroxyphenylpyruvate dioxygenase (EC 1.13.11.27)" | functionalHierarchy4=="Ferric uptake regulation protein FUR" | functionalHierarchy4=="Superoxide dismutase [Fe] (EC 1.15.1.1)" | functionalHierarchy4=="Peptide methionine sulfoxide reductase MsrA (EC 1.8.4.11)")

plot_bar(oxi, x="Silver_concentration", fill="functionalHierarchy4")

# For facets (divide into grids:)
plot_bar(oxi, x="Silver_concentration", fill="functionalHierarchy4") + facet_grid( ~ functionalHierarchy1  )

mem<-subset_taxa(ps_transformed, functionalHierarchy4=="Outer membrane lipoprotein carrier protein LolA" | functionalHierarchy4=="Flagellar hook-associated protein FliD" | functionalHierarchy4=="Flagellar basal-body rod protein FlgF" | functionalHierarchy4=="Apolipoprotein N-acyltransferase (EC 2.3.1.-)" | functionalHierarchy4=="Phosphocarrier protein, nitrogen regulation associated" | functionalHierarchy4=="Lipoprotein NlpD" | functionalHierarchy4=="Chemotaxis regulator - transmits chemoreceptor signals to flagelllar motor components CheY" | functionalHierarchy4=="Serine acetyltransferase (EC 2.3.1.30)" | functionalHierarchy4==" Membrane-bound lytic murein transglycosylase B precursor (EC 3.2.1.-)")

mis<-subset_taxa(ps_transformed, functionalHierarchy4=="Branched-chain amino acid aminotransferase (EC 2.6.1.42)" | functionalHierarchy4=="Methylglutaconyl-CoA hydratase (EC 4.2.1.18)" | functionalHierarchy4=="Glycolate dehydrogenase (EC 1.1.99.14), iron-sulfur subunit GlcF" | functionalHierarchy4=="UDP-glucose 4-epimerase (EC 5.1.3.2)" | functionalHierarchy4=="Xylose isomerase (EC 5.3.1.5)" | functionalHierarchy4=="protein YceG like" | functionalHierarchy4=="Cell division protein FtsK" | functionalHierarchy4=="Transmembrane regulator protein PrtR" | functionalHierarchy4=="Dihydropteroate synthase (EC 2.5.1.15)" | functionalHierarchy4=="GTP cyclohydrolase I (EC 3.5.4.16) type 2" | functionalHierarchy4=="Enoyl-CoA hydratase (EC 4.2.1.17)" | functionalHierarchy4=="Glutamate 5-kinase (EC 2.7.2.11)" | functionalHierarchy4=="Dihydroorotate dehydrogenase (EC 1.3.3.1)" | functionalHierarchy4=="OpgC protein" | functionalHierarchy4=="Alkanesulfonates-binding protein")









mainFactor="Group"

ps<-physeq
dds = phyloseq_to_deseq2(ps, as.formula(paste0("~",mainFactor)))

# Change factor levels ### CUSTOMIZE!!
# This shouldn't be necessary with new code added at start.
# colData(dds)[[mainFactor]] <- factor(colData(dds)[[mainFactor]], levels=myLevels)

# Geometric mean function
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

# Calculate geometric means
geoMeans = apply(counts(dds), 1, gm_mean)

# Run DESeq2
dds = estimateSizeFactors(dds)
dds = estimateSizeFactors(dds, geoMeans = geoMeans)

dds = DESeq(dds, fitType="local")
res = results(dds)
DESeq2::plotMA(dds)
# colData(dds)

# Create list of DESeq results for each contrast (factor)...
resList <- list()

for (i in 2:length(myLevels)) {
  message("Contrast index:")
  print(i)
  message("Contrast factor:")
  print(myLevels[i])
  message("List index:")
  print(i-1)
  message("Computing results...")
  ### This may need customization!
  ### Defaults to use mainFactor as the variable for contrast,
  ### Then the first level of that factor as control
  ### And each successive level as the contrasts
  myContrast <- c(mainFactor,
                  myLevels[i],
                  myLevels[1])
  print(myContrast)
  myResults <- results(dds,
                       contrast=myContrast,
                       cooksCutoff = FALSE)
  message("Done computing results. Adding to list.")
  resList[[i-1]] <- myResults
  message("Done adding to list.")
}


# Filter results table using adjusted p-value of alpha...
alpha = 0.01
sigtabList <- list()

for (i in 1:length(resList)) {
  print(i)
  sigTab <- resList[[i]]
  # Add taxonomy
  if (nrow(sigTab) == 0) {
    next
  } else {
    sigTab <- cbind(as(sigTab, "data.frame"),
                    as(tax_table(ps)[rownames(sigTab), ], "matrix"),
                    contrast=myLevels[[i+1]])
    sigTab <- sigTab[!is.na(sigTab$padj) & sigTab$padj < alpha, ]
    sigtabList[[i]] <- sigTab
  }
}
sigtabList <- sigtabList[!sapply(sigtabList, is.null)] 




summaryTable <- data.frame( baseMean=resList[[1]]$baseMean )
for (i in 1:length(resList)) {
  print(i)
  #print(myLevels[[i+1]])
  n <- myLevels[[i+1]]
  n <- paste0("log2FoldChange ",n)
  p <- paste0("padj ",n)
  message(n)
  summaryTable <- cbind(summaryTable, log2FoldChange=resList[[i]]$log2FoldChange, padj=resList[[i]]$padj)
  names(summaryTable)[[ncol(summaryTable)-1]] <- n
  names(summaryTable)[[ncol(summaryTable)]] <- p
}
summaryTable <- cbind(summaryTable, as(tax_table(ps)[rownames(summaryTable), ], "matrix"))





##############
# Plot results
theme_set(theme_bw())
#scale_fill_discrete <- function(palname = "Set1", ...) {
#  scale_fill_brewer(palette = palname, ...)
#}


# Write dataframe of all results
significantResults <- do.call(rbind, sigtabList)

#######################################
### Write results table from DESeq2
#######################################
write.table(significantResults, file="DESeq_output_significant.txt", quote=F, sep='\t', col.names=NA)
write.table(summaryTable, file="DESeq_output_all_genes.txt", quote=F, sep='\t', col.names=NA)
#######################################

# Create list of results tables where NA genera are removed
sigtableListNA.RM <-lapply(sigtabList, subset, !is.na(Genus))

# Create list of combined results
combined_sigtableListNA.RM <- do.call(rbind, sigtableListNA.RM)

attach(significantResults)
significantResultsOrdered <- significantResults[order(padj),]
significantResultsOrderedTop20 <- head(significantResultsOrdered, n=20)

topResults <- NULL
numResults=50
for (i in 2:length(levels(map$Group))) {
  topResults <- rbind(topResults,
                head(significantResultsOrdered %>% 
                dplyr::filter(contrast==levels(map$Group)[[i]]),
                n=numResults))
  }


## Genus plots, no N/A
ggplot(topResults,
       aes(x=functionalHierarchy4,
           y=log2FoldChange,
           color=functionalHierarchy4,
           size=padj)) +
  geom_point() +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) +
  facet_grid(~contrast, scales="free_x") +
  ggtitle(paste0("Genes with significantly different relative abundance, grouped by treatment")) +
  theme(legend.position = "none")


# Examine specific functions...
myPS<-ps
searchterm="Predicted cobalt transporter in Bacteroides_Porphyromonas"
num <- grep(searchterm, myPS@tax_table@.Data[,4])
myPS@otu_table@.Data[num,]
dds@assays[[1]][num,]
