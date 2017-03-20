library(RnBeads)
library(sva)
library(biomaRt)
library(Metrics)
library(PMCMR)
#String operations
library(stringr)

#Plotting
library(ggplot2)
library(extrafont)
library(Cairo)
library(heatmap.plus)
library(gplots) #for heatmap.2
library(RColorBrewer)
library(UpSetR)

#Reading and writing tables
library(readr)
library(openxlsx)

#library(matrixStats)
library(R.utils)

#Data arrangement
library(reshape2)
library(dplyr)
library(parallel)
library(data.table)

#Functional programming
library(magrittr)
library(purrr)
library(tidyr)
library(broom)
library(WGCNA)

rnb.execute.na.removal.fix <- function (rnb.set, threshold = rnb.getOption("filtering.missing.value.quantile")) {
    if (!inherits(rnb.set, "RnBSet")) {
        stop("invalid value for rnb.set")
    }
    if (!(is.double(threshold) && length(threshold) == 1 && (!is.na(threshold)))) {
        stop("invalid value for threshold")
    }
    if (!(0 <= threshold && threshold <= 1)) {
        stop("invalid value for threshold; expected a value between 0 and 1")
    }
    filterRes <- RnBeads:::rnb.execute.na.removal.internal(rnb.set, NULL,
        threshold)
    list(dataset.before = rnb.set, dataset = remove.sites(rnb.set,
        filterRes$filtered), filtered = filterRes$filtered, threshold = threshold,
        naCounts = filterRes$naCounts)
}

GetSizes <- function(dataset, pval.column, log.column, p.val) {
    dataset.sig <- filter_(dataset, str_c(pval.column, " < ", p.val))
    dataset.up <- filter_(dataset.sig, str_c(log.column, " > 0"))
    dataset.down <- filter_(dataset.sig, str_c(log.column, " < 0"))
    return(list(positive = dim(dataset.up)[1], negative = -(dim(dataset.down)[1])))
}

DecidePlot <- function(file.name, decide.plot, y.lab, bar.padding = 300, pos.hjust = -0.4, neg.hjust = 1.3) {
    y.max <- max(decide.plot$Num.Sites) + nchar(max(decide.plot$Num.Sites)) * bar.padding
    y.min = min(decide.plot$Num.Sites) - nchar(abs(min(decide.plot$Num.Sites))) * bar.padding

    p <- ggplot()
    p <- p + geom_bar(data = filter(decide.plot, Direction == "positive"),  aes(x = Comparison, y = Num.Sites), stat = "identity", colour = "black", fill = "red", position = "dodge")   
    p <- p + geom_text(data = filter(decide.plot, Direction == "positive"), stat = "identity", size = 4, aes(x = Comparison, y = Num.Sites, hjust = pos.hjust, label = Num.Sites), position = position_dodge(width = 1))
    p <- p + geom_bar(data = filter(decide.plot, Direction == "negative"),  aes(x = Comparison, y = Num.Sites), stat = "identity", colour = "black", fill = "green", position = "dodge") 
    p <- p + geom_text(data = filter(decide.plot, Direction == "negative"), stat = "identity", size = 4, aes(x = Comparison, y = Num.Sites, hjust = neg.hjust, label = abs(Num.Sites)), position = position_dodge(width = 1))
    p <- p + theme_bw() + coord_flip() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
    p <- p + theme(axis.title.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_text(hjust = 0), axis.text.x = element_blank(), axis.ticks.x = element_blank()) + ylab(y.lab)
    p <- p + theme(plot.background = element_blank()) + ylim(y.min, y.max) + facet_grid(Test + Total.Sites ~ .) 
    CairoPDF(file.name, width = 7, height = 8, bg = "transparent")
    print(p)
    dev.off()
}

Heatmap <- function(meths, ngenes, file.name, cluster.object) {
    sites.ordered <-  apply(meths, 1, mad) %>% order(decreasing = TRUE)
    sites.plot <- meths[sites.ordered[1:ngenes],]
    cluster.tree <- cluster.object[[7]]@result
    attr(cluster.tree, "class") <- "hclust"
    cluster.dendro <- as.dendrogram(cluster.tree)

    CairoPDF(file.name, width = 10, height = 10)
    heatmap.2(sites.plot, Rowv = TRUE, Colv = cluster.dendro, dendrogram = "both", scale = "none", trace = "none", labRow = FALSE)
    dev.off()
}

UniqueNames <- function(data.list) {
    comparison.name <- names(data.list)
    df.extract <- data.list[[1]]
    colnames(df.extract)[5] <- str_c(comparison.name, ".log2FC")
    colnames(df.extract)[6] <- str_c(comparison.name, ".p.val")
    colnames(df.extract)[7] <- str_c(comparison.name, ".p.val.adj")
    colnames(df.extract)[8] <- str_c(comparison.name, ".combinedRank")
    list(df.extract)
}

DMWorkbook <- function(dataset, filename) {
    pval.cols <- colnames(dataset) %>% str_detect("p.val$") %>% which
    adj.pval.cols <- colnames(dataset) %>% str_detect("p.val.adj") %>% which
    logfc.cols <- colnames(dataset) %>% str_detect("log2FC") %>% which
    description.cols <- colnames(dataset) %>% str_detect("Description") %>% which

    wb <- createWorkbook()
    addWorksheet(wb = wb, sheetName = "Sheet 1", gridLines = TRUE)
    writeDataTable(wb = wb, sheet = 1, x = dataset)
    sig.pvalues <- createStyle(fontColour = "red")
    conditionalFormatting(wb, 1, cols = pval.cols, rows = 1:nrow(dataset), rule = "<0.005", style = sig.pvalues)
    conditionalFormatting(wb, 1, cols = adj.pval.cols, rows = 1:nrow(dataset), rule = "<0.05", style = sig.pvalues)
    conditionalFormatting(wb, 1, cols = logfc.cols, rows = 1:nrow(dataset), style = c("#63BE7B", "white", "red"), type = "colourScale")
    setColWidths(wb, 1, cols = 1, widths = 20)
    setColWidths(wb, 1, cols = 2, widths = "auto")
    setColWidths(wb, 1, cols = 4:ncol(dataset), widths = "auto")
    setColWidths(wb, 1, cols = description.cols, widths = 45)
    pageSetup(wb, 1, orientation = "landscape", fitToWidth = TRUE)
    freezePane(wb, 1, firstRow = TRUE)
    showGridLines(wb, 1, showGridLines = TRUE)

    saveWorkbook(wb, filename, overwrite = TRUE) 
}

GetEnrichr <- function(comparison, submit.df, cutoff, enrichr.terms, site.type) {
    comparison.pval <- str_c(comparison, "p.val.adj", sep = ".")
    filter.df <- filter_(submit.df, str_c(comparison.pval, "<", cutoff)) %>% slice(1:1000)
    print(dim(filter.df))
    enrichr.data <- map(enrichr.terms, GetEnrichrData, filter.df, FALSE)
    enrichr.names <- enrichr.terms[!is.na(enrichr.data)]
    enrichr.data <- enrichr.data[!is.na(enrichr.data)]
    names(enrichr.data) <- enrichr.names
    trap1 <- map(names(enrichr.data), EnrichrWorkbook, enrichr.data, comparison, site.type)
}

FilterEnrichr <- function(enrichr.df, size = 100) {
    enrichr.df$Num.Genes <- map(enrichr.df$Genes, str_split, ",") %>% map(extract2, 1) %>% map_int(length)
    enrichr.filter <- filter(enrichr.df, Num.Genes > 4) %>% filter(P.value < 0.05)
    if (nrow(enrichr.df) > size) {
        enrichr.filter %<>% slice(1:size)
    }
    enrichr.filter
}

VolcanoPlot <- function(top.table, filename, maintitle, cutoff = 0.05, cutoff.column = "P.Value", log.column = "logFC", xlabel = "Log Fold Change", ylabel = "Log.Pvalue") {
    top.table$Significant <- factor(top.table[[cutoff.column]] < cutoff)
    top.table$Log.Pvalue <- -log10(top.table[[cutoff.column]])
    p <- ggplot(top.table, aes_string(x = log.column, y = "Log.Pvalue")) + geom_point(aes(color = Significant))
    p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    p <- p + theme(plot.title = element_text(hjust = 0.5)) + ggtitle(maintitle)
    p <- p + theme(legend.position = "none", plot.background = element_blank())
    p <- p + theme(panel.border = element_rect(color = "black", size = 1))
    p <- p + xlab(xlabel) + ylab(ylabel)
    CairoPDF(filename, width = 6, height = 6, bg = "transparent")
    print(p)
    dev.off()
}

#temporarily changed plot titles to prefix instead of full comparison
MapVolcanoPlots <- function(comparison, prefix, tab.list) {
    tab.reduce <- select(as.data.frame(tab.list[[comparison]]), mean.mean.quot.log2, comb.p.adj.fdr)
    comparison.title <- str_replace_all(comparison, "_", " ") %>% str_replace_all("m\\.", "m ") %>% str_replace_all("n\\.", "n ")
    colnames(tab.reduce) <- c("logFC", "P.Value")
    comparison.format <- str_replace(comparison, "/", "")
    VolcanoPlot(tab.reduce, str_c(prefix, comparison.format, sep = "_"), capitalize(prefix))
}

Top5Plot <- function(contrast, diffmeth.table, meth.matrix, pheno.df, prefix) {
    rank.column <- str_c(contrast, "combinedRank", sep = ".")
    top5.ensembl <- arrange_(diffmeth.table, rank.column)$Ensembl.ID[1:5]
    top5.symbol <- arrange_(diffmeth.table, rank.column)$Symbol[1:5]
    top5.beta <- t(meth.matrix[top5.ensembl,])
    colnames(top5.beta) <- top5.symbol
    top5.df <- data.frame(Combined = as.character(pheno.df$Combined.Factor), top5.beta) %>% gather(Gene, Methylation, -Combined)
    top5.df$Gene %<>% factor(levels = top5.symbol)
    top5.df$Combined %<>% str_replace_all("\\.", " ") %>% factor(levels = c("duodenum mucosa", "duodenum crypts", "colon mucosa", "colon crypts")) 

    contrast.title <- str_replace_all(contrast, "_", " ") %>% str_replace_all("m\\.", "m ") %>% str_replace_all("n\\.", "n ")
    p <- ggplot(top5.df, aes(x = Combined, y = Methylation, color = Combined)) + geom_jitter() + geom_boxplot() + theme_bw()
    p <- p + facet_wrap(~ Gene, ncol = 5, scales = "free") + theme(legend.position = "none")
    p <- p + theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_blank())
    p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.x = element_blank()) 
    p <- p + theme(plot.background = element_blank(), panel.border = element_rect(size = 1, color = "black"))
    p <- p + theme(plot.title = element_text(hjust = 0.5))
    p <- p + ggtitle(contrast.title) + ylim(0,1)
    CairoPDF(str_c("top5", prefix, contrast, sep = "_"), height = 4, width = 16, bg = "transparent")
    print(p)
    dev.off()
}

objects.size <- lapply(ls(), function(thing) print(object.size(get(thing)), units = 'auto')) 
names(objects.size) <- ls()
unlist(objects.size) %>% sort

source("../FRDA project/common_functions.R")
#sample.data <- read.xlsx("./111016 UNGC Web samples sheet Goldman lab.xlsx")
sample.data <- read_csv("./2016-9236 Sample Sheet.csv") %>% data.frame
sample.data$barcode <- str_c(sample.data$chip.ID, sample.data$stripe, sep = "_")

write_csv(sample.data, "./sample_annotations.csv")
idat.dir <- "./iScans"
sample.annotation <- "./sample_annotations.csv"
report.dir <- "./reports"

rnb.initialize.reports(report.dir)
logger.start(fname = NA)
parallel.setup(7)
#options(fftempdir = "~/tmp/Rtmp")
rnb.options(logging.disk = TRUE)
rnb.options(disk.dump.big.matrices=FALSE)

data.source <- c(idat.dir, sample.annotation)
#rnb.options(import.gender.prediction = FALSE)
result <- rnb.run.import(data.source = data.source, data.type = "infinium.idat.dir", dir.reports = report.dir)
rnb.set <- result$rnb.set
SaveRDSgz(rnb.set, "./save/rnb.set.rda")

rnb.g19.filter <- !grepl("^[1-8]$", pheno(rnb.set)$External.Sample.ID)
rnb.g19 <- remove.samples(rnb.set, rnb.g19.filter)

rnb.run.qc(rnb.g19, report.dir)
raw.beta <- meth(rnb.g19)
rownames(raw.beta) <- rownames(annotation(rnb.g19, type = "sites"))
write.csv(raw.beta, "./allsites_beta.csv")

biological.clock <- read_csv("./allsites_beta.output.csv")
g19.annotations <- read_csv("./g19_annotations.csv") %>% data.frame
g19.annotations$Cell.Type <- str_c(g19.annotations$Sort.Marker, "Div", g19.annotations$DIV, sep = "_") %>% factor

g19.plot <- data.frame(Cell.Type = g19.annotations$Cell.Type, DNAm.Age = biological.clock$DNAmAge)
p <- ggplot(g19.plot, aes(x = Cell.Type, y = DNAm.Age, color = Cell.Type)) + geom_point() + theme_bw()
p <- p + theme(legend.position = "none")
p <- p + theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_blank())
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.x = element_blank()) 
p <- p + theme(plot.background = element_blank(), panel.border = element_rect(size = 1, color = "black"))
p <- p + theme(plot.title = element_text(hjust = 0.5)) + ggtitle("DNA methylation age")
CairoPDF("biological_clock", height = 6, width = 6, bg = "transparent")
print(p)
dev.off()

rnb.g19 <- addPheno(rnb.g19, g19.annotations$Cell.Type, "Cell.Type")
rnb.filter <- rnb.execute.context.removal(rnb.g19)$dataset
rnb.filter <- rnb.execute.snp.removal(rnb.filter, snp = "any")$dataset
rnb.filter <- rnb.execute.sex.removal(rnb.filter)$dataset
rnb.greedy <- rnb.execute.greedycut(rnb.filter)
filter.sites <- rnb.greedy$sites
rnb.filter <- remove.sites(rnb.filter, filter.sites)

rnb.filter <- rnb.execute.na.removal.fix(rnb.filter, 0)$dataset
rnb.filter <- rnb.execute.variability.removal(rnb.filter, 0.005)$dataset
rnb.norm <- rnb.execute.normalization(rnb.filter, method = "bmiq", bgcorr.method = "methylumi.noob")
SaveRDSgz(rnb.norm, "./save/rnb.norm.rda")

dred.sites <- rnb.execute.dreduction(rnb.norm)
dred.promoters <- rnb.execute.dreduction(rnb.norm, target = "promoters")
dred.genes <- rnb.execute.dreduction(rnb.norm, target = "genes")

dred.plot <- data.frame(dred.promoters$mds$manhattan[,1:2])
colnames(dred.plot) <- c("PC1", "PC2")
dred.plot %<>% mutate(Cell.Type = pheno(rnb.norm)$Cell.Type)

#Convert this plot ggplot
p <- ggplot(data = dred.plot, aes(x = PC1, y = PC2, col = Cell.Type)) + geom_point()
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p <- p + theme(plot.background = element_blank(), legend.background = element_blank())
CairoPDF("promoters_celltype", width = 8, height = 6, bg = "transparent")
print(p)
dev.off()

dred.allsites <- data.frame(dred.sites$mds$manhattan[,1:2])
colnames(dred.allsites) <- c("PCA1", "PCA2")
dred.allsites %<>% mutate(Cell.Type = pheno(rnb.norm)$Cell.Type)

p <- ggplot(data = dred.allsites, aes(x = PCA1, y = PCA2, col = Cell.Type)) + geom_point()
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p <- p + theme(plot.background = element_blank(), legend.background = element_blank())
CairoPDF("allsites_celltype", width = 8, height = 6, bg = "transparent")
print(p)
dev.off()

dred.genes.plot <- data.frame(dred.genes$mds$manhattan[,1:2])
colnames(dred.genes.plot) <- c("PCA1", "PCA2")
dred.genes.plot %<>% mutate(Cell.Type = pheno(rnb.norm)$Cell.Type)

p <- ggplot(data = dred.genes.plot, aes(x = PCA1, y = PCA2, col = Cell.Type)) + geom_point()
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p <- p + theme(plot.background = element_blank(), legend.background = element_blank())
CairoPDF("genes_celltype", width = 8, height = 6, bg = "transparent")
print(p)
dev.off()

comp.cols <- "Cell.Type"
reg.types <- c("genes", "promoters")

diffmeth.adj <- rnb.execute.computeDiffMeth(rnb.norm, pheno.cols = comp.cols, pheno.cols.all.pairwise = comp.cols, region.types = reg.types)

comparisons <- get.comparisons(diffmeth.adj)
comparisons.format <- str_replace(comparisons, " \\(based on Cell.Type\\)", "") %>% str_replace_all(" ", "_") 

tab.sites <- map(comparisons, get.table, object = diffmeth.adj, region.type = "sites", return.data.frame = TRUE) 
names(tab.sites) <- comparisons.format
tab.promoters <- map(comparisons, get.table, object = diffmeth.adj, region.type = "promoters", return.data.frame = TRUE) 
names(tab.promoters) <- comparisons.format
tab.genes <- map(comparisons, get.table, object = diffmeth.adj, region.type = "genes", return.data.frame = TRUE) 
names(tab.genes) <- comparisons.format

SaveRDSgz(tab.sites, "./save/tab.sites.rda")

map(comparisons.format, MapVolcanoPlots, "promoters", tab.promoters)
map(comparisons.format, MapVolcanoPlots, "genes", tab.genes)
map(comparisons.format, MapVolcanoPlots, "sites", tab.sites)

#Get Thresholds
thresholds <- c(0.01, 0.005, 0.001)
thresholds.fdr <- c(0.05, 0.01, 0.005)

promoters.threshold <- map(thresholds, . %>% map(.x = tab.promoters, .f = GetSizes, pval.column = "comb.p.val", log.column = "mean.mean.quot.log2") %>% reduce(rbind) %>% data.frame) %>% 
    map2(.y = thresholds, .f = mutate) %>% 
    reduce(rbind) %>%
    mutate(Comparison = rep(comparisons, length(thresholds)))
promoters.threshold$Comparison %<>% str_replace(" \\(based on Combined.Factor\\)", "")
colnames(promoters.threshold)[3] <- "Test"
promoters.threshold$Test <- str_c("p < ", promoters.threshold$Test)

promoters.threshold.fdr <- map(thresholds.fdr, . %>% map(.x = tab.promoters, .f = GetSizes, pval.column = "comb.p.adj.fdr", log.column = "mean.mean.quot.log2") %>% reduce(rbind) %>% data.frame) %>% 
    map2(.y = thresholds.fdr, .f = mutate) %>% 
    reduce(rbind) %>%
    mutate(Comparison = rep(comparisons, length(thresholds.fdr)))
promoters.threshold.fdr$Comparison %<>% str_replace(" \\(based on Combined.Factor\\)", "")
colnames(promoters.threshold.fdr)[3] <- "Test"
promoters.threshold.fdr$Test <- str_c("FDR p < ", promoters.threshold.fdr$Test) 

promoters.threshold.combine <- rbind(promoters.threshold, promoters.threshold.fdr) 
promoters.threshold.combine$positive %<>% unlist
promoters.threshold.combine$negative %<>% unlist
promoters.threshold.combine$Test %<>% factor
promoters.threshold.genetotals <- select(promoters.threshold.combine, positive, negative, Test) %>% group_by(Test) %>% summarise(sum(sum(positive), sum(abs(negative))))
colnames(promoters.threshold.genetotals)[2] <- "Total.Sites"
promoters.threshold.genetotals$Total.Sites %<>% str_c(" Sites")
promoters.threshold.combine %<>% join(promoters.threshold.genetotals)
promoters.threshold.plot <- gather(promoters.threshold.combine, Direction, Num.Sites, -Test, -Comparison, -Total.Sites)

genes.threshold <- map(thresholds, . %>% map(.x = tab.genes, .f = GetSizes, pval.column = "comb.p.val", log.column = "mean.mean.quot.log2") %>% reduce(rbind) %>% data.frame) %>% 
    map2(.y = thresholds, .f = mutate) %>% 
    reduce(rbind) %>%
    mutate(Comparison = rep(comparisons, length(thresholds)))
genes.threshold$Comparison %<>% str_replace(" \\(based on Combined.Factor\\)", "")
colnames(genes.threshold)[3] <- "Test"
genes.threshold$Test <- str_c("p < ", genes.threshold$Test)

genes.threshold.fdr <- map(thresholds.fdr, . %>% map(.x = tab.genes, .f = GetSizes, pval.column = "comb.p.adj.fdr", log.column = "mean.mean.quot.log2") %>% reduce(rbind) %>% data.frame) %>% 
    map2(.y = thresholds.fdr, .f = mutate) %>% 
    reduce(rbind) %>%
    mutate(Comparison = rep(comparisons, length(thresholds.fdr)))
genes.threshold.fdr$Comparison %<>% str_replace(" \\(based on Combined.Factor\\)", "")
colnames(genes.threshold.fdr)[3] <- "Test"
genes.threshold.fdr$Test <- str_c("FDR p < ", genes.threshold.fdr$Test)

genes.threshold.combine <- rbind(genes.threshold, genes.threshold.fdr) 
genes.threshold.combine$positive %<>% unlist
genes.threshold.combine$negative %<>% unlist
genes.threshold.combine$Test %<>% factor
genes.threshold.genetotals <- select(genes.threshold.combine, positive, negative, Test) %>% group_by(Test) %>% summarise(sum(sum(positive), sum(abs(negative))))
colnames(genes.threshold.genetotals)[2] <- "Total.Sites"
genes.threshold.genetotals$Total.Sites %<>% str_c(" Sites")
genes.threshold.combine %<>% join(genes.threshold.genetotals)
genes.threshold.plot <- gather(genes.threshold.combine, Direction, Num.Sites, -Test, -Comparison, -Total.Sites)

sites.threshold <- map(thresholds, . %>% map(.x = tab.sites, .f = GetSizes, pval.column = "diffmeth.p.val", log.column = "mean.quot.log2") %>% reduce(rbind) %>% data.frame) %>% 
    map2(.y = thresholds, .f = mutate) %>% 
    reduce(rbind) %>%
    mutate(Comparison = rep(comparisons, length(thresholds)))
sites.threshold$Comparison %<>% str_replace(" \\(based on Combined.Factor\\)", "")
colnames(sites.threshold)[3] <- "Test"
sites.threshold$Test <- str_c("p < ", sites.threshold$Test)

sites.threshold.fdr <- map(thresholds.fdr, . %>% map(.x = tab.sites, .f = GetSizes, pval.column = "diffmeth.p.adj.fdr", log.column = "mean.quot.log2") %>% reduce(rbind) %>% data.frame) %>% 
    map2(.y = thresholds.fdr, .f = mutate) %>% 
    reduce(rbind) %>%
    mutate(Comparison = rep(comparisons, length(thresholds.fdr)))
sites.threshold.fdr$Comparison %<>% str_replace(" \\(based on Combined.Factor\\)", "")
colnames(sites.threshold.fdr)[3] <- "Test"
sites.threshold.fdr$Test <- str_c("FDR p < ", sites.threshold.fdr$Test)

sites.threshold.combine <- rbind(sites.threshold, sites.threshold.fdr) 
sites.threshold.combine$positive %<>% unlist
sites.threshold.combine$negative %<>% unlist
sites.threshold.combine$Test %<>% factor
sites.threshold.genetotals <- select(sites.threshold.combine, positive, negative, Test) %>% group_by(Test) %>% summarise(sum(sum(positive), sum(abs(negative))))
colnames(sites.threshold.genetotals)[2] <- "Total.Sites"
sites.threshold.genetotals$Total.Sites %<>% str_c(" Sites")
sites.threshold.combine %<>% join(sites.threshold.genetotals)
sites.threshold.plot <- gather(sites.threshold.combine, Direction, Num.Sites, -Test, -Comparison, -Total.Sites)

DecidePlot("promoter_threshold_selection", promoters.threshold.plot, "Differentially Methylated Promoters", bar.padding = 500)
DecidePlot("gene_threshold_selection", genes.threshold.plot, "Differentially Methylated Genes", bar.padding = 900)
DecidePlot("site_threshold_selection", sites.threshold.plot, "Differentially Methylated Sites", bar.padding = 10000)

#Gene Annotation
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
vega <- useMart("ENSEMBL_MART_VEGA", dataset = "hsapiens_gene_vega")

promoters.annotation <- annotation(rnb.norm, type = "promoters")
promoters.annotation$Ensembl.ID <- rownames(promoters.annotation)
SaveRDSgz(promoters.annotation, "./save/promoters_annotation.rda")

promoters.bm.table <- getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol', 'description', 'gene_biotype'), filters = 'ensembl_gene_id', values = as.character(promoters.annotation$Ensembl.ID), mart = ensembl)
promoters.bm.table$description %<>% str_replace_all(" \\[.*\\]$", "")
colnames(promoters.bm.table) <- c("Ensembl.ID", "Symbol", "Description", "Gene.Type")

promoters.bm.table.vega <- getBM(attributes = c('ens_gene', 'vega_gene_id', 'external_gene_name', 'description', 'gene_biotype'), filters = 'ens_gene', values = as.character(promoters.annotation$Ensembl.ID), mart = vega)
colnames(promoters.bm.table.vega) <- c("Ensembl.ID", "Vega.ID", "Gene.Name", "Description.Vega", "Gene.Type.Vega")

SaveRDSgz(promoters.bm.table, "./save/promoters.bm.table.rda")
SaveRDSgz(promoters.bm.table.vega, "./save/promoters.bm.table.vega.rda")

genes.annotation <- annotation(rnb.norm, type = "genes")
genes.annotation$Ensembl.ID <- rownames(genes.annotation)
SaveRDSgz(genes.annotation, "./save/genes_annotation.rda")

genes.bm.table <- getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol','description', 'gene_biotype'), filters = 'ensembl_gene_id', values = as.character(genes.annotation$Ensembl.ID), mart = ensembl)
genes.bm.table$description %<>% str_replace_all(" \\[.*\\]$", "")
colnames(genes.bm.table) <- c("Ensembl.ID", "Symbol", "Description", "Gene.Type")

genes.bm.table.vega <- getBM(attributes = c('ens_gene', 'vega_gene_id', 'external_gene_name', 'description', 'gene_biotype'), filters = 'ens_gene', values = as.character(genes.annotation$Ensembl.ID), mart = vega)
colnames(genes.bm.table.vega) <- c("Ensembl.ID", "Vega.ID", "Gene.Name", "Description.Vega", "Gene.Type.Vega")

SaveRDSgz(genes.bm.table, "./save/genes.bm.table.rda")
SaveRDSgz(genes.bm.table.vega, "./save/genes.bm.table.vega.rda")

tab.promoters.annot <- map(tab.promoters, mutate, Ensembl.ID = rownames(promoters.annotation)) %>% 
    map(join, promoters.annotation) %>%
    map(join, promoters.bm.table) %>%
    map(join, promoters.bm.table.vega) %>%
    map(dplyr::select, Ensembl.ID, Symbol:Gene.Type, mean.mean.quot.log2, comb.p.val, comb.p.adj.fdr, combinedRank, Chromosome:Strand, num.sites, CpG:G, Vega.ID:Gene.Type.Vega, entrezID, mean.mean.g1:mean.mean.diff, mean.num.na.g1:mean.nsamples.covg.thresh.g2) 
SaveRDSgz(tab.promoters.annot, "./save/tab.promoters.annot.rda")

tab.genes.annot <- map(tab.genes, mutate, Ensembl.ID = rownames(genes.annotation)) %>% 
    map(join, genes.annotation) %>%
    map(join, genes.bm.table) %>%
    map(join, genes.bm.table.vega) %>%
    map(dplyr::select, Ensembl.ID, Symbol:Gene.Type, mean.mean.quot.log2, comb.p.val, comb.p.adj.fdr, combinedRank, Chromosome:Strand, num.sites, CpG:G, Vega.ID:Gene.Type.Vega, entrezID, mean.mean.g1:mean.mean.diff, mean.num.na.g1:mean.nsamples.covg.thresh.g2) 
SaveRDSgz(tab.genes.annot, "./save/tab.genes.annot.rda")

names(tab.promoters.annot) <- comparisons.format
names(tab.genes.annot) <- comparisons.format

promoters.table.reduce <- map(tab.promoters.annot, dplyr::select, Ensembl.ID:Gene.Type.Vega) %>% 
    lmap(UniqueNames) %>% 
    reduce(join) %>%
    dplyr::select(Ensembl.ID:Gene.Type, dplyr::contains("log2FC"), dplyr::matches("p.val$"), dplyr::contains("p.val.adj"), dplyr::contains("combinedRank"), Chromosome:Gene.Type.Vega) 
promoters.table.reduce[,grepl("p.val|log2FC", colnames(promoters.table.reduce))] %<>% signif(3)
DMWorkbook(promoters.table.reduce, "./promoters.combined.xlsx")

genes.table.reduce <- map(tab.genes.annot, dplyr::select, Ensembl.ID:Gene.Type.Vega) %>% 
    lmap(UniqueNames) %>% 
    reduce(join) %>%
    dplyr::select(Ensembl.ID:Gene.Type, dplyr::contains("log2FC"), dplyr::matches("p.val$"), dplyr::contains("p.val.adj"), dplyr::contains("combinedRank"), Chromosome:Gene.Type.Vega) 
genes.table.reduce[,grepl("p.val|log2FC", colnames(genes.table.reduce))] %<>% signif(3)
DMWorkbook(genes.table.reduce, "./genes.combined.xlsx")

promoters.norm <- meth(rnb.norm, type = "promoters")
rownames(promoters.norm) <- rownames(annotation(rnb.norm, type = "promoters"))

rtkn.meth <- promoters.norm["ENSG00000114993",]
rtkn.df <- data.frame(Cell.Type = pheno(rnb.norm)$Cell.Type, RTKN = rtkn.meth)
