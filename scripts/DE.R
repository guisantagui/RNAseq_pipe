################################################################################
# Bulk RNAseq preprocessing pipeline: Perform differential expression analysis.#
################################################################################

if (!require("BiocManager", quietly = T))
        install.packages("BiocManager")
if (!require("DESeq2", quietly = T))
        BiocManager::install("DESeq2", update = F)
if (!require(ggplot2, quietly = T)) install.packages("ggplot2")
if (!require(ggtext, quietly = T)) install.packages("ggtext")
if (!require(ggrepel, quietly = T)) install.packages("ggrepel")
if(!require(factoextra, quietly = T)) install.packages("factoextra")
if(!require(ggpubr, quietly = T)) install.packages("ggpubr")
if (!require("argparser", quietly = T)) install.packages("argparser")

library(DESeq2)
library(ggplot2)
library(ggtext)
library(ggrepel)
library(factoextra)
library(ggpubr)
library(argparser)

# Terminal argument parser
################################################################################
parser <- arg_parser("Given a directory where aligmnent BAM files are and a gtf annotation file, obtain a counts matrix.")

parser <- add_argument(parser = parser,
                       arg = c("--gene_counts",
                               "--ercc_counts",
                               "--sample_info",
                               "--useERCCs",
                               "--outDir"),
                       help = c("CSV file with the gene counts.",
                                "CSV file with the ERCC counts.",
                                "CSV file with sample names in one column, and treatment info in another one.",
                                "If ERCC should be taken into account for normalization.",
                                "Output directory where placing the results."),
                       flag = c(F, F, F, T, F))

parsed <- parse_args(parser)

# Functions
################################################################################

# Create directory if it doesn't exist
createIfNot <- function(pth){
        if (!dir.exists(pth)){
                dir.create(pth, recursive = T)
        }
}

# Add a / if it's not at the end of a directory string
addSlashIfNot <- function(pth){
        lastChar <- substr(pth, nchar(pth), nchar(pth))
        if(lastChar != "/"){
                pth <- paste0(pth, "/")
        }
        return(pth)
}

# Obtaining the top N contributor variables
# to the specified components (in our case PC1 and PC2)
getTopContrib <- function(PC, topN = 12, x = "PC1", y = "PC2"){
        contrib <- facto_summarize(PC, 
                                   "var", 
                                   axes = c(as.numeric(gsub("PC", "", x)),
                                            as.numeric(gsub("PC", "", y))))
        contrib <- contrib[order(contrib$contrib, decreasing = T), 
                           c("name", "contrib")]
        topContrib <- as.character(contrib$name[1:topN])
        return(topContrib)
}

# Do and plot a PCA.
plotPCACust <- function(df, scale = F, x = "PC1", y = "PC2", ntop = NA,
                        labs = F, topNFeats = NULL, biplot = F){
        if (scale){
                df <- apply(df, 1, function(x) (x - mean(x))/sd(x))
        }else{
                df <- t(df)
        }
        rv <- colVars(df)
        sel <- order(rv, decreasing = TRUE)[seq_len(min(500, 
                                                        length(rv)))]
        
        if(scale){
                PC <- prcomp(df[, sel], scale. = F, center = F)
        }else{
                PC <- prcomp(df[, sel], scale. = F, center = T)
        }
        
        
        dat <- data.frame(obsnames=row.names(PC$x), PC$x)
        dat <- dat[, c("obsnames", x, y)]
        colDat <- data.frame(dds@colData)
        dat$treatment <- colDat$treatment[match(dat$obsnames, colDat$sample)]
        propVar <- summary(PC)$importance[2, c(x, y)]
        propX <- round(propVar[names(propVar) == x]*100, digits = 2)
        propY <- round(propVar[names(propVar) == y]*100, digits = 2)
        
        datapc <- data.frame(varnames=rownames(PC$rotation), 
                             PC$rotation)
        mult <- min(
                (max(dat[,y]) - min(dat[,y])/(max(datapc[,y])-min(datapc[,y]))),
                (max(dat[,x]) - min(dat[,x])/(max(datapc[,x])-min(datapc[,x])))
        )
        datapc <- transform(datapc,
                            v1 = .7 * mult * (get(x)),
                            v2 = .7 * mult * (get(y))
        )
        datapc$x0 <- rep(0, nrow(datapc))
        datapc$y0 <- rep(0, nrow(datapc))
        if(!is.null(topNFeats)){
                varPlotFilt <- getTopContrib(PC, topN = topNFeats, x = x, y = y)
                datapc <- datapc[datapc$varnames %in% varPlotFilt, ]
        }
        
        pcaPlt <- ggplot(dat, aes(PC1, PC2, color=treatment, label = obsnames)) +
                geom_point(size=3) +
                xlab(sprintf("PC1 (%s %%)", propX)) +
                ylab(sprintf("PC1 (%s %%)", propY)) + 
                coord_fixed() +
                theme(title = ggtext::element_markdown(),
                      axis.title.y = ggtext::element_markdown(),
                      panel.background = element_blank(),
                      panel.border = element_rect(colour = "black", fill=NA,
                                                  linewidth = 1),
                      panel.grid.major = element_line(colour = "#d4d4d4"),
                      legend.position = "right")
        if (labs){
                pcaPlt <- pcaPlt +
                        geom_text_repel()
        }
        if (biplot){
                pcaPlt <- pcaPlt +
                        geom_text_repel(data=datapc, 
                                        aes(x=v1, y=v2, label=varnames), 
                                        color = "black", 
                                        size = 3,
                                        max.overlaps = 100) + 
                        geom_segment(data = datapc, aes(x=x0, 
                                                        y=y0, 
                                                        xend=v1, 
                                                        yend=v2, 
                                                        label = varnames),
                                     arrow = arrow(length=unit(0.2,"cm"),
                                                   type = "closed",
                                                   angle = 20), 
                                     alpha=0.75, 
                                     color="black", 
                                     size = 0.5)
        }
        return(pcaPlt)
}

# Plot p-value histogram
plot_pval_hist <- function(res_df, what_pval = "pvalue",
                           bins = NULL,
                           ylims = NULL){
        valVec <- res_df[, what_pval]
        valVec <- valVec[!is.na(valVec)]
        df <- data.frame(value = valVec)
        hst <- ggplot(data = df, mapping = aes(x = value)) +
                geom_histogram(bins = bins) +
                xlab(what_pval) +
                theme(title = ggtext::element_markdown(),
                      axis.title.y = ggtext::element_markdown(),
                      panel.background = element_blank(),
                      panel.border = element_rect(colour = "black", fill=NA,
                                                  linewidth = 1),
                      panel.grid.major = element_line(colour = "#d4d4d4"),
                      legend.position = "right")
        if (!is.null(ylims)){
                hst <- hst +
                        ylim(ylims)
        }
        return(hst)
}

# Obtain a volcano plot
getVolcPlot <- function(fitObjct, alpha = 0.05,
                        useLog2FC4Sign = T,
                        log2FC_thshld = 0.5, outputSignGenes = F,
                        labsThrshld = NULL,
                        topNGenes = NULL,
                        topNGenes_by = NULL){
        #labsThrshld <- c(4, 25)
        #fitObjct <- results
        alphaMinLog10 <- -log10(alpha)
        #compVec <- strsplit(colnames(fitObjct$coefficients), split = "-")[[1]]
        modelResults <- fitObjct
        
        volcPlotDF <- data.frame(proteins = rownames(modelResults),
                                 adj.P.Val = modelResults$padj,
                                 adj.pVal_minLog10 = -log10(modelResults$padj),
                                 log2FC = modelResults$log2FoldChange)
        volcPlotDF <- volcPlotDF[!is.na(volcPlotDF$adj.pVal_minLog10), ]
        
        volcPlotDF$sign <- volcPlotDF$adj.pVal_minLog10 >= alphaMinLog10
        if(useLog2FC4Sign){
                log2FC_bool <- volcPlotDF$log2FC >= log2FC_thshld | 
                        volcPlotDF$log2FC <= -1 *log2FC_thshld
                volcPlotDF$sign <- volcPlotDF$sign & log2FC_bool
        }
        volcPlotDF$sign <- c("not_significant",
                             "significant")[as.numeric(volcPlotDF$sign) + 1]
        
        volcPlotDF$sign <- factor(volcPlotDF$sign)
        volcPlotDF$labs <- volcPlotDF$sign == "significant"
        print(sprintf("There are %s significant genes with significance level of %s.",
                      sum(volcPlotDF$labs),
                      alpha))
        if(!is.null(labsThrshld)){
                volcPlotDF$labs <- volcPlotDF$labs & ((volcPlotDF$log2FC <= -labsThrshld[1] | volcPlotDF$log2FC >= labsThrshld[1]) | volcPlotDF$adj.pVal_minLog10 >= labsThrshld[2])
        }
        if(!is.null(topNGenes) & !is.null(topNGenes_by)){
                signDF <- volcPlotDF[volcPlotDF$sign == "significant", ]
                topGenes <- signDF$proteins[order(abs(signDF[, topNGenes_by]),
                                                  decreasing = T)][1:topNGenes]
                volcPlotDF$labs <- volcPlotDF$labs & (volcPlotDF$proteins %in% topGenes)
        }
        volcPlotDF_labs <- volcPlotDF[volcPlotDF$labs, ]
        volcPlotDF_subSet <- volcPlotDF[volcPlotDF$sign == "significant", ]
        
        #x_label <- sprintf("log2 FC (%s vs %s)", compVec[1], compVec[2])
        #x_label <- gsub("sick", "patients", x_label)
        #x_label <- gsub("ctrl", "controls", x_label)
        
        volcPlot <- ggplot(volcPlotDF, mapping = aes(x = log2FC,
                                                     y = adj.pVal_minLog10,
                                                     color = sign)) +
                geom_point() +
                scale_color_manual(values = c("significant" = "red",
                                              "not_sigificant" = "black")) +
                geom_hline(yintercept = alphaMinLog10, linetype = "dashed") +
                labs(x = "log2 FC (TFGB vs controls)",
                     y ="-log10(adj. p-val)") +
                geom_text_repel(data = volcPlotDF_labs,
                                aes(x = log2FC,
                                    y = adj.pVal_minLog10,
                                    label = proteins)) +
                theme(axis.text.y = element_text(size=15),
                      axis.text.x = element_text(size=15),
                      panel.background = element_blank(),
                      panel.grid.major = element_line(colour = "gray"), 
                      panel.grid.minor = element_blank(),
                      axis.line = element_line(colour = "black"),
                      axis.line.y = element_line(colour = "black"),
                      panel.border = element_rect(colour = "black",
                                                  fill=NA, linewidth = 1),
                      legend.position = "none")
        if(useLog2FC4Sign){
                volcPlot <- volcPlot +
                        geom_vline(xintercept = log2FC_thshld,
                                   linetype = "dashed") +
                        geom_vline(xintercept = -1 * log2FC_thshld,
                                   linetype = "dashed")
        }
        if(outputSignGenes){
                signGenes <- volcPlotDF_subSet
                out <- list(significative = signGenes,
                            plot = volcPlot)
        }else{
                out <- volcPlot
        }
        return(out)
}

# Directory stuff
################################################################################
#gene_counts_file <- "/Users/guillemsantamariaaguilar/Documents/100cia/applications/carlos_simon_foundation/round2/analysis/results/preprocessing/counts/ref/counts.csv"
#ercc_counts_file <- "/Users/guillemsantamariaaguilar/Documents/100cia/applications/carlos_simon_foundation/round2/analysis/results/preprocessing/counts/ercc/counts.csv"
#samp_info_file <- "/Users/guillemsantamariaaguilar/Documents/100cia/applications/carlos_simon_foundation/round2/analysis/data/sample_info.csv"
#useERCCs <- T
#outDir <- "/Users/guillemsantamariaaguilar/Documents/100cia/applications/carlos_simon_foundation/round2/analysis/results/DE/w_ercc/"

gene_counts_file <- parsed$gene_counts
ercc_counts_file <- parsed$ercc_counts
samp_info_file <- parsed$sample_info
useERCCs <- parsed$useERCCs
outDir <- addSlashIfNot(parsed$outDir)

createIfNot(outDir)

# Load data
################################################################################
gene_counts <- read.csv(gene_counts_file, row.names = 1)
ercc_counts <- read.csv(ercc_counts_file, row.names = 1)
samp_info <- read.csv(samp_info_file, row.names = 1)

# Filter low expressed genes
gene_counts <- gene_counts[rowSums(gene_counts) != 0, ]
ercc_counts <- ercc_counts[rowSums(ercc_counts) != 0, ]

print(sprintf("%s genomic features don't have zeros accross all the samples",
              nrow(gene_counts)))
print(sprintf("%s ERCC features don't have zeros accross all the samples",
              nrow(ercc_counts)))

#gene_counts <- gene_counts[, colnames(gene_counts) != "SRR5223570"]
#ercc_counts <- ercc_counts[, colnames(ercc_counts) != "SRR5223570"]
#samp_info <- samp_info[samp_info$sample != "SRR5223570", ]

samp_info$treatment <- factor(samp_info$treatment, levels = c("ctrl", "TGFB"))

# Create DESeq object
################################################################################

counts_merged <- rbind.data.frame(gene_counts,
                                  ercc_counts)

if (useERCCs){
        dds <- DESeqDataSetFromMatrix(countData = counts_merged,
                                      colData = samp_info,
                                      design = ~treatment)
        dds <- estimateSizeFactors(dds,
                                   controlGenes = rownames(dds) %in%
                                           rownames(ercc_counts))
}else{
        dds <- DESeqDataSetFromMatrix(countData = gene_counts,
                                      colData = samp_info,
                                      design = ~treatment)
}

# Run DESeq
################################################################################
dds <- DESeq(dds)

# Obtain results dataframe
results <- results(dds)
# Shrink logFC, given that we have small sample size
resLFC <- lfcShrink(dds, coef="treatment_TGFB_vs_ctrl", type="apeglm")

# Save results dataframes
write.csv(as.data.frame(results),
          file = sprintf("%sDESeqRes.csv",
                         outDir))

write.csv(as.data.frame(resLFC),
          file = sprintf("%sDESeqRes_LFCShrunk.csv",
                         outDir))
# Plot diagnostic plots
pdf(file = sprintf("%sdispEsts.pdf", outDir), height = 5, width = 5)
plotDispEsts(dds[!rownames(dds) %in% rownames(ercc_counts), ],
             ylim = c(1e-06, 1e02))
dev.off()

pdf(file = sprintf("%sMAPlot.pdf", outDir), height = 5, width = 5)
DESeq2::plotMA(results)
dev.off()

pdf(file = sprintf("%sMAPlot_lfcShrunk.pdf", outDir), height = 5, width = 5)
DESeq2::plotMA(resLFC)
dev.off()

# Do PCA and plot it 
pcaPlot <- plotPCACust(assay(rlog(dds)),
                       scale = F,
                       labs = T,
                       biplot = T,
                       topNFeats = 5)

ggsave(filename = sprintf("%spca_rlog.pdf", outDir),
       width = 5,
       height = 5,
       plot = pcaPlot)

# Plot p-value histograms
pval_hist <- plot_pval_hist(results,
                            what_pval = "pvalue",
                            bins = 50,
                            ylims = c(0, 2500))
padj_hist <- plot_pval_hist(results,
                            what_pval = "padj",
                            bins = 50,
                            ylims = c(0, 2500))

pval_hists_comb <- ggarrange(pval_hist, padj_hist)

ggsave(filename = sprintf("%spval_hists.pdf", outDir),
       height = 4, width = 8,
       plot = pval_hists_comb)

# Do Volcano plots
volcPlot_top10Genes_logFC <- getVolcPlot(fitObjct = resLFC,
                                         topNGenes = 10,
                                         topNGenes_by = "log2FC",
                                         useLog2FC4Sign = F)
volcPlot_top10Genes_pAdj <- getVolcPlot(fitObjct = resLFC,
                                        topNGenes = 10,
                                        topNGenes_by = "adj.pVal_minLog10",
                                        useLog2FC4Sign = F)

volcPlotComb <- ggarrange(volcPlot_top10Genes_logFC,
                          volcPlot_top10Genes_pAdj)

ggsave(filename = sprintf("%svolcPlots.pdf", outDir),
       height = 4, width = 8,
       plot = volcPlotComb)

# Get FC of FOXP3 and plot it
pdf(file = sprintf("%sFOXP3_counts.pdf", outDir), height = 5, width = 5)
plotCounts(dds, "FOXP3", intgroup = "treatment")
dev.off()

capture.output(print(resLFC["FOXP3", ]), file = sprintf("%sFOXP3.txt", outDir))
