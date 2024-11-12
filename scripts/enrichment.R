################################################################################
# Bulk RNAseq preprocessing pipeline: Enrichment                               #
################################################################################

# Need to have installed datasets and dataformat command line tools
# from NCBI
if(!require("BiocManager", quietly = T)) install.packages("BiocManager")
if(!require("clusterProfiler", quietly = T)){
        BiocManager::install("clusterProfiler", update = F)
}
if(!require("org.Hs.eg.db", quietly = T)){
        BiocManager::install("org.Hs.eg.db", update = F)
}
if(!require(ggsignif, quietly = T)) install.packages("ggsignif")
if(!require(ggplot2, quietly = T)) install.packages("ggplot2")
if (!require("argparser", quietly = T)) install.packages("argparser")

library(clusterProfiler)
library(org.Hs.eg.db)
library(ggsignif)
library(ggplot2)
library(argparser)

# Terminal argument parser
################################################################################
parser <- arg_parser("Given a directory where aligmnent BAM files are and a gtf annotation file, obtain a counts matrix.")

parser <- add_argument(parser = parser,
                       arg = c("input",
                               "--alpha",
                               "--outDir"),
                       help = c("CSV DE file, obtained with DESeq2.",
                                "Alpha threshold for significance.",
                                "Output directory where placing the results."),
                       flag = c(F, F, F))

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

# Get dotplot for ORA
get_dotplot <- function(enrich, alph,
                        sort_by, topN = NULL){
        plot_df <- enrich@result[enrich@result$p.adjust <= alph, ]
        plot_df$gene_ratio <- sapply(plot_df$GeneRatio,
                                     function(x) as.numeric(strsplit(x, "/")[[1]][1])/as.numeric(strsplit(x, "/")[[1]][2]))
        if(sort_by == "p.adjust"){
                plot_df <- plot_df[order(plot_df$p.adjust), ]
        }else if(sort_by == "gene_ratio"){
                plot_df <- plot_df[order(plot_df$gene_ratio,
                                         decreasing = T), ]
        }
        if(!is.null(topN)){
                plot_df <- plot_df[1:topN, ]
        }
        plot_df$Description <- factor(plot_df$Description,
                                      levels = plot_df$Description[length(plot_df$Description):1])
        plt <- ggplot(plot_df, mapping = aes(x = Description, y = gene_ratio,
                                             color = p.adjust,
                                             size = Count)) +
                geom_point() + 
                scale_size(range = c(2, 10)) +
                labs(x = "GO name", y = "enrichment ratio",
                     color = "FDR",
                     size = "Count") +
                scale_color_gradient(low = "red", high = "blue") +
                guides(fill = guide_colorbar(reverse = TRUE)) +
                coord_flip() +
                scale_y_continuous(expand = expansion(mult = c(0.2, 0.2))) +
                theme(axis.text.y = element_text(size=15),
                      axis.text.x = element_text(size=10),
                      axis.title.x = element_text(size=18),
                      axis.title.y = element_blank(),
                      title = element_text(size=20),
                      panel.background = element_blank(),
                      panel.grid.major = element_line(colour = "gray"),
                      panel.grid.minor = element_blank(),
                      legend.text = element_text(size=12),
                      legend.title = element_text(size=13),
                      axis.line = element_line(colour = "black"),
                      axis.line.y = element_line(colour = "black"),
                      panel.border = element_rect(colour = "black",
                                                  fill=NA,
                                                  linewidth = 1))
        return(plt)
}

# Directory stuff
################################################################################
DE_file <- parsed$input
alph <- as.numeric(parsed$alpha)
outDir <- addSlashIfNot(parsed$outDir)

createIfNot(outDir)

# Load data
################################################################################
DE <- read.csv(DE_file,
               row.names = 1)

# Do ORA
################################################################################
DE_4ORA <- DE[!is.na(DE$padj), ]
DE_sign <- rownames(DE_4ORA)[DE_4ORA$padj <= alph]

ora_BP <- enrichGO(DE_sign,
                   org.Hs.eg.db,
                   keyType = "SYMBOL",
                   ont = "BP")

ora_CC <- enrichGO(DE_sign,
                   org.Hs.eg.db,
                   keyType = "SYMBOL",
                   ont = "CC")

ora_MF <- enrichGO(DE_sign,
                   org.Hs.eg.db,
                   keyType = "SYMBOL",
                   ont = "MF")

ora_BP_dtplt <- get_dotplot(ora_BP,
                            alph = alph,
                            sort_by = "gene_ratio",
                            topN = 10)
ora_CC_dtplt <- get_dotplot(ora_CC,
                            alph = alph,
                            sort_by = "gene_ratio",
                            topN = 10)
ora_MF_dtplt <- get_dotplot(ora_MF,
                            alph = alph,
                            sort_by = "gene_ratio",
                            topN = 10)

ggsave(sprintf("%sora_BP_dtplt.pdf", outDir), height = 4,
       width = 7, plot = ora_BP_dtplt)
ggsave(sprintf("%sora_CC_dtplt.pdf", outDir), height = 4,
       width = 7, plot = ora_CC_dtplt)
ggsave(sprintf("%sora_MF_dtplt.pdf", outDir), height = 4,
       width = 7, plot = ora_MF_dtplt)

print(sprintf("%s biological processes are significantly enriched (alpha <= %s).",
              sum(ora_BP@result$p.adjust <= alph),
              alph))
print(sprintf("%s cellular components are significantly enriched (alpha <= %s).",
              sum(ora_CC@result$p.adjust <= alph),
              alph))
print(sprintf("%s molecular functions are significantly enriched (alpha <= %s).",
              sum(ora_MF@result$p.adjust <= alph),
              alph))

write.csv(ora_BP@result, file = sprintf("%sora_BP.csv", outDir))
write.csv(ora_CC@result, file = sprintf("%sora_CC.csv", outDir))
write.csv(ora_MF@result, file = sprintf("%sora_MF.csv", outDir))