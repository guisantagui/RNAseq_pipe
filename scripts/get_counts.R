################################################################################
# Bulk RNAseq preprocessing pipeline: given a directory where aligmnent BAM    #
# files are and a gtf annotation file, obtain a counts matrix.                 #
################################################################################
if (!require("BiocManager", quietly = T))
        install.packages("BiocManager")
if (!require("Rsubread", quietly = T))
        BiocManager::install("Rsubread", update = F)
if (!require("argparser", quietly = T)) install.packages("argparser")

library(Rsubread)
library(argparser)

# Terminal argument parser
################################################################################
parser <- arg_parser("Given a directory where aligmnent BAM files are and a gtf annotation file, obtain a counts matrix.")

parser <- add_argument(parser = parser,
                       arg = c("input",
                               "--annotFile",
                               "--att",
                               "--strand",
                               "--outDir"),
                       help = c("Directory where aligned BAM files are. Each aligned BAM should be within a subdirectory, named according to the sample.",
                                "Annotation GTF file.",
                                "Attribute to aggregate counts.",
                                "Strandness. Can be 'none', 'forward' or 'reverse'.",
                                "Output directory where placing the results."),
                       flag = c(F, F, F, F, F))

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

# Directory stuff
################################################################################
alignDir <- addSlashIfNot(parsed$input)
annotFile <- parsed$annotFile
outDir <- addSlashIfNot(parsed$outDir)
att <- parsed$att
strand <- parsed$strand

createIfNot(outDir)

# Get counts
################################################################################
if (strand == "none"){
        std_int <- 0
}else if (strand == "forward"){
        std_int <- 1
}else if (strand == "reverse"){
        std_int <- 2
}

cntsLst <- list()
for(d in list.dirs(alignDir, recursive = F)){
        samp_name <- basename(d)
        bam_file <- list.files(d, full.names = T)
        bam_file <- bam_file[grep(".bam", bam_file)]
        cnts <- featureCounts(bam_file,
                              annot.ext = annotFile,
                              isPairedEnd = T,
                              isGTFAnnotationFile = T,
                              GTF.attrType = att,
                              GTF.featureType = "exon",
                              strandSpecific = std_int)
        cntsLst[[samp_name]] <- cnts
}

countsMat <- do.call(cbind, lapply(cntsLst, function(x) x$counts))
colnames(countsMat) <- names(cntsLst)

print(sprintf("The obtained counts matrix has %s features.",
              as.character(nrow(countsMat))))

write.csv(countsMat, file = sprintf("%scounts.csv", outDir))
print(sprintf("counts.csv saved at %s.", outDir))