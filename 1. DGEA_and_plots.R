### RNA-seq analysis: MEL - old and new samples ###

# Libraries #####
library(dplyr)
library(DESeq2)
library(edgeR)

# Preprocessing #####

# Load the files into R
files = list.files("Count_files/")
x = readDGE(paste0("Count_files/", files), 
            columns = c(1,2))

# Rename samples with more readable names
sample_names_prev = colnames(x)
samplenames_almost = substring(sample_names_prev, nchar("Count_files/") + 1,
                               nchar(sample_names_prev) - 7)

# Add prefix to indicate batch
samplenames = samplenames_almost
for (i in 1:ncol(x)) {
  if (i <= 8) {
    samplenames[i] = paste0("2023_", samplenames[i])
  } else {
    samplenames[i] = paste0("AUTH_", samplenames[i])
  }
}

colnames(x) = samplenames
rm(sample_names_prev, samplenames_almost, i, samplenames)

# Create batch variable
batch = c(rep("2023_exp", 8), rep("auth", 6)) # 2 batches

HMBA_ind = grepl("HMBA", colnames(x))
treatment = rep("control", 14)
treatment[HMBA_ind] = "HMBA" # 2 treatments

timepoint = c("24h", "24h", "48h", "48h", "72h", "72h", "48h", "48h",
              "24h", "48h", "72h", "24h", "48h", "72h")

# Add to the DGE object
x$samples$batch = factor(batch)
x$samples$batch = relevel(x$samples$batch, ref = "auth")
x$samples$treatment = factor(treatment)
x$samples$treatment = relevel(x$samples$treatment, ref = "control")
x$samples$timepoint = factor(timepoint)
x$samples$timepoint = relevel(x$samples$timepoint, ref = "24h")
rm(batch, treatment, timepoint, HMBA_ind)

# Map the genes to Entrez IDs so that we have at least two annotations available
library(org.Mm.eg.db)
official = org.Mm.egSYMBOL
mapped_genes_official = mappedkeys(official)
official_df = as.data.frame(official[mapped_genes_official])
official_df = official_df %>% dplyr::rename(EntrezGene.ID = gene_id, Gene.Symbol = symbol)
official_df$HGNC_Official = "Yes"
official_df = official_df[-which(duplicated(official_df$Gene.Symbol)), ]

rm(official, mapped_genes_official, files)

# Create a genes slot with gene annotation
gene_annot = as.data.frame(list(Gene.Symbol = rownames(x))) %>%
  left_join(official_df, by = "Gene.Symbol") %>%
  dplyr::select(-HGNC_Official)
x$genes = gene_annot

# DESeq2 objects preparation #####

# How many genes in this dataset are all zeros?
table(rowSums(x$counts==0)==14) # 14: number of samples

# FALSE  TRUE 
# 16520  8463

print(paste0(round(100*8463/(16520+8463), 2), "% of genes only have zero counts."))
# 33.88% of genes only have zero counts.

# Create a full design matrix for the filtering
design = model.matrix(~ x$samples$batch + x$samples$treatment + x$samples$timepoint +
                        x$samples$treatment:x$samples$timepoint)

# We need to filter out uninformative genes
keep.exprs = names(which(filterByExpr(x, design = design) == TRUE))
x_filt = x
x_filt$counts = x_filt$counts[keep.exprs,]
x_filt$genes = x_filt$genes[x_filt$genes$Gene.Symbol %in% keep.exprs,]
dim(x_filt)

# Create a TPM matrix
library(AnnotationHub)

# Create an AnnotationHub object
ah = AnnotationHub()

# Query for the latest Mus musculus EnsDb object
latest_mus_musculus = query(ah, c("EnsDb", "Mus musculus")) 

# Print the latest version
tail(latest_mus_musculus)

# ID of interest is "AH109650"
edb = ah[["AH109650"]]

# Get gene data
genes_data = genes(edb)
genes_map = data.frame(gene_id = genes_data$gene_id,
                       Gene.Symbol = genes_data$gene_name)

# Get transcript data
transcripts_data = transcripts(edb)
transcript_lengths = data.frame(transcript_id = transcripts_data$tx_id,
                                gene_id = transcripts_data$gene_id,
                                tx_length = width(transcripts_data)) %>%
  inner_join(genes_map, by = "gene_id") %>%
  dplyr::filter(!Gene.Symbol == "") %>%
  group_by(Gene.Symbol) %>%
  summarise(longest_transcript_length = max(tx_length, na.rm = TRUE))

# Calculate gene lengths in kilobases
transcript_lengths$transcript_length_kb = transcript_lengths$longest_transcript_length / 1000

# Map gene lengths to rownames of the counts matrix
rownames(transcript_lengths) = transcript_lengths$Gene.Symbol
mapped = transcript_lengths[transcript_lengths$Gene.Symbol %in% 
                              x_filt$genes$Gene.Symbol,]
mapped_lengths = mapped$transcript_length_kb

# Calculate counts per kilobase (CPK)
cpk = x_filt$counts[mapped$Gene.Symbol, ] / mapped_lengths
# Calculate the sum of CPK values for each sample
cpk_sum = colSums(cpk)

# Calculate TPM values
tpm = cpk / cpk_sum * 1e6
log2tpm = log2(tpm + 1)
rm(cpk, cpk_sum, keep.exprs, mapped_lengths, ah, edb)

removals = c(which(x_filt$genes$Gene.Symbol == "no_feature"),
             which(x_filt$genes$Gene.Symbol == "ambiguous"))
x_filt$genes = x_filt$genes[-removals, ]
x_filt$counts = x_filt$counts[-removals, ]

# DESeq2 object
dds = DESeqDataSetFromMatrix(countData = x_filt$counts,
                             colData = x_filt$samples,
                             design = ~ batch + treatment + timepoint +
                               treatment:timepoint)

featureData = as.data.frame(x_filt$genes)
mcols(dds) = DataFrame(mcols(dds), featureData)
mcols(dds)
rm(featureData)

# Density curve function
library(ggplot2)
library(tidyr)
create_density_curve <- function(matrix) {
  as.data.frame(matrix) %>%
    gather(.) %>%
    ggplot(aes(x = value, fill = "Density Curve")) + 
    geom_density(alpha = 0.5, position = "identity") +
    xlim(-1, NA) +
    labs(title = "Distribution of values by sample",
         x = "Value", y = "Density") +
    theme_classic()+
    theme(legend.position = "none")
}

# Per column density plot
create_density_plot_color <- function(matrix) {
  as.data.frame(matrix) %>%
    gather(., key = "sample", value = "value") %>%
    ggplot(aes(x = value, fill = sample)) + 
    geom_density(alpha = 0.5, position = "identity") +
    xlim(-1, NA) +
    labs(title = "Distribution of values by sample",
         x = "Value", y = "Density") +
    theme_classic()+
    theme(legend.position = "none")
}

# Create folder for distribution plots
if (!dir.exists("MEL_exploratory_plots")) {
  dir.create("MEL_exploratory_plots")
}

if (!dir.exists("MEL_exploratory_plots/MEL_distribution_plots")) {
  dir.create("MEL_exploratory_plots/MEL_distribution_plots")
}

# Density plots/curves for our data
create_density_curve(x_filt$counts)
ggsave(filename = "MEL_filtered_read_counts_density_curves.tiff",
       path = "MEL_exploratory_plots/MEL_distribution_plots",
       width = 1920, height = 1080, device = 'tiff', units = "px",
       dpi = 150, compression = "lzw")
dev.off()

# The previous plot contains count outliers (abnormally high counts) which make 
# it difficult to see the distribution of counts near zero. We replot with 
# xlim 200

create_density_curve(x_filt$counts)+
  scale_x_continuous(limits = c(-1, 200), breaks = seq(0, 200, 10))
ggsave(filename = "MEL_filtered_read_counts_density_curves_xlim200.tiff",
       path = "MEL_exploratory_plots/MEL_distribution_plots",
       width = 1920, height = 1080, device = 'tiff', units = "px",
       dpi = 150, compression = "lzw")
dev.off()

# Color density plots for counts
create_density_plot_color(x_filt$counts)
ggsave(filename = "MEL_filtered_read_counts_color_densities.tiff",
       path = "MEL_exploratory_plots/MEL_distribution_plots",
       width = 1920, height = 1080, device = 'tiff', units = "px",
       dpi = 150, compression = "lzw")
dev.off()

# With xlim 200
create_density_plot_color(x_filt$counts)+
  scale_x_continuous(limits = c(-1, 200), breaks = seq(0, 200, 10))
ggsave(filename = "MEL_filtered_read_counts_color_densities_xlim200.tiff",
       path = "MEL_exploratory_plots/MEL_distribution_plots",
       width = 1920, height = 1080, device = 'tiff', units = "px",
       dpi = 150, compression = "lzw")
dev.off()

# Now with TPM
create_density_curve(log2tpm)
ggsave(filename = "MEL_filtered_log2tpm_density_curves.tiff",
       path = "MEL_exploratory_plots/MEL_distribution_plots",
       width = 1920, height = 1080, device = 'tiff', units = "px",
       dpi = 150, compression = "lzw")
dev.off()

# Color density plots for counts
create_density_plot_color(log2tpm)
ggsave(filename = "MEL_filtered_log2tpm_color_densities.tiff",
       path = "MEL_exploratory_plots/MEL_distribution_plots",
       width = 1920, height = 1080, device = 'tiff', units = "px",
       dpi = 150, compression = "lzw")
dev.off()

# Further exploratory analysis #####
# VST and rlog transformations
vsd = vst(dds, blind = FALSE)
rld = rlog(dds, blind = FALSE)

# See https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#data-transformations-and-visualization
# for why blind = FALSE

dds = estimateSizeFactors(dds)

# Create folder for mean-variance plots and exploratory heatmaps
if (!dir.exists("MEL_exploratory_plots/Transformations_and_heatmaps")) {
  dir.create("MEL_exploratory_plots/Transformations_and_heatmaps")
}

# Produce a scatter plot of transformed counts for the first two samples
df = bind_rows(
  as.data.frame(log2(counts(dds, normalized=TRUE)[, 1:2]+1)) %>%
    mutate(transformation = "log2(x + 1)"),
  as.data.frame(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst"),
  as.data.frame(assay(rld)[, 1:2]) %>% mutate(transformation = "rlog"))

colnames(df)[1:2] <- c("x", "y")  

lvls = c("log2(x + 1)", "vst", "rlog")
df$transformation = factor(df$transformation, levels=lvls)

ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) +
  coord_fixed() + facet_grid( . ~ transformation)

# Produce mean - sd plots 
library(vsn)
library(ggpubr)
library(cowplot)

ntd = normTransform(dds) # log2(counts + 1) transformation
plot1 = meanSdPlot(assay(ntd))$gg + rremove("ylab") + rremove("xlab") + 
  ggtitle(bquote(bold("log2(counts + 1)"))) +
  theme(plot.title = element_text(hjust = 0.5, vjust = 0.5))
plot2 = meanSdPlot(assay(rld))$gg + rremove("ylab") + rremove("xlab") + 
  ggtitle(bquote(bold("rlog"))) +
  theme(plot.title = element_text(hjust = 0.5, vjust = 0.5))
plot3 = meanSdPlot(assay(vsd))$gg + rremove("ylab") + rremove("xlab") + 
  ggtitle(bquote(bold("vst"))) +
  theme(plot.title = element_text(hjust = 0.5, vjust = 0.5))
dev.off()

# Use ggarrange to align the three plots horizontally and set the x-axis and y-axis labels
require(grid)
fig = ggarrange(plot1, plot2, plot3, 
                ncol = 3, nrow = 1, 
                common.legend = TRUE, legend = "bottom",
                labels = c("A", "B", "C"),
                align = "hv",
                font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top"))
tiff("MEL_exploratory_plots/Transformations_and_heatmaps/meanSdPlot.tiff", res = 300, width = 4500, height = 1600, compression = "lzw")
fig
annotate_figure(fig, left = textGrob(bquote(bold("sd")), rot = 90, vjust = 1, gp = gpar(cex = 1)),
                bottom = textGrob(bquote(bold("rank (mean)")), gp = gpar(cex = 1)))
dev.off()

# Heatmap of the count matrix - top 50 genes
library(pheatmap)
library(rcartocolor)
select = order(rowMeans(counts(dds,normalized=TRUE)),
               decreasing=TRUE)[1:50]
df = as.data.frame(colData(dds)[,c("batch", "treatment", "timepoint")])

# Annotation colors
ann_colors = list(
  batch = c(`2023_exp` = "dodgerblue4", auth = "deeppink4"),
  treatment = c(HMBA = "purple4", control = "orange"),
  timepoint = c(`24h` = carto_pal(n = 7, "Emrld")[1], 
                `48h` = carto_pal(n = 7, "Emrld")[4], 
                `72h` = carto_pal(n = 7, "Emrld")[7])
)

# log2(counts + 1)
tiff("MEL_exploratory_plots/Transformations_and_heatmaps/top50_heatmap_log2(counts+1).tiff",
     res = 200, width = 1420, 
     height = 1080, compression = "lzw")
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df, main = "log2(counts + 1)",
         annotation_colors = ann_colors, fontsize = 6)
dev.off()

# rlog
tiff("MEL_exploratory_plots/Transformations_and_heatmaps/top50_heatmap_rlog.tiff",
     res = 200, width = 1420, 
     height = 1080, compression = "lzw")
pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df, main = "rlog",
         annotation_colors = ann_colors, fontsize = 6)
dev.off()

# vst
tiff("MEL_exploratory_plots/Transformations_and_heatmaps/top50_heatmap_vst.tiff", res = 200, width = 1420, 
     height = 1080, compression = "lzw")
pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df, main = "vst",
         annotation_colors = ann_colors, fontsize = 6)
dev.off()

# We will only use the vst-transformed data for the next visualisations

# Heatmap of sample distances
library(RColorBrewer)
library(tibble)

sampleVSTdists = dist(t(assay(vsd)))
sampleDistMatrix = as.matrix(sampleVSTdists)
#colnames(sampleDistMatrix) = NULL
colors = colorRampPalette(viridisLite::magma(10))(255)

tiff("MEL_exploratory_plots/Transformations_and_heatmaps/sample_distance_heatmap_vst.tiff", 
     res = 700, width = 6760, 
     height = 4760, compression = "lzw")
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleVSTdists,
         clustering_distance_cols = sampleVSTdists,
         col = colors,
         fontsize = 8,
         main = "Sample distance heatmap",
         annotation_col = df,
         annotation_colors = ann_colors)
dev.off()

# Heatmap of sample Poisson distances (more appropriate for count data)
library(PoiClaClu)
poisd = PoissonDistance(t(counts(dds)))
samplePoisDistMatrix = as.matrix(poisd$dd)
dimnames(samplePoisDistMatrix) = dimnames(sampleDistMatrix)
#colnames(samplePoisDistMatrix) = NULL

tiff("MEL_exploratory_plots/Transformations_and_heatmaps/sample_poisson_distance_heatmap.tiff", res = 700, width = 6760, 
     height = 4760, compression = "lzw")
pheatmap(samplePoisDistMatrix,
         clustering_distance_rows = poisd$dd,
         clustering_distance_cols = poisd$dd,
         col = colors,
         fontsize = 8,
         main = "Sample Poisson distance heatmap",
         annotation_col = df,
         annotation_colors = ann_colors)
dev.off()

# Dimensionality reduction #####
# Create dimensionality reduction plots folder 
if (!dir.exists("MEL_exploratory_plots/Dimensionality_reduction_plots")) {
  dir.create("MEL_exploratory_plots/Dimensionality_reduction_plots")
}

# Principal Component Analysis
library(ggrepel)
pcaData = plotPCA(vsd, intgroup = c("treatment", "batch"), 
                  returnData=TRUE) %>%
  dplyr::rename(Sample.ID = name)
time = as.data.frame(vsd@colData@listData) %>%
  rownames_to_column %>%
  dplyr::select(Sample.ID = rowname, timepoint)
pcaData = pcaData %>%
  inner_join(time, by = "Sample.ID")
percentVar = round(100 * attr(pcaData, "percentVar"))
pcaplot = ggplot(pcaData, aes(PC1, PC2, color = treatment, shape = batch)) +
  geom_point(size = 0.5, alpha = 0.65) +
  geom_text_repel(aes(label = timepoint), segment.linetype = 1,
                      segment.color = "black", segment.size = 0.1,
                  min.segment.length = 0,
                  size = 1, nudge_y = 0.03, nudge_x = 0.03, show.legend = FALSE) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  ggtitle("Principal Component Analysis: VST") +
  theme_classic() +
  scale_color_manual(name = "treatment", values = c("orange", "purple4"),
                     labels = c("control", "HMBA"))+
  theme(plot.title = element_text(size = 4, face = "bold", vjust = 0.5, hjust = 0.5),
        axis.text = element_text(size = 3, hjust = 0.5, vjust = 0.5, 
                                 color = "black"),
        axis.title = element_text(size = 4, face = "bold"),
        axis.ticks = element_line(linewidth = 0.1),
        axis.line = element_line(linewidth = 0.3),
        legend.position = "right",
        legend.key.size = unit(2, units = "mm"),
        legend.text = element_text(size = 3),
        legend.title = element_text(face = "bold", size = 3.5),
        legend.margin = ggplot2::margin(0, 0, 0, 0, unit = "mm"),
        legend.spacing.y = unit(0.5, units = "mm"),
        legend.spacing.x = unit(0.5, units = "mm"),
        legend.background = element_blank())

# Generalized PCA (more appropriate for count data)
library(glmpca)
gpca = glmpca(counts(dds), L=2)
gpca.dat = gpca$factors
gpca.dat$treatment = dds$treatment
gpca.dat$batch = dds$batch
gpca.dat$timepoint = dds$timepoint

gpcaplot = ggplot(gpca.dat, aes(x = dim1, y = dim2, color = treatment, shape = batch)) +
  geom_point(size = 0.5, alpha = 0.65) +
  geom_text_repel(aes(label = timepoint), segment.linetype = 1,
                      segment.color = "black", segment.size = 0.1,
                  min.segment.length = 0,
                  size = 1, nudge_y = 0.03, nudge_x = 0.03, show.legend = FALSE) +
  ggtitle("Generalized Principal Component Analysis (glmpca, counts)") +
  theme_classic() +
  scale_color_manual(name = "treatment", values = c("orange", "purple4"),
                     labels = c("control", "HMBA"))+
  theme(plot.title = element_text(size = 4, face = "bold", vjust = 0.5, hjust = 0.5),
        axis.text = element_text(size = 3, hjust = 0.5, vjust = 0.5, 
                                 color = "black"),
        axis.title = element_text(size = 4, face = "bold"),
        axis.ticks = element_line(linewidth = 0.1),
        axis.line = element_line(linewidth = 0.3),
        legend.position = "right",
        legend.key.size = unit(2, units = "mm"),
        legend.text = element_text(size = 2),
        legend.title = element_text(face = "bold", size = 3.5),
        legend.margin = ggplot2::margin(0, 0, 0, 0, unit = "mm"),
        legend.spacing.y = unit(0.5, units = "mm"),
        legend.spacing.x = unit(0.5, units = "mm"),
        legend.background = element_blank())

tiff("MEL_exploratory_plots/Dimensionality_reduction_plots/PCA.tiff", width = 2880, height = 1380, res = 700,
     compression = "lzw")
ggarrange(pcaplot, gpcaplot, 
          ncol = 2, nrow = 1, common.legend = TRUE, 
          labels = c("A", "B"), legend = "right",
          font.label = list(size = 4, color = "black", face = "bold"))
dev.off()

# Multidimensional scaling plots for VST and Poisson distances

# VST
mds <- as.data.frame(colData(vsd))  %>%
  cbind(cmdscale(sampleDistMatrix))
mdsplot = ggplot(mds, aes(x = `1`, y = `2`, color = treatment, 
                          shape = batch)) +
  geom_point(size = 0.5, alpha = 0.65) +
  geom_text_repel(aes(label = timepoint, segment.linetype = 1,
                      segment.color = "black", segment.size = 0.1),
                  min.segment.length = 0,
                  size = 1, nudge_y = 0.03, nudge_x = 0.03) +
  ggtitle("Multidimensional scaling: VST") +
  theme_classic() +
  scale_color_manual(name = "treatment", values = c("orange", "purple4"),
                     labels = c("control", "HMBA"))+
  theme(plot.title = element_text(size = 4, face = "bold", vjust = 0.5, hjust = 0.5),
        axis.text = element_text(size = 3, hjust = 0.5, vjust = 0.5, 
                                 color = "black"),
        axis.title = element_text(size = 4, face = "bold"),
        axis.ticks = element_line(linewidth = 0.1),
        axis.line = element_line(linewidth = 0.3),
        legend.position = "right",
        legend.key.size = unit(2, units = "mm"),
        legend.text = element_text(size = 3),
        legend.title = element_text(face = "bold", size = 3.5),
        legend.margin = ggplot2::margin(0, 0, 0, 0, unit = "mm"),
        legend.spacing.y = unit(0.5, units = "mm"),
        legend.spacing.x = unit(0.5, units = "mm"),
        legend.background = element_blank())+
  labs(x = "MDS1", y = "MDS2") +
  guides(alpha = "none")

# Poisson
mdsPois <- as.data.frame(colData(dds)) %>%
  cbind(cmdscale(samplePoisDistMatrix))
poisson_mdsplot = ggplot(mdsPois, aes(x = `1`, y = `2`, color = treatment, 
                                      shape = batch)) +
  geom_point(size = 0.5, alpha = 0.65) +
  geom_text_repel(aes(label = timepoint, segment.linetype = 1,
                      segment.color = "black", segment.size = 0.1),
                  min.segment.length = 0,
                  size = 1, nudge_y = 0.03, nudge_x = 0.03) +
  ggtitle("Multidimensional scaling: Poisson distance (counts)") +
  theme_classic() +
  scale_color_manual(name = "treatment", values = c("orange", "purple4"),
                     labels = c("control", "HMBA"))+
  theme(plot.title = element_text(size = 4, face = "bold", vjust = 0.5, hjust = 0.5),
        axis.text = element_text(size = 3, hjust = 0.5, vjust = 0.5, 
                                 color = "black"),
        axis.title = element_text(size = 4, face = "bold"),
        axis.ticks = element_line(linewidth = 0.1),
        axis.line = element_line(linewidth = 0.3),
        legend.position = "right",
        legend.key.size = unit(2, units = "mm"),
        legend.text = element_text(size = 3),
        legend.title = element_text(face = "bold", size = 3.5),
        legend.margin = ggplot2::margin(0, 0, 0, 0, unit = "mm"),
        legend.spacing.y = unit(0.5, units = "mm"),
        legend.spacing.x = unit(0.5, units = "mm"),
        legend.background = element_blank())+
  labs(x = "MDS1", y = "MDS2") +
  guides(alpha = "none")

fig = ggarrange(mdsplot +rremove("ylab"), poisson_mdsplot +rremove("ylab"), 
                ncol = 2, nrow = 1, 
                common.legend = TRUE, legend = "right",
                labels = c("A", "B"),
                align = "hv",
                font.label = list(size =4, color = "black", face = "bold"))
tiff("MEL_exploratory_plots/Dimensionality_reduction_plots/MDS.tiff", res = 700, width = 2880, height = 1080, compression = "lzw")
fig
annotate_figure(fig, left = textGrob(bquote(bold("MDS2")), rot = 90, vjust = 1, gp = gpar(cex = 0.45)))
dev.off()

# Differential Expression #####
# We are going to use a LRT as described in the DESeq2 vignette:
# https://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#time-series-experiments
# and the RNAseq workflow by Michael Love:
# http://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html#time-course-experiments

dds = DESeq(dds, test = "LRT", reduced = ~ batch + treatment + timepoint)
resultsNames(dds) # Names of the coefficients

contrasts = list(HMBA24h_vs_Control48h = list("treatment_HMBA_vs_control",
                 "timepoint_48h_vs_24h"),
                 HMBA48h_vs_Control48h = list(c("treatment_HMBA_vs_control", 
                                                "treatmentHMBA.timepoint48h")),
                 HMBA72h_vs_Control48h = list(c("treatment_HMBA_vs_control",
                                                "timepoint_72h_vs_24h",
                                                "treatmentHMBA.timepoint72h"),
                                              "timepoint_48h_vs_24h"),
                 HMBA48h_vs_HMBA24h = list(c("timepoint_48h_vs_24h",
                                              "treatmentHMBA.timepoint48h")),
                 HMBA72h_vs_HMBA24h = list(c("treatmentHMBA.timepoint72h",
                                             "timepoint_72h_vs_24h")),
                 HMBA72h_vs_HMBA48h = list(c("timepoint_48h_vs_24h",
                                             "treatmentHMBA.timepoint48h"),
                                           c("treatmentHMBA.timepoint72h",
                                             "timepoint_72h_vs_24h")))

# hmba24h vs control48h
res_hmba24h_vs_control48h <- results(dds, contrast=contrasts[["HMBA24h_vs_Control48h"]], 
                                     listValues=c(1, -1), test = "Wald")

# hmba48h vs control48h
res_hmba48h_vs_control48h <- results(dds, 
                                     contrast=contrasts[["HMBA48h_vs_Control48h"]], 
                                     listValues=c(1, -1), test = "Wald")

# hmba72h vs control48h
res_hmba72h_vs_control48h <- results(dds, contrast=contrasts[["HMBA72h_vs_Control48h"]],
                                     listValues=c(1, -1), test = "Wald")

# HMBA48h vs HMBA24h
res_HMBA48h_vs_HMBA24h <- results(dds, 
                                  contrast=contrasts[["HMBA48h_vs_HMBA24h"]],
                                  listValues=c(1, -1), test = "Wald")

# HMBA72h vs HMBA24h
res_HMBA72h_vs_HMBA24h <- results(dds, 
                                  contrast=contrasts[["HMBA72h_vs_HMBA24h"]], 
                                  listValues=c(1, -1), test = "Wald")

# HMBA72h vs HMBA48h
res_HMBA72h_vs_HMBA48h <- results(dds, contrast=contrasts[["HMBA72h_vs_HMBA48h"]],
                                  listValues=c(1, -1), test = "Wald")

results = list(res_hmba24h_vs_control48h, res_hmba48h_vs_control48h, 
               res_hmba72h_vs_control48h, res_HMBA48h_vs_HMBA24h,
               res_HMBA72h_vs_HMBA24h, res_HMBA72h_vs_HMBA48h)
names(results) = names(contrasts)

# Shrinking logFCs
# apeglm_shrunk_results = list() # apeglm does not work with contrasts
# normal_shrunk_results = list() # normal cannot be implemented with interactions designs
ashr_shrunk_results = list()

for (i in 1:length(results)) {
  # apeglm_shrunk_results[[i]] = lfcShrink(dds, contrast = contrasts[[i]], type = "apeglm")
  # normal_shrunk_results[[i]] = lfcShrink(dds, contrast = contrasts[[i]], type = "normal")
  ashr_shrunk_results[[i]] = lfcShrink(dds, res = results[[i]], type = "ashr")
}
names(ashr_shrunk_results) = names(results)

# Create a workbook and save data there to export later
library(openxlsx)
dgea_wb = createWorkbook()
dgea_wb_shrunk = createWorkbook()

# Create header and body styles for the output Excel Sheet
header_style = createStyle(textDecoration = "bold")
body_style = createStyle(halign = "left", valign = "center")

# Choose if you want the colored gene symbols or not (color coding is time-consuming)
color_handle = FALSE # Change that to FALSE if you don't want colored symbols

if (color_handle) {
  blue_style = createStyle(fontColour = "#279AE5", textDecoration = "bold")  # Blue and bold
  red_style = createStyle(fontColour = "#B91B60", textDecoration = "bold")  # Red and bold
  grey_style = createStyle(fontColour = "#818383", textDecoration = "bold")  # Grey and bold
  
  # Define a function to determine the color style based on the conditions
  getStyle <- function(log2foldchange, padj) {
    if (is.na(log2foldchange) | is.na(padj)) {
      return(grey_style)
    } else if (log2foldchange < 0 & padj < 0.05) {
      return(blue_style)
    } else if (log2foldchange > 0 & padj < 0.05) {
      return(red_style)
    } else {
      return(grey_style)
    }
  }
}

# Create data frames that will be exported into .xlsx files
for (i in 1:length(results)) {
  # Raw LFCs
  genes = rownames(results[[i]])
  res = results[[i]]@listData %>%
    as.data.frame() %>% 
    mutate(Gene.Symbol = genes) %>%
    dplyr::arrange(padj) %>%
    dplyr::select(Gene.Symbol, everything())
  
  addWorksheet(dgea_wb, names(results)[i])
  writeData(dgea_wb, names(results)[i], res)
  
  # Apply the styles
  addStyle(dgea_wb, names(results)[i], style = header_style, rows = 1, 
           cols = 1:ncol(res), gridExpand = TRUE) # first row only
  addStyle(dgea_wb, names(results)[i], style = body_style, rows = 2:nrow(res), 
           cols = 1:ncol(res), gridExpand = TRUE)
  
  if (color_handle) {
    color_styles = mapply(getStyle, res$log2FoldChange, res$padj)
    mapply(addStyle, MoreArgs = list(wb = dgea_wb, sheet = names(results)[i], 
                                     cols = which(colnames(res) == "Gene.Symbol")), 
           style = color_styles, rows = 2:(nrow(res)+1))
  }
  
  # ASHR LFCs
  genes = rownames(ashr_shrunk_results[[i]])
  res_ashr = ashr_shrunk_results[[i]]@listData %>%
    as.data.frame() %>% 
    mutate(Gene.Symbol = genes) %>%
    dplyr::arrange(padj) %>%
    dplyr::select(Gene.Symbol, everything())
  
  addWorksheet(dgea_wb_shrunk, names(ashr_shrunk_results)[i])
  writeData(dgea_wb_shrunk, names(ashr_shrunk_results)[i], res_ashr)
  
  # Apply the styles
  addStyle(dgea_wb_shrunk, names(ashr_shrunk_results)[i], style = header_style, 
           rows = 1, cols = 1:ncol(res_ashr), gridExpand = TRUE) # first row only
  addStyle(dgea_wb_shrunk, names(ashr_shrunk_results)[i], style = body_style, 
           rows = 2:nrow(res_ashr), cols = 1:ncol(res_ashr), gridExpand = TRUE)
  
  if (color_handle) {
    color_styles = mapply(getStyle, res_ashr$log2FoldChange, res_ashr$padj)
    mapply(addStyle, MoreArgs = list(wb = dgea_wb_shrunk, 
                                     sheet = names(ashr_shrunk_results)[i], 
                                     cols = which(colnames(res_ashr) == "Gene.Symbol")), 
           style = color_styles, rows = 2:(nrow(res_ashr)+1))
  }
}

# Create DGEA folder
if (!dir.exists("DGEA")) {
  dir.create("DGEA")
}

saveWorkbook(dgea_wb, "DGEA/DGEA_raw.xlsx",
             overwrite = TRUE)
saveWorkbook(dgea_wb_shrunk, "DGEA/DGEA_ashr.xlsx",
             overwrite = TRUE)
rm(genes, res, res_ashr, i); gc()

# MA-plots
if (!dir.exists("DGEA/MA_plots")) {
  dir.create("DGEA/MA_plots")
}

for (i in 1:length(results)) {
  tiff(paste0("DGEA/MA_plots/MA_plot_",
  names(results)[i], ".tiff"), res = 700,
       width = 8000, height = 4000, compression = "lzw")
  par(mfrow=c(1,2), mar=c(4,4,2,1))
  print(DESeq2::plotMA(results[[i]], main = "MA plot: raw LFCs",
                       xlab = substitute(paste(bold('Mean of normalized counts'))),
                       ylab = substitute(paste(bold('log fold change'))),
                       alpha = 0.05))
  print(DESeq2::plotMA(ashr_shrunk_results[[i]], main = "MA plot: ashr shrunk LFCs",
                       xlab = substitute(paste(bold('Mean of normalized counts'))),
                       ylab = substitute(paste(bold('log fold change'))),
                       alpha = 0.05))
  dev.off()
}

# Volcano plots
library(EnhancedVolcano)

if (!dir.exists("DGEA/Volcano_plots")) {
  dir.create("DGEA/Volcano_plots")
}

# Volcano loop (raw LFCs)
for (i in 1:length(results)) {
  # Construct the data frame
  genes = rownames(results[[i]])
  res = results[[i]]@listData %>%
    as.data.frame() %>% 
    mutate(Gene.Symbol = genes) %>%
    dplyr::arrange(padj)
  
  # create custom key-value pairs for stat. sig genes (p.adj < 0.05) and n.s genes
  keyvals.colour <- ifelse(
    res$log2FoldChange < -1 & res$padj < 0.05, 'royalblue',
    ifelse(res$log2FoldChange > 1 & res$padj < 0.05, 'red4',
           ifelse(abs(res$log2FoldChange) < 1 & res$padj < 0.05, 'pink',
                  'grey')))
  
  # keyvals.colour[is.na(keyvals.colour)] <- 'black'
  names(keyvals.colour)[keyvals.colour == 'royalblue'] <- 'Down-regulated'
  names(keyvals.colour)[keyvals.colour == 'red4'] <- 'Up-regulated'
  names(keyvals.colour)[keyvals.colour == 'pink'] <- '|log2FC| < 1'
  names(keyvals.colour)[keyvals.colour == 'grey'] <- 'p.adj > 0.05'
  
  res$aes = keyvals.colour
  
  volcano = EnhancedVolcano(res,
                            lab = res[, "Gene.Symbol"],
                            caption = NULL,
                            x = 'log2FoldChange',
                            y = 'padj',
                            pCutoff = 0.05,
                            FCcutoff = 1,
                            cutoffLineType = "dashed",
                            cutoffLineWidth = 0.3,
                            cutoffLineCol = "black",
                            colCustom = keyvals.colour,
                            colAlpha = 0.7,
                            ylab = bquote(bold(-log[10]("BH adj. p-value"))),
                            xlab = bquote(bold(log[2]("Fold Change"))),
                            pointSize = 1,
                            axisLabSize = 7,
                            subtitle = NULL,
                            labSize = 2,
                            selectLab = res[1:20, "Gene.Symbol"],
                            legendLabSize = 8,
                            legendIconSize = 4,
                            title = paste0("Volcano plot: contrast = ",
                            names(results)[i]),
                            labFace = "bold",
                            boxedLabels = TRUE,
                            drawConnectors = TRUE,
                            typeConnectors = "closed",
                            arrowheads = FALSE,
                            widthConnectors = 0.3,
                            max.overlaps = Inf) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(linewidth = 0.4),
          plot.title = element_text(face = "bold", size = 10, hjust = 0.5, vjust = 0.5),
          axis.title = element_text(face = "bold", size = 8),
          axis.line = element_line(colour = "black", linewidth = 0.4),
          axis.ticks = element_line(colour = "black", linewidth = 0.4),
          axis.ticks.length = unit(1, units = "mm"),
          legend.position = "right",
          legend.title = element_blank(),
          legend.margin = ggplot2::margin(1, 1, 1, 1, unit = "mm"),
          legend.spacing.y = unit(1, units = "mm"),
          legend.spacing.x = unit(1, units = "mm"),
          legend.background = element_blank(),
          legend.key.size = unit(0.5, units = "mm"),
          legend.text = element_text(size = 6),
          legend.box.background = element_rect(colour = "black"))
  
  volcano
  ggsave(filename = paste0("Volcano_raw_LFCs_", names(results)[i], ".tiff"),
         path = "DGEA/Volcano_plots", 
         width = 152, height = 142, device = 'tiff', units = "mm",
         dpi = 700, compression = "lzw")
}
rm(genes, i, res, keyvals.colour, volcano)

# Volcano loop (ashr-shrunk LFCs)
for (i in 1:length(ashr_shrunk_results)) {
  # Construct the data frame
  genes = rownames(ashr_shrunk_results[[i]])
  res = ashr_shrunk_results[[i]]@listData %>%
    as.data.frame() %>% 
    mutate(Gene.Symbol = genes) %>%
    dplyr::arrange(padj)
  
  # create custom key-value pairs for stat. sig genes (p.adj < 0.05) and n.s genes
  keyvals.colour <- ifelse(
    res$log2FoldChange < -1 & res$padj < 0.05, 'royalblue',
    ifelse(res$log2FoldChange > 1 & res$padj < 0.05, 'red4',
           ifelse(abs(res$log2FoldChange) < 1 & res$padj < 0.05, 'pink',
                  'grey')))
  
  # keyvals.colour[is.na(keyvals.colour)] <- 'black'
  names(keyvals.colour)[keyvals.colour == 'royalblue'] <- 'Down-regulated'
  names(keyvals.colour)[keyvals.colour == 'red4'] <- 'Up-regulated'
  names(keyvals.colour)[keyvals.colour == 'pink'] <- '|log2FC| < 1'
  names(keyvals.colour)[keyvals.colour == 'grey'] <- 'p.adj > 0.05'
  
  res$aes = keyvals.colour
  
  volcano = EnhancedVolcano(res,
                            lab = res[, "Gene.Symbol"],
                            caption = NULL,
                            x = 'log2FoldChange',
                            y = 'padj',
                            pCutoff = 0.05,
                            FCcutoff = 1,
                            cutoffLineType = "dashed",
                            cutoffLineWidth = 0.3,
                            cutoffLineCol = "black",
                            colCustom = keyvals.colour,
                            colAlpha = 0.7,
                            ylab = bquote(bold(-log[10]("BH adj. p-value"))),
                            xlab = bquote(bold(log[2]("Fold Change"))),
                            pointSize = 1,
                            axisLabSize = 7,
                            subtitle = NULL,
                            labSize = 2,
                            selectLab = res[1:20, "Gene.Symbol"],
                            legendLabSize = 8,
                            legendIconSize = 4,
                            title = paste0("Volcano plot: contrast = ",
                                           names(results)[i]),
                            labFace = "bold",
                            boxedLabels = TRUE,
                            drawConnectors = TRUE,
                            typeConnectors = "closed",
                            arrowheads = FALSE,
                            widthConnectors = 0.3,
                            max.overlaps = Inf) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(linewidth = 0.4),
          plot.title = element_text(face = "bold", size = 10, hjust = 0.5, vjust = 0.5),
          axis.title = element_text(face = "bold", size = 8),
          axis.line = element_line(colour = "black", linewidth = 0.4),
          axis.ticks = element_line(colour = "black", linewidth = 0.4),
          axis.ticks.length = unit(1, units = "mm"),
          legend.position = "right",
          legend.title = element_blank(),
          legend.margin = ggplot2::margin(1, 1, 1, 1, unit = "mm"),
          legend.spacing.y = unit(1, units = "mm"),
          legend.spacing.x = unit(1, units = "mm"),
          legend.background = element_blank(),
          legend.key.size = unit(0.5, units = "mm"),
          legend.text = element_text(size = 6),
          legend.box.background = element_rect(colour = "black"))
  
  volcano
  ggsave(filename = paste0("Volcano_ashr_shrunk_LFCs_", names(ashr_shrunk_results)[i], ".tiff"),
         path = "DGEA/Volcano_plots", 
         width = 152, height = 142, device = 'tiff', units = "mm",
         dpi = 700, compression = "lzw")
}
rm(genes, i, res, keyvals.colour, volcano)

# Gene clustering with mean-centered values from the VST data
library(genefilter)

# Create output directories
if (!dir.exists("DGEA/Post_DGEA_heatmaps")) {
  dir.create("DGEA/Post_DGEA_heatmaps")
}

# Subdirectories
if (!dir.exists("DGEA/Post_DGEA_heatmaps/raw_lfcs")) {
  dir.create("DGEA/Post_DGEA_heatmaps/raw_lfcs")
}

if (!dir.exists("DGEA/Post_DGEA_heatmaps/ashr_shrunk")) {
  dir.create("DGEA/Post_DGEA_heatmaps/ashr_shrunk")
}

# Number of features for plots
n_degs = 20

# Top n_degs variable genes fromt the VST object
topVarGenes = head(order(rowVars(assay(vsd)), decreasing = TRUE), n_degs)
matvar  = assay(vsd)[topVarGenes, ]
matvar  = matvar - rowMeans(matvar)
colnames(matvar) = colnames(assay(vsd))
anno = as.data.frame(colData(vsd)[, c("batch", "treatment", "timepoint")])

# Create a custom palette with 20*14 = 280 color breaks (one for each tile)
colors2 = colorRampPalette(c(viridisLite::mako(280)[140:280],
                             rev(viridisLite::rocket(280)[140:280])))(280)

tiff(paste0("DGEA/Post_DGEA_heatmaps/VST_top20vargenes_samples_heatmap.tiff"),
     res = 700, width = 6760, 
     height = 4760, compression = "lzw")
pheatmap(matvar, annotation_col = anno,
         col = colors2,
         fontsize = 8,
         cutree_col = 2,
         main = "Top 20 variable VST genes (centered) heatmap",
         annotation_colors = ann_colors,
         cluster_rows = FALSE)
dev.off()

# Heatmaps for raw LFCs
for (i in 1:length(results)) {
  # Construct the data frame
  genes = rownames(results[[i]])
  res = results[[i]]@listData %>%
    as.data.frame() %>% 
    mutate(Gene.Symbol = genes) %>%
    dplyr::arrange(padj)
  
  topDEGs = res[which(res$padj < 0.05),]
  
  # Gene clustering with raw counts and Poisson distance
  topDEGs = topDEGs[1:n_degs, ]
  matcount = assay(dds)[topDEGs$Gene.Symbol,]
  DEGpoisd = PoissonDistance(t(counts(dds[topDEGs$Gene.Symbol,])))
  DEGsamplePoisDistMatrix = as.matrix(DEGpoisd$dd)
  rownames(DEGsamplePoisDistMatrix) = dds$cell_line
  colnames(DEGsamplePoisDistMatrix) = NULL
  
  # Same for genes
  DEGpoisd_genes = PoissonDistance(counts(dds[topDEGs$Gene.Symbol,]))
  DEGgenesPoisDistMatrix = as.matrix(DEGpoisd_genes$dd)
  rownames(DEGgenesPoisDistMatrix) = topDEGs$Gene.Symbol
  colnames(DEGgenesPoisDistMatrix) = NULL
  
  # Log-transformed matrix
  logmat = log2(matcount + 1)
  
  # The same matrix but order samples by metastatic status and genes by logFC (raw)
  genes_order = topDEGs %>% dplyr::filter(padj < 0.05) %>%
    arrange(log2FoldChange)
  gaps_row = length(which(genes_order$log2FoldChange < 0))
  genes_order = genes_order$Gene.Symbol
  
  sample_order = anno %>% dplyr::arrange(treatment)
  gaps_col = length(which(sample_order$treatment == "control"))
  sample_order = rownames(sample_order)
  
  matcount_ordered = matcount[genes_order, sample_order]
  
  # Mean-centering
  logmat_znorm = (logmat - rowMeans(logmat))/sqrt(rowVars(logmat))
  logmat_ordered_znorm = logmat_znorm[genes_order, sample_order]
  
  # In the plot below the cluster distances are not affected by the log-transformation
  # of the matrix. It is only used to make a smoother color scale
  tiff(paste0("DGEA/Post_DGEA_heatmaps/raw_lfcs/",
              names(results)[i], "_DEG_Poisson_clustered_samples_and_genes.tiff"),
       res = 700, width = 7260, 
       height = 5760, compression = "lzw")
  pheatmap(logmat_znorm, annotation_col = anno,
           breaks = seq(-4, 4, 8/length(colors2)),
           annotation_colors = ann_colors,
           clustering_distance_rows = DEGpoisd_genes$dd,
           clustering_distance_cols = DEGpoisd$dd,
           col = colors2,
           cutree_cols = 2,
           cutree_rows = 2,
           fontsize = 8,
           fontsize_row = 8,
           fontsize_col = 6,
           main = paste0(names(results)[i], 
                         " - DEGs vs. samples: Poisson heatmap, cluster cut = 2"),
           legend_breaks = c(-4, 4),
           legend_labels = c("lower expression", "higher expression"))
  dev.off()
  
  
  # Ordered genes and samples
  tiff(paste0("DGEA/Post_DGEA_heatmaps/raw_lfcs/",
              names(results)[i], "_DEG_samples_ordered.tiff"),
       res = 700, width = 7260, 
       height = 5760, compression = "lzw")
  pheatmap(logmat_ordered_znorm, annotation_col = anno,
           breaks = seq(-4, 4, 8/length(colors2)),
           annotation_colors = ann_colors,
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           col = colors2,
           gaps_row = gaps_row,
           gaps_col = gaps_col,
           cutree_cols = 2,
           fontsize = 8,
           fontsize_row = 8,
           fontsize_col = 6,
           border_color = "grey60",
           main = paste0(names(results)[i], 
                         " - DEGs vs. samples heatmap (ordered)"),
           legend_breaks = c(-4, 4),
           legend_labels = c("lower expression", "higher expression"))
  dev.off()
}
rm(genes, i, res, topDEGs, matcount, matcount_ordered, DEGgenesPoisDistMatrix,
   DEGpoisd, DEGpoisd_genes, DEGsamplePoisDistMatrix, logmat, logmat_ordered_znorm,
   logmat_znorm, sample_order, genes_order, gaps_col)

# Heatmaps for ashr-shrunk LFCs
for (i in 1:length(ashr_shrunk_results)) {
  # Construct the data frame
  genes = rownames(ashr_shrunk_results[[i]])
  res = ashr_shrunk_results[[i]]@listData %>%
    as.data.frame() %>% 
    mutate(Gene.Symbol = genes) %>%
    dplyr::arrange(padj)
  
  topDEGs = res[which(res$padj < 0.05),]
  
  # Gene clustering with raw counts and Poisson distance
  topDEGs = topDEGs[1:n_degs, ]
  matcount = assay(dds)[topDEGs$Gene.Symbol,]
  DEGpoisd = PoissonDistance(t(counts(dds[topDEGs$Gene.Symbol,])))
  DEGsamplePoisDistMatrix = as.matrix(DEGpoisd$dd)
  rownames(DEGsamplePoisDistMatrix) = dds$cell_line
  colnames(DEGsamplePoisDistMatrix) = NULL
  
  # Same for genes
  DEGpoisd_genes = PoissonDistance(counts(dds[topDEGs$Gene.Symbol,]))
  DEGgenesPoisDistMatrix = as.matrix(DEGpoisd_genes$dd)
  rownames(DEGgenesPoisDistMatrix) = topDEGs$Gene.Symbol
  colnames(DEGgenesPoisDistMatrix) = NULL
  
  # Log-transformed matrix
  logmat = log2(matcount + 1)
  
  # The same matrix but order samples by metastatic status and genes by logFC (raw)
  genes_order = topDEGs %>% dplyr::filter(padj < 0.05) %>%
    arrange(log2FoldChange)
  gaps_row = length(which(genes_order$log2FoldChange < 0))
  genes_order = genes_order$Gene.Symbol
  
  sample_order = anno %>% dplyr::arrange(treatment)
  gaps_col = length(which(sample_order$treatment == "control"))
  sample_order = rownames(sample_order)
  
  matcount_ordered = matcount[genes_order, sample_order]
  
  # Mean-centering
  logmat_znorm = (logmat - rowMeans(logmat))/sqrt(rowVars(logmat))
  logmat_ordered_znorm = logmat_znorm[genes_order, sample_order]
  
  # In the plot below the cluster distances are not affected by the log-transformation
  # of the matrix. It is only used to make a smoother color scale
  tiff(paste0("DGEA/Post_DGEA_heatmaps/ashr_shrunk/",
  names(ashr_shrunk_results)[i], "_DEG_Poisson_clustered_samples_and_genes.tiff"),
       res = 700, width = 7260, 
       height = 5760, compression = "lzw")
  pheatmap(logmat_znorm, annotation_col = anno,
           breaks = seq(-4, 4, 8/length(colors2)),
           annotation_colors = ann_colors,
           clustering_distance_rows = DEGpoisd_genes$dd,
           clustering_distance_cols = DEGpoisd$dd,
           col = colors2,
           cutree_cols = 2,
           cutree_rows = 2,
           fontsize = 8,
           fontsize_row = 8,
           fontsize_col = 6,
           main = paste0(names(ashr_shrunk_results)[i], 
                         " - DEGs vs. samples: Poisson heatmap, cluster cut = 2"),
           legend_breaks = c(-4, 4),
           legend_labels = c("lower expression", "higher expression"))
  dev.off()
  
  
  # Ordered genes and samples
  tiff(paste0("DGEA/Post_DGEA_heatmaps/ashr_shrunk/",
              names(ashr_shrunk_results)[i], "_DEG_samples_ordered.tiff"),
       res = 700, width = 7260, 
       height = 5760, compression = "lzw")
  pheatmap(logmat_ordered_znorm, annotation_col = anno,
           breaks = seq(-4, 4, 8/length(colors2)),
           annotation_colors = ann_colors,
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           col = colors2,
           gaps_row = gaps_row,
           gaps_col = gaps_col,
           cutree_cols = 2,
           fontsize = 8,
           fontsize_row = 8,
           fontsize_col = 6,
           border_color = "grey60",
           main = paste0(names(ashr_shrunk_results)[i], 
                         " - DEGs vs. samples heatmap (ordered)"),
           legend_breaks = c(-4, 4),
           legend_labels = c("lower expression", "higher expression"))
  dev.off()
}

rm(genes, i, res, topDEGs, matcount, matcount_ordered, DEGgenesPoisDistMatrix,
   DEGpoisd, DEGpoisd_genes, DEGsamplePoisDistMatrix, logmat, logmat_ordered_znorm,
   logmat_znorm, sample_order, genes_order, gaps_col)

gc()

# Session info
sessionInfo()