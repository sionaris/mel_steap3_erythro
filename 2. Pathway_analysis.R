# This script uses the lists produced by DGEA for the purpose of ORA
# (Over-representation analysis) and GSEA
# (Gene Set Enrichment Analysis)

library(dplyr)
library(ggplot2)
library(clusterProfiler)
library(openxlsx)
library(org.Mm.eg.db)
library(ReactomePA)

# Import data (6 lists from DGEA)
ashr_dgea = list()
for (i in 1:6) {
  ashr_dgea[[i]] = openxlsx::read.xlsx("DGEA/DGEA_ashr.xlsx", sheet = i)
}
rm(i); gc()
names(ashr_dgea) = c("HMBA24h_vs_Control48h",
                     "HMBA48h_vs_Control48h",
                     "HMBA72h_vs_Control48h",
                     "HMBA48h_vs_HMBA24h",
                     "HMBA72h_vs_HMBA24h",
                     "HMBA72h_vs_HMBA48h")

# Create output folders
if (!dir.exists("Pathways")) {
  dir.create("Pathways")
}

for (i in 1:length(ashr_dgea)) {
  if (!dir.exists(paste0("Pathways/", names(ashr_dgea)[i]))) {
    dir.create(paste0("Pathways/", names(ashr_dgea)[i]))
  }
}

pathway_approaches = c("ORA", "GSEA", "pathfindR")

for (i in 1:length(pathway_approaches)) {
  for (j in 1:length(ashr_dgea)) {
    if (!dir.exists(paste0("Pathways/", names(ashr_dgea)[j], "/",
                           pathway_approaches[[i]]))) {
      dir.create(paste0("Pathways/", names(ashr_dgea)[j], "/",
                        pathway_approaches[[i]]))
    }
  }
}

rm(i, j); gc()

# Create lists of workbooks
ora_wb_list = list()
gsea_wb_list = list()
pathfindR_wb_list = list()

for (i in 1:length(ashr_dgea)) {
  ora_wb_list[[i]] = createWorkbook()
  gsea_wb_list[[i]] = createWorkbook()
  pathfindR_wb_list[[i]] = createWorkbook()
}
names(ora_wb_list) = names(gsea_wb_list) = names(pathfindR_wb_list) =
  names(ashr_dgea)
rm(i); gc()

# Create lists of output
ora_results = list()
gsea_results = list()
pathfindR_results = list()

# Set download option to libcurl for clusterProfiler
R.utils::setOption("clusterProfiler.download.method","libcurl")

# Objects to keep after every loop
keepers = c("ashr_dgea", "gsea_wb_list", "ora_wb_list", "pathfindR_wb_list",
            "gsea_results", "ora_results", "pathfindR_results", "pathway_approaches",
            "i")

# ORA and GSEA #####
for (i in 1:length(ashr_dgea)) {
  
  # Define the cutoff for significant genes
  sig_genes = ashr_dgea[[i]]$padj < 0.05
  
  # Extract the log2 fold changes and gene names of the significant genes
  log2fc = ashr_dgea[[i]]$log2FoldChange[sig_genes]
  genes = ashr_dgea[[i]]$Gene.Symbol[sig_genes]
  gene_list = log2fc
  names(gene_list) = genes
  rm(log2fc, genes)
  gene_list = sort(gene_list, decreasing = TRUE) # symbols
  gene_list_conv = suppressWarnings(bitr(names(gene_list), fromType = "SYMBOL",
                          toType = "ENTREZID", OrgDb = org.Mm.eg.db))
  gene_list_entrez = gene_list
  names(gene_list_entrez) = as.character(gene_list_conv$ENTREZID)
  gene_list_entrez = sort(gene_list_entrez, decreasing = TRUE) # Entrez ids
  entrez_universe = suppressWarnings(bitr(ashr_dgea[[i]]$Gene.Symbol, fromType = "SYMBOL",
                                          toType = "ENTREZID", OrgDb = org.Mm.eg.db))$ENTREZID
  
  # Gene Set Enrichment Analysis (GSEA) #####
  
  # GSEA temporary list
  gsea_output = list()
  
  # Perform GO gene set enrichment analysis using clusterProfiler
  RNGversion("4.2.2")
  set.seed(123)
  gseaGO = suppressWarnings(gseGO(geneList = gene_list,
                 ont = "ALL",
                 OrgDb = "org.Mm.eg.db",
                 keyType = "SYMBOL",
                 pvalueCutoff = 0.05,
                 pAdjustMethod = "BH",
                 minGSSize = 3,
                 maxGSSize = 800,
                 seed = TRUE))
  
  if (exists("gseaGO")) {
    if (nrow(gseaGO) > 0) {
      # gseaGO = clusterProfiler::setReadable(gseaGO, 'org.Mm.eg.db')
      gsea_output[["Gene Ontology"]] = gseaGO@result
      
      addWorksheet(gsea_wb_list[[i]], "GSEA_GO")
      writeData(gsea_wb_list[[i]], "GSEA_GO", as.data.frame(gseaGO@result))
    }
  }
  
  # KEGG
  RNGversion("4.2.2")
  set.seed(123)
  gseaKEGG = suppressWarnings(gseKEGG(geneList = gene_list_entrez,
                     organism = "mmu",
                     keyType = "ncbi-geneid",
                     pvalueCutoff = 0.05,
                     pAdjustMethod = "BH",
                     minGSSize = 3,
                     maxGSSize = 800,
                     seed = TRUE))
  
  if (exists("gseaKEGG")) {
    if (nrow(gseaKEGG) > 0) {
      gseaKEGG = clusterProfiler::setReadable(gseaKEGG, 'org.Mm.eg.db', keyType = "ENTREZID")
      gsea_output[["KEGG"]] = gseaKEGG@result
      
      addWorksheet(gsea_wb_list[[i]], "GSEA_KEGG")
      writeData(gsea_wb_list[[i]], "GSEA_KEGG", as.data.frame(gseaKEGG@result))
    }
  }
  
  # Reactome
  RNGversion("4.2.2")
  set.seed(123)
  gseaReactome = suppressWarnings(gsePathway(geneList = gene_list_entrez,
                            organism = "mouse",
                            pvalueCutoff = 0.05,
                            pAdjustMethod = "BH",
                            minGSSize = 3,
                            maxGSSize = 800,
                            seed = TRUE))
  
  if (exists("gseaReactome")) {
    if (nrow(gseaReactome) > 0) {
      gseaReactome = clusterProfiler::setReadable(gseaReactome, 'org.Mm.eg.db',
                                                  keyType = "ENTREZID")
      gsea_output[["Reactome"]] = gseaReactome@result
      
      addWorksheet(gsea_wb_list[[i]], "GSEA_Reactome")
      writeData(gsea_wb_list[[i]], "GSEA_Reactome", as.data.frame(gseaReactome@result))
    }
  }
  
  # WikiPathways
  RNGversion("4.2.2")
  set.seed(123)
  gseaWP = suppressWarnings(gseWP(geneList = gene_list_entrez,
                 organism = "Mus musculus",
                 pvalueCutoff = 0.05,
                 pAdjustMethod = "BH",
                 minGSSize = 3,
                 maxGSSize = 800,
                 seed = TRUE))
  
  if (exists("gseaWP")) {
    if (nrow(gseaWP) > 0) {
      gseaWP = clusterProfiler::setReadable(gseaWP, 'org.Mm.eg.db',
                                            keyType = "ENTREZID")
      gsea_output[["WikiPathways"]] = gseaWP@result
      
      addWorksheet(gsea_wb_list[[i]], "GSEA_WP")
      writeData(gsea_wb_list[[i]], "GSEA_WP", as.data.frame(gseaWP@result))
    }
  }
  
  # DO analysis (Disease Enrichment)
  RNGversion("4.2.2")
  set.seed(123)
  gseaDO = suppressWarnings(DOSE::gseDO(geneList = gene_list_entrez,
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "BH",
                       minGSSize = 3,
                       maxGSSize = 800,
                       seed = TRUE))
  
  if (exists("gseaDO")) {
    if (nrow(gseaDO) > 0) {
      gseaDO = clusterProfiler::setReadable(gseaDO, 'org.Mm.eg.db',
                                            keyType = "ENTREZID")
      gsea_output[["Disease Ontology"]] = gseaDO@result
      
      addWorksheet(gsea_wb_list[[i]], "GSEA_DO")
      writeData(gsea_wb_list[[i]], "GSEA_DO", as.data.frame(gseaDO@result))
    }
  }
  
  # DO in the Network Cancer Gene
  # http://ncg.kcl.ac.uk/
  # RNGversion("4.2.2")
  # set.seed(123)
  # gseaNCG = suppressWarnings(DOSE::gseNCG(geneList = gene_list_entrez,
  #                       pvalueCutoff = 0.05,
  #                       pAdjustMethod = "BH",
  #                       minGSSize = 3,
  #                       maxGSSize = 800,
  #                       seed = TRUE))
  #
  # if (exists("gseaNCG")) {
  #  if (nrow(gseaNCG) > 0) {
  #    gseaNCG = clusterProfiler::setReadable(gseaNCG, 'org.Mm.eg.db',
  #                                           keyType = "ENTREZID")
  #    gsea_output[["Network Cancer Gene (NCG)"]] = gseaNCG@result
  #    
  #    addWorksheet(gsea_wb_list[[i]], "GSEA_NCG")
  #    writeData(gsea_wb_list[[i]], "GSEA_NCG", as.data.frame(gseaNCG@result))
  #  }
  #}
  
  # DO in the Disease Gene Network (DisGeNET)
  # http://disgenet.org/
  # RNGversion("4.2.2")
  # set.seed(123)
  # gseaDGN = suppressWarnings(DOSE::gseDGN(geneList = gene_list_entrez,
  #                       pvalueCutoff = 0.05,
  #                       pAdjustMethod = "BH",
  #                       minGSSize = 3,
  #                       maxGSSize = 800,
  #                       seed = TRUE))
  
  # if (exists("gseaDGN")) {
  #  if (nrow(gseaDGN) > 0) {
  #    gseaDGN = clusterProfiler::setReadable(gseaDGN, 'org.Mm.eg.db',
  #                                           keyType = "ENTREZID")
  #    gsea_output[["Disease Gene Network (DisGeNET)"]] = gseaDGN@result
  #    
  #    addWorksheet(gsea_wb_list[[i]], "GSEA_DGN")
  #    writeData(gsea_wb_list[[i]], "GSEA_DGN", as.data.frame(gseaDGN@result))
  #  }
  # }
  
  gsea_results[[i]] = gsea_output
  
  # Over-representation analysis (ORA) #####
  
  # ORA temporary list
  ora_output = list()
  
  # GO
  oraGO = suppressWarnings(enrichGO(gene = names(gene_list),
                   ont = "ALL",
                   universe = ashr_dgea[[i]]$Gene.Symbol,
                   OrgDb = "org.Mm.eg.db",
                   keyType = "SYMBOL",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.1,
                   pAdjustMethod = "BH",
                   minGSSize = 3,
                   maxGSSize = 800,
                   pool = TRUE))
  
  if (exists("oraGO")) {
    if (nrow(oraGO) > 0) {
      oraGO_q0.1 = oraGO@result %>% dplyr::filter(qvalue <= 0.1)
      ora_output[["Gene Ontology"]] = oraGO@result
      
      addWorksheet(ora_wb_list[[i]], "ORA_GO")
      writeData(ora_wb_list[[i]], "ORA_GO", as.data.frame(oraGO_q0.1))
    }
  }
  
  # KEGG
  oraKEGG = suppressWarnings(enrichKEGG(gene = names(gene_list_entrez),
                       universe = entrez_universe,
                       organism = "mmu",
                       keyType = "ncbi-geneid",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.1,
                       pAdjustMethod = "BH",
                       minGSSize = 3,
                       maxGSSize = 800))
  
  if (exists("oraKEGG")) {
    if (nrow(oraKEGG) > 0) {
      oraKEGG = clusterProfiler::setReadable(oraKEGG, 'org.Mm.eg.db', keyType = "ENTREZID")
      oraKEGG_q0.1 = oraKEGG@result %>% dplyr::filter(qvalue <= 0.1)
      ora_output[["KEGG"]] = oraKEGG@result
      
      addWorksheet(ora_wb_list[[i]], "ORA_KEGG")
      writeData(ora_wb_list[[i]], "ORA_KEGG", as.data.frame(oraKEGG_q0.1))
    }
  }
  
  # WikiPathways
  oraWP = suppressWarnings(enrichWP(gene = names(gene_list_entrez),
                   universe = entrez_universe,
                   organism = "Mus musculus",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.1,
                   pAdjustMethod = "BH",
                   minGSSize = 3,
                   maxGSSize = 800))
  
  if (exists("oraWP")) {
    if (nrow(oraWP) > 0) {
      oraWP = clusterProfiler::setReadable(oraWP, 'org.Mm.eg.db', keyType = "ENTREZID")
      oraWP_q0.1 = oraWP@result %>% dplyr::filter(qvalue <= 0.1)
      ora_output[["WikiPathways"]] = oraWP@result
      
      addWorksheet(ora_wb_list[[i]], "ORA_WP")
      writeData(ora_wb_list[[i]], "ORA_WP", as.data.frame(oraWP_q0.1))
    }
  }
  
  # Reactome
  oraReactome = suppressWarnings(enrichPathway(gene = names(gene_list_entrez),
                              universe = entrez_universe,
                              organism = "mouse",
                              pvalueCutoff = 0.05,
                              qvalueCutoff = 0.1,
                              pAdjustMethod = "BH",
                              minGSSize = 3,
                              maxGSSize = 800,
                              readable = TRUE))
  
  if (exists("oraReactome")) {
    if (nrow(oraReactome) > 0) {
      oraReactome_q0.1 = oraReactome@result %>% dplyr::filter(qvalue <= 0.1)
      ora_output[["Reactome"]] = oraReactome@result
      
      addWorksheet(ora_wb_list[[i]], "ORA_Reactome")
      writeData(ora_wb_list[[i]], "ORA_Reactome", as.data.frame(oraReactome_q0.1))
    }
  }
  
  # DO
  # oraDO = suppressWarnings(DOSE::enrichDO(gene = names(gene_list_entrez),
  #                       universe = entrez_universe,
  #                       pvalueCutoff = 0.05,
  #                       qvalueCutoff = 0.1,
  #                       pAdjustMethod = "BH",
  #                       minGSSize = 3,
  #                       maxGSSize = 800,
  #                       readable = TRUE))
  
  #if (exists("oraDO")) {
  #  if (nrow(oraDO) > 0) {
  #    oraDO_q0.1 = oraDO@result %>% dplyr::filter(qvalue <= 0.1)
  #    ora_output[["Disease Ontology"]] = oraDO@result
  #    
  #    addWorksheet(ora_wb_list[[i]], "ORA_DO")
  #    writeData(ora_wb_list[[i]], "ORA_DO", as.data.frame(oraDO_q0.1))
  #  }
  #}
  
  # DO: NCG
  #oraNCG = suppressWarnings(DOSE::enrichNCG(gene = names(gene_list_entrez),
  #                         universe = entrez_universe,
  #                         pvalueCutoff = 0.05,
  #                         qvalueCutoff = 0.1,
  #                         pAdjustMethod = "BH",
  #                         minGSSize = 3,
  #                         maxGSSize = 800,
  #                         readable = TRUE))
#  
#  if (exists("oraNCG")) {
#    if (nrow(oraNCG) > 0) {
#      oraNCG_q0.1 = oraNCG@result %>% dplyr::filter(qvalue <= 0.1)
#      ora_output[["Network Cancer Gene (NCG)"]] = oraNCG@result
#      
#      addWorksheet(ora_wb_list[[i]], "ORA_NCG")
#      writeData(ora_wb_list[[i]], "ORA_NCG", as.data.frame(oraNCG_q0.1))
#    }
#  }
  
  # DO: DGN
  #oraDGN = suppressWarnings(DOSE::enrichDGN(gene = names(gene_list_entrez),
  #                         universe = entrez_universe,
  #                         pvalueCutoff = 0.05,
  #                         qvalueCutoff = 0.1,
  #                         pAdjustMethod = "BH",
  #                         minGSSize = 3,
  #                         maxGSSize = 800,
  #                         readable = TRUE))
#  
  #if (exists("oraDGN")) {
  #  if (nrow(oraDGN) > 0) {
  #    oraDGN_q0.1 = oraDGN@result %>% dplyr::filter(qvalue <= 0.1)
  #    ora_output[["Disease Gene Network (DisGeNET)"]] = oraDGN@result
  #    
  #    addWorksheet(ora_wb_list[[i]], "ORA_DGN")
  #    writeData(ora_wb_list[[i]], "ORA_DGN", as.data.frame(oraDGN_q0.1))
  #  }
  #}
  
  ora_results[[i]] = ora_output
  
  # Export the workbooks
  saveWorkbook(ora_wb_list[[i]], paste0("Pathways/", names(ashr_dgea)[i],
                                        "/ORA/ORA_output_", names(ashr_dgea)[i]),
               overwrite = TRUE)
  saveWorkbook(gsea_wb_list[[i]], paste0("Pathways/", names(ashr_dgea)[i],
                                        "/GSEA/GSEA_output_", names(ashr_dgea)[i]),
               overwrite = TRUE)
  
  # Keep track of the loop's progress
  cat(paste("Done with", names(ashr_dgea)[i], "\n"))
  
  # Remove garbage
  rm(list=setdiff(ls(), keepers))
}