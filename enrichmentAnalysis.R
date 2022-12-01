library(clusterProfiler)  #4.0.5
library(org.Hs.eg.db)     #3.13.0
library(BiocParallel)     #1.26.2
library(KEGGREST)
library(xlsx)
library(futile.logger)
library(VennDiagram)
library(ggplot2)          #3.3.5
library(tidyr)            #1.1.4
library(reshape2)


GSEA_paramP <- function(data, pathGO, save_in){
  # To convert NCBI ids to human entrez ids. This is needed to run the package. There are ways to adapt it for nfur only, but for now I do everything based on human orthologs
  hSymbols = read.table("C:/Users/tehil/Dropbox/Projects/NCBI-Human-orthologs.txt", head = T, sep = "\t")
  
  # This input to GSEA is a ranked list of genes. Read the input ranked list.
  head(data)
  
  # Get human ortholog symbols based on the BLAST results file using org.Hs.eg.db package
  # Some Ids will fail to map and will be ignored
  dataH = merge(hSymbols, data, by.x = "ncbi", by.y = "Gene") 
  entrezIds = bitr(as.character(dataH[,2]), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db") # Get entrez ids for annotation
  dataHE = merge(dataH, entrezIds, by.x = "human", by.y = "SYMBOL") # Get human symbols
  head(dataHE)
  
  # There can be duplicate values because of paralogs, I take average of those for quantitative score
  unique = aggregate(dataHE[,3], list(dataHE$human), mean)
  dataHEU = merge(unique, entrezIds, by.x = "Group.1", by.y = "SYMBOL") 
  colnames(dataHEU) = c("human", "mlog10QvalxFC", "entrez")
  head(dataHEU)
  
  geneList = dataHEU[,2]  # gene list for GO 
  names(geneList) = as.character(dataHEU[,1]) # with entrez ids as names
  
  geneListKegg = geneList # gene list for KEGG
  names(geneListKegg) = as.character(dataHEU[,3]) #  with humna symbols as names
  
  # *** Sort the gene list based on quantitative score in decreasing order. This is critical for GSEA  
  geneList = sort(geneList, decreasing = TRUE)
  geneListKegg = sort(geneListKegg, decreasing = TRUE)
  
  head(geneList)
  tail(geneList)
  
  head(geneListKegg)
  tail(geneListKegg)
  
  # *****************  Now do different enrichment analyses *****************************
  
  # Gene Ontology ------------------------------------------------------------------------------------------------------------------------------------
  ego3 <- gseGO(geneList     = geneList,
                OrgDb        = org.Hs.eg.db,
                keyType      = 'SYMBOL',
                ont          = c("ALL"),
                pvalueCutoff = 1)
  #head(ego3)
  #dotplot(ego3)
  #gseaplot(ego3, geneSetID = 1, by="runningScore",title = edo2$Description[1])
  save(ego3, file= paste0(paste(pathGO, save_in, sep = '/'), "_GOGSEA.Rdata"))
  write.table(ego3, paste0(paste(pathGO, save_in, sep = '/'), "_GOGSEA.csv"), sep = ",", quote = T, row.names = F)
  
  # KEGG Gene Set Enrichment Analysis ----------------------------------------------------------------------------------------------------------------------------
  #kk2 <- gseKEGG(geneList     = geneListKegg,
  #               organism     = 'hsa',
  #               pvalueCutoff = 1)
  #kk2 = setReadable(kk2, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
  
  #save(kk2, file= paste0(paste(pathGO, save_in, sep = '/'), "_KEGG.Rdata"))
  #write.table(kk2, paste0(paste(pathGO, save_in, sep = '/'), "_KEGG.csv"), sep = ",", quote = T, row.names = F)
}

enrichmentClusterProfile <- function(geneListHuman, folder, fileName, backgroundGeneList = 'none') {
  
  # geneListHuman: human gene list by Symbol annotation
  # folder, fileName:  folder name and fileName for saving the files.
  # backgroundGeneList: (optional) background list. if provide use it to the universe.
  
  if (length(backgroundGeneList) == 1) {
    backgroundGeneList <- NULL
    entrezIdBG <- NULL
  } else
    entrezIdBG = bitr(as.character(backgroundGeneList), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")$ENTREZID # Get entrez ids for annotation
  
  # go 
  ego3 <- enrichGO(gene        = geneListHuman,
                   OrgDb        = org.Hs.eg.db,
                   universe     = backgroundGeneList,
                   keyType      = 'SYMBOL',
                   ont          = c("ALL"),
                   pvalueCutoff = 1)
 
  if (length(ego3) > 0)
    if (nrow(ego3) > 0)
      write.table(ego3, paste0(folder, fileName, "_GO.csv"), sep = ",", quote = T, row.names = F)
}

ont <- 'BP'
DBs <- c('GO') 


extract_sigGenes <- function(path, group, ngroup, ont= '', pvalCutoff=0.05) {
  # extract list of significant genes and enrichment table from certain database, kind of test, 
  #kind of ontology(optional) and condition. there is option to choose p-value cutoff
  
  file_name <-  paste0(path, group[ngroup], '_GO.csv') #
  #file_name <-  paste0(path, group[ngroup], '_GOGSEA.csv')
  print(file_name)
  Extable <- read.csv(file_name)
  Extable <- Extable[!is.na(Extable$qvalue) & Extable$ONTOLOGY == ont ,]
  #genes <- Extable[Extable$qvalues < pvalCutoff & Extable$ONTOLOGY == ont ,]$ID # p.adjust
  genes <- Extable[Extable$qvalue < pvalCutoff & Extable$ONTOLOGY == ont ,]$ID # p.adjust
  return(list(table=Extable, genes=genes))
}

threeGroupsCompare <- function(path, group, DB, nameG, ont= '', pvalCutoff=0.05){
  # create table with all the significant enriched pathway in at least one of the given comparisons,
  # drawing Venn diagram with the number of significant pathways in each group (less than 5 groups, if more- without venn diagram)
  # 
  # return table with the go ID, description and for each group:
  # 1. NES (normalized enrichment score) 2.adjust p.value 3.column that contain the core genes in each pathway
  
  tables <- list()
  genes <- list()
  for (i in 1:length(group)){
    GG <- extract_sigGenes(path, group, i, ont)
    tables[[i]] <- GG$table
    genes[[i]] <- GG$genes
  }
  
  names(genes) <- names(group)
  #https://www.rdocumentation.org/packages/VennDiagram/versions/1.6.20/topics/venn.diagram
  if (length(genes) < 6){
    v <- venn.diagram(genes, category.names = names(group), main= paste(DB, ont),
                      filename = NULL)
    ggsave(v, file=paste0(path, 'MergedAnalysis/', nameG, '_venn.pdf'), device = 'pdf', width=5.5, height = 6)
  }
  # one column contain which group have significant enrichment in the pathway
  vennGroups <- get.venn.partitions(genes)
  vennGroups$..set.. <- gsub("\\).*", '', vennGroups$..set..)
  vennGroups$..set.. <- gsub("\\(", '', vennGroups$..set..)
  vennGroups$..set.. <- gsub("n", '#', vennGroups$..set..)
  vennGroups$..set.. <- as.character(lapply(vennGroups$..set.., as.name))
  vennGroups$..set.. <- gsub("n", '+', vennGroups$..set..)
  vennGroups$..set.. <- gsub("#", 'n', vennGroups$..set..)
  
  tableGroups <- data.frame(group=rep('NAN', length(unique(unlist(genes)))), ID=unique(unlist(genes)))
  
  for (gr in 1:nrow(vennGroups)){
    if (length(vennGroups$..values..[[gr]]))
      tableGroups[tableGroups$ID %in% vennGroups$..values..[[gr]],]$group <- vennGroups$..set..[gr]
  }
  print(table(tableGroups$group))
  
  #col names for overlapping tables
  #colNamesGroup <- c('ID', 'NES', 'qvalues', 'core_enrichment')
  colNamesGroup <- c('ID', 'GeneRatio', 'BgRatio', 'qvalue', 'geneID')   # 'p.adjust''Count', 
  N <- length(colNamesGroup) - 1
  
  #col names for describe table
  colNamesGeneral <- c('ONTOLOGY', 'ID', 'Description')
  
  
  tablePvalNES <- merge(tables[[1]][colNamesGroup], tables[[2]][colNamesGroup], all=T, by='ID', suffixes=paste0('.', names(group)[1:2]))
  tablePathwaysNames <- merge(tables[[1]], tables[[2]], by=colNamesGeneral,all=T, suffixes=paste0('.', names(group)[1:2]))
  
  if (length(group) > 2){ # more than 2 groups in the comparisons
    for (k in 3:length(group)){
      tablePvalNES <- merge(tablePvalNES, tables[[k]][colNamesGroup], all=T, by='ID')
      
      colnames(tablePvalNES) <- c(colnames(tablePvalNES)[1:(1+(k-1)*N)], paste0(colnames(tablePvalNES)[(2+(k-1)*N):(1+k*N)], '.', names(group[k])))
      tablePathwaysNames <- merge(tablePathwaysNames, tables[[k]], by=colNamesGeneral, all=T)
    }
  }
  
  # merge the column with the groups are significant in each pathways with the pathway id and description
  tableDescriptionValues <- merge(tableGroups, tablePathwaysNames[colNamesGeneral], by='ID', all.x=T)
  
  #merge everything together
  merge(tableDescriptionValues, tablePvalNES, by='ID', all.x=T)
  
}

sendToCompareSaveXlsx <- function(path, group, DBs, nameG, test, ont=''){
  xlsx_file_name <- paste0(path, 'MergedAnalysis/', DBs, ont,'_', nameG, '.xlsx')
  print(xlsx_file_name)
  #write.xlsx(t(data.frame(DB=DBs, ontology=ont, ageUp=c('young-positive enrichment score'), genotypeUp=c('Het-positive enricment score'))),
  #          file = xlsx_file_name, sheetName='summary', append=FALSE)
  write.xlsx(t(data.frame(DB=DBs, ontology=ont)),
             file = xlsx_file_name, sheetName='summary', append=FALSE)
  
  test_com <- threeGroupsCompare(path, group, DBs, nameG, 'BP', pvalCutoff=0.05) #)
  if(nrow(test_com) > 0)
    write.xlsx(test_com, file = xlsx_file_name, sheetName=nameG, append=TRUE)
  
}

######gonads DND and WT#######
#GSEA
pathGO <- 'C:/Users/tehil/Dropbox/Projects/single-cell gonads/differential_expression/'
files <- gsub('.csv', '', list.files(path=paste0(pathGO, 'GeneSets/'), pattern="*.csv"))
for (i in 1:length(files)){
  data = read.csv(paste0(pathGO, 'GeneSets/', files[i], '.csv'))
  data = data[data$maxEx > 0,]
  data$mlog10QvalxFC = -log10(data$pvals_adj+10^-6) * data$logFC #
  data$Gene = data$NCBI
  ego3 <- GSEA_paramP(data[c('Gene', 'mlog10QvalxFC')], pathGO, paste0('Results/', files[i]))
}

#GO
for (i in 1:length(files)){
  data = read.csv(paste0(pathGO, 'GeneSets/', files[i], '.csv'))
  data = data[data$maxEx > 0,]
  sigGenesP = data[data$pvals_adj < 0.1 & data$logFC > log2(1.5), ]$Human
  sigGenesN = data[data$pvals_adj < 0.1 & data$logFC < -log2(1.5), ]$Human
  sigGenes = data[data$pvals_adj < 0.1 & abs(data$logFC) > log2(1.5), ]$Human
  background = unique(data$Human)
  #enrichmentClusterProfile(unique(sigGenesN),paste0(pathGO, 'GO/neg_') , files[i], background)
  #enrichmentClusterProfile(unique(sigGenesP),paste0(pathGO, 'GO/pos_') , files[i], background)
  #enrichmentClusterProfile(unique(sigGenes),paste0(pathGO, 'GO/') , files[i], background)
  
  print(files[i])
  print(dim(data))
  print(length(sigGenesP))
  print(length(sigGenesN))
}


###combined results to one file####
names(files) <- files

sendToCompareSaveXlsx(paste0(pathGO, 'Results/'), files[1:3], DBs,'F_GSEA_0.1_', ont)
sendToCompareSaveXlsx(paste0(pathGO, 'Results/'), files[4:7], DBs,'M_GSEA_0.1_', ont)

sendToCompareSaveXlsx(paste0(pathGO, 'GO/'), files[1:4], DBs,'F_GO', ont)
sendToCompareSaveXlsx(paste0(pathGO, 'GO/'), files[5:9], DBs,'M_GO', ont)
sendToCompareSaveXlsx(paste0(pathGO, 'GO/'), files, DBs,'ALL_GO', ont)


sendToCompareSaveXlsx(paste0(pathGO, 'GO/pos/'), files[1:3], DBs,'F_pos', ont)
sendToCompareSaveXlsx(paste0(pathGO, 'GO/pos/'), files[4:7], DBs,'M_pos', ont)

sendToCompareSaveXlsx(paste0(pathGO, 'GO/neg/'), files[1:3], DBs,'F_neg', ont)
sendToCompareSaveXlsx(paste0(pathGO, 'GO/neg/'), files[4:7], DBs,'M_neg', ont)


#combined up and down
files <- gsub('_GO.csv', '', list.files(path=paste0(pathGO, 'GO/'), pattern="*.csv"))
names(files) <- files

sendToCompareSaveXlsx(paste0(pathGO, 'GO/'), files[c(1,10,18)], DBs,'F_Fibroblast', ont)
sendToCompareSaveXlsx(paste0(pathGO, 'GO/'), files[c(1,10,18)+1], DBs,'F_Netrophils', ont)
sendToCompareSaveXlsx(paste0(pathGO, 'GO/'), files[c(1,10,18)+2], DBs,'F_OvarianEpithelium', ont)
sendToCompareSaveXlsx(paste0(pathGO, 'GO/'), files[c(1,10,18)+3], DBs,'F_SMC', ont)
sendToCompareSaveXlsx(paste0(pathGO, 'GO/'), files[c(1,10,18)+4], DBs,'M_Fibroblast', ont)
sendToCompareSaveXlsx(paste0(pathGO, 'GO/'), files[c(1,10,18)+5], DBs,'M_Leydig', ont)
sendToCompareSaveXlsx(paste0(pathGO, 'GO/'), files[c(7,24)], DBs,'M_Netrophils', ont)
sendToCompareSaveXlsx(paste0(pathGO, 'GO/'), files[c(7,15,24)+1], DBs,'M_Sertoli', ont)
sendToCompareSaveXlsx(paste0(pathGO, 'GO/'), files[c(7,15,24)+2], DBs,'M_SMC', ont)


sendToCompareSaveXlsx(paste0(pathGO, 'GO/'), files, DBs,'ALLcomparsions', ont)
sendToCompareSaveXlsx(paste0(pathGO, 'GO/'), files[10:26], DBs,'ALL', ont)
sendToCompareSaveXlsx(paste0(pathGO, 'GO/'), files[10:13], DBs,'F_neg', ont)
sendToCompareSaveXlsx(paste0(pathGO, 'GO/'), files[14:18], DBs,'M_neg', ont)
sendToCompareSaveXlsx(paste0(pathGO, 'GO/'), files[19:21], DBs,'F_pos', ont)
sendToCompareSaveXlsx(paste0(pathGO, 'GO/'), files[22:26], DBs,'M_pos', ont)

## GO visualization----

GoCompareEnrichmentGroup = read.xlsx('C:/Users/tehil/Dropbox/Projects/single-cell gonads/differential_expression/GO/MergedAnalysis/GO_ALL_orginal.xlsx', sheetName = 'ALL')
chosenPathways <- read.csv("C:/Users/tehil/Dropbox/Projects/single-cell gonads/differential_expression/pathways2.csv", row.names=1)

GOdots <- function(GoCompareEnrichmentGroup, chosenPathways, group, title="GO", grouplevel=c(), labelss=c(), cirles=c()){
  # extract the pathways and the group 
  GoCompareEnrichmentGroup <- GoCompareEnrichmentGroup[GoCompareEnrichmentGroup$ID %in% chosenPathways$ID, -grep('geneID.', colnames(GoCompareEnrichmentGroup))]
  cnames <- c(colnames(GoCompareEnrichmentGroup)[1:5], colnames(GoCompareEnrichmentGroup)[grep(paste(group, collapse = '|'), colnames(GoCompareEnrichmentGroup))])
  GoCompareEnrichmentGroup <- GoCompareEnrichmentGroup[,cnames]
  
  #order the pathways according FDR or the original order
  pathwaysOrder <- chosenPathways$Description
  #
  
  #select the pathways with significant value at least in one group
  if (length(grep('qvalue', colnames(GoCompareEnrichmentGroup))) > 1){
    checkTable <- GoCompareEnrichmentGroup[, grep('qvalue', colnames(GoCompareEnrichmentGroup))]
    row.names(checkTable)
    checkTable[is.na(checkTable)] = 1
    GoCompareEnrichmentGroup <- GoCompareEnrichmentGroup[rowSums(checkTable < 0.05) > 0,]
    #colnames(GoCompareEnrichmentGroup) <- gsub('p.adjust', 'padjust', colnames(GoCompareEnrichmentGroup))
  }
  
  
  
  meltEnrich <- melt(GoCompareEnrichmentGroup, id.vars = c('ID', 'Description'))
  meltEnrich <- meltEnrich[!meltEnrich$variable %in% c('ONTOLOGY', 'NA.', 'group'),]
  meltEnrich <- separate(data = meltEnrich, col = variable, into = c("measure", "Group"), sep = "\\.")
  meltEnrich$Group <- factor(meltEnrich$Group, levels = unique(meltEnrich$Group))
  #meltEnrich$value <- as.numeric(meltEnrich$value)
  GeneRatio <- meltEnrich[meltEnrich$measure == 'GeneRatio', c('ID', 'Group', 'Description', 'value')]
  GeneRatio <- separate(data = GeneRatio, col = value, into = c("numerator", "denominator"), sep = "\\/")
  GeneRatio$numerator = as.numeric(GeneRatio$numerator)
  GeneRatio$denominator = as.numeric(GeneRatio$denominator)
  GeneRatio$GeneRatio = GeneRatio$numerator / GeneRatio$denominator
  pval <- meltEnrich[meltEnrich$measure == 'qvalue', c('ID', 'Group', 'Description', 'value')]
  colnames(GeneRatio) <- gsub('value', 'GeneRatio', colnames(GeneRatio))
  colnames(pval) <- gsub('value', 'FDR', colnames(pval))
  goVisDF <- merge(pval, GeneRatio, by=c('ID', 'Group', 'Description'))
  
  #goVisDF[is.na(goVisDF)] = 0
  goVisDF$FDR = as.numeric(goVisDF$FDR)
  
  
  #goVisDF$p_adjust[goVisDF$p_adjust > 0.05] = 0
  goVisDF$sig <- c(goVisDF$FDR < 0.05)
  goVisDF$sig <- replace(goVisDF$sig, goVisDF$sig==T, 'Sig')
  goVisDF$sig <- replace(goVisDF$sig, goVisDF$sig==F, 'NS')
  goVisDF$sig <- factor(goVisDF$sig, levels = c('Sig', 'NS'))
  
  goVisDF$Description <- paste0(toupper(substr(goVisDF$Description, 1, 1)), substr(goVisDF$Description, 2, nchar(goVisDF$Description)))
  pathwaysOrder <- paste0(toupper(substr(pathwaysOrder, 1, 1)), substr(pathwaysOrder, 2, nchar(pathwaysOrder)))
  goVisDF$Description <- factor(goVisDF$Description, levels = rev(pathwaysOrder))
  
  #print(goVisDF)
  goVisDF$Group = factor(goVisDF$Group, levels = grouplevel)
  
  one <- ggplot(data = goVisDF, aes(x = '1', y = Description, color= FDR)) + 
  #  geom_point(aes(shape=sig, size = GeneRatio)) +  #-log10()
  #one <- ggplot(data = goVisDF, aes(x = GeneRatio, y = Description, color= FDR )) + 
    geom_point(aes(shape=sig, size = numerator)) +  #-log10()
  #  scale_color_gradient(low = "blue", high = "red") +
    scale_color_gradient(low = "red", high = "blue", limits = c(0,0.05)) +
    scale_shape_manual(values = c(19, 21) , breaks = waiver())+
    theme_bw() +
    ylab("") + 
    xlab("") + 
    #theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
    facet_grid(.~Group) + 
    ggtitle(title)
  
  if (length(labelss) == 0 | length(cirles) == 0)
    return(one)
  
  two <- ggplot(data = goVisDF, aes(x = '1', y = Description, color= FDR, shape = sig, size = GeneRatio)) + 
    geom_point() +  
    scale_color_gradient(low = "blue", high = "red") +
    scale_size_continuous(name='GeneRatio', breaks=labelss, labels=labelss)+
    scale_shape_manual(values = c(19, 21)) +
    guides(size = guide_legend(override.aes = list(shape = cirles)))+
    theme_bw() +
    ylab("") + 
    xlab("") + 
    theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
    facet_grid(.~Group) + 
    ggtitle(title)
  
  return(two)
}

grouplevel = c('neg_M_SMC_WTvsKO', 'neg_M_Fibroblast_WTvsKO', 'neg_M_Sertoli_WTvsKO', 'neg_M_Leydig_WTvsKO',
               'neg_F_SMC_WTvsKO', 'neg_F_Fibroblast_WTvsKO', 'neg_F_Netrophils_WTvsKO', 'neg_F_OvarianEpithelium_WTvsKO')
godots <- GOdots(GoCompareEnrichmentGroup, chosenPathways, c('neg_'), "neg", grouplevel = grouplevel)
pdf('C:/Users/tehil/Dropbox/Projects/single-cell gonads/differential_expression/GOvis.pdf', width = 7, height = 4)
plot(godots)
dev.off()

chosenPathways <- read.csv("C:/Users/tehil/Dropbox/Projects/single-cell gonads/differential_expression/pathwaysPos.csv", row.names=1)
grouplevelP = c('pos_M_SMC_WTvsKO', 'pos_M_Fibroblast_WTvsKO', 'pos_M_Netrophils_WTvsKO', 'pos_M_Sertoli_WTvsKO', 'pos_M_Leydig_WTvsKO',
                'pos_F_SMC_WTvsKO', 'pos_F_Fibroblast_WTvsKO', 'pos_F_Netrophils_WTvsKO', 'pos_F_OvarianEpithelium_WTvsKO')
godotsP <- GOdots(GoCompareEnrichmentGroup, chosenPathways, c('pos_'), "pos", grouplevel = grouplevelP)
pdf('C:/Users/tehil/Dropbox/Projects/single-cell gonads/differential_expression/GOvisPos.pdf', width = 7, height = 3)
plot(godotsP)
dev.off()

godotsPN <- GOdots(GoCompareEnrichmentGroup, chosenPathways, c('_WTvsKO'), "both", grouplevel = c(grouplevel, grouplevelP))


