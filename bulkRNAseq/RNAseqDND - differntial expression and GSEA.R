library(edgeR)            #3.34.1
library(ggfortify)        #0.4.13
library(ggplot2)          #3.3.5
library(stats)            #4.1.1
library(factoextra)       #1.0.7
library(ggConvexHull)     #0.1.0
library(DOSE)             #3.18.3
library(clusterProfiler)  #4.0.5
library(org.Hs.eg.db)     #3.13.0
library(BiocParallel)     #1.26.2
library(gridExtra)
library(xlsx)             #0.6.5
library(sva)
library(qusage)


# sources----
path <- 'RNAseq/'
expressionPath <- 'matrix_expression_RNAseq_star2306.csv'
sampleInfoPath <- 'infoTable.csv'


#function
createDEGobject <- function(path, expPath, infoPath){
  # take two table- 1. RNA seq with genes on the columns and samples in row 2. information table- age, genotype, feed and sex for each sample
  # create DEG object, order the samples(first - young WT male full)
  # return the DE-object
  
  # read the expression (counts) table and order the names of genes
  exMatrixT <- read.csv(paste0(path, expPath), row.names = 1)
  exMatrix <- data.frame(t(exMatrixT))
  colnames(exMatrix) <- gsub('smp', '', row.names(exMatrixT))

  # read the samples information table, make the parameter as factor and re-factor them(wt, male, full, young as default)
  infoTable <- read.csv(paste0(path, infoPath), row.names = 1, stringsAsFactors=T)
  infoTable$genotype <- relevel(infoTable$genotype, ref = "WT")
  infoTable$sex <- relevel(infoTable$sex, ref = "male")
  infoTable$age <- relevel(infoTable$age, ref = "young")
  infoTable$batch <- factor(infoTable$batch)
  infoTable$sample <- row.names(infoTable)
  
  infoTable <- infoTable[, names(infoTable) %in% c("age", "genotype", "sex", "sample", "batch")]
  infoTable <- infoTable[colnames(exMatrix),] # very important line!! order the info according the RNA-seq data; the first column in count table will be corresponded to the first row in sample table
  
  # create DGE object
  DEobj <- DGEList(exMatrix, samples = infoTable, group=paste(infoTable$age, infoTable$genotype, infoTable$sex),
                    genes = row.names(exMatrix))
  
  DEobj <- DEobj[,order(DEobj$samples$genotype)]
  DEobj <- DEobj[,order(DEobj$samples$age)]
  DEobj <- DEobj[,order(DEobj$samples$sex)]
  DEobj$samples$group <- factor(DEobj$samples$group, levels = unique(DEobj$samples$group))
  return(DEobj)
}
filterNormDEGobject <-function(DEobj){
  # filtering low expressed genes
  # normalized each sample by weight average of non-DE genes using TMM
  
  # filtering out lower expressed genes
  keep <- filterByExpr(DEobj, model.matrix(~genotype+age+sex, data=DEobj$samples))
  print(table(keep))
  DEobj <- DEobj[keep,, keep.lib.sizes=FALSE]
  
  #Normalize data using TMM
  # TMM- 
  DEobj <- calcNormFactors(DEobj, method = 'TMM')
  par(mfrow = c(1, 2))
  boxplot(cpm(DEobj, normalized.lib.sizes = F, log = T), outline=F, col='white', main='Before')
  boxplot(cpm(DEobj, normalized.lib.sizes = T, log = T), outline=F, col='white', main='After')
  return(DEobj)
}
dataExplorer <- function(DEobj, titleG='', logTag = T, clusterSamples = F){
  # dataExplorer accept DGE object
  # return PCA (PC1/2 and PC1/3) with and without labels, and figure of variances against the number of dimensions
  
  #using log2(CPM)
  cpm.data <- log2(cpm(DEobj, normalized.lib.sizes = T)+1)
  cpm.data <- cpm.data[apply(cpm.data, 1, sd) != 0,]
  infoTable <- DEobj$samples[,c('age', 'genotype', 'sex', 'sample', 'batch', 'group')]
  
  #batch correction
  #cpm.data <- removeBatchEffect(cpm.data, batch=infoTable$batch)
  
  
  #*****************PCA************
  #https://cran.r-project.org/web/packages/ggfortify/vignettes/plot_pca.html
  pca.de <- prcomp(t(cpm.data), scale. = T)
  fz <-fviz_eig(pca.de)
  # PCA with labels PC1/2 and PC1/3         #'group'  infoTable$group infoTable$group
  pcaL12 <- autoplot(pca.de, data = infoTable, colour = 'group', x=1, y=2,label = TRUE, title=titleG)+
    geom_convexhull(aes(fill = infoTable$group, color = infoTable$group), alpha=0.2, show.legend = F)+
    ggtitle(titleG)+
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), aspect.ratio=1)
  
  pcaL13 <- autoplot(pca.de, data = infoTable, colour = 'group', x=1, y=3,label = TRUE, title=titleG)+
    geom_convexhull(aes(fill = infoTable$group, color = infoTable$group), alpha=0.2, show.legend = F)+
    ggtitle(titleG)+
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), aspect.ratio=1)
  
  #PCA with dots
  pcaD12 <- autoplot(pca.de, data = infoTable, colour = 'group', shape='age', x=1, y=2,label =F,
                     size = 2.5, title=titleG)+#, shape = 19
    geom_convexhull(aes(fill = infoTable$group, color = infoTable$group, alpha=infoTable$group),
                    show.legend = T)+
    #scale_color_manual(name='Group', values = c('#2E3192', '#A00716', '#1D71BB', '#E30613'))+
    #scale_fill_manual(name='Group', values = c('#2E3192', '#A00716', '#1D71BB', '#E30613'))+
    scale_color_manual(name='Group', values = c('#2E3192', '#1D71BB', '#2E3192','#1D71BB', '#A00716',  '#E30613', '#A00716', '#E30613'))+
    scale_fill_manual(name='Group', values = c('#2E3192', '#1D71BB', '#2E3192','#1D71BB', '#A00716',  '#E30613', '#A00716', '#E30613'))+
    scale_alpha_manual(name='Group', values = c(0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2))+
    scale_shape_manual(name='Group', values = c(19,17))+
    ggtitle(titleG)+
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), aspect.ratio=1)
  
  pcaD13 <- autoplot(pca.de, data = infoTable, colour = 'group', shape='age', x=1, y=3,label =F,
                     size = 2.5, title=titleG)+#, shape = 19
    geom_convexhull(aes(fill = infoTable$group, color = infoTable$group, alpha=infoTable$group),
                    show.legend = T)+
    scale_color_manual(name='Group', values = c('#2E3192', '#1D71BB', '#2E3192','#1D71BB', '#A00716',  '#E30613', '#A00716', '#E30613'))+
    scale_fill_manual(name='Group', values = c('#2E3192', '#1D71BB', '#2E3192','#1D71BB', '#A00716',  '#E30613', '#A00716', '#E30613'))+
    scale_alpha_manual(name='Group', values = c(0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2))+
    scale_shape_manual(name='Group', values = c(19,17))+
    ggtitle(titleG)+
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), aspect.ratio=1)
  
  return(list(knee= fz, PCAlabels12 = pcaL12, PCAlabels13 = pcaL13, PCAdots12 = pcaD12, PCAdots13 = pcaD13))
}


# create DGE object and separate to three datasets- male fully fed, male fasted and female fasted.
DEobj <- createDEGobject(path, expressionPath, sampleInfoPath)

#batch
combatDEobj = DGEList(counts = ComBat_seq(DEobj$counts, batch=DEobj$samples$batch, group=DEobj$samples$group, full_mod=TRUE),
                       samples = DEobj$samples)


combatDEobjY <- filterNormDEGobject(combatDEobj[,combatDEobj$samples$age == 'young'])
combatDEobjO <- filterNormDEGobject(combatDEobj[,combatDEobj$samples$age == 'old'])
combatDEobj <- filterNormDEGobject(combatDEobj)

cpmAll <- log2(cpm(combatDEobj, normalized.lib.sizes = T)+1)
hm <- actualHeatmap(cpmAll, combatDEobj$samples, rownames(cpmAll), 'all') 
pdf(paste0(path, '/figures/heatmap all.pdf'), width=20, height=18) 
draw(hm, column_title = "heatmap")
dev.off()

combatDEobj$samples$group = paste(combatDEobj$samples$age, combatDEobj$samples$genotype, combatDEobj$samples$sex)
AllGraphsCOM <- dataExplorer(combatDEobj, 'all_combat')

pdf(paste0(path, 'figures/PCA RNAseq all13.pdf'), width=10, height=5)
grid.arrange(AllGraphsCOM$PCAlabels13, AllGraphsCOM$PCAdots13, nrow = 1)
dev.off()

youngGraphs <- dataExplorer(DEobjY, 'young')
oldGraphs <- dataExplorer(DEobjO, 'old')
youngGraphsB <- dataExplorer(combatDEobjY, 'youngB')
oldGraphsB <- dataExplorer(combatDEobjO, 'oldB')


#save the PCA graphs
pdf(paste0(path, '/figures/PCA RNAseq.pdf'), width=20, height=18) 
grid.arrange(youngGraphs$PCAlabels12,  oldGraphs$PCAlabels12,  youngGraphsB$PCAlabels12, oldGraphsB$PCAlabels12,
             youngGraphs$PCAlabels13,  oldGraphs$PCAlabels13,  youngGraphsB$PCAlabels13, oldGraphsB$PCAlabels13,
             youngGraphs$knee,         oldGraphs$knee,         youngGraphsB$knee,        oldGraphsB$knee,  nrow = 3)
dev.off()

# To convert NCBI ids to human entrez ids. There are ways to adapt it for Nfur only, but for now I do everything based on human orthologs
humanKfConversion = read.table(paste0(path, "inputFiles/NCBI-Human-orthologs.txt"), head = T, sep = "\t")

# -------DE------
# Differential expression analysis between ages, between genotypes and in the interaction between them in the three datasets

DErankingSaving <- function(DEobj, design, con, test, pathGL, saveRanking=F){
  DEobj <- estimateDisp(DEobj, design)
  fit <- glmQLFit(DEobj, design, robust = T)
  FC <- data.frame(row.names = row.names(DEobj))
  pval <- data.frame(row.names = row.names(DEobj))
  FDR <- data.frame(row.names = row.names(DEobj))
  for (g in 1:length(con)){
    qlf <- glmQLFTest(fit, contrast = con[[g]])
    top <- topTags(qlf, n=Inf)$table
    print(names(con[g]))
    print(paste('+', dim(top[top$FDR < 0.05 & top$logFC > 0,])))
    print(paste('-', dim(top[top$FDR < 0.05 & top$logFC < 0,])))
    print(head(top[top$FDR < 0.05 & top$logFC < 0,]))
    top$mlog10QvalxFC <- -log10(top$FDR) * top$logFC
    FC[[names(con)[g]]] <- top[rownames(FC),]$logFC
    pval[[names(con)[g]]] <- top[rownames(pval),]$PValue
    FDR[[names(con)[g]]] <- top[rownames(FDR),]$FDR
    
    #saving the ranking list
    ranking_genes <- data.frame(Gene = row.names(top), mlog10QvalxFC = top$mlog10QvalxFC, FC=top$logFC, pval=top$PValue, FDR=top$FDR)
    if (saveRanking)
      write.csv(ranking_genes, paste0(pathGL, 'ranking_', names(con[g]), ' ', test[g], '.csv'), row.names = F)
  }
  return(list(fc=FC, fdr=FDR, pval=pval))
}

# all, without batch correction
# "young WT male"  "young dnd male",  "old WT male", "old dnd male", "young WT female","young dnd female","old WT female","old dnd female"
con <- list(WT_M = c(-1,0,1,0,0,0,0,0),
            DND_M = c(0,-1,0,1,0,0,0,0),
            WT_F = c(0,0,0,0,-1,0,1,0),
            DND_F = c(0,0,0,0,0,-1,0,1),
            M_Young = c(-1,1,0,0,0,0,0,0),
            M_Old = c(0,0,-1,1,0,0,0,0),
            F_Young = c(0,0,0,0,-1,1,0,0),
            F_Old = c(0,0,0,0,0,0,-1,1),
            WT_Young = c(-1,0,0,0,1,0,0,0),
            DND_Young = c(0,-1,0,0,0,1,0,0),
            WT_Old = c(0,0,-1,0,0,0,1,0),
            DND_Old = c(0,0,0,-1,0,0,0,1))

tests <- c('age', 'age','age', 'age',
           'genotype', 'genotype', 'genotype', 'genotype',
           'sex', 'sex', 'sex', 'sex')

experimentGroups <- list(ALL = DEobj)

# batch correction
# "young WT male"    "young dnd male"   "young WT female"  "young dnd female"
con <- list(Male = c(-1,1,0,0),
            Female = c(0,0,-1,1),
            c(1,-1,-1,1))

tests <- c('genotype', 'genotype', 'sex+genotype')

experimentGroups <- list(young = combatDEobjY, 
                         old = combatDEobjO)

FCs <- list()
pvals <- list()
FDRs <- list()

for (i in 1:length(experimentGroups)){
  y <- experimentGroups[[i]]
  design <- model.matrix(~0+group, data=y$samples)
  colnames(design) <- gsub('group', '', colnames(design))
  print(colnames(design))
  
  contrast_new <- con
  names(contrast_new) <- paste(names(experimentGroups)[i], names(contrast_new), sep = '_')
  
  de <- DErankingSaving(y, design, contrast_new, tests, paste0(path, 'GO/GeneSets/'), F)
  FCs[[names(experimentGroups)[i]]] = de$fc
  pvals[[names(experimentGroups)[i]]] = de$pval
  FDRs[[names(experimentGroups)[i]]] = de$fdr
}

# batch as co-factor
# "young WT male" "young dnd male" "young WT female" "young dnd female" "batch2" "batch3" "batch4" "batch5" 
con <- list(Male = c(-1,1,0,0,0,0,0),
            Female = c(0,0,-1,1,0,0,0),
            c(1,-1,-1,1,0,0,0))

con <- list(Male = c(-1,1,0,0,0),
            Female = c(0,0,-1,1,0),
            c(1,-1,-1,1,0))

tests <- c('genotype', 'genotype', 'sex+genotype')

experimentGroups <- list(young = DEobjY, 
                         old = DEobjO)

FCs <- list()
pvals <- list()
FDRs <- list()

for (i in 1:length(experimentGroups)){
  y <- experimentGroups[[i]]
  y$samples$batch <- factor(y$samples$batch)
  design <- model.matrix(~0+group+batch, data=y$samples)
  colnames(design) <- gsub('group', '', colnames(design))
  print(colnames(design))
  
  contrast_new <- con
  names(contrast_new) <- paste(names(experimentGroups)[i], names(contrast_new), sep = '_')
  
  de <- DErankingSaving(y, design, contrast_new, tests, paste0(path, 'GO/GeneSets/B_'), T)
  FCs[[names(experimentGroups)[i]]] = de$fc
  pvals[[names(experimentGroups)[i]]] = de$pval
  FDRs[[names(experimentGroups)[i]]] = de$fdr
}

#----GSEA ----
#drug datasets
drugDS = read.gmt('DrugDB/DSigDB_All.gmt')
drugCTD = read.csv('DrugDB/CTD_chem_gene_ixns.csv', skip = 27, header = TRUE)[c('ChemicalName', 'GeneSymbol')] #https://ctdbase.org

GSEAparam <- function(data, savePath, path){ # thanks for Param Priya Singh
  # This input to GSEA is a ranked list of genes. Read the input ranked list.
  
  # Get human ortholog symbols based on the BLAST results file using org.Hs.eg.db package
  # Some Ids will fail to map and will be ignored
  dataH = merge(humanKfConversion, data, by.x = "ncbi", by.y = "Gene") 
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
  #ego3 <- gseGO(geneList     = geneList,
  #              OrgDb        = org.Hs.eg.db,
  #              keyType      = 'SYMBOL',
  #              ont          = c("ALL"),
  #              pvalueCutoff = 1)
  #write.table(ego3, paste0(path, 'GO/Results/', savePath, "_GOGSEA.csv"), sep = ",", quote = T, row.names = F)
  
  #mgsig <- GSEA(geneList     = geneList,
  #              TERM2GENE = msigdbr_t2g, 
  #              pvalueCutoff = 1)
  #write.table(mgsig, paste0(path, 'GO/Results/', savePath, "_OTHER.csv"), sep = ",", quote = T, row.names = F)
  
  drug <- GSEA(geneList      = geneList,
               pvalueCutoff  = 1,
               TERM2GENE     = drugCTD)
  write.table(drug, paste0(path, 'GO/Results/', savePath, "_DRUGCTD.csv"), sep = ",", quote = T, row.names = F)
  
}


# run GSEA on all the ranking lists from the previous step
files <- list.files(path=paste0(path, 'GO/GeneSets/'), pattern="*.csv")
for (i in 1:length(files)){ #c(6,7,8,9) 
  f <- gsub('.csv|ranking_', '', files[i])
  data = read.csv(paste0(path, 'GO/GeneSets/', files[i]))
  ego3 <- GSEAparam(data[,c('Gene', 'mlog10QvalxFC')], f, path)
}

# calculate FC matrix for y/o and wt/het----
cpmAll <- log2(cpm(DEobj, normalized.lib.sizes = T)+1)
cpmY <- log2(cpm(combatDEobjY, normalized.lib.sizes = T)+1)
cpmO <- log2(cpm(combatDEobjO, normalized.lib.sizes = T)+1)

FCall <- merge(FCs$young, FCs$old, by=0)
rownames(FCall) = FCall$Row.names
FCall = FCall[2:length(FCall)]

pvalall <- merge(pvals$young, pvals$old, by=0)
rownames(pvalall) = pvalall$Row.names
pvalall = pvalall[2:length(pvalall)]

FDRall <- merge(FDRs$young, FDRs$old, by=0)
rownames(FDRall) = FDRall$Row.names
FDRall = FDRall[2:length(FDRall)]

# saving the object to the visualization script----
save(DEobj, combatDEobjY, combatDEobjO, cpmY, cpmO, FCall, pvalall, FDRall, file = paste0(path,'outputFiles/data.RDataRev'))

write.csv(FCall, paste0(path, 'outputFiles/differential analysis_FC.csv'))
write.csv(pvalall, paste0(path, 'outputFiles/differential analysis_Pval.csv'))
write.csv(FDRall, paste0(path, 'outputFiles/differential analysis_FDR.csv'))

#####################Mouse data###################
# compare to mouse data: https://doi.org/10.1016/j.cmet.2019.06.018
library(msigdbr)
subC = c('CP:BIOCARTA', 'CP:KEGG', 'CP:REACTOME', 'GO:BP', 'GO:MF')
all_gene_sets = msigdbr(species = "human")
msigdbr_t2g  = all_gene_sets[all_gene_sets$gs_subcat %in% subC, c('gs_name', 'gene_symbol')]
CGP = all_gene_sets[all_gene_sets$gs_subcat %in% c('CGP'), c('gs_name', 'gene_symbol')]

GoGenotype = read.xlsx(paste0(path, '/GO/Results/MergedAnalysis/GO_compare_genotype in rev.xlsx'), sheetName = 'genotype')
GoGenotype$ID = gsub('_', ' ', GoGenotype$ID)
GoGenotype$ID = gsub('GOBP', 'GO', GoGenotype$ID)
GoGenotype$ID = gsub('GOMF', 'GO', GoGenotype$ID)
GoGenotype = GoGenotype[,colnames(GoGenotype) %in% c("ID", colnames(GoGenotype)[grep('NES|qvalue', colnames(GoGenotype))])]
colnames(GoGenotype) = gsub('Female.', 'F', gsub('Male.', 'M', colnames(GoGenotype)))
GoGenotypeM = GoGenotype[,-grep('Female', colnames(GoGenotype))]
GoGenotypeF = GoGenotype[,-grep('Male', colnames(GoGenotype))]


miceData = read.csv(paste0(path, 'inputFiles/mice_data.csv'),row.names = 1)
miceDataFish = merge(miceData, GoGenotype, by.x=0, by.y='ID', all.x=T)
rownames(miceDataFish) = miceDataFish$Row.names


mm = round(max(abs(miceDataFish[rownames(miceDataFish) %in% chosenPathways$ID,])))
col_fun = circlize::colorRamp2(c(-4, 0, 4), c("blue", "white", "red"))

miceDataFish = miceDataFish[,grep('NES', colnames(miceDataFish))]
#miceDataFishM = miceDataFish[,grep('.M', colnames(miceDataFish))]
#miceDataFishF = miceDataFish[,grep('.F', colnames(miceDataFish))]

miceDataFish =  miceDataFish[,-grep('old', colnames(miceDataFish))] #miceDataFishM  #
colnames(miceDataFish)
newOrderCol = c("NES.estradiol.F", "NES.estradiol.M", "NES.young_M", "NES.young_F",  "NES.Acarbose.F", "NES.MR.M", "NES.CR.F", 
                "NES.CR.M", "NES.Rapamycin.M", "NES.Protandim.M", "NES.Protandim.F", "NES.GHRKO.M", 
                "NES.Snell.dwarf.M", "NES.Rapamycin.F", "NES.Acarbose.M")
miceDataFish = miceDataFish[,newOrderCol]

infoMice = data.frame(group=colnames(miceDataFish), sex=sub('.*(?=.$)', '', colnames(miceDataFish), perl=T))
infoMice$sex = factor(gsub('F', 'female', gsub('M', 'male', infoMice$sex)), levels = c('male', 'female'))  
hrorder = actualHeatmap(miceDataFish, infoMice, row.names(miceDataFish), 'mice', convertNames = F, clusterSamples = T, scaleT = F, show_row_names = F, col=col_fun)
hrorderSelect = actualHeatmap(miceDataFish, infoMice, chosenPathways$ID, 'miceSelected', geneOrder = chosenPathways$ID, convertNames = F, clusterRows = F, col=col_fun,  clusterSamples = F, scaleT = F)
write.csv(hrorder, paste0(path, 'outputfiles/hr order heatmap mice.csv'))

pdf(paste0(path, 'figures/heatmap young male.pdf'), width = 6.5, height = 8)
pdf(paste0(path, 'figures/heatmap young male long.pdf'), width = 6.5, height = 30)
show(hrorder)
dev.off()


pdf(paste0(path, 'figures/heatmap selected2.pdf'), width = 6.5, height = 5)
show(hrorderSelect)
dev.off()

plot(miceDataFish$NES.young_F, miceDataFish$NES.young_M)
aaa = miceDataFish[,c('NES.young_F', 'NES.young_M')]
aaa = miceDataFish[,c('NES.young_F', 'NES.estradiol.F')]
aaa = miceDataFish[,c('NES.young_M', 'NES.estradiol.M')]
rownames(aaa) = rownames(miceDataFish)
aaa$dif = abs(aaa[[1]] - aaa[[2]])
aaa$color = 'black'
aaa$color[aaa$dif < 0.3] = 'red'
plot(aaa[[1]], aaa[[2]], col=aaa$color)
rownames(aaa[aaa$color == 'red',])
