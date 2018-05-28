library("topGO")

geneID2GO <- readMappings(file = "E:/Data/Laupala/RNAseq_Transcriptome/v057/annotation_nodup_noiso_nopart_GOuniverse.txt")
geneUniverse <- names(geneID2GO) # all genes with GO term

# target genes/transcripts of interest, in this case all 1-LOD scaffolds, see QTL.r for details
QTL_genes<-read.delim("E:/Data/Laupala/CerEukCross/cereuk_QTL_transcripts_annotated_nodup_noiso_nopart.txt") 

 # specify for each gene in the gene universe whether it is present in the list of target genes
geneList <- factor(as.integer(geneUniverse %in% as.vector(QTL_genes$transcript)))
names(geneList) <- geneUniverse

# create GOdata object for Biological Processes
myGOdata_BP <- new("topGOdata", description="CerEuk_QTL", ontology="BP",allGenes=geneList,  annot = annFUN.gene2GO, gene2GO = geneID2GO)

# calculate fisher statistics using parentchild (Grossman et al 2007) algortihm: assess the annotation of the parental nodes of the significant term to correct for the inheritance problem
resultFisher_BP_pc <- runTest(myGOdata_BP, algorithm="parentchild", statistic="fisher")

# calculate FDR for all GO terms with at least one gene in the outlier set
allRes_BP_pc <- GenTable(myGOdata_BP, pcFisher = resultFisher_BP_pc, topNodes=4177)
allRes_BP_pc <- allRes_BP_pc[which(allRes_BP_pc$Significant > 0),]
allRes_BP_pc$fdr_pvalue<-p.adjust(allRes_BP_pc$pcFisher,"fdr",nrow(allRes_BP_pc))

# per linkage group:

pb <- txtProgressBar(min = 0, max = length(levels(QTL_genes$LG)), style = 3) # progress bar

for(i in 1:length(levels(QTL_genes$LG))) {
	
	tempQTL_genes<-QTL_genes[QTL_genes[,"LG"]==levels(QTL_genes$LG)[i],]
	tempgeneList <- factor(as.integer(geneUniverse %in% as.vector(tempQTL_genes$transcript)))
	names(tempgeneList) <- geneUniverse
	tempmyGOdata_BP <- new("topGOdata", description=paste0("CerEuk_QTL_",i), ontology="BP",allGenes=tempgeneList,  annot = annFUN.gene2GO, gene2GO = geneID2GO)
	assign(paste0("LOD1_resultFisher_BP_pc_",i),runTest(tempmyGOdata_BP, algorithm="parentchild", statistic="fisher"))
	assign(paste0("LOD1_allRes_BP_pc_",i),GenTable(tempmyGOdata_BP, pcFisher = get(paste0("LOD1_resultFisher_BP_pc_",i)), topNodes=100))
	setTxtProgressBar(pb, i) # update progress bar
	}
	
close(pb) # close progress bar




