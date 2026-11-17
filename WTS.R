# WTS analysis. Set the working directory to your private input folder.
# Required inputs retain their original filenames below; they are not public data.
# The original ENSG2SYM function from Toolbox.R is included in this file.
# Main analysis: the four comparisons in sections [07]-[10].
# Later exploratory sections are retained inside if (FALSE) for manual use.
# Install the packages listed below before running the script.

#┌──────────────────────────────────┐
#│ [01] library install & load      │----
#└──────────────────────────────────┘
rm(list=ls())
gc()
#Useless
# .libPaths("path/to/R/library")
# ip = as.data.frame(installed.packages())
# head(ip)
# ip = subset(ip, !grepl("MRO", ip$LibPath))
# ip = ip[!(ip[,"Priority"] %in% c("base", "recommended")),]
# path.lib = unique(ip$LibPath)
# pkgs.to.remove = ip[,1]
# head(pkgs.to.remove)
# sapply(pkgs.to.remove, remove.packages, lib = path.lib)

# install.packages("MASS")
# install.packages("BiocManager")

# BiocManager::install("dplyr")
# BiocManager::install("ggplot2")
# BiocManager::install("ggrepel")
# BiocManager::install("stringr")
# BiocManager::install("biomaRt")
# BiocManager::install("scales")
# BiocManager::install("affy")
# BiocManager::install("S4Vectors")
# BiocManager::install("DESeq2")
# BiocManager::install("clusterProfiler")
# BiocManager::install("reshape2")
# BiocManager::install("viridis")
# BiocManager::install("ggsignif")
# BiocManager::install("enrichplot")
# BiocManager::install("pathview")
# BiocManager::install("psych")
# BiocManager::install("vegan")
# BiocManager::install("msigdbr",force=T)
# BiocManager::install("ggpubr",force=T)
####
library(S4Vectors)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(stringr)
library(biomaRt)
library(scales)
library(affy)
library(DESeq2)
library(matrixStats)
library(clusterProfiler)
library(reshape2)
library(viridis)
library(ggsignif)
library(enrichplot)
library(pathview)
library(psych)
library(vegan)
library(msigdbr)
library(pheatmap)
library(ggpubr)
#library(adespatial)
#library(ade4)
#┌──────────────────────────────────┐
#│ [02] Set seed & options          │----
#└──────────────────────────────────┘
set.seed(48575)
# Working directory: folder containing the private input files.
R.utils::setOption("clusterProfiler.download.method","auto")
enrichR.DBs=c("KEGG_2021_Human","GO_Molecular_Function_2021","MSigDB_Hallmark_2020")
range=3
col2=colorRampPalette(c("dodgerblue", "white", "coral"))
#┌──────────────────────────────────┐
#│ [03] input data load             │----
#└──────────────────────────────────┘
# Nomalized expression level data load
TPM = read.table("H.pylori-Gastric_MALToma-WTS-TPM.tsv",header = T,check.names = F,sep = "\t")
colnames(TPM) = gsub("-", ".", colnames(TPM)) # "-"가 특정 tool에서 colname으로 받아들여지지 않을 때가 있어, 이를 피하기 위해
rownames(TPM) = TPM$gene_id # gene_id를 rowname으로 넣고 불필요한 column 삭제
TPM = TPM[,-c(1,2)]

# Expected count data load
CNT = read.table("H.pylori-Gastric_MALToma-WTS-COUNT.tsv",header = T,check.names = F,sep = "\t")
colnames(CNT) = gsub("-",".",colnames(CNT))
rownames(CNT) = CNT$gene_id
CNT=CNT[,-c(1,2)]

# reprocess colnames
tmp = str_split(colnames(TPM),"\\.")%>%as.data.frame()%>%t()%>%as.data.frame()
colnames(TPM)=paste0(tmp$V3,tmp$V5,"_",tmp$V1,tmp$V2)
colnames(CNT)=paste0(tmp$V3,tmp$V5,"_",tmp$V1,tmp$V2)
rm(tmp)

# Metadata
#Group=read.table("grouptable.tsv",header=T) # Early
#Group=read.table("grouptable2.tsv",header=T) # previous sample metadata
#Group=read.table("grouptable3.tsv",header=T) # previous sample metadata
Group=read.table("grouptable_new.tsv",header=T) # original sample metadata

Group=Group%>%
   mutate(Type=case_when(
      (sample=="Normal" & HP=="positive") ~ "PNormal",
      (sample=="Tumor" & HP=="positive") ~ "PTumor",
      (sample=="Eradicated" & HP=="positive") ~ "PFollowup",
      (sample=="Normal" & HP=="negative") ~ "NNormal",
      (sample=="Tumor" & HP=="negative") ~ "NTumor",
      (sample=="Eradicated" & HP=="negative") ~ "NFollowup"
   )
   )%>%
   mutate(ID=substr(Name,4,5))
Group$Type=factor(Group$Type,levels = c("PNormal","NNormal","PTumor","NTumor","PFollowup","NFollowup"))
Group=Group[order(Group$Type),]
Group$Type=as.character(Group$Type)

# 10 excluding
TPM = TPM[,!substr(colnames(TPM),4,5)%in%c("10")]
CNT = CNT[,!substr(colnames(CNT),4,5)%in%c("10")]
Group = Group%>%filter(ID!="10")

# Exclude followup
Group=Group%>%
   filter(!grepl("Followup",Group$Type))
TPM=TPM[,colnames(TPM)%in%(Group$Name)]
CNT=CNT[,colnames(CNT)%in%(Group$Name)]

table(Group$sample,Group$HP)

col_ann=data.frame(Group[order(Group$Type),]%>%dplyr::select("Name","Type"))
rownames(col_ann) = col_ann$Name
col_ann=col_ann%>%dplyr::select("Type")
# GSEA_Ddatabase
# for (A in c("C1","C2","C3","C4","C5","C6","C7","C8","H")){
#    assign(paste0(A),readRDS(paste0('GSEA/',tolower(A),'.all.v7.5.1.symbols.RDS')))}
# GSEAdb = mget(c("C1","C2","C3","C4","C5","C6","C7","C8","H"))
# rm(list = c("C1","C2","C3","C4","C5","C6","C7","C8","H"),A)
msigdb=msigdbr(species = "Homo sapiens")

#┌──────────────────────────────────┐
#│ [04] ENSEMBL ID -> HGNC SYMBOL   │----
#└──────────────────────────────────┘
# Original mapping function from Toolbox.R
ENSG2SYM <- function(ENSG,
                     method,
                     integrated=FALSE){

   if (integrated==TRUE){
      library(biomaRt)
      library(stringr)
      library(AnnotationDbi)
      library(org.Hs.eg.db)

      if (any(grepl("\\.",ENSG))) {
         ensg_id <- ENSG%>%
            str_replace(pattern = "\\.",replacement = "," )%>%
            str_split(pattern = ",",simplify = TRUE) %>%
            as.data.frame()
         ENSG <- ensg_id$V1
      }
      SYM<-left_join(ENSG2SYM(ENSG,"AnnotationDbi"),
                     ENSG2SYM(ENSG,"biomaRt"),
                     by="ensembl_gene_id",multiple="all")
      SYM$hgnc_symbol_bio[which(SYM$hgnc_symbol_bio=="")]<-NA
      SYM<-SYM%>%
         mutate(integrated=ifelse(is.na(hgnc_symbol_bio)&is.na(hgnc_symbol_Ann),ensembl_gene_id,
                                  ifelse(!is.na(hgnc_symbol_bio)&is.na(hgnc_symbol_Ann),hgnc_symbol_bio,
                                         ifelse(is.na(hgnc_symbol_bio)&!is.na(hgnc_symbol_Ann),hgnc_symbol_Ann,
                                                ifelse(hgnc_symbol_Ann==hgnc_symbol_bio,hgnc_symbol_bio,hgnc_symbol_bio)))))
      SYM<-SYM[!duplicated(SYM$ensembl_gene_id),]
      SYM<-SYM %>%
         mutate(integrated=ifelse(duplicated(SYM[,4]),integrated,ensembl_gene_id))
      colnames(SYM)<-c("ensembl_gene_id", "hgnc_symbol_AnnotationDbi", "hgnc_symbol_biomaRt", "integrated")
      return(SYM)
   }
   if (method=="biomaRt") {
      library(biomaRt)
      library(stringr)

      if (any(grepl("\\.",ENSG))) {
         ensg_id <- ENSG%>%
            str_replace(pattern = "\\.",replacement = "," )%>%
            str_split(pattern = ",",simplify = TRUE) %>%
            as.data.frame()
         ENSG <- ensg_id$V1
      }
      SYM <- getBM(filters="ensembl_gene_id",
                   attributes=c("ensembl_gene_id","hgnc_symbol"),
                   values=ENSG,
                   mart= mart, # loaded biomaRt data
                   useCache = FALSE)
      colnames(SYM)<-c("ensembl_gene_id","hgnc_symbol_biomaRt")
      return(SYM)
   }
   if (method=="AnnotationDbi"){
      library(AnnotationDbi)
      library(org.Hs.eg.db)

      SYM<-data.frame(ensembl_gene_id=ENSG,
                      hgnc_symbol_Ann=mapIds(org.Hs.eg.db,
                                             keys=ENSG,
                                             column="SYMBOL",
                                             keytype="ENSEMBL",
                                             multiVals="first"),row.names = NULL)
      return(SYM)
   }
}


mart=useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl",mirror="useast")
TRANSFORM=ENSG2SYM(rownames(TPM),integrated = TRUE)

tmp=cbind(TPM,symbol=TRANSFORM[,3])%>%
   filter(!is.na(symbol))%>%
   filter(!duplicated(symbol))
rownames(tmp)=tmp$symbol
TPM=tmp[,-c(length(tmp))]

tmp=cbind(CNT,symbol=TRANSFORM[,3])%>%
   filter(!is.na(symbol))%>%
   filter(!duplicated(symbol))
rownames(tmp)=tmp$symbol
CNT=tmp[,-c(length(tmp))]

TPM_log = mutate_all(TPM, function(x) log2(x+1))

#┌──────────────────────────────────┐
#│ [05] PCA Plot                    │----
#└──────────────────────────────────┘
pcr=prcomp(t(TPM))
obj=pcr$x%>%
   as.data.frame()
obj=obj%>%
   mutate(ID=Group$ID[match(rownames(obj),Group$Name)],
          Name=Group$Name[match(rownames(obj),Group$Name)],
          Type=Group$Type[match(rownames(obj),Group$Name)])
tmp=summary(pcr)
tmp=tmp$importance%>%as.data.frame

obj%>%
   ggplot(aes(PC1,PC2,color=Type))+
   geom_point()+
   theme_minimal()+
   xlab(paste0("PC1(",tmp$PC1[2]*100,"%)"))+
   ylab(paste0("PC2(",tmp$PC2[2]*100,"%)"))+
   stat_ellipse()
obj%>%
   ggplot(aes(PC1,PC2,color=Type))+
   geom_text(aes(label=Name))+
   theme_minimal()+
   xlab(paste0("PC1(",tmp$PC1[2]*100,"%)"))+
   ylab(paste0("PC2(",tmp$PC2[2]*100,"%)"))+
   stat_ellipse()

rm(pcr,tmp,obj)
#png("Plot/PCA_TPM_point.png",width = 600,height = 500)
#dev.off()

#┌──────────────────────────────────┐
#│ [06] Function call               │----
#└──────────────────────────────────┘
dir.create("Table/Final2",recursive=T,showWarnings=F)
dir.create("Plot/Final2",recursive=T,showWarnings=F)

#  [function] DEG GSEA
degsea= function(Count_table=NULL,
                 Metadata=NULL,
                 Contrast_vector=NULL,
                 p.adjust.methods="BH",
                 cols=c("up" = "#ffad73", "down" = "#26b3ff", "ns" = "grey"),
                 sizes=c("up" = 2, "down" = 2, "ns" = 1),
                 alphas=c("up" = 1, "down" = 1, "ns" = 0.5),
                 Explanation="",
                 GSEA_GO_ont="BP",
                 paired=FALSE,
                 PATH=NULL,
                 H=TRUE
){
   pre=Count_table[,c(colnames(Count_table)%in%Metadata$Name)]
   pre=mutate_all(pre, function(x) as.integer(x))
   # pre = pre[-which(rowSums(pre) == 0), ]
   pre=pre[matrixStats::rowMedians(pre%>%as.matrix()) >= 10,,drop=F]

   #aaaaaaaaaaaaaaaaaaaaaaa
   Metadata=Metadata[order(Metadata$Name),]
   pre=pre[,order(colnames(pre))]
   #aaaaaaaaaaaaaaaaaaaaaaa

   if (paired) {
      med=DESeqDataSetFromMatrix(countData = pre, colData = Metadata, design = ~ ID + Type, tidy = F)
   } else {
      med=DESeqDataSetFromMatrix(countData = pre, colData = Metadata, design = ~ Type, tidy = F)
   }


   obj=DESeq(med,quiet = T)

   res=results(obj,contrast = c("Type",Contrast_vector),pAdjustMethod = "fdr")
   res = res[order(res$padj),]

   wrt=res%>%as.data.frame()
   wrt=wrt%>%mutate(GENE=rownames(wrt),.before=baseMean)
   #write.table(wrt,paste0("Table/",PATH,"/DEG_",Explanation,"_",Contrast_vector[1],"_",Contrast_vector[2],".DEG.results.tsv"),sep="\t",quote = F,row.names=F)
   writexl::write_xlsx(wrt,paste0("Table/",PATH,"/DEG_",Explanation,"_",Contrast_vector[1],"_",Contrast_vector[2],".DEG.results.xlsx"))

   cut=res%>%as.data.frame()%>%
      mutate(gene_type = case_when(log2FoldChange >= 1 & pvalue <= 0.05 ~ "up",
                                   log2FoldChange <= -1 & pvalue <= 0.05 ~ "down",
                                   TRUE ~ "ns"))
   UP=cut%>%filter(gene_type=="up")
   UP=UP%>%mutate(GENE=rownames(UP),.before=baseMean)
   DOWN=cut%>%filter(gene_type=="down")
   DOWN=DOWN%>%mutate(GENE=rownames(DOWN),.before=baseMean)

   # Volcano plots
   cut=cut%>%
      mutate(gene=rownames(cut))
   cut2=cut%>%
      filter(!is.na(padj))
   tmp=cut2%>%
      ggplot(aes(x = log2FoldChange,
                 y = -log10(padj),
                 fill = gene_type,
                 size = gene_type,
                 alpha = gene_type)) +
      geom_point(shape = 21,
                 colour = "black") + # Specify shape and colour as fixed local parameters
      geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
      geom_vline(xintercept = c(log2(0.5), log2(2)), linetype = "dashed") +
      theme_bw() +
      geom_text_repel(data=head(cut2, 100), aes(label=gene),max.overlaps = Inf)+
      scale_fill_manual(values = cols) + # Modify point colour
      scale_size_manual(values = sizes) + # Modify point size
      scale_alpha_manual(values = alphas) # Modify point transparency
   tmp
   ggsave(plot = tmp,filename = paste0("Plot/",PATH,"/Volcano_",Explanation,"_",Contrast_vector[1],"_",Contrast_vector[2],".png"),
          device = "png",
          dpi = 640,
          width = 12,
          height = 12)

   # GSEA(Gene Set Enrichment Analysis) : Ontology
   gnl = cut$log2FoldChange
   names(gnl) = rownames(cut)
   gnl=na.omit(gnl)
   gnl = sort(gnl, decreasing = TRUE)
   gnl2=data.frame(gnl)%>%mutate(rank=rank(gnl,ties.method = "random"))%>%arrange()
   gnl=as.vector(gnl2$rank)
   names(gnl)=rownames(gnl2)
   gnl = sort(gnl, decreasing = TRUE)

   gse = gseGO(geneList=gnl,
               ont = "BP", #ALL
               keyType = "SYMBOL",
               #nPerm = 10000,
               minGSSize = 3,
               #maxGSSize = 800,
               pvalueCutoff = 0.05,
               verbose = TRUE,
               OrgDb = org.Hs.eg.db,
               pAdjustMethod = p.adjust.methods,
               scoreType="pos",
               eps=0)
   enrichplot::dotplot(gse, showCategory=15, split=".sign") + facet_grid(.~.sign)
   ggsave(paste0("Plot/",PATH,"/GSEA_GO_","BP","_",Explanation,"_",Contrast_vector[1],"_",Contrast_vector[2],".png"),
          device = "png",
          dpi = "retina",
          width = 12,
          height = 17)
   gse = gseGO(geneList=gnl,
               ont = "ALL", #ALL
               keyType = "SYMBOL",
               #nPerm = 10000,
               minGSSize = 3,
               #maxGSSize = 800,
               pvalueCutoff = 0.05,
               verbose = TRUE,
               OrgDb = org.Hs.eg.db,
               pAdjustMethod = p.adjust.methods,
               scoreType="pos",
               eps=0)
   enrichplot::dotplot(gse, showCategory=15, split=".sign") + facet_grid(.~.sign)
   ggsave(paste0("Plot/",PATH,"/GSEA_GO_","ALL","_",Explanation,"_",Contrast_vector[1],"_",Contrast_vector[2],".png"),
          device = "png",
          dpi = "retina",
          width = 12,
          height = 17)
   #write.table(gse@result, file=paste0("Table/",PATH,"/GSEA_GO","_",GSEA_GO_ont,"_",Explanation,"_",Contrast_vector[1],"_",Contrast_vector[2],".tsv"),sep="\t",row.names = F)
   writexl::write_xlsx(gse@result,paste0("Table/",PATH,"/GSEA_GO","_",GSEA_GO_ont,"_",Explanation,"_",Contrast_vector[1],"_",Contrast_vector[2],".xlsx"))

   # GSEA(Gene Set Enrichment Analysis) : KEGG pathway
   cut2=cut[-which(is.na(match(rownames(cut),mappedkeys(org.Hs.egSYMBOL2EG)))),]
   ids=bitr(rownames(cut2), fromType = "SYMBOL", toType = "ENTREZID", OrgDb=org.Hs.eg.db)
   #ids=select(org.Hs.eg.db,keys=rownames(cut2),columns = c("ENTREZID","SYMBOL"),keytype = "SYMBOL")
   dedup_ids = ids[!duplicated(ids[c("SYMBOL")]),]
   kegg_gene_list = cut$log2FoldChange
   kegg_gene_list = kegg_gene_list[match(dedup_ids$SYMBOL,rownames(cut))]
   names(kegg_gene_list) = dedup_ids$ENTREZID
   kegg_gene_list=na.omit(kegg_gene_list)
   kegg_gene_list=sort(kegg_gene_list, decreasing = TRUE)

   kk2 = gseKEGG(geneList     = kegg_gene_list,
                 organism     = "hsa",
                 #nPerm        = 10000,
                 minGSSize    = 3, # 3
                 #maxGSSize    = 800,
                 pvalueCutoff = 0.05,
                 verbose = TRUE,
                 pAdjustMethod= p.adjust.methods,
                 keyType      = "ncbi-geneid")
   #write.table(kk2@result, file=paste0("Table/",PATH,"/GSEA_KEGG","_",Explanation,"_",Contrast_vector[1],"_",Contrast_vector[2],".tsv"),sep="\t",row.names = F)
   writexl::write_xlsx(kk2@result,paste0("Table/",PATH,"/GSEA_KEGG","_",Explanation,"_",Contrast_vector[1],"_",Contrast_vector[2],".xlsx"))

   enrichplot::dotplot(kk2, showCategory = 15, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)
   ggsave(paste0("Plot/",PATH,"/GSEA_KEGG_",Explanation,"_",Contrast_vector[1],"_",Contrast_vector[2],".png"),
          device = "png",
          dpi = "retina",
          width = 12,
          height = 17)
   emapplot(pairwise_termsim(kk2))
   ggsave(paste0("Plot/",PATH,"/GSEA_emapplot_KEGG_",Explanation,"_",Contrast_vector[1],"_",Contrast_vector[2],".png"),
          device = "png",
          dpi = "retina",
          width = 12,
          height = 12)

   if (H) {
      # GSEA(Gene Set Enrichment Analysis) : MsigDB Hallmark genes
      msigdbr.H=msigdb%>%
         as.data.frame()%>%
         filter(gs_cat=="H")%>%
         dplyr::select(gs_name,human_gene_symbol)
      msd=GSEA(gnl,TERM2GENE = msigdbr.H,scoreType="pos",eps=0)
      enrichplot::dotplot(msd, showCategory = 15, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)
      ggsave(paste0("Plot/",PATH,"/GSEA_MSigDB.H_",Explanation,"_",Contrast_vector[1],"_",Contrast_vector[2],".png"),
             device = "png",
             dpi = "retina",
             width = 12,
             height = 17)
      #write.table(msd@result, file=paste0("Table/",PATH,"/GSEA_MSigDB.H_",Explanation,"_",Contrast_vector[1],"_",Contrast_vector[2],".tsv"),sep="\t",row.names = F)
      writexl::write_xlsx(msd@result,paste0("Table/",PATH,"/GSEA_MSigDB.H_",Explanation,"_",Contrast_vector[1],"_",Contrast_vector[2],".xlsx"))

      assign(paste0("Result_",Explanation,"_",Contrast_vector[1],"_",Contrast_vector[2]),mget(c("pre","res","obj","med","wrt","cut","UP","DOWN","gse","kk2","msd")),envir = environment())
   } else {
      assign(paste0("Result_",Explanation,"_",Contrast_vector[1],"_",Contrast_vector[2]),mget(c("pre","res","obj","med","wrt","cut","UP","DOWN","gse","kk2")),envir = environment())
   }

   return(eval(parse(text=paste0("Result_",Explanation,"_",Contrast_vector[1],"_",Contrast_vector[2]))))
}
pheatcustom=function(
      degsea_result,
      all_samples=F,
      gene_count=40
){
   Genelist=degsea_result$cut$gene[1:gene_count]
   if (all_samples) {
      tmp=colnames(TPM_log)
      obj=as.data.frame(t(apply(TPM_log,1,scale))); colnames(obj)=tmp

      pheatmap(obj[Genelist,rownames(col_ann)],
               cluster_rows = T,
               cluster_cols = F,
               annotation_col=col_ann,
               legend_breaks = c(-range:range),
               breaks=seq(-range,range,length.out=11),
               color=col2(11),
               scale = "none")
   } else {
      tmp=colnames(TPM_log)
      Genelist=degsea_result$cut$gene[1:gene_count]
      obj=TPM_log[,Metadata$Name]
      tmp=colnames(obj)
      obj=as.data.frame(t(apply(obj,1,scale))); colnames(obj)=tmp
      anno=col_ann%>%filter(Type%in%Metadata$Type)

      pheatmap(obj[Genelist,rownames(anno)],
               cluster_rows = T,
               cluster_cols = F,
               annotation_col=anno,
               legend_breaks = c(-range:range),
               breaks=seq(-range,range,length.out=11),
               color=col2(11),
               scale = "none",
               cellwidth=10,
               cellheight=10)
   }
}


#┌──────────────────────────────────┐
#│ [07] PNormal vs NNormal          │----
#└──────────────────────────────────┘
Metadata=Group%>%filter(Type%in%c("PNormal","NNormal"))
Metadata
table(Metadata$HP,Metadata$Type)
PNormal__NNormal=degsea(Count_table = CNT, Metadata = Metadata,
                        Contrast_vector = c("PNormal","NNormal"),
                        paired = F,
                        GSEA_GO_ont = "ALL",
                        PATH = "Final2")
pheatcustom(PNormal__NNormal,all_samples = F)
pheatcustom(PNormal__NNormal,all_samples = T)

# tmp<-getwd()
# setwd("Plot/Final2/pathview/PTumor_PNormal/")
# pathlist <- c("hsa04640")
# pathview(gene.data = PTumor__PNormal$kk2@geneList,
#          pathway.id = pathlist,
#          species = "hsa")
# for (i in pathlist) {
#    file.remove(paste0(i,".png"))
#    file.remove(paste0(i,".xml"))
# }
# setwd(tmp)

#┌──────────────────────────────────┐
#│ [08] PTumor vs PNormal           │----
#└──────────────────────────────────┘
Metadata=Group%>%filter(Type%in%c("PTumor","PNormal"))
Metadata
table(Metadata$HP,Metadata$Type)
PTumor__PNormal=degsea(Count_table = CNT, Metadata = Metadata,
                       Contrast_vector = c("PTumor","PNormal"),paired = T,
                       GSEA_GO_ont = "ALL",
                       PATH = "Final2")
pheatcustom(PTumor__PNormal,all_samples = F)
pheatcustom(PTumor__PNormal,all_samples = T)

#┌──────────────────────────────────┐
#│ [09] PTumor vs NTumor            │----
#└──────────────────────────────────┘
Metadata=Group%>%filter(Type%in%c("PTumor","NTumor"))
Metadata
table(Metadata$HP,Metadata$Type)
PTumor__NTumor=degsea(Count_table = CNT, Metadata = Metadata,
                      Contrast_vector = c("PTumor","NTumor"),
                      paired = F,
                      GSEA_GO_ont = "ALL",
                      PATH = "Final2")
pheatcustom(PTumor__NTumor,all_samples = F)
pheatcustom(PTumor__NTumor,all_samples = T)

## opposite
# Metadata=Group%>%filter(Type%in%c("NTumor","PTumor"))
# Metadata
# table(Metadata$HP,Metadata$Type)
# NTumor__PTumor=degsea(Count_table = CNT, Metadata = Metadata,
#                       Contrast_vector = c("NTumor","PTumor"),
#                       paired = F,
#                       GSEA_GO_ont = "ALL",
#                       PATH = "Final2")
# pheatcustom(NTumor__PTumor,all_samples = F)
# pheatcustom(NTumor__PTumor,all_samples = T)


#┌──────────────────────────────────┐
#│ [10] NTumor vs NNormal           │----
#└──────────────────────────────────┘
Metadata=Group%>%filter(Type%in%c("NTumor","NNormal"))
Metadata
table(Metadata$HP,Metadata$Type)
NTumor__NNormal=degsea(Count_table = CNT, Metadata = Metadata,
                       Contrast_vector = c("NTumor","NNormal"),
                       paired = F,
                       GSEA_GO_ont = "ALL",
                       PATH = "Final2")
pheatcustom(NTumor__NNormal,all_samples = F)
pheatcustom(NTumor__NNormal,all_samples = T)

# Original exploratory code below requires its own inputs and group settings.
if (FALSE) {
# Gene fusion vs non
Metadata=Group%>%dplyr::filter(Type=="NTumor")
Metadata$Type=c("noGF","noGF","GF","noGF","GF","noGF")
GF_noGF=degsea(Count_table = CNT, Metadata = Metadata,
               Contrast_vector = c("GF","noGF"),
               paired = F,
               GSEA_GO_ont = "ALL",
               PATH = "Final2")


tmp<-getwd()
setwd("Plot/Final2/pathview/GF_noGF/")
pathlist <- c("hsa05160")
pathview(gene.data = GF_noGF$kk2@geneList,
         pathway.id = pathlist,
         species = "hsa")
for (i in pathlist) {
   file.remove(paste0(i,".png"))
   file.remove(paste0(i,".xml"))
}
setwd(tmp)

#┌──────────────────────────────────┐
#│ [11] Gene fusion visualization   │----
#└──────────────────────────────────┘
library(chimeraviz)

defuse833ke <- system.file(
   "extdata",
   "defuse_833ke_results.filtered.tsv",
   package="chimeraviz")
fusions <- import_defuse(defuse833ke, "hg19",1)


#┌──────────────────────────────────┐
#│ [12] DEG_GENE_identification     │----
#└──────────────────────────────────┘
pheat_col_ann=col_ann%>%
   mutate(sample_type=case_when(
      Type=="NNormal" ~ "Normal",
      Type=="PNormal" ~ "Normal",
      Type=="NTumor" ~ "Tumor",
      Type=="PTumor" ~ "Tumor"
   ))%>%
   mutate("H. pylori"=case_when(
      substr(Type,1,1)=="N" ~ "Negative",
      substr(Type,1,1)=="P" ~ "Positive"
   ))%>%
   dplyr::select("sample_type","H. pylori")
# Heatmap
DEG_gene_list=c(
   head(PNormal__NNormal$cut$gene,10),
   head(PTumor__PNormal$cut$gene,10),
   #head(PTumor__NTumor$cut$gene,10),
   head(NTumor__NNormal$cut$gene,10),
   c("TNFAIP3","EZH2","MYC","EED","CDT1","RPA1")
)
DEG_gene_list=DEG_gene_list[! DEG_gene_list %in% c("IGLV7-46","GPR153","CLCA1","ADA","AQP1")]
DEG_gene_list=DEG_gene_list[!DEG_gene_list%>%duplicated()]

tmp=colnames(TPM_log)
obj=as.data.frame(t(apply(TPM_log,1,scale))); colnames(obj)=tmp
pheatmap(obj[DEG_gene_list,rownames(pheat_col_ann)],
         cluster_rows = T,
         cluster_cols = F,
         annotation_col=pheat_col_ann,
         legend_breaks = c(-range:range),
         breaks=seq(-range,range,length.out=14),
         color=col2(13),
         scale = "none",
         cellwidth=20,
         cellheight=20,
         angle_col=45
)

# chronic inflammatory response
DEG_gene_list=c("CXCL13","SCN11A","ADORA2B","CYP19A1","UNC13D","IL10","IDO1","LTA","FOXP3","S100A8","S100A9","CCL11","NFKBIZ","THBS1","TNF","TNFAIP3","VCAM1","VNN1")
DEG_gene_list=c("AHCY","BDNF","CCL11","CXCL13","GJA1","IL1RN","IDO1","IL10","S100A9","VCAM1","TNF","THBS1","VNN1","TNFAIP3","NF1","EZH2","MYC","EED")
pheat_col_ann2=pheat_col_ann%>%filter(!(`sample_type`!="Normal" & `H. pylori`!="Positive"))
tmp=colnames(TPM_log)
obj=as.data.frame(t(apply(TPM_log,1,scale))); colnames(obj)=tmp
pheatmap(obj[DEG_gene_list,rownames(pheat_col_ann2)],
         cluster_rows = T,
         cluster_cols = F,
         annotation_col=pheat_col_ann2,
         legend_breaks = c(-range:range),
         breaks=seq(-range,range,length.out=11),
         color=col2(11),
         scale = "none",
         cellwidth=20,
         cellheight=20,
         angle_col=45
)

# Intersection
intersect(
   PTumor__NTumor$cut%>%filter(gene_type!="ns")%>%dplyr::select(gene),
   NTumor__NNormal$cut%>%filter(gene_type!="ns")%>%dplyr::select(gene)
)
intersect(
   PTumor__NTumor$gse@result%>%dplyr::select(Description),
   PNormal__NNormal$gse@result%>%dplyr::select(Description)
)

linear_obj=obj[c("NKX6-3","CDT1","RPA1"),]%>%t()%>%as.data.frame()
lm_fit=lm("NKX6-3" ~ "CDT1", obj=linear_obj)
#┌──────────────────────────────────┐
#│ [13] Microbiome detection        │----
#└──────────────────────────────────┘
# diversity
alpha=read.table("bracken/0.05_alpha_diversity.txt",sep="\t",
                 col.names = c("sample",
                               "Berger-parker",
                               "Shannon",
                               "Fisher_index",
                               "Simpson_index",
                               "Simpson_reciprocal_index"))
alpha=alpha[,c(1,3,5)]
alpha=alpha%>%
   mutate(Type=Group$Type[match(alpha$sample,Group$Name)])%>%
   filter(Type%in%c("PNormal","NNormal","PTumor","NTumor"))
alpha=alpha%>%melt(id.vars=c("sample","Type"))
alpha=alpha%>%mutate(Type2=ifelse(substr(Type,1,1)=="P","pos","neg"))


my_comparisons=list(
   c("NNormal","PNormal"),
   c("NTumor","PTumor"),
   c("NTumor","NNormal"),
   c("PTumor","PNormal")
)
alpha%>%
   filter(variable=="Shannon")%>%
   ggplot(aes(x=Type,y=value,group=Type,color=Type))+
   geom_boxplot()+
   theme_bw()+
   xlab("Sample type") + ylab("Shannon")+
   geom_jitter()

alpha%>%
   filter(variable=="Simpson_index")%>%
   ggplot(aes(x=Type,y=value,group=Type,color=Type))+
   geom_boxplot()+
   theme_bw()+
   xlab("Sample type") + ylab("Simpson_index")+
   geom_jitter()

alpha%>%
   filter(variable=="Simpson_index")%>%
   ggplot(aes(x=Type2,y=value,group=Type2,color=Type2))+
   geom_boxplot()+
   theme_bw()+
   xlab("Sample type") + ylab("Simpson_index")+
   geom_jitter()
alpha%>%
   filter(variable=="Shannon")%>%
   ggplot(aes(x=Type2,y=value,group=Type2,color=Type2))+
   geom_boxplot()+
   theme_bw()+
   xlab("Sample type") + ylab("Shannon")+
   geom_jitter()

#┌──────────────────────────────────┐
#│ [14] Junk                        │----
#└──────────────────────────────────┘

# differentially expressed microbiome
cts=read.table("bracken_species_all",row.names=1, sep="\t",header=T,check.names=F)
cts=cts[,grepl("*_num", names(cts))]
colnames(cts)=gsub("^Report_|\\.species.bracken_num$|\\.phylum.bracken_num$|.genus.bracken_num$|.bracken_num$","",colnames(cts))
files_h = list.files(path="k2report/", pattern="*_*", full.names=TRUE, recursive=FALSE)
transcriptome_size=c()
for (i in files_h){
   t=read.table(i,sep="\t", quote= "")
   total_reads=t[1,2] + t[2,2] #get total reads abundance of a sample
   fn=gsub("\\.k2report","",basename(i)) #grep the sample name
   transcriptome_size=c(transcriptome_size,setNames(total_reads,fn)) #add sample name with value to the list lh
}
transcriptome_size=data.frame(sample_name=names(transcriptome_size),Size=transcriptome_size)
coldata=merge(Group2,transcriptome_size,by.x="Name",by.y="sample_name")
tmp=t(cts)%>%as.data.frame()
tmp=tmp%>%mutate(`Helicobacter pylori`=`Helicobacter pylori_BU`+`Helicobacter pylori_C`)%>%t()%>%as.data.frame()
cts=tmp[-match(c("Helicobacter pylori_BU","Helicobacter pylori_C"),rownames(tmp)),]


var=c("H.P_neg_GML","H.P_neg_Nor")
coldata2=coldata[coldata$Type%in%var,]
cts2=cts[,colnames(cts)%in%coldata2$Name]


dds = DESeqDataSetFromMatrix(countData = cts2,
                             colData = coldata2,
                             design= ~ Type)
obj = DESeq(dds)

res=results(obj,contrast=c("Type",var[1],var[2]))
res = res[order(res$padj),]
wrt=res%>%as.data.frame()
wrt=wrt%>%mutate(GENE=rownames(wrt),.before=baseMean)
cut=res%>%as.data.frame()%>%
   mutate(gene_type = case_when(log2FoldChange >= 1 & padj <= 0.05 ~ "up",
                                log2FoldChange <= -1 & padj <= 0.05 ~ "down",
                                TRUE ~ "ns"))
UP=cut%>%filter(gene_type=="up")
UP=UP%>%mutate(GENE=rownames(UP),.before=baseMean)
DOWN=cut%>%filter(gene_type=="down")
DOWN=DOWN%>%mutate(GENE=rownames(DOWN),.before=baseMean)

cut=cut%>%
   mutate(gene=rownames(cut))
cut%>%
   ggplot(aes(x = log2FoldChange,
              y = -log10(padj),
              fill = gene_type,
              size = gene_type,
              alpha = gene_type)) +
   geom_point(shape = 21,
              colour = "black") + # Specify shape and colour as fixed local parameters
   geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
   geom_vline(xintercept = c(log2(0.5), log2(2)), linetype = "dashed") +
   theme_bw() +
   geom_text_repel(data=head(cut, 40), aes(label=gene))+
   scale_fill_manual(values = cols) + # Modify point colour
   scale_size_manual(values = sizes) + # Modify point size
   scale_alpha_manual(values = alphas) # Modify point transparency
ggsave(paste0("Plot/microbiome_test",var[1],var[2],"_.png"),
       device = "png",
       dpi = "retina",
       width = 7,
       height = 7)






#stacked bar
library(forcats)

genus_list=read.table("genuslist.txt")
genus_list=genus_list[!duplicated(genus_list),]
genus_list=genus_list[!grepl("\\.",genus_list)]

cts=read.table("bracken/0.05_genus.tsv",row.names=1, sep="\t",header=T,check.names=F)
cts=cts%>%filter(!grepl("_|MGYG|-|[0-9]",rownames(cts)))
cts=cts[,-c(1:6,length(cts))]

#colnames(cts)=gsub(".cladeReads$|.cladeReads..$|.cladeReads.z.score$|.cladeReads...z.score$","",colnames(cts))
colname_pavian_obj=colnames(cts)
colname_pavian_obj_parsed=paste(substr(colname_pavian_obj,1,9),
                                rep(c("reads","reads_pct","Zscore","Zscore_pct"),37),
                                sep = "_" )
colnames(cts)=colname_pavian_obj_parsed
rm(colname_pavian_obj,colname_pavian_obj_parsed)
Reads=cts[,c(grep("reads$",colnames(cts)))]
Zscore=cts[,c(grep("Zscore$",colnames(cts)))]
MELT=Reads%>%t()%>%as.data.frame()
MELT=MELT%>%
   mutate(sample=rownames(MELT))%>%
   melt(value.name = "value")
MELT$sample=substr(MELT$sample,1,9)

MELT=MELT%>%filter(MELT$variable%in%genus_list)
#MELT=MELT[!grepl("[NT][23].*",MELT$sample),]

MELT=MELT%>%
   mutate(Type=Group$Type[match(MELT$sample,Group$Name)])


MELT=MELT%>%
   mutate(Type=case_when(
      Type=="PNormal" ~ "H.P+ Normal",
      Type=="NNormal" ~ "H.P- Normal",
      Type=="PTumor" ~ "H.P+ Tumor",
      Type=="NTumor" ~ "H.P- Tumor"
      # Type=="Eneg" ~ "H.P- Followup",
      # Type=="Epos" ~ "H.P+ Followup"
   )
   )
MELT2=MELT[order(MELT$value),]
MELT2$variable=as.factor(MELT2$variable)
MELT2%>%
   ggplot(aes(x=sample,y=value)) +
   geom_bar(aes(fill=variable),position="fill",stat="identity") +
   facet_grid(~Type,scales="free_x",space="free_x")+
   theme_bw()+
   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#linear regression?
alpha=read.table("bracken/0.05_alpha_diversity.txt",sep="\t",
                 col.names = c("sample",
                               "Berger-parker",
                               "Shannon",
                               "Fisher_index",
                               "Simpson_index",
                               "Simpson_reciprocal_index"))
alpha=alpha[,c(1,3,5)]
alpha=alpha%>%
   mutate(Type=Group3$Type[match(alpha$sample,Group$Name)])%>%
   filter(Type%in%c("Npos","Nneg","Tpos","Tneg","Epos","Eneg"))
alpha=alpha%>%mutate(Type2=ifelse(substr(Type,2,4)=="pos","pos","neg"))

cts=read.table("bracken/0.05_genus.tsv",row.names=1, sep="\t",header=T,check.names=F)
cts=cts%>%filter(!grepl("_|MGYG|-|[0-9]",rownames(cts)))
cts=cts[,-c(1:6,length(cts))]

#colnames(cts)=gsub(".cladeReads$|.cladeReads..$|.cladeReads.z.score$|.cladeReads...z.score$","",colnames(cts))
colname_pavian_obj=colnames(cts)
colname_pavian_obj_parsed=paste(substr(colname_pavian_obj,1,9),
                                rep(c("reads","reads_pct","Zscore","Zscore_pct"),37),
                                sep = "_" )
colnames(cts)=colname_pavian_obj_parsed
rm(colname_pavian_obj,colname_pavian_obj_parsed)
reads_pct=cts[,c(grep("reads_pct$",colnames(cts)))]%>%t()%>%as.data.frame()
rownames(reads_pct)=substr(rownames(reads_pct),1,9)
reads_pct=reads_pct%>%filter(rownames(reads_pct)%in%alpha$sample)
df=cbind(reads_pct,alpha[,c(2,3)])

df2=df%>%mutate(sample=rownames(df))

df2%>%ggplot(aes(x=Helicobacter,y=Simpson_index))+
   geom_point()+geom_smooth(method='lm')+ theme_bw()
res=lm(Simpson_index ~ Helicobacter,data=df2)
summary(res)

df2%>%ggplot(aes(x=Helicobacter,y=Shannon))+
   geom_point()+geom_smooth(method='lm')+ theme_bw()
res=lm(Shannon ~ Helicobacter,data=df2)
summary(res)

df2%>%ggplot(aes(x=Rothia,y=Simpson_index))+
   geom_point()+geom_smooth(method='lm')+ theme_bw()
res=lm(Simpson_index ~ Rothia,data=df2)
summary(res)

df2%>%ggplot(aes(x=Rothia,y=Shannon))+
   geom_point()+geom_smooth(method='lm')+ theme_bw()
res=lm(Shannon ~ Rothia,data=df2)
summary(res)

res=lm(Shannon ~ Helicobacter + Rothia + Veillonella + Porphyromonas ,data=df2)
summary(res)
library(jtools)
plot_summs(res)




df3=df2[,c(1,3,159,160,161)]%>%mutate(Type=Group3$Type[match(df2$sample,Group$Name)])
df3[order(df3$Shannon),]%>%filter(Type=="Tpos")
df3[order(df3$Shannon),]%>%filter(Type=="Tneg")
df3[order(df3$Shannon),]%>%filter(Type=="Npos")
df3[order(df3$Shannon),]%>%filter(Type=="Tneg")


library(progeny)


tmp=PNormal__NNormal$gse@result$core_enrichment[grep("B cell",PNormal__NNormal$gse@result$Description)]
a=c()
for (i in length(tmp)) {
   a=c(a,strsplit(tmp[i],"/")%>%unlist())
}
B_cell_genelist=a[!duplicated(a)]

tmp=colnames(TPM_log)
obj=as.data.frame(t(apply(TPM_log,1,scale))); colnames(obj)=tmp
pheatmap(obj[B_cell_genelist,rownames(pheat_col_ann)],
         cluster_rows = T,
         cluster_cols = F,
         annotation_col=pheat_col_ann,
         legend_breaks = c(-range:range),
         breaks=seq(-range,range,length.out=11),
         color=col2(11),
         scale = "none",
         cellwidth=10,
         cellheight=10,
         angle_col=45
)

#└──────────────────────────────────┘----



















gse = gseGO(geneList=PTumor__PNormal$gse@geneList,
            ont = "BP", #ALL
            keyType = "SYMBOL",
            #nPerm = 10000,
            minGSSize = 3,
            #maxGSSize = 800,
            pvalueCutoff = 0.05,
            verbose = TRUE,
            OrgDb = org.Hs.eg.db,
            pAdjustMethod = p.adjust.methods,
            scoreType="pos",
            eps=0
)
enrichplot::dotplot(gse, showCategory=30, split=".sign") + facet_grid(.~.sign)
}
