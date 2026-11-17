# WTS analysis, retaining ENSG2SYM, degsea and pheatcustom from 250326.R.
# Run from the repository directory. See README.md for input columns and settings.

#┌──────────────────────────────────┐
#│ [01] library install & load      │----
#└──────────────────────────────────┘
library(S4Vectors)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(stringr)
library(DESeq2)
library(matrixStats)
library(reshape2)
library(pheatmap)
library(writexl)

#┌──────────────────────────────────┐
#│ [02] Set seed & options          │----
#└──────────────────────────────────┘
set.seed(48575)
Input.dir="data/WTS"
Run.ID=format(Sys.time(),"%Y%m%d-%H%M%S")
OUT=file.path("results/WTS",Run.ID)
Count.file=file.path(Input.dir,"H.pylori-Gastric_MALToma-WTS-COUNT.tsv")
TPM.file=file.path(Input.dir,"H.pylori-Gastric_MALToma-WTS-TPM.tsv")
Group.file=file.path(Input.dir,"grouptable_new.tsv")
Annotation.file=file.path(Input.dir,"gene_symbols.tsv")
Annotation.label="User-supplied gene-symbol snapshot; record its origin/release here"
Query.Ensembl=FALSE # Explicit fallback when no local snapshot is available.
Ensembl.version=NULL # NULL uses the current service, not a historical release.
Exclude.ID=c("10") # Exclude source ID 10 from DEG analysis.
Count.mode="truncate" # Original as.integer(expected count); not a tximport analysis.
Run.GSEA=TRUE
Run.Hallmark=FALSE # Optional source analysis; needs an identified local TERM2GENE file.
Hallmark.file=file.path(Input.dir,"hallmark_term2gene.tsv")
Reference.dir=file.path(Input.dir,"GSEA_references")
Use.reference.cache=FALSE # TRUE reuses GO_BP.rds / GO_ALL.rds / KEGG.rds.
Make.plots=TRUE
Table.dir=file.path(OUT,"Table")
Plot.dir=file.path(OUT,"Plot")
range=3
col2=colorRampPalette(c("dodgerblue","white","coral"))
options(clusterProfiler.download.method="auto")
if (dir.exists(OUT)) stop("Run directory exists. Change Run.ID.")
if (!Count.mode%in%c("truncate","integer")) stop("Count.mode must be truncate or integer.")
if (Run.Hallmark && !Run.GSEA) stop("Hallmark requires Run.GSEA=TRUE.")
if (!all(file.exists(c(Count.file,TPM.file,Group.file)))) stop("Missing WTS input. See README.md.")
if (!file.exists(Annotation.file) && !Query.Ensembl) stop("Provide gene_symbols.tsv, or explicitly enable Query.Ensembl.")
if (Run.GSEA) {
   library(clusterProfiler)
   library(AnnotationDbi)
   library(org.Hs.eg.db)
   library(enrichplot)
   BiocParallel::register(BiocParallel::SerialParam())
}
dir.create(OUT,recursive=T)
dir.create(Table.dir,recursive=T)
dir.create(Plot.dir,recursive=T)

#┌──────────────────────────────────┐
#│ [03] input data load             │----
#└──────────────────────────────────┘
TPM=read.table(TPM.file,header=T,check.names=F,sep="\t",quote="",comment.char="")
CNT=read.table(Count.file,header=T,check.names=F,sep="\t",quote="",comment.char="")
if (names(CNT)[1]!="gene_id" || !identical(CNT$gene_id,TPM$gene_id) || !identical(names(CNT),names(TPM))) stop("COUNT/TPM gene IDs, sample names or order differ.")
if (anyNA(CNT$gene_id) || anyDuplicated(sub("[.].*$","",CNT$gene_id))) stop("Missing or duplicated gene_id after version-suffix removal.")
# The source files contain gene_id, one empty annotation column, then samples.
if (names(CNT)[2]!="" || any(!is.na(CNT[[2]]) & nzchar(as.character(CNT[[2]]))) || any(!is.na(TPM[[2]]) & nzchar(as.character(TPM[[2]])))) stop("Expected one empty second annotation column.")
rownames(TPM)=TPM$gene_id
rownames(CNT)=CNT$gene_id
TPM=TPM[,-c(1,2),drop=F]
CNT=CNT[,-c(1,2),drop=F]
if (!all(vapply(CNT,is.numeric,logical(1))) || any(!is.finite(as.matrix(CNT))) || any(as.matrix(CNT)<0) || any(as.matrix(CNT)>.Machine$integer.max)) stop("Invalid nonnegative COUNT values.")
if (!all(vapply(TPM,is.numeric,logical(1))) || any(!is.finite(as.matrix(TPM))) || any(as.matrix(TPM)<0)) stop("Invalid TPM values.")

# Original sample-name transformation, now checked rather than inferred silently.
tmp=strsplit(gsub("-",".",colnames(CNT),fixed=T),".",fixed=T)
if (any(lengths(tmp)!=5L)) stop("Expected five source sample-name components.")
Short.names=vapply(tmp,function(x)paste0(x[3],x[5],"_",x[1],x[2]),character(1))
if (anyDuplicated(Short.names)) stop("Sample-name transformation is not unique.")
colnames(TPM)=Short.names
colnames(CNT)=Short.names
Group=read.table(Group.file,header=T,sep="\t",stringsAsFactors=F,check.names=F,quote="",comment.char="")
if (!all(c("Name","sample","HP")%in%names(Group)) || anyNA(Group[,c("Name","sample","HP")]) || anyDuplicated(Group$Name)) stop("Metadata requires unique Name, sample and HP.")
if (!setequal(Group$Name,Short.names)) stop("Metadata and COUNT/TPM samples differ.")
if (any(!Group$sample%in%c("Normal","Tumor","Eradicated")) || any(!Group$HP%in%c("positive","negative"))) stop("Unknown sample/HP group.")
Group=Group%>%mutate(Type=case_when(
   sample=="Normal" & HP=="positive" ~ "PNormal",
   sample=="Tumor" & HP=="positive" ~ "PTumor",
   sample=="Eradicated" & HP=="positive" ~ "PFollowup",
   sample=="Normal" & HP=="negative" ~ "NNormal",
   sample=="Tumor" & HP=="negative" ~ "NTumor",
   sample=="Eradicated" & HP=="negative" ~ "NFollowup"),
   ID=substr(Name,4,5),excluded_ID=ID%in%Exclude.ID,excluded_followup=sample=="Eradicated")
Group$included=!Group$excluded_ID & !Group$excluded_followup
write.table(Group,file.path(OUT,"all_sample_metadata.tsv"),sep="\t",quote=F,row.names=F)
Group=Group%>%filter(included)
Group$Type=factor(Group$Type,levels=c("PNormal","NNormal","PTumor","NTumor","PFollowup","NFollowup"))
Group=Group[order(Group$Type),]
Group$Type=as.character(Group$Type)
TPM=TPM[,Group$Name,drop=F]
CNT=CNT[,Group$Name,drop=F]
col_ann=data.frame(Type=Group$Type,row.names=Group$Name)

#┌──────────────────────────────────┐
#│ [04] ENSEMBL ID -> HGNC SYMBOL   │----
#└──────────────────────────────────┘
ENSG2SYM <- function(ENSG,
                     method,
                     integrated=FALSE){
   if (integrated && file.exists(Annotation.file)) {
      SYM<-read.table(Annotation.file,header=T,sep="\t",quote="",comment.char="",stringsAsFactors=F,check.names=F)
      if (!all(c("gene_id","symbol")%in%names(SYM)) || anyNA(SYM$gene_id)) stop("gene_symbols.tsv requires gene_id and symbol columns.")
      # Retain the first mapping per ENSG, as in the original function.
      SYM=SYM[!duplicated(sub("[.].*$","",SYM$gene_id)),]
      SYM<-SYM[match(sub("[.].*$","",ENSG),sub("[.].*$","",SYM$gene_id)),]
      return(data.frame(ensembl_gene_id=ENSG,hgnc_symbol_AnnotationDbi=NA_character_,
                        hgnc_symbol_biomaRt=SYM$symbol,integrated=SYM$symbol,stringsAsFactors=F))
   }


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


if (!file.exists(Annotation.file)) {
   library(biomaRt)
   mart=useEnsembl(biomart="genes",dataset="hsapiens_gene_ensembl",mirror="useast",version=Ensembl.version)
}
TRANSFORM=ENSG2SYM(rownames(TPM),integrated=TRUE)
TRANSFORM=TRANSFORM[match(sub("[.].*$","",rownames(TPM)),sub("[.].*$","",TRANSFORM$ensembl_gene_id)),]
TRANSFORM$ensembl_gene_id=rownames(TPM)
TRANSFORM$hgnc_symbol_biomaRt[TRANSFORM$hgnc_symbol_biomaRt==""]=NA_character_
write.table(TRANSFORM,file.path(OUT,"gene_mapping_used.tsv"),sep="\t",quote=F,row.names=F,na="NA")
write.table(data.frame(gene_id=TRANSFORM$ensembl_gene_id,symbol=TRANSFORM$hgnc_symbol_biomaRt),file.path(OUT,"gene_symbols_used.tsv"),sep="\t",quote=F,row.names=F,na="NA")
tmp=cbind(TPM,symbol=TRANSFORM$hgnc_symbol_biomaRt)%>%filter(!is.na(symbol))%>%filter(!duplicated(symbol))
rownames(tmp)=tmp$symbol
TPM=tmp[,-ncol(tmp),drop=F]
tmp=cbind(CNT,symbol=TRANSFORM$hgnc_symbol_biomaRt)%>%filter(!is.na(symbol))%>%filter(!duplicated(symbol))
rownames(tmp)=tmp$symbol
CNT=tmp[,-ncol(tmp),drop=F]
if (!nrow(CNT) || !identical(rownames(CNT),rownames(TPM))) stop("No usable or consistent gene symbols.")
TPM_log=mutate_all(TPM,function(x) log2(x+1))
if (Run.Hallmark) {
   if (!file.exists(Hallmark.file)) stop("Provide an identified local Hallmark TERM2GENE snapshot.")
   msigdbr.H=read.table(Hallmark.file,header=T,sep="\t",quote="",comment.char="",stringsAsFactors=F)
   if (!all(c("gs_name","human_gene_symbol")%in%names(msigdbr.H))) stop("Hallmark columns: gs_name, human_gene_symbol.")
}
Files=c(Count.file,TPM.file,Group.file,if(file.exists(Annotation.file)) Annotation.file,if(Run.Hallmark) Hallmark.file)
write.table(data.frame(file=Files,md5=unname(tools::md5sum(Files))),file.path(OUT,"input_hashes.tsv"),sep="\t",quote=F,row.names=F)
dput(mget(c("Exclude.ID","Count.mode","Run.GSEA","Run.Hallmark","Annotation.label","Query.Ensembl","Ensembl.version","Use.reference.cache")),file.path(OUT,"settings_used.R"))

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

PCA.plot=obj%>%
   ggplot(aes(PC1,PC2,color=Type))+
   geom_point()+
   theme_minimal()+
   xlab(paste0("PC1(",tmp$PC1[2]*100,"%)"))+
   ylab(paste0("PC2(",tmp$PC2[2]*100,"%)"))+
   stat_ellipse()
if (Make.plots) ggsave(file.path(Plot.dir,"PCA_TPM.pdf"),plot=PCA.plot,width=7,height=6)
PCA.labels=obj%>%
   ggplot(aes(PC1,PC2,color=Type))+
   geom_text(aes(label=Name))+
   theme_minimal()+
   xlab(paste0("PC1(",tmp$PC1[2]*100,"%)"))+
   ylab(paste0("PC2(",tmp$PC2[2]*100,"%)"))+
   stat_ellipse()

if (Make.plots) ggsave(file.path(Plot.dir,"PCA_TPM_labels.pdf"),plot=PCA.labels,width=7,height=6)
rm(pcr,tmp,obj)
#png("Plot/PCA_TPM_point.png",width = 600,height = 500)
#dev.off()

#┌──────────────────────────────────┐
#│ [06] Function call               │----
#└──────────────────────────────────┘
dir.create(file.path(Table.dir,"Final2"),recursive=T,showWarnings=F)
dir.create(file.path(Plot.dir,"Final2"),recursive=T,showWarnings=F)

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
   if (!all(Metadata$Name%in%colnames(Count_table)) || anyDuplicated(Metadata$Name)) stop("Invalid comparison sample metadata.")
   Metadata=Metadata[order(Metadata$Name),]
   pre=Count_table[,Metadata$Name,drop=F]
   if (Count.mode=="integer" && any(as.matrix(pre)!=floor(as.matrix(pre)))) stop("Fractional counts found in integer mode.")
   # Retain the original expected-count truncation when Count.mode is truncate.
   pre=mutate_all(pre,function(x) as.integer(x))
   pre=pre[matrixStats::rowMedians(as.matrix(pre))>=10,,drop=F]
   if (!nrow(pre)) stop("No genes pass the original median-count filter.")
   if (!setequal(unique(as.character(Metadata$Type)),Contrast_vector)) stop("Comparison groups differ from Contrast_vector.")
   Metadata$Type=factor(Metadata$Type)
   if (paired && anyDuplicated(paste(Metadata$ID,Metadata$Type))) stop("Paired design cannot include duplicate ID/type observations.")
   design_matrix=model.matrix(if(paired) ~ ID + Type else ~ Type,data=Metadata)
   if (qr(design_matrix)$rank<ncol(design_matrix)) stop("Design matrix is not full rank.")
   Prefix=paste0(Explanation,"_",Contrast_vector[1],"_",Contrast_vector[2])
   write.table(Metadata,paste0(Table.dir,"/",PATH,"/metadata",Prefix,".tsv"),sep="\t",quote=F,row.names=F)
   write.table(design_matrix,paste0(Table.dir,"/",PATH,"/design",Prefix,".tsv"),sep="\t",quote=F,col.names=NA)

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
   #write.table(wrt,paste0(Table.dir,"/",PATH,"/DEG_",Explanation,"_",Contrast_vector[1],"_",Contrast_vector[2],".DEG.results.tsv"),sep="\t",quote = F,row.names=F)
   writexl::write_xlsx(wrt,paste0(Table.dir,"/",PATH,"/DEG_",Explanation,"_",Contrast_vector[1],"_",Contrast_vector[2],".DEG.results.xlsx"))

   cut=res%>%as.data.frame()%>%
      mutate(gene_type = case_when(log2FoldChange >= 1 & padj <= 0.05 ~ "up",
                                   log2FoldChange <= -1 & padj <= 0.05 ~ "down",
                                   TRUE ~ "ns"))
   UP=cut%>%filter(gene_type=="up")
   UP=UP%>%mutate(GENE=rownames(UP),.before=baseMean)
   DOWN=cut%>%filter(gene_type=="down")
   DOWN=DOWN%>%mutate(GENE=rownames(DOWN),.before=baseMean)
   write.table(wrt,paste0(Table.dir,"/",PATH,"/DEG",Prefix,".tsv"),sep="\t",quote=F,row.names=F,na="NA")
   write.table(UP,paste0(Table.dir,"/",PATH,"/DEG_up",Prefix,".tsv"),sep="\t",quote=F,row.names=F,na="NA")
   write.table(DOWN,paste0(Table.dir,"/",PATH,"/DEG_down",Prefix,".tsv"),sep="\t",quote=F,row.names=F,na="NA")


   # Volcano plots
   cut=cut%>%
      mutate(gene=rownames(cut))
   cut2=cut%>%
      filter(!is.na(padj))
   if (Make.plots) {
   tmp=cut2%>%
      ggplot(aes(x = log2FoldChange,
                 y = -log10(pmax(padj,.Machine$double.xmin)),
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
   ggsave(plot = tmp,filename = paste0(Plot.dir,"/",PATH,"/Volcano_",Explanation,"_",Contrast_vector[1],"_",Contrast_vector[2],".png"),
          device = "png",
          dpi = 640,
          width = 12,
          height = 12)

   }

   gse=NULL
   kk2=NULL
   msd=NULL
   if (Run.GSEA) {
      # Source GO/Hallmark ranking: positive ranks of LFC, with seeded random ties.
      # Reset the source seed per comparison, so optional plotting cannot alter GSEA draws.
      set.seed(48575)
      saveRDS(.Random.seed,paste0(Table.dir,"/",PATH,"/rng_before_GO",Prefix,".rds"))
      gnl=cut$log2FoldChange
      names(gnl)=rownames(cut)
      gnl=sort(gnl[is.finite(gnl)],decreasing=TRUE)
      gnl=rank(gnl,ties.method="random")
      gnl=sort(gnl,decreasing=TRUE)
      write.table(data.frame(gene=names(gnl),rank=unname(gnl)),paste0(Table.dir,"/",PATH,"/GO_ranking",Prefix,".tsv"),sep="\t",quote=F,row.names=F)
      versions=vapply(c("clusterProfiler","DOSE","fgsea","org.Hs.eg.db","GO.db"),function(x)as.character(packageVersion(x)),character(1))
      dir.create(file.path(OUT,"GSEA_references"),showWarnings=F)
      for (ont in c("BP","ALL")) {
         cache=file.path(Reference.dir,paste0("GO_",ont,Prefix,".rds"))
         if (!file.exists(cache)) cache=file.path(Reference.dir,paste0("GO_",ont,".rds"))
         if (Use.reference.cache) {
            if (!file.exists(cache)) stop("Missing GO reference cache: ",cache)
            ref=readRDS(cache)
            if (!is.null(ref$input_genes) && !setequal(ref$input_genes,names(gnl))) stop("GO cache was captured for another tested gene universe.")
            if (!identical(ref$package_versions,versions)) stop("Reference cache package versions differ; restore its recorded environment.")
            gse=clusterProfiler::GSEA(gnl,TERM2GENE=ref$term2gene,TERM2NAME=ref$term2name,
                                    minGSSize=3,pvalueCutoff=0.05,pAdjustMethod=p.adjust.methods,scoreType="pos",eps=0,verbose=FALSE)
         } else {
            gse=gseGO(geneList=gnl,ont=ont,keyType="SYMBOL",minGSSize=3,pvalueCutoff=0.05,
                      verbose=FALSE,OrgDb=org.Hs.eg.db,pAdjustMethod=p.adjust.methods,scoreType="pos",eps=0)
         }
         saveRDS(gse,paste0(Table.dir,"/",PATH,"/GSEA_GO_",ont,Prefix,".rds"))
         gse.table=if(is.null(gse)) data.frame(ID=character(),Description=character()) else as.data.frame(gse)
         writexl::write_xlsx(gse.table,paste0(Table.dir,"/",PATH,"/GSEA_GO_",ont,Prefix,".xlsx"))
         if (Make.plots && nrow(gse.table)) {
            tmp=enrichplot::dotplot(gse,showCategory=15,split=".sign")+facet_grid(.~.sign)
            ggsave(paste0(Plot.dir,"/",PATH,"/GSEA_GO_",ont,Prefix,".png"),plot=tmp,width=12,height=17,dpi="retina")
         }
         if (Use.reference.cache) {
            file.copy(cache,file.path(OUT,"GSEA_references",basename(cache)),overwrite=TRUE)
            ff=file.path(OUT,"reference_hashes.tsv")
            write.table(data.frame(file=cache,md5=unname(tools::md5sum(cache))),ff,sep="\t",quote=F,row.names=F,col.names=!file.exists(ff),append=file.exists(ff))
         } else if (!is.null(gse) && length(gse@geneSets)) {
            tt=stack(gse@geneSets)[,c("ind","values")]
            names(tt)=c("term","gene")
            labels=data.frame(term=names(gse@geneSets),name=names(gse@geneSets))
            ii=match(gse.table$ID,labels$term)
            labels$name[ii[!is.na(ii)]]=gse.table$Description[!is.na(ii)]
            ref=list(term2gene=tt,term2name=labels,input_genes=names(gnl),package_versions=versions,
                     eps=0,created_utc=format(Sys.time(),tz="UTC"),origin="Captured from this run's tested GO gene sets; not a historical study reference")
            saveRDS(ref,file.path(OUT,"GSEA_references",paste0("GO_",ont,Prefix,".rds")))
         }
      }
      # Preserve the original signed LFC ranking for KEGG.
      ids=bitr(rownames(cut),fromType="SYMBOL",toType="ENTREZID",OrgDb=org.Hs.eg.db)
      dedup_ids=ids[!duplicated(ids$SYMBOL),]
      kegg_gene_list=cut$log2FoldChange[match(dedup_ids$SYMBOL,rownames(cut))]
      names(kegg_gene_list)=dedup_ids$ENTREZID
      kegg_gene_list=sort(kegg_gene_list[is.finite(kegg_gene_list)],decreasing=TRUE)
      write.table(data.frame(gene=names(kegg_gene_list),log2FoldChange=unname(kegg_gene_list)),paste0(Table.dir,"/",PATH,"/KEGG_ranking",Prefix,".tsv"),sep="\t",quote=F,row.names=F)
      saveRDS(.Random.seed,paste0(Table.dir,"/",PATH,"/rng_before_KEGG",Prefix,".rds"))
      cache=file.path(Reference.dir,paste0("KEGG",Prefix,".rds"))
      if (!file.exists(cache)) cache=file.path(Reference.dir,"KEGG.rds")
      if (Use.reference.cache) {
         if (!file.exists(cache)) stop("Missing KEGG reference cache.")
         ref=readRDS(cache)
         if (!is.null(ref$input_genes) && !setequal(ref$input_genes,names(kegg_gene_list))) stop("KEGG cache was captured for another tested gene universe.")
         if (!identical(ref$package_versions,versions)) stop("Reference cache package versions differ; restore its recorded environment.")
         kk2=clusterProfiler::GSEA(kegg_gene_list,TERM2GENE=ref$term2gene,TERM2NAME=ref$term2name,
                                  minGSSize=3,pvalueCutoff=0.05,pAdjustMethod=p.adjust.methods,eps=ref$eps,verbose=FALSE)
         file.copy(cache,file.path(OUT,"GSEA_references",basename(cache)),overwrite=TRUE)
         ff=file.path(OUT,"reference_hashes.tsv")
         write.table(data.frame(file=cache,md5=unname(tools::md5sum(cache))),ff,sep="\t",quote=F,row.names=F,col.names=!file.exists(ff),append=file.exists(ff))
      } else {
         kk2=gseKEGG(geneList=kegg_gene_list,organism="hsa",minGSSize=3,pvalueCutoff=0.05,
                     verbose=FALSE,pAdjustMethod=p.adjust.methods,keyType="ncbi-geneid",eps=1e-10)
      }
      saveRDS(kk2,paste0(Table.dir,"/",PATH,"/GSEA_KEGG",Prefix,".rds"))
      kk2.table=if(is.null(kk2)) data.frame(ID=character(),Description=character()) else as.data.frame(kk2)
      writexl::write_xlsx(kk2.table,paste0(Table.dir,"/",PATH,"/GSEA_KEGG",Prefix,".xlsx"))
      if (Make.plots && nrow(kk2.table)) {
         tmp=enrichplot::dotplot(kk2,showCategory=15,split=".sign")+facet_grid(.~.sign)
         ggsave(paste0(Plot.dir,"/",PATH,"/GSEA_KEGG",Prefix,".png"),plot=tmp,width=12,height=17,dpi="retina")
      }
      if (!Use.reference.cache && !is.null(kk2) && length(kk2@geneSets)) {
         tt=stack(kk2@geneSets)[,c("ind","values")]
         names(tt)=c("term","gene")
         labels=data.frame(term=names(kk2@geneSets),name=names(kk2@geneSets))
         ii=match(kk2.table$ID,labels$term)
         labels$name[ii[!is.na(ii)]]=kk2.table$Description[!is.na(ii)]
         ref=list(term2gene=tt,term2name=labels,input_genes=names(kegg_gene_list),package_versions=versions,
                  eps=1e-10,created_utc=format(Sys.time(),tz="UTC"),origin="Captured from this run's tested KEGG gene sets; not a historical study reference")
         saveRDS(ref,file.path(OUT,"GSEA_references",paste0("KEGG",Prefix,".rds")))
      }
      capture.output(AnnotationDbi::metadata(org.Hs.eg.db),file=file.path(OUT,"OrgDb_metadata.txt"))
   }

   if (H && Run.GSEA) {
      # GSEA(Gene Set Enrichment Analysis) : MsigDB Hallmark genes
      msd=GSEA(gnl,TERM2GENE = msigdbr.H,scoreType="pos",eps=0)
      if (Make.plots && !is.null(msd) && nrow(as.data.frame(msd))) {
         tmp=enrichplot::dotplot(msd,showCategory=15,split=".sign")+facet_grid(.~.sign)
         ggsave(paste0(Plot.dir,"/",PATH,"/GSEA_MSigDB.H",Prefix,".png"),plot=tmp,width=12,height=17,dpi="retina")
      }
      #write.table(msd@result, file=paste0(Table.dir,"/",PATH,"/GSEA_MSigDB.H_",Explanation,"_",Contrast_vector[1],"_",Contrast_vector[2],".tsv"),sep="\t",row.names = F)
      writexl::write_xlsx(if(is.null(msd)) data.frame(ID=character(),Description=character()) else as.data.frame(msd),paste0(Table.dir,"/",PATH,"/GSEA_MSigDB.H_",Explanation,"_",Contrast_vector[1],"_",Contrast_vector[2],".xlsx"))

   }
   Result=mget(c("pre","res","obj","med","wrt","cut","UP","DOWN","gse","kk2"))
   if (H) Result$msd=msd
   attr(Result,"comparison")=paste(Contrast_vector,collapse="_")
   return(Result)
}

pheatcustom=function(
      degsea_result,
      all_samples=F,
      gene_count=40
){
   if (!Make.plots) return(invisible(NULL))
   Genelist=head(degsea_result$cut$gene,gene_count)
   if (!length(Genelist)) return(invisible(NULL))
   Metadata=as.data.frame(SummarizedExperiment::colData(degsea_result$med))
   Samples=if(all_samples) rownames(col_ann) else Metadata$Name
   obj=TPM_log[Genelist,Samples,drop=F]
   obj=t(scale(t(as.matrix(obj))))
   obj=obj[apply(obj,1,function(x)all(is.finite(x))),,drop=F]
   if (!nrow(obj)) return(invisible(NULL))
   anno=col_ann[Samples,,drop=F]
   filename=paste0(Plot.dir,"/Final2/Heatmap_",attr(degsea_result,"comparison"),if(all_samples) "_all" else "",".pdf")
   pheatmap(obj,cluster_rows=nrow(obj)>1,cluster_cols=F,annotation_col=anno,
            legend_breaks=c(-range:range),breaks=seq(-range,range,length.out=11),
            color=col2(11),scale="none",filename=filename,cellwidth=10,cellheight=10)
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
                        PATH = "Final2",H=Run.Hallmark)
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
                       PATH = "Final2",H=Run.Hallmark)
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
                      PATH = "Final2",H=Run.Hallmark)
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
#                       PATH = "Final2",H=Run.Hallmark)
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
                       PATH = "Final2",H=Run.Hallmark)
pheatcustom(NTumor__NNormal,all_samples = F)
pheatcustom(NTumor__NNormal,all_samples = T)


capture.output(sessionInfo(),file=file.path(OUT,"sessionInfo.txt"))
