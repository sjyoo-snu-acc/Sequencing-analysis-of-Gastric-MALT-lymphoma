# WES analysis. Set the working directory to your private input folder.
# Required: the ANNOVAR/VEP files named below and sample_info.tsv.
# sample_info.tsv: sample_column, case_number, include_wes (TRUE/FALSE), hp_status.
# It stores the original specimen-to-case map and exclusions; do not publish it.
# Optional declared_fusions.tsv: case_number, gene_symbol (and provenance).
# Those rows are manually declared plotting annotations, not WES fusion calls.

#┌──────────────────────────────────┐
#│ [01] library install & load      │----
#└──────────────────────────────────┘
library(plyr)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(reshape2)
library(maftools)

Sample<-read.table("sample_info.tsv",header=T,sep="\t",check.names=F,stringsAsFactors=F)
if (!all(c("sample_column","case_number","include_wes","hp_status")%in%colnames(Sample))) stop("Check sample_info.tsv columns.")
ANV <- read.table("R.GML_H.pylori_n10_Mutect2_paired_merged.AnnoVar.hg19_multianno.txt",header=T,check.names = F,sep="\t")
VEP <- read.table("R.GML_H.pylori_n10_Mutect2_paired_merged.VEP.result.txt",header=T,check.names = F,sep="\t")
colnames(ANV)[which(duplicated(colnames(ANV)))]<-paste0(colnames(ANV)[duplicated(colnames(ANV))],"_2")

#┌──────────────────────────────────┐
#│ [02] MERGE & simple trimming     │----
#└──────────────────────────────────┘

#GML<-full_join(ANV,VEP,by="ID")
GML<-left_join(ANV,VEP,by="ID")
colnames(GML)


table(GML$Func.refGene)
F1.GML<- GML%>%
   filter(Func.refGene!="intronic" ,
          Func.refGene!="intergenic" ,
          Func.refGene!="downstream" ,
          Func.refGene!="upstream" ,
          Func.refGene!="upstream;downstream" ,
          Func.refGene!="UTR3" ,
          Func.refGene!="UTR5" ,
          Func.refGene!="ncRNA_exonic" ,
          Func.refGene!="ncRNA_exonic;splicing" ,
          Func.refGene!="ncRNA_intronic" ,
          Func.refGene!="ncRNA_splicing")

table(F1.GML$Func.refGene)

mergeF<-F1.GML[,c(1:73,94:110,74:93)]
MELT<-melt(mergeF,id.vars=colnames(mergeF)[1:90],colnames(mergeF[91:110]))
MELTs<-subset(MELT,!grepl("./.:.:.:.:.",value))

tmp<-MELTs%>%
   mutate(index=c(1:dim(MELTs)[1]))%>%
   dplyr::select("INFO","FORMAT","value","index")
tmp1<-tmp%>%filter(FORMAT==names(table(tmp$FORMAT))[1]) # GT:AD:AF:DP:F1R2:F2R1:PGT:PID:PS:SB
tmp11<-ldply(tmp1$value%>%str_split(pattern = ":"),rbind)
tmp12<-tmp1$FORMAT[1]%>%str_split(pattern = ":")
colnames(tmp11)<-tmp12[[1]]
colnames(tmp11)[3]<-"AF2"
tmp13<-cbind(tmp11[,c("GT","AD","AF2","DP","SB","F1R2","F2R1","PGT","PID","PS")],index=tmp1$index)

tmp2<-tmp%>%filter(FORMAT==names(table(tmp$FORMAT))[2]) # GT:AD:AF:DP:F1R2:F2R1:SB
tmp21<-ldply(tmp2$value%>%str_split(pattern = ":"),rbind)
tmp22<-tmp2$FORMAT[1]%>%str_split(pattern = ":")
colnames(tmp21)<-tmp22[[1]]
colnames(tmp21)[3]<-"AF2"
tmp23<-cbind(tmp21[,c("GT","AD","AF2","DP","SB","F1R2","F2R1")],index=tmp2$index)

tmp3<-tmp%>%filter(FORMAT==names(table(tmp$FORMAT))[3]) # GT:AD:AF:DP:F1R2:F2R1:SB:PGT:PID:PS
tmp31<-ldply(tmp3$value%>%str_split(pattern = ":"),rbind)
tmp32<-tmp3$FORMAT[1]%>%str_split(pattern = ":")
colnames(tmp31)<-tmp32[[1]]
colnames(tmp31)[3]<-"AF2"
tmp33<-cbind(tmp31[,c("GT","AD","AF2","DP","SB","F1R2","F2R1","PGT","PID","PS")],index=tmp3$index)

tmp4<-rbind.fill(tmp13,tmp23,tmp33)
tmp5<-tmp4[order(tmp4$index),]

tmp5$AD<-gsub(",\\.,",",",tmp5$AD)
tmp5$AD<-gsub(",\\.$","",tmp5$AD)

tmp6 <- data.frame(str_split_fixed(tmp5$AD, pattern = ",",n = Inf))
tmp66 <- data.frame(do.call('rbind', strsplit(as.character(tmp5$SB), split = ',', fixed = TRUE)))
colnames(tmp6)<-c("AD_ref","AD_alt")
colnames(tmp66)<-c("ref_fw","ref_rv","alt_fw","alt_rv")

tmp7<-cbind(tmp5,tmp6,tmp66)
tmp8<-tmp7[,c("GT","AD","AD_ref","AD_alt","AF2","DP","SB","ref_fw","ref_rv","alt_fw","alt_rv","F1R2","F2R1","PGT","PID","PS","index")]

tmp9<-cbind(MELTs,tmp8)

tmp9[,which(colnames(tmp9)%in%c("DP","AD_ref","AD_alt","ref_fw","ref_rv","alt_fw","alt_rv"))]<-
   sapply(tmp9[,which(colnames(tmp9)%in%c("DP","AD_ref","AD_alt","ref_fw","ref_rv","alt_fw","alt_rv"))], FUN = "as.numeric")

BEF<-tmp9%>%
   mutate(VAF=(AD_alt/(AD_ref+AD_alt)),
          SBratio=(ref_fw+alt_fw)/(ref_rv+alt_rv),
          COSMIC=ifelse(cosmic96_coding!=".","In COSMIC","Not in COSMIC"),
          varID=ID,
          samID=variable)

rm(list=ls(pattern = "tmp"))
RAW<-BEF
#hist(BEF$VAF,breaks = 200,xlim = c(0.1,0.9),ylim = c(0,100))
#hist(as.numeric(BEF$DP),breaks = 1000,xlim = c(0,200))
#hist(BEF$SBratio)

BEF<-RAW%>%filter(
   # AD_ref>=5,
   # AD_alt>=5,
   # alt_fw>=5,
   # alt_rv>=5,
   # DP>=20,
   # SBratio<=10,
   # SBratio>=0.1,
   #VAF>0.01,
   # VAF<0.99,
   #AF_eas<0.01,
   #as.numeric(VAF)>0.00001,
   #CANONICAL=="YES",
   grepl("1",RAW$GT),
   #BIOTYPE%in%c("protein_coding"),
   #ExonicFunc.refGene%in%c("nonsynonymous SNV","frameshift insertion","frameshift deletion","stopgain","startloss")
   #!(ExonicFunc.refGene%in%c(".","unknown"))
)

#write_tsv(BEF, file="BeforeFiltering.tsv")

#┌──────────────────────────────────┐
#│ [03] Filtering                   │----
#└──────────────────────────────────┘

FIN<-BEF%>%filter(
   DP>=50, #50?
   #AD_ref>=40, #50?
   AD_alt>=5,
   #alt_fw>=10,
   #alt_rv>=10,
   SBratio<=4,
   SBratio>=0.25,
   CANONICAL=="YES",
   BIOTYPE%in%c("protein_coding"),
   grepl("haematopoietic_and_lymphoid_tissue",BEF$cosmic96_coding) |   grepl("stomach",BEF$cosmic96_coding),
   #COSMIC=="In COSMIC",
   #grepl("haematopoietic_and_lymphoid_tissue",BEF$cosmic96_coding),
   ExonicFunc.ensGene%in%c("nonsynonymous SNV","frameshift insertion","frameshift deletion","stopgain","startloss","stoploss"),
   !grepl("normal_artifact",BEF$FILTER),
   #FILTER=="PASS",
   #(AF_eas<0.01 | AF_eas=="."),
   AF_eas<0.1,
   # ExAC_EAS<0.1
   VAF>=0.05,
   VAF<0.40
)%>%mutate(ID2=paste0(varID,".",samID))

A<-table(FIN[,c("samID","varID")])%>%as.data.frame()
A<-A%>%mutate(ID2=paste0(varID,".",samID))
FIN<-FIN[match(A$ID2,FIN$ID2),]
FIN<-FIN[!is.na(FIN$ID2),]
vc_cols = RColorBrewer::brewer.pal(n = 9, name = 'Paired')
names(vc_cols) = c('Frame_Shift_Del','Missense_Mutation','Nonsense_Mutation','Multi_Hit','Frame_Shift_Ins','In_Frame_Ins','Splice_Site','In_Frame_Del','Gene_Fusion')



tmp<-FIN[,c("variable","Gene.refGene","ExonicFunc.ensGene","Func.ensGene","ID","Protein_position", # 1 - 6
            "Amino_acids","cosmic96_coding","COSMIC","INFO","CLNSIG_220425","BIOTYPE","FILTER", # 7 - 12
            "value","FORMAT","Chr","Start","End","Ref","Alt", # 13 - 19
            "Func.refGene","GeneDetail.refGene","AAChange.refGene","ExonicFunc.refGene")] # 20 - 23
colnames(tmp)[1] <- "Tumor_Sample_Barcode"

tmp<-tmp%>%
   filter(!Tumor_Sample_Barcode%in%Sample$sample_column[!Sample$include_wes])

tmp$Tumor_Sample_Barcode<-Sample$case_number[match(tmp$Tumor_Sample_Barcode,Sample$sample_column)]
if (anyNA(tmp$Tumor_Sample_Barcode)) stop("A sample is missing from sample_info.tsv.")

# Original manual plotting placeholders; genomic fields are not called variants.
Fusion<-data.frame(case_number=character(),gene_symbol=character())
if (file.exists("declared_fusions.tsv")) {
   Fusion<-read.table("declared_fusions.tsv",header=T,sep="\t",check.names=F,stringsAsFactors=F)
   for (i in seq_len(nrow(Fusion))) {
      tmp<-rbind(tmp,c(Fusion$case_number[i],Fusion$gene_symbol[i],"Gene_Fusion","exonic","chr1_117150564_T_C","1","T/A","ID=","In COSMIC","-",".","protein_coding",".",".",".","chr1","1","1","T","C","exonic",".",".","Gene_Fusion"))
   }
}

write_tsv(tmp, file="BeforeFiltering.tsv")
MAF<-annovarToMaf("BeforeFiltering.tsv",tsbCol="Tumor_Sample_Barcode")
MAF$Variant_Classification[MAF$Hugo_Symbol%in%Fusion$gene_symbol]<-"Gene_Fusion"
before<-read.maf(maf = MAF,clinicalData = NULL,vc_nonSyn = c("Missense_Mutation","Nonsense_Mutation","Translation_Start_Site","Gene_Fusion"))
before@clinical.data<-before@clinical.data%>%mutate(
   H.pylori=Sample$hp_status[match(Tumor_Sample_Barcode,Sample$case_number)])


oncoplot(maf=before,
         cohortSize = 9,
         top = 2000,
         colors=vc_cols,
         gene_mar = 10,
         clinicalFeatures = 'H.pylori',
         sortByAnnotation = T,showTumorSampleBarcodes = T,showPct = FALSE)
