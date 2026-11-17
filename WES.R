# WES analysis, retaining the analysis order in Trimmedddd.R.
# Run from the repository directory after preparing data/WES (README.md).

#┌──────────────────────────────────┐
#│ [01] library install & load      │----
#└──────────────────────────────────┘
library(dplyr)
library(stringr)
library(ggplot2)
library(reshape2)
library(maftools)

Input.dir<-"data/WES"
Run.ID<-format(Sys.time(),"%Y%m%d-%H%M%S")
OUT<-file.path("results/WES",Run.ID)
Annovar.file<-file.path(Input.dir,"R.GML_H.pylori_n10_Mutect2_paired_merged.AnnoVar.hg19_multianno.txt")
VEP.file<-file.path(Input.dir,"R.GML_H.pylori_n10_Mutect2_paired_merged.VEP.result.txt")
Sample.file<-file.path(Input.dir,"sample_info.tsv")
Fusion.file<-file.path(Input.dir,"declared_fusions.tsv")
AF.column<-"AF_eas" # First ANNOVAR block; AF_eas_2 is the second block.
AF.cutoff<-0.1      # Allele fraction (10%), retaining the source threshold.
Keep.missing.AF<-TRUE # The source retained missing population annotations ('.').

if (dir.exists(OUT)) stop("Run directory exists. Change Run.ID.")
if (!all(file.exists(c(Annovar.file,VEP.file,Sample.file)))) stop("Missing WES input. See README.md.")
Sample<-read.table(Sample.file,header=T,sep="\t",check.names=F,stringsAsFactors=F,quote="",comment.char="")
if (!all(c("sample_column","case_number","include_wes","hp_status")%in%colnames(Sample))) stop("Check sample_info.tsv columns.")
for (i in intersect(c("include_wes","include_plot"),names(Sample))) {
   if (!all(as.character(Sample[[i]])%in%c("TRUE","FALSE"))) stop("Inclusion flags must be TRUE/FALSE.")
   Sample[[i]]<-as.logical(Sample[[i]])
}
if (!"include_plot"%in%names(Sample)) Sample$include_plot<-Sample$include_wes
if (anyNA(Sample$sample_column) || anyNA(Sample$case_number)) stop("Use empty text for unmapped controls, not NA.")
if (anyDuplicated(Sample$sample_column[nzchar(Sample$sample_column)])) stop("Duplicated sample_column.")
if (any(!nzchar(Sample$case_number[Sample$include_wes | Sample$include_plot]))) stop("Included cases require case_number.")
if (anyDuplicated(Sample$case_number[Sample$include_plot])) stop("Plot cases must be unique.")
if (any(!Sample$hp_status[Sample$include_plot]%in%c("positive","negative","unknown"))) stop("Invalid hp_status.")
ANV<-read.table(Annovar.file,header=T,check.names=F,sep="\t",quote="",comment.char="",stringsAsFactors=F)
VEP<-read.table(VEP.file,header=T,check.names=F,sep="\t",quote="",comment.char="#",stringsAsFactors=F)
colnames(ANV)[which(duplicated(colnames(ANV)))]<-paste0(colnames(ANV)[duplicated(colnames(ANV))],"_2")
if (anyDuplicated(names(ANV)) || anyDuplicated(names(VEP))) stop("Unresolved duplicated annotation headers.")
if (!AF.column%in%names(ANV)) stop("AF.column is absent from ANNOVAR.")
dir.create(OUT,recursive=T)
Files<-c(Annovar.file,VEP.file,Sample.file,if(file.exists(Fusion.file)) Fusion.file)
write.table(data.frame(file=Files,md5=unname(tools::md5sum(Files))),file.path(OUT,"input_hashes.tsv"),sep="\t",row.names=F,quote=F)
write.table(Sample,file.path(OUT,"sample_info_used.tsv"),sep="\t",row.names=F,quote=F)
dput(list(AF.column=AF.column,AF.cutoff=AF.cutoff,Keep.missing.AF=Keep.missing.AF),file.path(OUT,"settings_used.R"))

#┌──────────────────────────────────┐
#│ [02] MERGE & simple trimming     │----
#└──────────────────────────────────┘

#GML<-full_join(ANV,VEP,by="ID")
GML<-left_join(ANV,VEP,by="ID",relationship="many-to-many")
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

# Select genotype columns by metadata, independently of annotation column order.
Sample.columns<-Sample$sample_column[nzchar(Sample$sample_column)]
if (!all(Sample.columns%in%names(F1.GML))) stop("A sample column is absent from the annotation table.")
mergeF<-F1.GML
MELT<-reshape2::melt(mergeF,id.vars=setdiff(names(mergeF),Sample.columns),measure.vars=Sample.columns)
MELTs<-MELT[!is.na(MELT$value) & !grepl("^\\./\\.|^\\.\\|\\.",MELT$value),]
if (!nrow(MELTs)) stop("No nonmissing genotype records.")

# FORMAT tags can appear in different orders, with or without optional tags.
tmp<-vector("list",length(unique(MELTs$FORMAT)))
for (i in seq_along(tmp)) {
   idx<-which(MELTs$FORMAT==unique(MELTs$FORMAT)[i])
   tags<-strsplit(as.character(MELTs$FORMAT[idx[1]]),":",fixed=T)[[1]]
   if (anyDuplicated(tags) || !all(c("GT","AD","DP","SB")%in%tags)) stop("FORMAT requires unique GT, AD, DP and SB tags.")
   values<-strsplit(as.character(MELTs$value[idx]),":",fixed=T)
   if (any(lengths(values)!=length(tags))) stop("Genotype field length differs from FORMAT.")
   obj<-as.data.frame(do.call(rbind,values),stringsAsFactors=F)
   names(obj)<-tags
   names(obj)[names(obj)=="AF"]<-"AF2"
   for (j in setdiff(c("GT","AD","AF2","DP","SB","F1R2","F2R1","PGT","PID","PS"),names(obj))) obj[[j]]<-NA_character_
   tmp[[i]]<-cbind(obj[,c("GT","AD","AF2","DP","SB","F1R2","F2R1","PGT","PID","PS")],index=idx)
}
tmp5<-bind_rows(tmp)
tmp5<-tmp5[order(tmp5$index),]
# Records without allele 1 do not enter the source analysis. Filter before AD parsing.
idx<-which(grepl("(^|[/|])1([/|]|$)",tmp5$GT))
MELTs<-MELTs[idx,,drop=F]
tmp5<-tmp5[idx,,drop=F]
if (nrow(tmp5)) {
tmp5$AD<-gsub(",\\.,",",",tmp5$AD)
tmp5$AD<-gsub(",\\.$","",tmp5$AD)
tmp6<-strsplit(tmp5$AD,",",fixed=T)
tmp66<-strsplit(tmp5$SB,",",fixed=T)
if (any(lengths(tmp6)!=2L) || any(lengths(tmp66)!=4L)) stop("Expected biallelic AD and four SB counts; normalize multiallelic records first.")
tmp6<-as.data.frame(do.call(rbind,tmp6),stringsAsFactors=F)
tmp66<-as.data.frame(do.call(rbind,tmp66),stringsAsFactors=F)
colnames(tmp6)<-c("AD_ref","AD_alt")
colnames(tmp66)<-c("ref_fw","ref_rv","alt_fw","alt_rv")
tmp9<-cbind(MELTs,tmp5,tmp6,tmp66)
} else {
   tmp9<-cbind(MELTs,tmp5,data.frame(AD_ref=numeric(),AD_alt=numeric(),ref_fw=numeric(),ref_rv=numeric(),alt_fw=numeric(),alt_rv=numeric()))
}

for (i in c("DP","AD_ref","AD_alt","ref_fw","ref_rv","alt_fw","alt_rv")) {
   x<-as.character(tmp9[[i]])
   y<-suppressWarnings(as.numeric(x))
   if (any(is.na(y) & !is.na(x) & !x%in%c(".","")) || any(!is.finite(y[!is.na(y)])) || any(y<0,na.rm=T)) stop("Invalid genotype field: ",i)
   tmp9[[i]]<-y
}

# Numeric comparison also handles scientific notation; missing annotations stay NA.
x<-as.character(tmp9[[AF.column]])
tmp9$Population_AF<-suppressWarnings(as.numeric(x))
if (any(is.na(tmp9$Population_AF) & !is.na(x) & !x%in%c(".",""))) stop("Non-numeric population AF; check AF.column.")
if (any(!is.finite(tmp9$Population_AF[!is.na(tmp9$Population_AF)])) || any(tmp9$Population_AF<0 | tmp9$Population_AF>1,na.rm=T)) stop("Population AF must be a fraction in [0,1].")

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
   grepl("(^|[/|])1([/|]|$)",RAW$GT),
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
   (!is.na(Population_AF) & Population_AF<AF.cutoff) | (Keep.missing.AF & is.na(Population_AF)),
   # ExAC_EAS<0.1
   VAF>=0.05,
   VAF<0.40
)%>%mutate(ID2=paste0(varID,".",samID))

FIN<-FIN[!duplicated(FIN$ID2),]
FIN<-FIN%>%filter(samID%in%Sample$sample_column[Sample$include_wes])
FIN$case_number<-Sample$case_number[match(FIN$samID,Sample$sample_column)]
write.table(FIN,file.path(OUT,"selected_variants.tsv"),sep="\t",quote=F,row.names=F,na="NA")
Counts<-data.frame(case_number=Sample$case_number[Sample$include_wes])
Counts$variants<-as.integer(table(factor(FIN$case_number,levels=Counts$case_number)))
write.table(Counts,file.path(OUT,"variant_counts.tsv"),sep="\t",quote=F,row.names=F)
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

# Fusions are externally declared WTS findings, separate from called SNVs/indels.
# No dummy genomic variants are added to the MAF.
Fusion<-data.frame(case_number=character(),gene_symbol=character())
Clinical<-data.frame(Tumor_Sample_Barcode=Sample$case_number[Sample$include_plot],H.pylori=Sample$hp_status[Sample$include_plot])
if (!nrow(Clinical)) stop("At least one include_plot case is required.")
if (file.exists(Fusion.file)) {
   Fusion<-read.table(Fusion.file,header=T,sep="\t",check.names=F,stringsAsFactors=F,quote="",comment.char="")
   if (!all(c("case_number","gene_symbol")%in%names(Fusion)) || anyNA(Fusion[,c("case_number","gene_symbol")])) stop("Check declared_fusions.tsv.")
   if (any(!Fusion$case_number%in%Clinical$Tumor_Sample_Barcode) || any(!nzchar(Fusion$gene_symbol))) stop("Fusion cases must belong to the plotting cohort.")
}
write.table(Fusion,file.path(OUT,"declared_fusions_used.tsv"),sep="\t",quote=F,row.names=F)
write.table(tmp,file.path(OUT,"BeforeFiltering.tsv"),sep="\t",quote=F,row.names=F,na="NA")
if (nrow(tmp)) {
   MAF<-annovarToMaf(file.path(OUT,"BeforeFiltering.tsv"),tsbCol="Tumor_Sample_Barcode")
   write.table(MAF,file.path(OUT,"mutations.maf.tsv"),sep="\t",quote=F,row.names=F)
   before<-read.maf(maf=MAF,clinicalData=Clinical,
                    vc_nonSyn=c("Missense_Mutation","Nonsense_Mutation","Translation_Start_Site","Nonstop_Mutation","Frame_Shift_Ins","Frame_Shift_Del"))
   # read.maf keeps mutation-bearing samples; retain the declared negative cases for plotting.
   before@clinical.data<-data.table::as.data.table(Clinical)
   before@variant.classification.summary$Tumor_Sample_Barcode<-factor(before@variant.classification.summary$Tumor_Sample_Barcode,levels=Clinical$Tumor_Sample_Barcode)
   pdf(file.path(OUT,"oncoplot_snv_indel.pdf"),width=12,height=10)
   oncoplot(maf=before,cohortSize=nrow(Clinical),top=2000,colors=vc_cols,gene_mar=10,
            clinicalFeatures=if(nrow(Clinical)>1) 'H.pylori' else NULL,sortByAnnotation=nrow(Clinical)>1,showTumorSampleBarcodes=T,showPct=FALSE,
            sampleOrder=Clinical$Tumor_Sample_Barcode,removeNonMutated=FALSE)
   dev.off()
}
Events<-data.frame(case_number=FIN$case_number,gene=FIN$Gene.refGene,event=FIN$ExonicFunc.ensGene)
if (nrow(Fusion)) Events<-bind_rows(Events,data.frame(case_number=Fusion$case_number,gene=Fusion$gene_symbol,event="Gene_Fusion"))
Events<-unique(Events)
write.table(Events,file.path(OUT,"oncoplot_events.tsv"),sep="\t",quote=F,row.names=F)
if (nrow(Events)) {
   Events<-Events%>%filter(case_number%in%Clinical$Tumor_Sample_Barcode)%>%
      group_by(case_number,gene)%>%summarise(event=if(n()>1) "Multi_Hit" else first(event),.groups="drop")
   Events$case_number<-factor(Events$case_number,levels=Clinical$Tumor_Sample_Barcode)
   obj<-ggplot(Events,aes(x=case_number,y=gene,fill=event))+geom_tile(colour="white")+
      scale_x_discrete(drop=F)+theme_bw()+xlab("Case")+ylab("Gene")
   ggsave(file.path(OUT,"oncoplot_declared_events.pdf"),plot=obj,width=12,height=10)
}
capture.output(sessionInfo(),file=file.path(OUT,"sessionInfo.txt"))
