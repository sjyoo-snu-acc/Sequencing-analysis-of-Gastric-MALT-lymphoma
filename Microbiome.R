# Microbiome analysis from the genus/diversity/lm blocks in 250326.R.
# Run independently from WTS.R. Input formats are described in README.md.

#┌──────────────────────────────────┐
#│ [01] library & input settings    │----
#└──────────────────────────────────┘
library(dplyr)
library(ggplot2)
library(reshape2)

set.seed(48575)
Input.dir="data/Microbiome"
Run.ID=format(Sys.time(),"%Y%m%d-%H%M%S")
OUT=file.path("results/Microbiome",Run.ID)
Genus.file=file.path(Input.dir,"genus_counts.tsv")
Relative.file=file.path(Input.dir,"genus_relative_percent.tsv")
Alpha.file=file.path(Input.dir,"alpha_diversity.tsv")
Group.file=file.path(Input.dir,"grouptable_new.tsv")
Predictor.metric="relative_percent" # Source lm used reads_pct; set reads only for a deliberate reanalysis.
Alpha.mode="import" # Preserve supplied diversity indices; calculate is an optional new calculation.
Relative.mode="provided" # Optional sum_genus_counts uses the sum of supplied genus reads as denominator.
Genus.list=c("Helicobacter","Rothia","Veillonella","Porphyromonas")
# Record the actual tool/database/reference and denominator when preparing inputs.
Input.description="User-supplied genus tables and alpha indices; upstream versions/denominator must be documented"

if (dir.exists(OUT)) stop("Run directory exists. Change Run.ID.")
if (!Alpha.mode%in%c("import","calculate") || !Relative.mode%in%c("provided","sum_genus_counts")) stop("Invalid diversity or relative-abundance mode.")
if (!Predictor.metric%in%c("relative_percent","reads")) stop("Predictor.metric must be relative_percent or reads.")
Files=c(Genus.file,Group.file,if(Relative.mode=="provided") Relative.file,if(Alpha.mode=="import") Alpha.file)
if (!all(file.exists(Files))) stop("Missing Microbiome input. See README.md.")

#┌──────────────────────────────────┐
#│ [02] genus matrix & metadata     │----
#└──────────────────────────────────┘
cts=read.table(Genus.file,header=T,sep="\t",check.names=F,quote="",comment.char="",stringsAsFactors=F)
if (names(cts)[1]!="genus" || ncol(cts)<3 || anyNA(cts$genus) || any(!nzchar(cts$genus)) || anyDuplicated(cts$genus)) stop("genus_counts.tsv requires unique genus rows and sample columns.")
rownames(cts)=cts$genus
cts=cts[,-1,drop=F]
if (anyDuplicated(names(cts)) || !all(vapply(cts,is.numeric,logical(1))) || any(!is.finite(as.matrix(cts))) || any(as.matrix(cts)<0)) stop("Invalid genus read matrix.")
if (any(as.matrix(cts)!=floor(as.matrix(cts)))) stop("Use untransformed read counts in genus_counts.tsv.")
if (any(!Genus.list%in%rownames(cts))) stop("Supply all four named genus rows, including explicit zero counts where applicable.")
Samples=names(cts)
Group=read.table(Group.file,header=T,sep="\t",check.names=F,quote="",comment.char="",stringsAsFactors=F)
if (!all(c("Name","sample","HP")%in%names(Group)) || anyNA(Group[,c("Name","sample","HP")]) || anyDuplicated(Group$Name)) stop("Metadata requires unique Name, sample and HP.")
if (!all(Samples%in%Group$Name)) stop("A genus sample is missing from metadata.")
if (any(!Group$sample%in%c("Normal","Tumor","Eradicated")) || any(!Group$HP%in%c("positive","negative"))) stop("Unknown tissue/HP group.")
Group=Group%>%mutate(Type=case_when(
   sample=="Normal" & HP=="positive" ~ "PNormal",
   sample=="Normal" & HP=="negative" ~ "NNormal",
   sample=="Tumor" & HP=="positive" ~ "PTumor",
   sample=="Tumor" & HP=="negative" ~ "NTumor",
   sample=="Eradicated" & HP=="positive" ~ "PFollowup",
   sample=="Eradicated" & HP=="negative" ~ "NFollowup"),
   included_microbiome=Name%in%Samples)
# Include follow-up samples; the separate WTS exclusions are not applied here.
dir.create(OUT,recursive=T)
write.table(Group,file.path(OUT,"all_sample_metadata.tsv"),sep="\t",row.names=F,quote=F)
Group=Group[match(Samples,Group$Name),]
Group$Type=factor(Group$Type,levels=c("PNormal","NNormal","PTumor","NTumor","PFollowup","NFollowup"))
write.table(data.frame(file=Files,md5=unname(tools::md5sum(Files))),file.path(OUT,"input_hashes.tsv"),sep="\t",row.names=F,quote=F)
dput(mget(c("Predictor.metric","Alpha.mode","Relative.mode","Input.description","Genus.list")),file.path(OUT,"settings_used.R"))

#┌──────────────────────────────────┐
#│ [03] relative abundance          │----
#└──────────────────────────────────┘
if (Relative.mode=="provided") {
   reads_pct=read.table(Relative.file,header=T,sep="\t",check.names=F,quote="",comment.char="",stringsAsFactors=F)
   if (names(reads_pct)[1]!="genus" || anyDuplicated(reads_pct$genus) || anyNA(reads_pct$genus)) stop("Relative abundance requires unique genus rows.")
   rownames(reads_pct)=reads_pct$genus
   reads_pct=reads_pct[,-1,drop=F]
   if (anyDuplicated(names(reads_pct)) || !setequal(names(reads_pct),Samples) || !setequal(rownames(reads_pct),rownames(cts))) stop("Read/relative matrices differ in sample or genus IDs.")
   reads_pct=reads_pct[rownames(cts),Samples,drop=F]
   if (!all(vapply(reads_pct,is.numeric,logical(1))) || any(!is.finite(as.matrix(reads_pct))) || any(as.matrix(reads_pct)<0 | as.matrix(reads_pct)>100)) stop("Relative abundance must be a percentage in [0,100].")
   if (any(colSums(reads_pct)>100+1e-6)) stop("Summed genus percentages exceed 100; check rank, units and denominator.")
   if (any(as.matrix(cts)==0 & as.matrix(reads_pct)>0)) stop("Zero reads cannot have positive relative abundance.")
} else {
   if (any(colSums(cts)==0)) stop("Cannot normalize a sample with no genus reads.")
   reads_pct=as.data.frame(sweep(as.matrix(cts),2,colSums(cts),"/")*100,check.names=F)
}
write.table(cbind(genus=rownames(reads_pct),reads_pct),file.path(OUT,"genus_relative_percent_used.tsv"),sep="\t",quote=F,row.names=F)
MELT=reshape2::melt(as.matrix(reads_pct),varnames=c("genus","sample"),value.name="value")
MELT$Type=Group$Type[match(MELT$sample,Group$Name)]
MELT$sample=factor(MELT$sample,levels=Group$Name[order(Group$Type,Group$Name)])
obj=ggplot(MELT,aes(x=sample,y=value,fill=genus))+geom_bar(stat="identity")+
   facet_grid(~Type,scales="free_x",space="free_x")+theme_bw()+
   theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1))+
   xlab("Sample")+ylab("Relative abundance (%)")
ggsave(file.path(OUT,"Figure3A_genus_abundance.pdf"),plot=obj,width=14,height=7)

#┌──────────────────────────────────┐
#│ [04] alpha diversity             │----
#└──────────────────────────────────┘
if (Alpha.mode=="import") {
   alpha=read.table(Alpha.file,header=T,sep="\t",check.names=F,quote="",comment.char="",stringsAsFactors=F)
   if (!all(c("sample","Shannon","Simpson_index")%in%names(alpha)) || anyNA(alpha$sample) || anyDuplicated(alpha$sample)) stop("Alpha columns: sample, Shannon, Simpson_index.")
   if (!setequal(alpha$sample,Samples)) stop("Alpha-diversity and genus sample IDs differ.")
   alpha=alpha[match(Samples,alpha$sample),c("sample","Shannon","Simpson_index")]
} else {
   # Optional recalculation at genus rank: Shannon uses natural log; Simpson is 1-sum(p^2).
   # This is not substituted for the original alpha file without an explicit choice.
   if (any(colSums(cts)==0)) stop("Alpha diversity is undefined for an empty sample.")
   p=sweep(as.matrix(cts),2,colSums(cts),"/")
   plogp=p
   plogp[p>0]=p[p>0]*log(p[p>0])
   plogp[p==0]=0
   alpha=data.frame(sample=Samples,Shannon=-colSums(plogp),Simpson_index=1-colSums(p^2),row.names=NULL)
}
if (!all(vapply(alpha[,c("Shannon","Simpson_index")],is.numeric,logical(1))) || any(!is.finite(as.matrix(alpha[,c("Shannon","Simpson_index")])))) stop("Non-finite alpha diversity values.")
if (any(alpha$Shannon<0) || any(alpha$Simpson_index<0 | alpha$Simpson_index>1)) stop("Invalid diversity range; Simpson_index must be the 0-1 index, not its reciprocal.")
alpha$Type=Group$Type[match(alpha$sample,Group$Name)]
alpha$HP=Group$HP[match(alpha$sample,Group$Name)]
write.table(alpha,file.path(OUT,"alpha_diversity_used.tsv"),sep="\t",quote=F,row.names=F)
obj=reshape2::melt(alpha,id.vars=c("sample","Type","HP"))%>%
   ggplot(aes(x=HP,y=value,color=HP))+geom_boxplot()+geom_jitter(width=0.15)+
   facet_wrap(~variable,scales="free_y")+theme_bw()+xlab("H. pylori status")+ylab("Alpha diversity")
ggsave(file.path(OUT,"alpha_diversity_by_HP.pdf"),plot=obj,width=8,height=5)

#┌──────────────────────────────────┐
#│ [05] linear regression           │----
#└──────────────────────────────────┘
# The source used reads_pct, not the row order of the imported alpha table.
Predictors=if(Predictor.metric=="relative_percent") reads_pct else cts
df2=as.data.frame(t(Predictors[Genus.list,Samples,drop=F]))
df2$sample=rownames(df2)
df2$Shannon=alpha$Shannon[match(df2$sample,alpha$sample)]
df2$Simpson_index=alpha$Simpson_index[match(df2$sample,alpha$sample)]
if (anyNA(df2) || nrow(df2)<=length(Genus.list)+1) stop("Missing values or insufficient samples for the four-genus model.")
Models=list()
Estimates=list()
for (outcome in c("Shannon","Simpson_index")) {
   for (genus in Genus.list) {
      formula=reformulate(genus,response=outcome)
      X=model.matrix(formula,data=df2)
      if (qr(X)$rank<ncol(X)) stop("A regression predictor has no variation: ",genus)
      res=lm(formula,data=df2)
      nm=paste(outcome,genus,sep="_")
      Models[[nm]]=res
      coef_table=coef(summary(res))
      ci=confint(res)
      Estimates[[nm]]=data.frame(model="univariate",outcome=outcome,term=rownames(coef_table),
         estimate=coef_table[,1],std_error=coef_table[,2],statistic=coef_table[,3],pvalue=coef_table[,4],
         lower95=ci[,1],upper95=ci[,2],n=nobs(res),row.names=NULL)
   }
}
# Original multivariable response and four genus predictors are retained.
formula=Shannon ~ Helicobacter + Rothia + Veillonella + Porphyromonas
X=model.matrix(formula,data=df2)
if (qr(X)$rank<ncol(X)) stop("The four-genus regression is singular; review the supplied samples/predictors.")
res=lm(Shannon ~ Helicobacter + Rothia + Veillonella + Porphyromonas,data=df2)
Models$Shannon_multivariable=res
coef_table=coef(summary(res))
ci=confint(res)
Estimates$Shannon_multivariable=data.frame(model="multivariable",outcome="Shannon",term=rownames(coef_table),
   estimate=coef_table[,1],std_error=coef_table[,2],statistic=coef_table[,3],pvalue=coef_table[,4],
   lower95=ci[,1],upper95=ci[,2],n=nobs(res),row.names=NULL)
Estimates=bind_rows(Estimates)
Estimates$predictor_metric=Predictor.metric
write.table(df2,file.path(OUT,"regression_data.tsv"),sep="\t",quote=F,row.names=F)
write.table(Estimates,file.path(OUT,"regression_coefficients.tsv"),sep="\t",quote=F,row.names=F)
write.table(X,file.path(OUT,"multivariable_design.tsv"),sep="\t",quote=F,col.names=NA)
saveRDS(Models,file.path(OUT,"regression_models.rds"))
capture.output(lapply(Models,summary),file=file.path(OUT,"regression_summaries.txt"))
obj=reshape2::melt(df2,id.vars=c(Genus.list,"sample"),measure.vars=c("Shannon","Simpson_index"))%>%
   ggplot(aes(x=Helicobacter,y=value))+geom_point()+geom_smooth(method="lm",formula=y~x)+
   facet_wrap(~variable,scales="free_y")+theme_bw()+
   xlab(paste("Helicobacter (",Predictor.metric,")",sep=""))+ylab("Alpha diversity")
ggsave(file.path(OUT,"Figure3B_Helicobacter_diversity.pdf"),plot=obj,width=8,height=5)
obj=Estimates%>%filter(model=="multivariable",term!="(Intercept)")%>%
   ggplot(aes(x=term,y=estimate,ymin=lower95,ymax=upper95))+
   geom_pointrange()+geom_hline(yintercept=0,linetype="dashed")+coord_flip()+theme_bw()+
   xlab("Genus")+ylab("Shannon regression coefficient (95% CI)")
ggsave(file.path(OUT,"Figure3C_multivariable_Shannon.pdf"),plot=obj,width=8,height=5)
capture.output(sessionInfo(),file=file.path(OUT,"sessionInfo.txt"))
