data<-read.csv("I:/00_GROUP EXPERIMENTS/151125_allHela_chabbels_chambers/151126_only_4good_exp_HeLa_chambers/151126_group_150514_151116_151028_hela_fucci.csv")
a=3
title="151016_151028_hela_fucci_rosco20_stat_analyzed_151127"
dir="I:/00_GROUP EXPERIMENTS/151127_group_hela_rsco_151016_151028/make sub-dataframes/"
data$dvdt<-data$dv_tot/data$dt_tot

## make subdataframe for each phase (manually because I'm tired today) ####
# 240,5,1.6 ####
df<-data
t1 = c(240,5,240) ## enter the sequence of cycl_stage, in the order of appearance for one cell, from birth to division
t2= c(5,1.6,1.6)
names = c("birth","begining_g1_s","mitosis") ## enter all the timepoint you want to analyze for volume, time
names_phase = c("g1","s_g2_m","tot") ## enter all the phases you want to analyze for dv or dt, should be the name of the coumn in df1 -"dv" or "dt"
t_phase=c(240,5,240)
names_t_phase = c("birth","begining_g1_s","birth")


for (l in c(1:length(t1))) {
  K=which(colnames(df)==paste("v_",t1[l],sep=""))
  L=which(colnames(df)==paste("dv_",names_phase[l],sep=""))
  M=which(colnames(df)==paste("dt_",names_phase[l],sep=""))
  N=which(colnames(df)==paste("v_",t2[l],sep=""))
  d<-data.frame(df[,K],df[,L],df[,M],df[,N])
  colnames(d)=c(paste("v_",t1[l],sep=""),paste("dv_",names_phase[l],sep=""),paste("dt_",names_phase[l],sep=""),paste("v_",t2[l],sep=""))
  d$cond_date=df$cond_date
  d$cond_cell_type=df$cond_cell_type
  write.csv(d,paste(dir,title,"_sub_df_",names_phase[l],".csv",sep=""))
}

## chose df ####
df<-read.csv("I:/00_GROUP EXPERIMENTS/151127_group_hela_rsco_151016_151028/make sub-dataframes/151016_151028_hela_fucci_rosco20_stat_analyzed_151127_sub_df_tot.csv")
df2<-df
df2$cond_date=0 # to have statistiqcs on the entire dataset
df<-rbind(df,df2)
K= which(colnames(df)=="cond_date")
cond=unique(df$cond_date)
name_cond="_rosco20_subdf_g1_date"

library(robust)
library(ggplot2)
library(grid)
library(gridExtra)
library(MASS)
library(xtable)


####### bin by experiment #####
#alpha of vf vs vi plot ####
mat=matrix(NA,nrow=length(cond),ncol=11)
j=0
for (i in cond){
  j=j+1
  x=df$vi[which(df[K]==i & is.na(df$dv)==FALSE)]
  y=df$vf[which(df[K]==i & is.na(df$dv)==FALSE)]
  testa=lm(y~x)
  testb=lmrob(y~x)
  alpha=signif(coef(testb)[2],a)
  beta=signif(coef(testb)[1],a)
  mat[j,1]=i
  mat[j,2]=mean(y)
  mat[j,3]=sd(y)  
  mat[j,5]=signif(coef(testb)[2],a)
  mat[j,6]=signif(coef(testb)[1],a)
#   mat[j,7]=signif(anova(testa)$'Pr(>F)',a)
  mat[j,8]=signif(summary(testa)$r.squared,a)
  mat[j,9]=length(x) 
  mat[j,10]=signif(cor(x,y),a)
  mat[j,11]=signif(cor.test(x,y)$p.value,a)
}
colnames(mat)=c("cond","mean","sd","sem","alpha_vf_vi","beta_vf_vi","p_lm","R2_lm","n","pearson","p_pearson")
# table=tableGrob(mat,equal.width = F)
# grid.newpage()
# h <- grobHeight(table)
# w <- grobWidth(table)
# title_temp <- textGrob(paste(title,"\nsvf_vi_coeff_sorting_by ",name_cond), y=unit(0.5,"npc") + 0.6*h, 
#                        vjust=0, gp=gpar(fontsize=15))
# gt <- gTree(children=gList(table, title_temp))
# tiff(paste(dir,title,"_vf_vi_coeff_table_sorting_",name_cond,".tiff",sep=""), width=1000,height=200)
# grid.draw(gt)
# dev.off()
write.csv(as.data.frame(mat),paste(dir,title,"_vf_vi_coeff_table_sorting_",name_cond,".csv",sep=""))

#corr of dt vs vi plot ####
mat=matrix(NA,nrow=length(cond),ncol=21)
cc=NULL
j=0
for (i in cond){
  j=j+1
  x=df$dt[which(df[K]==i & is.na(df$dv)==FALSE)]
  y=df$dv[which(df[K]==i & is.na(df$dv)==FALSE)]
  z=df$vf[which(df[K]==i & is.na(df$dv)==FALSE)]
  zz=df$vi[which(df[K]==i & is.na(df$dv)==FALSE)]
  w=df$dvdt[which(df[K]==i & is.na(df$dv)==FALSE)]
  c=unique(as.character(df$cond_measure[which(df[K]==i)]))
  cc=c(cc,c)
  testa=lm(x~zz)
  testb=lmrob(x~zz)
  testaa=lm(y~zz)
  testbb=lmrob(y~zz)
  mat[j,1]=i
  mat[j,2]=NA
  mat[j,3]=mean(y)
  mat[j,4]=sd(y) 
  mat[j,5]=mean(x)
  mat[j,6]=sd(x) 
  mat[j,7]=mean(z)
  mat[j,8]=sd(z) 
  mat[j,9]=mean(zz)
  mat[j,10]=sd(zz)
  mat[j,11]=mean(w)
  mat[j,12]=sd(w)
  mat[j,13]=signif(cor(zz,x),a)
  mat[j,14]=signif(cor.test(zz,x)$p.value,a)
  mat[j,15]=cor.test(zz,x)$conf.int[1]
  mat[j,16]=cor.test(zz,x)$conf.int[2]
  mat[j,17]=signif(cor(zz,y),a)
  mat[j,18]=signif(cor.test(zz,y)$p.value,a)
  mat[j,19]=cor.test(zz,y)$conf.int[1]
  mat[j,20]=cor.test(zz,y)$conf.int[2]
  mat[j,21]=length(w)
}
colnames(mat)=c("cond_date","cond_measure","mean_dv","sd_dv",
                "mean_dt","sd_dt","mean_vf","sd_vf","mean_vi","sd_vi",
                "mean_dvdt","sd_dvdt",
                "Person_vi_dt","Pearson_pval_vi_dt","Pearson_ci_vi_dt_1","Pearson_ci_vi_dt_2",
                "Person_vi_dv","Pearson_pval_vi_dv","Pearson_ci_vi_dv_1","Pearson_ci_vi_dv_2",
                "n")

mat<-as.data.frame(mat)
mat$cond_measure<-cc
write.csv(as.data.frame(mat),paste(dir,title,"_dt_dv_vi_pearson_table_sorting_",name_cond,".csv",sep=""))

# table=tableGrob(mat,equal.width = F)
# grid.newpage()
# h <- grobHeight(table)
# w <- grobWidth(table)
# title_temp <- textGrob(paste(title,"\nsvf_vi_coeff_sorting_by ",name_cond), y=unit(0.5,"npc") + 0.6*h, 
#                        vjust=0, gp=gpar(fontsize=15))
# gt <- gTree(children=gList(table, title_temp))
# tiff(paste(dir,title,"_vf_vi_coeff_table_sorting_",name_cond,".tiff",sep=""), width=1000,height=200)
# grid.draw(gt)
# dev.off()

######## bin by mean dvdt

####### bin by dvdt #####
# enter parameters here + change line 5 in the loop , to add sorting-variable in df #####
df<-data[which(data$cond_measure=="chamber"),]
n = 10
sort = df$dvdt[which(is.na(df$dvdt)==F)]
sort2 = df$dvdt
bin = (max(sort)-min(sort))/n
df$sort_dvdt_all_cond = NA
bin_sort_dvdt_all_cond=NULL
name_tri="_bin_by_GR"
name_cond="_chamber"
cond_date=unique(df$cond_date)

for (i in c(1:nrow(df))){
  for (j in c(1:n)){
    if (is.na(sort2[i])==F
        && (sort2[i] <= min(sort)+bin*j)==TRUE && (sort2[i] > min(sort)-1+bin*(j-1))==TRUE) {
      df$sort_dvdt_all_cond[i] = j ## change this line !
    }else{}
  }
} 

# save a vector with bining values, make and save a table with informations about the sorting variable ####
bin_sort_dvdt_all_cond=NULL
for (j in c(0:n)){
  bin_sort_dvdt_all_cond=c(bin_sort_dvdt_all_cond,signif((min(sort)+bin*j),a))
}
res=matrix(ncol=3+length(cond),nrow=n)
res=as.data.frame(res)
for (i in c(1:n)){
  res[i,1]=signif(min(sort)-1+bin*(i-1),a)
  res[i,2]=signif(min(sort)+bin*i,a)
  res[i,3]=length(df$sort_dvdt_all_cond[which(df$sort_dvdt_all_cond==i)])
  for (j in c(4:(3+length(cond_date)))){
    if (length(df$sort_dvdt_all_cond[which(df$sort_dvdt_all_cond==i & df$cond_date==cond_date[j-3])])>0){
      res[i,j]=signif(length(df$sort_dvdt_all_cond[which(df$sort_dvdt_all_cond==i & df$cond_date==cond_date[j-3])])*100/res[i,3],a)
    }else {
      res[i,j]=0
    }
    colnames(res)[j]=paste("%", cond_date[j-3])
  }
}
colnames(res)[1:3]=c('min','max','n_tot')
table=tableGrob(res,equal.width = F)
grid.newpage()
h <- grobHeight(table)
w <- grobWidth(table)
title_temp <- textGrob(paste(title,name_cond,"\n",name_tri), y=unit(0.5,"npc") + 0.6*h, 
                       vjust=0, gp=gpar(fontsize=15))
gt <- gTree(children=gList(table, title_temp))
grid.draw(gt)
tiff(paste(dir,title,name_cond,"table",name_tri,".tiff",sep=""), width=500,height=300)
grid.draw(gt)
dev.off()
write.csv(as.data.frame(res),paste(dir,title,name_cond,"table_",name_tri,".csv",sep=""))

# corr of dt vs vi plot ####
K= which(colnames(df)=="sort_dvdt_all_cond")
cond=c(1:n)
name_cond="_subdf_g1"
name_tri="_date"
mat=matrix(NA,nrow=length(cond),ncol=27)
cc=NULL
j=0
for (i in cond){
  j=j+1
  x=df$dt[which(df[K]==i & is.na(df$dv)==FALSE)]
  y=df$dv[which(df[K]==i & is.na(df$dv)==FALSE)]
  z=df$vf[which(df[K]==i & is.na(df$dv)==FALSE)]
  zz=df$vi[which(df[K]==i & is.na(df$dv)==FALSE)]
  w=df$dvdt[which(df[K]==i & is.na(df$dv)==FALSE)]
  c=unique(as.character(df$cond_measure[which(df[K]==i)]))
  cc=c(cc,c)
  testa=lm(z~zz)
  testb=lmrob(z~zz)
  mat[j,1]=i
  mat[j,2]=NA
  mat[j,3]=mean(y)
  mat[j,4]=sd(y) 
  mat[j,5]=mean(x)
  mat[j,6]=sd(x) 
  mat[j,7]=mean(z)
  mat[j,8]=sd(z) 
  mat[j,9]=mean(zz)
  mat[j,10]=sd(zz)
  mat[j,11]=mean(w)
  mat[j,12]=sd(w)
  if (length(w)>3){
    mat[j,13]=signif(cor(zz,x),a)
    mat[j,14]=signif(cor.test(zz,x)$p.value,a)
    mat[j,15]=cor.test(zz,x)$conf.int[1]
    mat[j,16]=cor.test(zz,x)$conf.int[2]
    mat[j,17]=signif(cor(zz,y),a)
    mat[j,18]=signif(cor.test(zz,y)$p.value,a)
    mat[j,19]=cor.test(zz,y)$conf.int[1]
    mat[j,20]=cor.test(zz,y)$conf.int[2]
    mat[j,21]=signif(coef(testb)[2],a)
    mat[j,22]=signif(coef(testb)[1],a)
    mat[j,23]=signif(anova(testa)$'Pr(>F)'[1],a)
    mat[j,24]=signif(summary(testa)$r.squared,a)
    mat[j,25]=signif(cor(x,y),a)
    mat[j,26]=signif(cor.test(x,y)$p.value,a)    
  } else {}  
  mat[j,27]=length(w)
}
colnames(mat)=c("cond_date","cond_measure","mean_dv","sd_dv",
                "mean_dt","sd_dt","mean_vf","sd_vf","mean_vi","sd_vi",
                "mean_dvdt","sd_dvdt",
                "Person_vi_dt","Pearson_pval_vi_dt","Pearson_ci_vi_dt_1","Pearson_ci_vi_dt_2",
                "Person_vi_dv","Pearson_pval_vi_dv","Pearson_ci_vi_dv_1","Pearson_ci_vi_dv_2",
                "alpha_vf_vi","beta_vf_vi","p_lm_vf_vi","R2_lm_vf_vi","Person_vf_vi","Pearson_pval_vf_vi",
                "n")

mat<-as.data.frame(mat)
mat$cond_measure<-cc
write.csv(as.data.frame(mat),paste(dir,title,"_dt_dv_vi_alpha_pearson_table_sorting_",name_cond,name_tri,".csv",sep=""))

##### fucci analysis use subdataframe geenrated above ####
# choose dataframe of the phase to analyse ####
df<-read.csv("I:/00_GROUP EXPERIMENTS/151125_allHela_chabbels_chambers/151126_only_4good_exp_HeLa_chambers/group_all_3_hela_fucci_analyzed_151126_sub_df_tot.csv")
df<-df[,-1]
K= which(colnames(df)=="cond_date")
cond=unique(df$cond_date)
name_cond="_rosco20_subdf_tot"
name_tri="_date"

# corr of dt vs vi plot ####
mat=matrix(NA,nrow=length(cond),ncol=32)
cc=NULL
j=0
for (i in cond){
  j=j+1
#   x=df[,3][which(df[K]==i & is.na(df[2])==FALSE & is.na(df[3])==FALSE)]#dt
#   y=df[,2][which(df[K]==i & is.na(df[2])==FALSE & is.na(df[3])==FALSE)]#dv
#   z=df[,4][which(df[K]==i & is.na(df[2])==FALSE & is.na(df[3])==FALSE)]#vf
#   zz=df[,1][which(df[K]==i & is.na(df[2])==FALSE & is.na(df[3])==FALSE)]#vi
x=df$dt[which(df[K]==i & is.na(df$dv)==FALSE)]
y=df$dv[which(df[K]==i & is.na(df$dv)==FALSE)]
z=df$vf[which(df[K]==i & is.na(df$dv)==FALSE)]
zz=df$vi[which(df[K]==i & is.na(df$dv)==FALSE)]
  w=y/x
  testa=lm(z~zz)
  testb=lmrob(z~zz)
  mat[j,1]=i
  mat[j,2]=NA
  mat[j,3]=mean(y)
  mat[j,4]=sd(y)
  mat[j,5]=sd(y)/sqrt(length(y[which(is.na(y)==FALSE)]))
  mat[j,6]=mean(x)
  mat[j,7]=sd(x) 
  mat[j,8]=sd(x)/sqrt(length(x[which(is.na(x)==FALSE)]))
  mat[j,9]=mean(z)
  mat[j,10]=sd(z)
  mat[j,11]=sd(z)/sqrt(length(z[which(is.na(z)==FALSE)]))
  mat[j,12]=mean(zz)
  mat[j,13]=sd(zz)
  mat[j,14]=sd(zz)/sqrt(length(zz[which(is.na(zz)==FALSE)]))
  mat[j,15]=mean(w)
  mat[j,16]=sd(w)
  mat[j,17]=sd(w)/sqrt(length(w[which(is.na(w)==FALSE)]))
  if (length(w)>3){
    mat[j,18]=signif(cor(zz,x),a)
    mat[j,19]=signif(cor.test(zz,x)$p.value,a)
    mat[j,20]=cor.test(zz,x)$conf.int[1]
    mat[j,21]=cor.test(zz,x)$conf.int[2]
    mat[j,22]=signif(cor(zz,y),a)
    mat[j,23]=signif(cor.test(zz,y)$p.value,a)
    mat[j,24]=cor.test(zz,y)$conf.int[1]
    mat[j,25]=cor.test(zz,y)$conf.int[2]
    mat[j,26]=signif(coef(testb)[2],a)
    mat[j,27]=signif(coef(testb)[1],a)
    mat[j,28]=signif(anova(testa)$'Pr(>F)'[1],a)
    mat[j,29]=signif(summary(testa)$r.squared,a)
    mat[j,30]=signif(cor(x,y),a)
    mat[j,31]=signif(cor.test(x,y)$p.value,a)    
  } else {}  
  mat[j,32]=length(w [which(is.na(w)==FALSE)])
}
colnames(mat)=c("cond_date","cond_measure","mean_dv","sd_dv","sem_dv",
                "mean_dt","sd_dt","sem_dt","mean_vf","sd_vf","sem_vf","mean_vi","sd_vi","sem_vi",
                "mean_dvdt","sd_dvdt","sem_dvdt",
                "Person_vi_dt","Pearson_pval_vi_dt","Pearson_ci_vi_dt_1","Pearson_ci_vi_dt_2",
                "Person_vi_dv","Pearson_pval_vi_dv","Pearson_ci_vi_dv_1","Pearson_ci_vi_dv_2",
                "alpha_vf_vi","beta_vf_vi","p_lm_vf_vi","R2_lm_vf_vi","Person_vf_vi","Pearson_pval_vf_vi",
                "n")

mat<-as.data.frame(mat)
# mat$cond_measure<-cc
write.csv(as.data.frame(mat),paste(dir,title,"_dt_dv_vi_alpha_pearson_table_sorting",name_cond,name_tri,".csv",sep=""))
