#############   LINEAGE ANALYSIS VERSION 2 ##################################

#cette version permet:
#-d'avoir les données rangées dans un ordre aléatoire car elle n'utilise QUE les numéros d'identité des cellules
#- d'analyser et de comparer plusieurs sets de données

##import des données (modèle de dataframe obtenu à partir du programme growth_in_lineage_v2) ####
data1a <- read.csv("H:/150423_Raji_bilan/150423_combine Raji_131122_131220//131122_Raji_analyzed_150422_data_frame.csv")
data2a <- read.csv("H:/150423_Raji_bilan/150423_combine Raji_131122_131220//131220_Raji_analysis_150403 _data_frame.csv")
data2a = data2a[-which(data2a$lineage==2),]
data1=rbind(data1a,data2a)
write.csv(data1,paste(dir,title,"df.csv",sep=""))

data1b <-read.csv("H:/150423_Raji_bilan/150423_combine Raji_131122_131220//131122_Raji_analysis_150425_df_abs_values.csv")
data2b <-read.csv("H:/150423_Raji_bilan/150423_combine Raji_131122_131220//131220_Raji_analysis_wo_l2_150425_df_abs_values.csv")
data1b$cond_cell_type = "Raji"
data2b$cond_cell_type = "Raji"
data1b$cond_chamber = "vol_15.5"
data2b$cond_chamber = "vol_15.5"
data1b$cond_date = "131122"
data2b$cond_date = "131220"
data2=rbind(data1b,data2b)
write.csv(data2,paste(dir,title,"df_abs-val.csv",sep=""))

## remove data points that were found not normal from analysis below and re-save data frma ewith name = dataset2 for ex ####
title= paste(title,"_dataset2",sep="")
lineage=c(3,10,10,7)
number = c(6,20,21,5)
cond="131220"
for (i in c(1:length(lineage))){
  data1=data1[-which(data1$lineage==lineage[i] & data1$number_in_lineage==number[i] & data1$cond_date==cond),]
}

for (i in c(1:length(lineage))){
  data2b=data2b[-which(data2b$lineage==lineage[i] & data2b$number_in_lineage==number[i]),]
}
data2=rbind(data2b,data1b)
write.csv(data2,paste(dir,title,"df_abs-val.csv",sep=""))
write.csv(data1,paste(dir,title,"df.csv",sep=""))

## parameters to enter ####
dir=("H:/150423_Raji_bilan/150423_combine Raji_131122_131220/")
title="Raji_131122_131220_analyzed_150425"
a = 2 #number of significative digits after coma
cond1 = "150213"
cond2 = "150417"
cond3= "150417"
cond=c("150213","150417")

## call packages ####
library(ggplot2)
library(grid)
library(gridExtra)
library(MASS)
library(xtable)

###################### ANALYSIS ############################################

######### controls for ccl and volume ########## 
#### 1) CCL #####
##### choose parameters for analysis ####
df<-data1

m = 1.6 # timepoint chosen for end of cell cycle
b = 240 # timepoint chosen for bgining of cell cycle
e = 3 # timepoitn where cell is lost alive (e for exit)
# end = max(df$Time)/time_resolution # final frame of the movie
a=2
tri=df$lineage_complete
name_tri='lineages'

##### tri ####
ccl_complete = NULL
ccl_uncomplete = NULL
cell_number_complete = NULL
cell_number_uncomplete = NULL
generation_uncomplete=NULL
generation_complete = NULL
lineage_complete = NULL
lineage_uncomplete = NULL
ti_complete = NULL
tf_complete = NULL
cond_complete = NULL
cond_uncomplete= NULL

for (k in cond){
  end = max(df$Time[which(df$cond_date==k)])
  for (i in unique(df$lineage[which(df$cond_date==k)])){
    for (j in unique(c(df$number_in_lineage[which(df$lineage==i & df$number_in_lineage>1 & df$cond_date==k)]))){
      if (m %in% df$cycle_stage[which(df$lineage==i & df$number_in_lineage==j & df$cond_date==k)]
          && b %in% df$cycle_stage[which(df$lineage==i & df$number_in_lineage==j & df$cond_date==k)]){
        t1 = (df$Time[which(df$lineage==i & df$number_in_lineage==j & df$cycle_stage==m & df$cond_date==k)]-
          df$Time[which(df$lineage==i & df$number_in_lineage==j & df$cycle_stage==b & df$cond_date==k)])/60
        ccl_complete = c(ccl_complete,t1)
        cell_number_complete = c(cell_number_complete,j)
        generation_complete = c(generation_complete, df$generation[which(df$lineage==i & df$number_in_lineage==j & df$cycle_stage==b & df$cond_date==k)])
        lineage_complete = c(lineage_complete, i)
        ti_complete=c(ti_complete,df$Time[which(df$lineage==i & df$number_in_lineage==j & df$cycle_stage==b & df$cond_date==k)])
        tf_complete = c(tf_complete,df$Time[which(df$lineage==i & df$number_in_lineage==j & df$cycle_stage==m & df$cond_date==k)])
        cond_complete=c(cond_complete,k)
        rm (t1) 
      }else if (b %in% df$cycle_stage[which(df$lineage==i & df$number_in_lineage==j & df$cond_date==k)]){
        if (e %in% df$cycle_stage[which(df$lineage==i & df$number_in_lineage==j & df$cond_date==k)]) {
          t2 = (df$Time[which(df$lineage==i & df$number_in_lineage==j & df$cycle_stage==e & df$cond_date==k)] -
                  df$Time[which(df$lineage==i & df$number_in_lineage==j & df$cycle_stage==b & df$cond_date==k)])/60
          ccl_uncomplete = c(ccl_uncomplete,t2)
          cell_number_uncomplete = c(cell_number_uncomplete,j)
          generation_uncomplete = c(generation_uncomplete, df$generation[which(df$lineage==i & df$number_in_lineage==j & df$cycle_stage==b & df$cond_date==k)])
          lineage_uncomplete = c(lineage_uncomplete,i)
          cond_uncomplete=c(cond_uncomplete,k)
          rm (t2)
        }else{
          t2 = (end -
            df$Time[which(df$lineage==i & df$number_in_lineage==j & df$cycle_stage==b & df$cond_date==k)])/60
          ccl_uncomplete = c(ccl_uncomplete,t2)
          cell_number_uncomplete = c(cell_number_uncomplete,j)
          generation_uncomplete = c(generation_uncomplete, df$generation[which(df$lineage==i & df$number_in_lineage==j & df$cycle_stage==b & df$cond_date==k)])
          lineage_uncomplete = c(lineage_uncomplete,i)
          cond_uncomplete=c(cond_uncomplete,k)
          rm (t2)           
        }
      }
    }
  }
}

##### plots ####
type <- c(rep('complete',length(ccl_complete)),rep('uncomplete',length(ccl_uncomplete)))
cell_number = c(cell_number_complete,cell_number_uncomplete)
lineage <- c(lineage_complete, lineage_uncomplete)

###### hist  CCL complètes/non complètes ####
df <- data.frame(ccl_complete,cond_complete,cell_number_complete,lineage_complete)

xlab=" "
for (i in cond){
  xlab=paste(xlab,"\n",i,": mean=",signif(mean(df$ccl_complete[which(df$cond_complete==i)]),a), 
             ", n=",length(df$ccl_complete[which(df$cond_complete==i)]),sep="")
}
hist1 <- ggplot(df,aes(x=ccl_complete,fill=as.factor(cond_complete)))+
  geom_histogram(position="dodge")+
  theme_bw()+
  labs(title=paste(title, "complete ccl"),
       fill="cond",
       x= paste('DT [hours]',xlab,"\n","all conditions: mean=",signif(mean(ccl_complete),a), ", n=",length(ccl_complete),sep=""))
rm(df)

df <- data.frame(ccl_uncomplete,cond_uncomplete,cell_number_uncomplete,lineage_uncomplete)
xlab=" "
for (i in cond){
  xlab=paste(xlab,"\n",i,": n=",length(df$ccl_uncomplete[which(df$cond_uncomplete==i)]),sep="")
}
hist2 <- ggplot(df,aes(ccl_uncomplete, fill=as.factor(cond_uncomplete)))+
  geom_histogram(position="dodge")+
  theme_bw()+
  labs(title=paste(title," uncomplete ccl"),
       fill="cond",
       x= paste('DT [hours]',xlab,"\n","all conditions: n=", length(ccl_uncomplete),sep=""))
rm(df)

###### plot  CCL time ####
df <- data.frame(ccl_complete,generation_complete,ti_complete/60,tf_complete/60,lineage_complete,cond_complete)
xlab=" "
for (i in cond){
  x=df$ti_complete.60[which(df$cond_complete==i)]
  y=df$ccl_complete[which(df$cond_complete==i)]
  test=lm(y~x)
  xlab=paste(xlab,"\n",i,", lm: R2=",signif(summary(test)$r.squared,a),", p=",signif(anova(test)$'Pr(>F)',a), ", n= ",length(x),sep="")
  rm(x,y,test)
}

x=df$ti_complete.60
y=df$ccl_complete
test=lm(y~x)
testb=rlm(y~x)
pa <- ggplot(df,aes(x=ti_complete.60,y=ccl_complete,color=factor(cond_complete))) +
  geom_point() +
  labs(title= title, y="DT [hrs]", 
       x=paste("Time_i [hrs]","\n lm:","R2=",signif(summary(test)$r.squared,a),", p=",signif(anova(test)$'Pr(>F)',a),
               "a=",signif(coef(test)[2],a), ", b=", signif(coef(test)[1],2),
               "\n rlm: a=",signif(coef(testb)[2],a), ", b=", signif(coef(testb)[1],a),xlab,sep=""))+
  scale_colour_discrete (name="cond") +
#      geom_abline(intercept = coef(testb)[1], slope = coef(testb)[2])+
  #   stat_smooth(method='lm',se=FALSE)+
  theme_bw()
sink(file=paste(dir,title,"ctrl_dt_ti_lm_rlm.txt",sep=""))
summary(test)
summary(testb)
sink()
rm(test,testb,y,x)

xlab=" "
for (i in cond){
  x=df$tf_complete.60[which(df$cond_complete==i)]
  y=df$ccl_complete[which(df$cond_complete==i)]
  test=lm(y~x)
  xlab=paste(xlab,"\n",i,", lm: R2=",signif(summary(test)$r.squared,a),", p=",signif(anova(test)$'Pr(>F)',a), ", n= ",length(x),sep="")
  rm(x,y,test)
}
x=df$tf_complete.60
y=df$ccl_complete
test=lm(y~x)
testb=rlm(y~x)
pb <- ggplot(df,aes(x=tf_complete.60,y=ccl_complete,color=factor(cond_complete))) +
  geom_point() +
  labs(title= title, y="DT [hrs]", 
       x=paste("Time_f [hrs]","\n lm:","R2=",signif(summary(test)$r.squared,a),", p=",signif(anova(test)$'Pr(>F)',a),
               "a=",signif(coef(test)[2],a), ", b=", signif(coef(test)[1],2),
               "\n rlm: a=",signif(coef(testb)[2],a), ", b=", signif(coef(testb)[1],a),xlab,sep=""))+
  scale_colour_discrete (name="cond") +
  #   geom_abline(intercept = coef(testb)[1], slope = coef(testb)[2])+
  #   stat_smooth(method='lm',se=FALSE)+
  theme_bw()
sink(file=paste(dir,title,"ctrl_dt_tf_lm_rlm.txt",sep=""))
summary(test)
summary(testb)
sink()
rm(test,testb,y,x)

###### make figure ####
tiff(paste(dir,title,"_hist_plot_ccl_ctrl.tiff",sep=""), width=1000,height=500)
grid.arrange(hist1,pa,hist2,pb,ncol=2)
dev.off()
rm(df,tf_complete,ti_complete,ccl_complete,ccl_uncomplete,generation_complete,generation_uncomplete,lineage_complete,lineage_uncomplete,cell_number,ccl,cell_number_complete,cell_number_uncomplete)
rm(hist1,hist2,pa,pb)

#### 2) VF, Vi, DV #####
##### choose parameters for analysis ######
df <- data2
a=2
cond=c("130306", "140801")

#### hist ####
xlab=" "
for (i in cond){
  xlab=paste(xlab,"\n",i,": mean=",signif(mean(df$vf[which(df$cond_date==i & is.na(df$vf)==F)]),a), ", n=",length(df$vf[which(df$cond_date==i & is.na(df$vf)==F)]),sep="")
}
hist3 <- ggplot(df,aes(x=vf,fill=as.factor(cond_date)))+
  geom_histogram(position="dodge")+
  theme_bw()+
  labs(title=title,
       fill="cond",
       x= paste('Vf [µm3]','\n all conditions: ', 'mean=',signif(mean(df$vf,na.rm=TRUE),a*2)," µm3",
                ", n=",length(df$vf[which(is.na(df$vf)==F)]),xlab,sep=""))

xlab=" "
for (i in cond){
  xlab=paste(xlab,"\n",i,": mean=",signif(mean(df$vi[which(df$cond_date==i & is.na(df$vi)==F)]),a), ", n=",length(df$vi[which(df$cond_date==i & is.na(df$vi)==F)]),sep="")
}
hist4 <- ggplot(df,aes(x=vi,fill=as.factor(cond_date)))+
  geom_histogram(position="dodge")+
  theme_bw()+
  labs(title=title,
       fill="generation",
       x= paste('Vi [µm3]','\nall conditions: ', 'mean=',signif(mean(df$vi,na.rm=TRUE),a*2),
                ", n=",length(df$vi[which(is.na(df$vi)==F)]),xlab,sep=""))

xlab=" "
for (i in cond){
  xlab=paste(xlab,"\n",i,": mean=",signif(mean(df$dv[which(df$cond_date==i & is.na(df$dv)==F)]),a), ", n=",length(df$dv[which(df$cond_date==i & is.na(df$dv)==F)]),sep="")
}
hist5 <- ggplot(df,aes(x=dv,fill=as.factor(cond_date)))+
  geom_histogram(position="dodge")+
  theme_bw()+
  labs(title=title,
       fill="cond",
       x= paste('DV [µm3]','\nall conditions: ', 'mean=',signif(mean(df$dv,na.rm=TRUE),a*2),
                ', n=',length(df$dv[which(is.na(df$dv)==F)]),xlab,sep=""))

xlab=" "
for (i in cond){
  xlab=paste(xlab,"\n",i,": mean=",signif(mean(df$dvdt[which(df$cond_date==i & is.na(df$dvdt)==F)]),a), ", n=",length(df$dvdt[which(df$cond_date==i & is.na(df$dvdt)==F)]),sep="")
}
hist6 <- ggplot(df,aes(x=dvdt,fill=as.factor(cond_date)))+
  geom_histogram(position="dodge")+
  theme_bw()+
  labs(title=title,
       fill="cond",
       x= paste('DV/DT [µm3/hrs]','\nall conditions: ', 'mean=',signif(mean(df$dvdt,na.rm=TRUE),a*2),
                ', n=',length(df$dvdt[which(is.na(df$dvdt)==F)]),xlab,sep=""))

xlab=" "
for (i in cond){
  xlab=paste(xlab,"\n",i,": mean=",signif(mean(df$dt[which(df$cond_date==i & is.na(df$dt)==F)]),a), ", n=",length(df$dt[which(df$cond_date==i & is.na(df$dt)==F)]),sep="")
}
hist7 <- ggplot(df,aes(x=dt,fill=as.factor(cond_date)))+
  geom_histogram(position="dodge")+
  theme_bw()+
  labs(title=title,
       fill="cond",
       x= paste('DT [hrs]','\nall conditions: ', 'mean=',signif(mean(df$dt,na.rm=TRUE),a*2),
                ', n=',length(df$dt[which(is.na(df$dt)==F)]),xlab,sep=""))

#### plot ####
xlab=" "
for (i in cond){
  x=df$ti[which(df$cond_date==i)]
  y=df$vf[which(df$cond_date==i)]
  test=lm(y~x)
  xlab=paste(xlab,"\n",i,", lm: R2=",signif(summary(test)$r.squared,a),", p=",signif(anova(test)$'Pr(>F)',a), ", n= ",length(x),sep="")
  rm(x,y,test)
}
test=lm(df$vf~df$ti)
p1 <- ggplot(df,aes(x=ti,y=vf, color=factor(cond_date))) +
  geom_point() +
  labs(title= title, y="Vf [µm3]", x=paste("Time_i [hours]","\n lm: ","R2=",signif(summary(test)$r.squared,2),
                                           ", p=",signif(anova(test)$'Pr(>F)',2),
                                           ", n=",length(df$vf[which(is.na(df$vf)==F)]),xlab)) +
  scale_colour_discrete (name="date") +
  #stat_smooth(method="lm")+
  theme_bw()
sink(file=paste(dir,title,"ctrl_lm_vf_ti.txt",sep=""))
summary(test)
sink() 
rm(test)

xlab=" "
for (i in cond){
  x=df$ti[which(df$cond_date==i)]
  y=df$dv[which(df$cond_date==i)]
  test=lm(y~x)
  xlab=paste(xlab,"\n",i,", lm: R2=",signif(summary(test)$r.squared,a),", p=",signif(anova(test)$'Pr(>F)',a), ", n= ",length(x),sep="")
  rm(x,y,test)
}
test = lm(df$dv~df$ti)
p2 <- ggplot(df,aes(x=ti,y=dv,color=factor(cond_date))) +
  geom_point() +
  labs(title= title, y="DV [µm3]", x=paste("Time_i [hours]",
                                           "\n lm:","R2=",signif(summary(test)$r.squared,2),", p=",signif(anova(test)$'Pr(>F)',2),", n=",length(df$dv[which(is.na(df$dv)==F)]),
                                           xlab,sep="")) +
  scale_colour_discrete (name="date") +
  #   stat_smooth(method="lm",se=FALSE)+
  theme_bw()
sink(file=paste(dir,title,"ctrl_lm_dv_ti.txt",sep=""))
summary(test)
sink() 
rm(test)

xlab=" "
for (i in cond){
  x=df$ti[which(df$cond_date==i)]
  y=df$dvdt[which(df$cond_date==i)]
  test=lm(y~x)
  xlab=paste(xlab,"\n",i,", lm: R2=",signif(summary(test)$r.squared,a),", p=",signif(anova(test)$'Pr(>F)',a), ", n= ",length(x),sep="")
  rm(x,y,test)
}
y=df$dvdt
x=df$ti
test=lm(df$dvdt~df$ti)
p3 <- ggplot(df,aes(x=ti,y=dvdt,color=factor(cond_date))) +
  geom_point() +
  labs(title= title, y="DV/DT [µm3/hrs]", 
       x=paste("Time_i [hours]",
               "\n","R2=",signif(summary(test)$r.squared,2),", p=",signif(anova(test)$'Pr(>F)',2),", n=",length(df$dvdt[which(is.na(df$dv)==F)]),
               xlab,sep="")) +
  scale_colour_discrete (name="date") +
  theme_bw()
sink(file=paste(dir,title,"ctrl_lm_dvdt_ti.txt",sep=""))
summary(test)
sink() 
rm(test)

#### make figure ####
tiff(paste(dir,title,"_hist_vol_ctrl.tiff",sep=""), width=500,height=750)
grid.arrange(hist3,hist4,hist5)
dev.off()

tiff(paste(dir,title,"_plot_vol_ctrl.tiff",sep=""), width=500,height=750)
grid.arrange(p1,p2,p3, ncol=1)
dev.off()

tiff(paste(dir,title,"_hist_plot_vol_ctrl.tiff",sep=""), width=1000,height=1000)
grid.arrange(hist3,p1,hist4,p2,hist5,p3,hist6, ncol=2)
dev.off()

rm(df, p1,p2,p3,hist3,hist4,hist5,hist6,hist7)

########## absolute values dv, vf, vi, dt, dv/dt ########## 
#### choose df ####
df=data2
#### add parameter for sorting ####
### for time_i ####
# enter parameters here + change line 5 in the loop , to add sorting-variable in df #####
n = 4
sort = df$ti[which(is.na(df$ti)==F)]
sort2 = df$ti
bin = (max(sort)-min(sort))/n
df$sort_ti_all_cond = NA
bin_sort_ti_all_cond=NULL
name_tri="time_i"
cond=c("150213","150417")

for (i in c(1:nrow(df))){
  for (j in c(1:n)){
    if (is.na(sort2[i])==F
        && (sort2[i] <= min(sort)+bin*j)==TRUE && (sort2[i] > min(sort)-1+bin*(j-1))==TRUE) {
      df$sort_ti_all_cond[i] = j ## change this line !
    }else{}
  }
} 

# make vector and save a table with bining informations about the sorting variable ####
bin_sort_ti_all_cond=NULL
for (j in c(0:n)){
  bin_sort_ti_all_cond=c(bin_sort_ti_all_cond,signif((min(sort)+bin*j)/60,a))
}
res=matrix(ncol=3+length(cond),nrow=n)
res=as.data.frame(res)
for (i in c(1:n)){
  res[i,1]=min(sort)-1+bin*(i-1)
  res[i,2]=min(sort)+bin*i
  res[i,3]=length(df$sort_ti_all_cond[which(df$sort_ti_all_cond==i)])
  for (j in c(4:(3+length(cond)))){
    res[i,j]=signif(length(df$sort_ti_all_cond[which(df$sort_ti_all_cond==i & df$cond_date==cond[j-3])])*100/res[i,3],a)
    colnames(res)[j]=paste("%", cond[j-3])
  }
}
colnames(res)[1:3]=c('min','max','n_tot')
table=tableGrob(res,equal.width = F)
grid.newpage()
h <- grobHeight(table)
w <- grobWidth(table)
title_temp <- textGrob(paste(title,"\nsorting by ",name_tri), y=unit(0.5,"npc") + 0.6*h, 
                       vjust=0, gp=gpar(fontsize=15))
gt <- gTree(children=gList(table, title_temp))
grid.draw(gt)
tiff(paste(dir,title,"_abs_val_table_sorting_",name_tri,".tiff",sep=""), width=500,height=300)
grid.draw(gt)
dev.off()

# write.csv(df, paste(dir,title,"_df_abs_values.csv",sep=""))

### for average_gr ####
# enter parameters here + change line 5 in the loop , to add sorting-variable in df #####
n = 4
sort = df$dvdt[which(is.na(df$dvdt)==F)]
sort2 = df$dvdt
bin = (max(sort)-min(sort))/n
df$sort_dvdt_all_cond = NA
bin_sort_dvdt_all_cond=NULL
name_tri="average_GR"

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
  for (j in c(4:(3+length(cond)))){
    res[i,j]=signif(length(df$sort_dvdt_all_cond[which(df$sort_dvdt_all_cond==i & df$cond_date==cond[j-3])])*100/res[i,3],a)
    colnames(res)[j]=paste("%", cond[j-3])
  }
}
colnames(res)[1:3]=c('min','max','n_tot')
table=tableGrob(res,equal.width = F)
grid.newpage()
h <- grobHeight(table)
w <- grobWidth(table)
title_temp <- textGrob(paste(title,"\nsorting by ",name_tri), y=unit(0.5,"npc") + 0.6*h, 
                       vjust=0, gp=gpar(fontsize=15))
gt <- gTree(children=gList(table, title_temp))
grid.draw(gt)
tiff(paste(dir,title,"_abs_val_table_sorting_",name_tri,".tiff",sep=""), width=500,height=300)
grid.draw(gt)
dev.off()


# write.csv(df, paste(dir,title,"_df_abs_values.csv",sep=""))

### for Vi ####
# enter parameters here + change line 5 in the loop , to add sorting-variable in df #####
n = 4
sort = df$vi[which(is.na(df$vi)==F)]
sort2 = df$vi
bin = (max(sort)-min(sort))/n
df$sort_vi_all_cond = NA
bin_sort_vi_all_cond=NULL
name_tri="vi"

for (i in c(1:nrow(df))){
  for (j in c(1:n)){
    if (is.na(sort2[i])==F
        && (sort2[i] <= min(sort)+bin*j)==TRUE && (sort2[i] > min(sort)-1+bin*(j-1))==TRUE) {
      df$sort_vi_all_cond[i] = j ## change this line !
    }else{}
  }
} 

# save a vector with bining values, make and save a table with informations about the sorting variable ####
bin_sort_vi_all_cond=NULL
for (j in c(0:n)){
  bin_sort_vi_all_cond=c(bin_sort_vi_all_cond,signif((min(sort)+bin*j),a))
}
res=matrix(ncol=3+length(cond),nrow=n)
res=as.data.frame(res)
for (i in c(1:n)){
  res[i,1]=signif(min(sort)-1+bin*(i-1),a)
  res[i,2]=signif(min(sort)+bin*i,a)
  res[i,3]=length(df$sort_vi_all_cond[which(df$sort_vi_all_cond==i)])
  for (j in c(4:(3+length(cond)))){
    res[i,j]=signif(length(df$sort_vi_all_cond[which(df$sort_vi_all_cond==i & df$cond_date==cond[j-3])])*100/res[i,3],a)
    colnames(res)[j]=paste("%", cond[j-3])
  }
}
colnames(res)[1:3]=c('min','max','n_tot')
table=tableGrob(res,equal.width = F)
grid.newpage()
h <- grobHeight(table)
w <- grobWidth(table)
title_temp <- textGrob(paste(title,"\nsorting by ",name_tri), y=unit(0.5,"npc") + 0.6*h, 
                       vjust=0, gp=gpar(fontsize=15))
gt <- gTree(children=gList(table, title_temp))
grid.draw(gt)
tiff(paste(dir,title,"_abs_val_table_sorting_",name_tri,".tiff",sep=""), width=500,height=300)
grid.draw(gt)
dev.off()


# write.csv(df, paste(dir,title,"_df_abs_values.csv",sep=""))


#### choose parameters for analysis ####
tri = df$sort_ti_all_cond # column in dataframe with the sorting variable
K= which(colnames(df)=="sort_ti_all_cond") # number of the column with the sorting variable
name_tri = 'time_i' # title of the sorting variable
cond2 = bin_sort_ti_all_cond # list of the different values taken by the sorting variable
cond=c(1,2,3,4)
a=2 # number of significant digits

#### plots ####
## p1 dv/dt vs vi ####
xlab=" "
for (i in cond){
  x=df$vi[which(df[K]==i & is.na(df$vi)==F & is.na(df$dvdt)==F)]
  y=df$dvdt[which(df[K]==i & is.na(df$vi)==F & is.na(df$dvdt)==F)]
  test=lm(y~x)
  testb=rlm(y~x)
  xlab=paste(xlab,"\n",cond[i],", lm: R2=",signif(summary(test)$r.squared,a),", p=",signif(anova(test)$'Pr(>F)',a),",rlm: a=",signif(coef(testb)[2],a), ", b=", signif(coef(testb)[1],a), ", n= ",length(x),sep="")
  rm(x,y,test,testb)
}
x=df$vi[which(is.na(df$vi)==FALSE & is.na(df$dvdt)==FALSE)]
y=df$dvdt[which(is.na(df$vi)==FALSE & is.na(df$dvdt)==FALSE)]
test = lm(y~x)
testb = rlm(y~x)
p1 <- ggplot(df,aes(x=vi,y=dvdt, color=factor(tri))) +
  geom_point() +
  labs(title= title, y="DV/DT [µm3/hrs]", 
       x=paste("Vol_i [µm3]","\nall cond: lm:","R2=",signif(summary(test)$r.squared,a),", p=",signif(anova(test)$'Pr(>F)',a),
               ", n=",length(df$dvdt[which(is.na(df$vi)==FALSE & is.na(df$dvdt)==FALSE)]), 
               ", rlm: a=",signif(coef(testb)[2],a), ", b=", signif(coef(testb)[1],a),
               "\n Pearson: r=", signif(cor(x,y),a), ", p=",signif(cor.test(x,y)$p.value,a),xlab))+
  scale_colour_discrete (name=name_tri) +
  geom_abline(intercept = coef(testb)[1], slope = coef(testb)[2])+
  stat_smooth(method='lm',se=FALSE)+
  theme_bw()
sink(file=paste(dir,title,"_lm_rlm_coeff-corr_dvdt_vi.txt",sep=""))
summary(test)
summary(testb)
cor.test(x,y)
sink()
rm(test,testb,y,x)

## p2 vf/vi vs vi ####
xlab=" "
for (i in cond){
  x=df$vi[which(df[K]==i & is.na(df$vi)==F & is.na(df$vfvi)==F)]
  y=df$vfvi[which(df[K]==i & is.na(df$vi)==F & is.na(df$vfvi)==F)]
  test=lm(y~x)
  testb=rlm(y~x)
  xlab=paste(xlab,"\n",i,", lm: R2=",signif(summary(test)$r.squared,a),", p=",signif(anova(test)$'Pr(>F)',a),",rlm: a=",signif(coef(testb)[2],a), ", b=", signif(coef(testb)[1],a), ", n= ",length(x),sep="")
  print(i)
  rm(x,y,test,testb)
}
x=df$vi[which(is.na(df$vi)==FALSE & is.na(df$vfvi)==FALSE)]
y=df$vfvi[which(is.na(df$vi)==FALSE & is.na(df$vfvi)==FALSE)]
test = lm(y~x)
testb = rlm(y~x)
p2 <- ggplot(df,aes(x=vi,y=vfvi, color=factor(tri))) +
  geom_point() +
  labs(title= title, y="vf/vi %", 
       x=paste("Vol_i [µm3]","\n lm: ","R2=",signif(summary(test)$r.squared,2),", p=",signif(anova(test)$'Pr(>F)',2),
               ", n=",length(df$vfvi[which(is.na(df$vi)==FALSE & is.na(df$vfvi)==FALSE)]),
               ", rlm: a=",signif(coef(testb)[2],a), ", b=", signif(coef(testb)[1],a),
               "\n Pearson: r=", signif(cor(x,y),a), "p=",signif(cor.test(x,y)$p.value,a),xlab))+
  geom_abline(intercept = coef(testb)[1], slope = coef(testb)[2])+
  scale_colour_discrete (guide=FALSE)+#name=name_tri) +
  stat_smooth(method='lm',se=FALSE)+
  theme_bw()
sink(file=paste(dir,title,"_lm_rlm_coeff-corr_vfvi_vi.txt",sep=""))
summary(test)
summary(testb)
cor.test(x,y)
sink()
rm(test,testb,y,x)

## p3 dv vs vi ####
xlab=" "
for (i in cond){
  x=df$vi[which(df[K]==i & is.na(df$vi)==F & is.na(df$dv)==F)]
  y=df$dv[which(df[K]==i & is.na(df$vi)==F & is.na(df$dv)==F)]
  test=lm(y~x)
  testb=rlm(y~x)
  xlab=paste(xlab,"\n",i,", lm: R2=",signif(summary(test)$r.squared,a),", p=",signif(anova(test)$'Pr(>F)',a),",rlm: a=",signif(coef(testb)[2],a), ", b=", signif(coef(testb)[1],a), ", n= ",length(x),sep="")
  print(i)
  rm(x,y,test,testb)
}
x=df$vi[which(is.na(df$vi)==FALSE & is.na(df$dv)==FALSE)]
y=df$dv[which(is.na(df$vi)==FALSE & is.na(df$dv)==FALSE)]
test=lm(y~x)
testb=rlm(y~x)
p3 <- ggplot(df,aes(x=vi,y=dv, color=factor(tri))) +
  geom_point() +
  labs(title= title, y="DV [µm3]", 
       x=paste("Vol_i [µm3]","\n","lm: R2=",signif(summary(test)$r.squared,2),", p=",signif(anova(test)$'Pr(>F)',2),
               ", n=",length(df$dv[which(is.na(df$vi)==FALSE & is.na(df$dv)==FALSE)]),
               ", rlm: a=",signif(coef(testb)[2],2), ", b=", signif(coef(testb)[1],2),
               "\n Pearson: r=", signif(cor(x,y),a), "p=",signif(cor.test(x,y)$p.value,a),xlab))+
  scale_colour_discrete (guide=FALSE)+#name="name_tri") +
  geom_abline(intercept = coef(testb)[1], slope = coef(testb)[2])+
  stat_smooth(method='lm',se=FALSE)+
  theme_bw()
sink(file=paste(dir,title,"_lm_rlm_coeff-corr_dv_vi.txt",sep=""))
summary(test)
summary(testb)
cor.test(x,y)
sink()
rm(test,testb,y,x)

## p4 dt vs vi ####
xlab=" "
for (i in cond){
  x=df$vi[which(df[K]==i & is.na(df$vi)==FALSE & is.na(df$dt)==FALSE)]
  y=df$dt[which(df[K]==i & is.na(df$vi)==FALSE & is.na(df$dt)==FALSE)]
  test=lm(y~x)
  testb=rlm(y~x)
  xlab=paste(xlab,"\n",i,", lm: R2=",signif(summary(test)$r.squared,a),", p=",signif(anova(test)$'Pr(>F)',a),",rlm: a=",signif(coef(testb)[2],a), ", b=", signif(coef(testb)[1],a), ", n= ",length(x),sep="")
  print(i)
  rm(x,y,test,testb)
}
x=df$vi[which(is.na(df$vi)==FALSE & is.na(df$dt)==FALSE)]
y=df$dt[which(is.na(df$vi)==FALSE & is.na(df$dt)==FALSE)]
test=lm(df$dt~df$vi)
testb=rlm(y~x)
p4 <- ggplot(df,aes(x=vi,y=dt, color=factor(tri))) +
  geom_point() +
  labs(title= title, y="DT [hrs]", 
       x=paste("Vol_i [µm3]","\n","lm: R2=",signif(summary(test)$r.squared,a),", p=",signif(anova(test)$'Pr(>F)',a),
               ", n=",length(df$dt[which(is.na(df$vi)==FALSE & is.na(df$dt)==FALSE)]),
               ", rlm: a=",signif(coef(testb)[2],2), ", b=", signif(coef(testb)[1],2),
               "\n Pearson: r=", signif(cor(x,y),a), "p=",signif(cor.test(x,y)$p.value,a),xlab))+
  scale_colour_discrete (guide=FALSE)+#name="name_tri") +
  geom_abline(intercept = coef(testb)[1], slope = coef(testb)[2])+
  stat_smooth(method='lm',se=FALSE)+
  theme_bw()
sink(file=paste(dir,title,"_lm_rlm_coeff-corr_dt_vi.txt",sep=""))
summary(test)
summary(testb)
cor.test(x,y)
sink()
rm(test,testb,y,x)

## p5 vf vs vi ####
xlab=" "
for (i in cond){
  x=df$vi[which(df[K]==i & is.na(df$dv)==FALSE)]
  y=df$vf[which(df[K]==i & is.na(df$dv)==FALSE)]
  z=df$dv[which(df[K]==i & is.na(df$dv)==FALSE)]
  test=lm(y~x)
  testb=rlm(y~x)
  xlab=paste(xlab,"\n",i,", lm: R2=",signif(summary(test)$r.squared,a),", p=",signif(anova(test)$'Pr(>F)',a),
             ",rlm: a=",signif(coef(testb)[2],a), ", b=", signif(coef(testb)[1],a),
             ", m=", signif(mean(z),a),
             ", n= ",length(x),sep="")
  rm(x,y,test,testb)
}
x=df$vi[which(is.na(df$vi)==FALSE & is.na(df$vf)==FALSE)]
y=df$vf[which(is.na(df$vi)==FALSE & is.na(df$vf)==FALSE)]
test=lm(df$vf~df$vi)
testb=rlm(y~x)
p5 <- ggplot(df,aes(x=vi,y=vf, color=factor(tri))) +
  geom_point() +
  labs(title= title, y="Vf [µm3]", 
       x=paste("Vol_i [µm3]","\n","lm: R2=",signif(summary(test)$r.squared,a),", p=",signif(anova(test)$'Pr(>F)',a),
               ", n=",length(x),
               "a=",signif(coef(test)[2],a), ", b=", signif(coef(test)[1],a),
               ", rlm: a=",signif(coef(testb)[2],2), ", b=", signif(coef(testb)[1],2),
               "\n Pearson: r=", signif(cor(x,y),a), "p=",signif(cor.test(x,y)$p.value,a),
               ", mean DV=",signif(mean(df$dv,na.rm=TRUE),a),xlab))+
  scale_colour_discrete (guide=FALSE)+#name="name_tri") +
  geom_abline(intercept = coef(testb)[1], slope = coef(testb)[2])+
  stat_smooth(method='lm',se=FALSE)+
  theme_bw()
sink(file=paste(dir,title,"_lm_rlm_coeff-corr_vf_vi.txt",sep=""))
summary(test)
summary(testb)
cor.test(x,y)
sink()
rm(test,testb,y,x)

## p6 dv vs dt ####
xlab=" "
for (i in cond){
  x=df$dt[which(df[K]==i & is.na(df$dt)==FALSE & is.na(df$dv)==FALSE)]
  y=df$dv[which(df[K]==i & is.na(df$dt)==FALSE & is.na(df$dv)==FALSE)]
  test=lm(y~x)
  testb=rlm(y~x)
  xlab=paste(xlab,"\n",i,", lm: R2=",signif(summary(test)$r.squared,a),", p=",signif(anova(test)$'Pr(>F)',a),",rlm: a=",signif(coef(testb)[2],a), ", b=", signif(coef(testb)[1],a), ", n= ",length(x),sep="")
  print(i)
  rm(x,y,test,testb)
}
x=df$dt[which(is.na(df$dt)==FALSE & is.na(df$dv)==FALSE)]
y=df$dv[which(is.na(df$dt)==FALSE & is.na(df$dv)==FALSE)]
test=lm(y~x)
testb=rlm(y~x)
p6 <- ggplot(df,aes(x=dt,y=dv, color=factor(tri))) +
  geom_point() +
  labs(title= title, y="DV [µm3]", 
       x=paste("DT [hrs]","\n","R2=",signif(summary(test)$r.squared,a),", p=",signif(anova(test)$'Pr(>F)',a),
               ", n=",length(df$dt[which(is.na(df$dt)==FALSE & is.na(df$dv)==FALSE)]),
               ", rlm: a=",signif(coef(testb)[2],2), ", b=", signif(coef(testb)[1],2),
               "\n Pearson: r=", signif(cor(x,y),a), "p=",signif(cor.test(x,y)$p.value,a),xlab))+
  scale_colour_discrete (guide=FALSE)+#name="name_tri") +
  #   geom_abline(intercept = coef(testb)[1], slope = coef(testb)[2])+
  stat_smooth(method='lm',se=FALSE)+
  theme_bw()
sink(file=paste(dir,title,"_lm_rlm_coeff-corr_dv_dt.txt",sep=""))
summary(test)
summary(testb)
cor.test(x,y)
sink()
rm(test,testb,y,x)

## p7 dv/v vs dt ####
xlab=" "
for (i in cond){
  x=df$dt[which(df[K]==i & is.na(df$dt)==FALSE & is.na(df$dvvi)==FALSE)]
  y=df$dvvi[which(df[K]==i & is.na(df$dt)==FALSE & is.na(df$dvvi)==FALSE)]
  test=lm(y~x)
  testb=rlm(y~x)
  xlab=paste(xlab,"\n",i,", lm: R2=",signif(summary(test)$r.squared,a),", p=",signif(anova(test)$'Pr(>F)',a),",rlm: a=",signif(coef(testb)[2],a), ", b=", signif(coef(testb)[1],a), ", n= ",length(x),sep="")
  print(i)
  rm(x,y,test,testb)
}
x=df$dt[which(is.na(df$dt)==FALSE & is.na(df$dvvi)==FALSE)]
y=df$dvvi[which(is.na(df$dt)==FALSE & is.na(df$dvvi)==FALSE)]
test=lm(y~x)
testb=rlm(y~x)
p7 <- ggplot(df,aes(x=dt,y=dvvi,color=factor(tri))) +
  geom_point() +
  labs(title= title, y="DV/Vi %", 
       x=paste("DT [hrs]","\n","R2=",signif(summary(test)$r.squared,a),", p=",signif(anova(test)$'Pr(>F)',a),
               ", n=",length(x),
               ", rlm: a=",signif(coef(testb)[2],2), ", b=", signif(coef(testb)[1],2),
               "\n Pearson: r=", signif(cor(x,y),a), "p=",signif(cor.test(x,y)$p.value,a),xlab))+
  scale_colour_discrete (guide=FALSE)+#name="name_tri") +
  geom_abline(intercept = coef(testb)[1], slope = coef(testb)[2])+
  stat_smooth(method="lm", se=FALSE)+
  theme_bw()
sink(file=paste(dir,title,"_lm_rlm_coeff-corr_dvvi_dt.txt",sep=""))
summary(test)
summary(testb)
cor.test(x,y)
sink()
rm(test,testb,y,x)

## p8 dv/v vs vi ####
xlab=" "
for (i in cond){
  x=df$vi[which(df[K]==i & is.na(df$vi)==FALSE & is.na(df$dvvi)==FALSE)]
  y=df$dvvi[which(df[K]==i & is.na(df$vi)==FALSE & is.na(df$dvvi)==FALSE)]
  test=lm(y~x)
  testb=rlm(y~x)
  xlab=paste(xlab,"\n",i,", lm: R2=",signif(summary(test)$r.squared,a),", p=",signif(anova(test)$'Pr(>F)',a),",rlm: a=",signif(coef(testb)[2],a), ", b=", signif(coef(testb)[1],a), ", n= ",length(x),sep="")
  print(i)
  rm(x,y,test,testb)
}
x=df$vi[which(is.na(df$vi)==FALSE & is.na(df$dvvi)==FALSE)]
y=df$dvvi[which(is.na(df$vi)==FALSE & is.na(df$dvvi)==FALSE)]
test=lm(y~x)
testb=rlm(y~x)
p8 <- ggplot(df,aes(x=vi,y=dvvi, color=factor(tri))) +
  geom_point() +
  labs(title= title, y="DV/Vi %", 
       x=paste("vi [µm3]","\n","R2=",signif(summary(test)$r.squared,a),", p=",signif(anova(test)$'Pr(>F)',a),
               ", n=",length(x),
               ", rlm: a=",signif(coef(testb)[2],2), ", b=", signif(coef(testb)[1],2),
               "\n Pearson: r=", signif(cor(x,y),a), "p=",signif(cor.test(x,y)$p.value,a),xlab))+
  scale_colour_discrete (guide=FALSE) +
  geom_abline(intercept = coef(testb)[1], slope = coef(testb)[2])+
  stat_smooth(method="lm", se=FALSE)+
  theme_bw()
sink(file=paste(dir,title,"_lm_rlm_coeff-corr_dvvi_vi.txt",sep=""))
summary(test)
summary(testb)
cor.test(x,y)
sink()
rm(test,testb,y,x)

## p9 dt vs dv/dt ####
xlab=" "
for (i in cond){
  x=df$dvdt[which(df[K]==i & is.na(df$dvdt)==FALSE & is.na(df$dt)==FALSE)]
  y=df$dt[which(df[K]==i & is.na(df$dvdt)==FALSE & is.na(df$dt)==FALSE)]
  test=lm(y~x)
  testb=rlm(y~x)
  xlab=paste(xlab,"\n",i,", lm: R2=",signif(summary(test)$r.squared,a),", p=",signif(anova(test)$'Pr(>F)',a),",rlm: a=",signif(coef(testb)[2],a), ", b=", signif(coef(testb)[1],a), ", n= ",length(x),sep="")
  print(i)
  rm(x,y,test,testb)
}
x=df$dvdt[which(is.na(df$dvdt)==FALSE & is.na(df$dt)==FALSE)]
y=df$dt[which(is.na(df$dvdt)==FALSE & is.na(df$dt)==FALSE)]
test=lm(y~x)
testb=rlm(y~x)
p9 <- ggplot(df,aes(x=dvdt,y=dt, color=factor(tri))) +
  geom_point() +
  labs(title= title, y="DT [hrs]", 
       x=paste("DV/DT [µm3/hrs]","\n","R2=",signif(summary(test)$r.squared,a),", p=",signif(anova(test)$'Pr(>F)',a),
               ", n=",length(x),
               ", rlm: a=",signif(coef(testb)[2],2), ", b=", signif(coef(testb)[1],2),
               "\n Pearson: r=", signif(cor(x,y),a), "p=",signif(cor.test(x,y)$p.value,a),xlab))+
  scale_colour_discrete (guide=FALSE) +
  geom_abline(intercept = coef(testb)[1], slope = coef(testb)[2])+
  stat_smooth(method="lm", se=FALSE)+
  theme_bw()
sink(file=paste(dir,title,"_lm_rlm_coeff-corr_dt_dvdt.txt",sep=""))
summary(test)
summary(testb)
cor.test(x,y)
sink()
rm(test,testb,y,x)

## make figure #####
tiff(paste(dir,title,"_abs_val_",name_tri,".tiff",sep=""), width=1000,height=1200)
grid.arrange(p2,p1,p3,p4,p5,p6,p7,p8,ncol=2)
dev.off()

tiff(paste(dir,title,"_abs_val_",name_tri,".tiff",sep=""), width=1000,height=1300)
grid.arrange(p2,p1,p3,p4,p5,p6,p7,p8,p9,gt,ncol=2)
dev.off()
# rm(df)

####### comparison in lineages, correlations between sisters/mothers/cousins ######
#### 1) ratios between sisters for vf, vi, dv, dt ####
m = 1.6
b = 240
cond=c("150213","150417")
df <- data1
### tri ####
vf_1 = NULL
vf_2 = NULL
vi_1 = NULL
vi_2 = NULL
dt_1 = NULL
dt_2 = NULL
dvdt_1 = NULL
dvdt_2 = NULL
lineage = NULL
number_in_lineage = NULL
cond_date=NULL
for (k in cond){
  for (i in unique(df$lineage[which(df$cond_date==k)])){
    for (j in unique(df$number_in_lineage[which(df$lineage==i & df$number_in_lineage %%2==F & df$cond_date==k)])){
      if (m %in% df$cycle_stage[which(df$lineage==i & df$number_in_lineage==j & df$cond_date==k)]
          && m %in% df$cycle_stage[which(df$lineage==i & df$number_in_lineage==j+1 & df$cond_date==k)]
          && b %in% df$cycle_stage[which(df$lineage==i & df$number_in_lineage==j & df$cond_date==k)]
          && b %in% df$cycle_stage[which(df$lineage==i & df$number_in_lineage==j+1 & df$cond_date==k)]
          && is.na(df$vol[which(df$lineage==i & df$number_in_lineage==j & df$cycle_stage==m & df$cond_date==k)])==F
          && is.na(df$vol[which(df$lineage==i & df$number_in_lineage==j+1 & df$cycle_stage==m & df$cond_date==k)])==F
          && is.na(df$vol[which(df$lineage==i & df$number_in_lineage==j & df$cycle_stage==b & df$cond_date==k)])==F
          && is.na(df$vol[which(df$lineage==i & df$number_in_lineage==j+1 & df$cycle_stage==b & df$cond_date==k)])==F){
        vf_1 = c(vf_1,df$vol[which(df$lineage==i & df$number_in_lineage==j & df$cycle_stage==m & df$cond_date==k)])
        vf_2 = c(vf_2,df$vol[which(df$lineage==i & df$number_in_lineage==j+1 & df$cycle_stage==m & df$cond_date==k)])
        vi_1 = c(vi_1,df$vol[which(df$lineage==i & df$number_in_lineage==j & df$cycle_stage==b & df$cond_date==k)])
        vi_2 = c(vi_2,df$vol[which(df$lineage==i & df$number_in_lineage==j+1 & df$cycle_stage==b & df$cond_date==k)])
        dt_1 = c(dt_1,df$Time[which(df$lineage==i & df$number_in_lineage==j & df$cycle_stage==m & df$cond_date==k)]-
                   df$Time[which(df$lineage==i & df$number_in_lineage==j & df$cycle_stage==b & df$cond_date==k)])
        dt_2 = c(dt_2,df$Time[which(df$lineage==i & df$number_in_lineage==j+1 & df$cycle_stage==m & df$cond_date==k)]
                 - df$Time[which(df$lineage==i & df$number_in_lineage==j+1 & df$cycle_stage==b & df$cond_date==k)])
        dvdt_1 = c(dvdt_1, ((df$vol[which(df$lineage==i & df$number_in_lineage==j & df$cycle_stage==m & df$cond_date==k)]-
                              df$vol[which(df$lineage==i & df$number_in_lineage==j & df$cycle_stage==b & df$cond_date==k)])/
                              ((df$Time[which(df$lineage==i & df$number_in_lineage==j & df$cycle_stage==m & df$cond_date==k)]-
                                 df$Time[which(df$lineage==i & df$number_in_lineage==j & df$cycle_stage==b & df$cond_date==k)])/60)))
        dvdt_2=c(dvdt_2, ((df$vol[which(df$lineage==i & df$number_in_lineage==j+1 & df$cycle_stage==m & df$cond_date==k)]-
                             df$vol[which(df$lineage==i & df$number_in_lineage==j+1 & df$cycle_stage==b & df$cond_date==k)])/
                            ((df$Time[which(df$lineage==i & df$number_in_lineage==j+1 & df$cycle_stage==m & df$cond_date==k)]
                               - df$Time[which(df$lineage==i & df$number_in_lineage==j+1 & df$cycle_stage==b & df$cond_date==k)])/60)))
        lineage = c(lineage,i)
        number_in_lineage = c(number_in_lineage,j)
        cond_date=c(cond_date,k)
      }
    }
  }
}

### plots ####
df <- data.frame(vf_1,vf_2,vi_1,vi_2,dt_1,dt_2, lineage, number_in_lineage)
df$vf_ratio <- df$vf_2/df$vf_1 *100
df$vi_ratio <- df$vi_2/df$vi_1 *100
df$dv_ratio <- (df$vf_2-df$vi_2)/(df$vf_1-df$vi_1) *100
df$dt_ratio <- df$dt_2/df$dt_1 *100
rm(vf_1,vf_2,vi_1,vi_2,dt_1,dt_2, lineage, number_in_lineage)
# write.csv(df,paste(dir,title,"_sisters_df.csv"))

test=lm(df$vf_ratio~df$vi_ratio)
testb=rlm(df$vf_ratio~df$vi_ratio)
p1 <- ggplot(df,aes(x=vi_ratio,y=vf_ratio)) +
  geom_point() +
  labs(title= title, y="Vf ratio d1/d2 %", 
       x=paste("Vi ratio d1/d2 %","\n","R2=",signif(summary(test)$r.squared,a),", p=",signif(anova(test)$'Pr(>F)',a),
               ", n=",length(df$vf_ratio),
               ", rlm: a=",signif(coef(testb)[2],2), ", b=", signif(coef(testb)[1],2)
#                ,
#                "\n Pearson: r=", signif(cor(x,y),a), "p=",signif(cor.test(x,y)$p.value,a),xlab
               ))+
  geom_abline(intercept = coef(testb)[1], slope = coef(testb)[2])+
  theme_bw()
sink(file=paste(dir,title,"_ratio_sisters_lm_vf_vi.txt",sep=""))
summary(test)
summary(testb)
sink()
rm(test)

test=lm(df$dt_ratio~df$vi_ratio)
testb=rlm(df$dt_ratio~df$vi_ratio)
p2 <- ggplot(df,aes(x=vi_ratio,y=dt_ratio)) +
  geom_point() +
  labs(title= title, y="DT ratio d1/d2 %", 
       x=paste("Vi ratio d1/d2 %","\n","R2=",signif(summary(test)$r.squared,a),", p=",signif(anova(test)$'Pr(>F)',a),
               ", n=",length(df$vi_ratio),
               ", rlm: a=",signif(coef(testb)[2],2), ", b=", signif(coef(testb)[1],2)
#               , "\n Pearson: r=", signif(cor(x,y),a), "p=",signif(cor.test(x,y)$p.value,a),xlab
               ))+
  #scale_colour_discrete (name="lineage", labels=c("daughter 1","daughter 2")) +
  theme_bw()
sink(file=paste(dir,title,"_ratio_sisters_lm_dt_vi.txt",sep=""))
summary(test)
summary(testb)
sink()
rm(test)

test=lm(df$dv_ratio~df$vi_ratio)
testb=rlm(df$dv_ratio~df$vi_ratio)
p3 <- ggplot(df,aes(x=vi_ratio,y=dv_ratio)) +
  geom_point() +
  labs(title= title, y="DV ratio d1/d2 %", 
       x=paste("Vi ratio d1/d2 %","\n","R2=",signif(summary(test)$r.squared,a),", p=",signif(anova(test)$'Pr(>F)',a),
               ", n=",length(df$vi_ratio),
               ", rlm: a=",signif(coef(testb)[2],2), ", b=", signif(coef(testb)[1],2)
#                ,
#                "\n Pearson: r=", signif(cor(x,y),a), "p=",signif(cor.test(x,y)$p.value,a),xlab
               ))+
  #scale_colour_discrete (name="lineage", labels=c("daughter 1","daughter 2")) +
  theme_bw()
sink(file=paste(dir,title,"_ratio_sisters_lm_dv_vi.txt",sep=""))
summary(test)
summary(testb)
sink()
rm(test)

test=lm(df$dv_ratio~df$dt_ratio)
testb=rlm(df$dv_ratio~df$dt_ratio)
p4 <- ggplot(df,aes(x=dt_ratio,y=dv_ratio)) +
  geom_point() +
  labs(title= title, y="DV ratio d1/d2 %", 
       x=paste("DT ratio d1/d2 %","\n","R2=",signif(summary(test)$r.squared,a),", p=",signif(anova(test)$'Pr(>F)',a),
               ", n=",length(df$dv_ratio),
               ", rlm: a=",signif(coef(testb)[2],2), ", b=", signif(coef(testb)[1],2)
#                ,
#                "\n Pearson: r=", signif(cor(x,y),a), "p=",signif(cor.test(x,y)$p.value,a),xlab
               ))+
  #scale_colour_discrete (name="lineage", labels=c("daughter 1","daughter 2")) +
  theme_bw()
sink(file=paste(dir,title,"_ratio_sisters_lm_dv_dt.txt",sep=""))
summary(test)
summary(testb)
sink()
rm(test)

hist3 <- ggplot(df,aes(x=vi_ratio))+#,fill=as.factor(cond_date)))+
  geom_histogram()+
  theme_bw()+
  labs(title=title)

### save figure ####
tiff(paste(dir,title,"_ratios_sisters_dv_dt_vi_vf.tiff",sep=""), width=750,height=500)
grid.arrange(p1,p2,p3,p4,ncol=2)
dev.off()
rm(p1,p2,p3,p4)

tiff(paste(dir,title,"_ratios_sisters_vi_distrib.tiff",sep=""), width=175,height=250)
grid.arrange(hist3)
dev.off()
#### 2) comparison big/small ratios decrease asymmetry ###########################
### start with df from just above and make data-frame with big and small sisters ####

vi_big <- NULL
vi_small <-NULL
vf_big <- NULL
vf_small <- NULL
DT_big <- NULL
DT_small <- NULL
DV_big <- NULL
DV_small <- NULL

colnames(df)[1:4]=c("d1_vol_f","d2_vol_f","d1_vol_i","d2_vol_i")

for (i in c(1:nrow(df))){
  if (is.na(df$d1_vol_i[i])==F
      && is.na(df$d1_vol_f[i])==F
      && is.na(df$d2_vol_i[i])==F
      && is.na(df$d2_vol_f[i])==F) {
    if (df$d1_vol_i[i] <= df$d2_vol_i[i]){
      vi_big <- c(vi_big,df$d2_vol_i[i])
      vi_small <- c(vi_small,df$d1_vol_i[i])
      vf_big <- c(vf_big,df$d2_vol_f[i])
      vf_small <- c(vf_small,df$d1_vol_f[i])
    } else {
      vi_big <- c(vi_big,df$d1_vol_i[i])
      vi_small <- c(vi_small,df$d2_vol_i[i])
      vf_big <- c(vf_big,df$d1_vol_f[i])
      vf_small <- c(vf_small,df$d2_vol_f[i])
    }
  } else {}
}

ratio_i <- vi_big/vi_small *100
ratio_f <- vf_big/vf_small *100
initial <- rep(1,length(ratio_i))
final <- rep(2,length(ratio_f))
d <- data.frame(c(ratio_i,ratio_f),c(initial,final)) 
names(d)=c("ratio","stage")
write.csv(data.frame(ratio_i,ratio_f),paste(dir,title,"_big_small_df.csv"))

p <- ggplot(d,aes(factor(d$stage),d$ratio))+
  geom_boxplot()+
  labs(title=title,
       x=paste("t test, p=",signif(t.test(ratio_i,ratio_f)$p.value,a),", n=",length(ratio_i),
               "\nratio_i, m=", signif(mean(ratio_i),a),",sd=",signif(sd(d$ratio[which(d$stage==1)]),a),
               "\nratio_f, m=",signif(mean(ratio_f),a),",sd=",signif(sd(d$ratio[which(d$stage==2)]),a)),
       y="ratio volume big/small sister %")+
  scale_x_discrete (labels=c("initial","final"))+
  theme_bw()+
  theme(axis.text.x=element_text(size=12),
        text=element_text(size=12))

d <- data.frame(c(df$vf_ratio,df$vi_ratio),c(rep(2,length(df$vf_ratio)),rep(1,length(df$vi_ratio))))
colnames(d)=c("ratio","stage")
p1 <- ggplot(df,aes(factor(d$stage),d$ratio))+
  geom_boxplot()+
  labs(title=title,
       x=paste("t test, p=",signif(t.test(ratio_i,ratio_f)$p.value,a),", n=",length(ratio_i),
               "\nratio_i, m=", signif(mean(ratio_i),a),",sd=",signif(sd(d$ratio[which(d$stage==1)]),a),
                "\nratio_f, m=",signif(mean(ratio_f),a),",sd=",signif(sd(d$ratio[which(d$stage==2)]),a)),
       y="ratio volume in pair of sisters %")+
  scale_x_discrete (labels=c("initial","final"))+
  theme_bw()+
  theme(axis.text.x=element_text(size=12),
        text=element_text(size=12))

tiff(paste(dir,title,"_ratio_decrease_sister_boxplot.tiff",sep=""), width=500,height=250)
grid.arrange(p,p1,ncol=2)
dev.off()

#### 3) lineages for vf and dt : comparison mother/daughter sister/sister cousin/cousin ####
### enter the following parameters : ####
m=1.6
b=240
cond=c("131122","131220")

## ENORMOUS loop to sort the cells  ####
K= which(colnames(df)=="cond_date")
cond=c("131122","131220")
vf_md_m = NULL
vf_md_d = NULL
dt_md_m = NULL
dt_md_d = NULL
dvdt_md_m = NULL
dvdt_md_d = NULL
vf_mgd_m = NULL
vf_mgd_gd = NULL
dt_mgd_m = NULL
dt_mgd_gd = NULL
dvdt_mgd_m = NULL
dvdt_mgd_gd = NULL
vf_mggd_m = NULL
vf_mggd_ggd = NULL
dt_mggd_m = NULL
dt_mggd_ggd = NULL
dvdt_mggd_m = NULL
dvdt_mggd_ggd = NULL


for (c in cond){
  print(c)
  for (l in (unique(df$lineage[which(df[K]==c)]))){
    print("l")
    print(l)
    for (g in c(0:(max(df$generation[which(df[K]==c & df$lineage==l)])-1))){
      print("g")
      print(g)
      for (i in unique(df$number_in_lineage[which(df[K]==c & df$lineage==l & df$generation==g)])){
        print("i")
        print(i)
        for (j in unique(df$number_in_lineage[which(df[K]==c & df$lineage==l & df$generation==(g+1))])){ ## mother/daughter
          if (j %in% c((i*2):(i*2+1))
              && m %in% df$cycle_stage[which(df[K]==c & df$lineage==l & df$number_in_lineage==j)]
              && m %in% df$cycle_stage[which(df[K]==c & df$lineage==l & df$number_in_lineage==i)]){ ## mother and daughter exist and finished a cell cycle
            if (is.na(df$vol[which(df[K]==c & df$lineage==l & df$number_in_lineage==j & df$cycle_stage == m)])==F
                && is.na(df$vol[which(df[K]==c & df$lineage==l & df$number_in_lineage==i & df$cycle_stage == m)])==F){
              vf_md_m = c(vf_md_m,df$vol[which(df[K]==c & df$lineage==l & df$number_in_lineage==i & df$cycle_stage==m)])
              vf_md_d = c(vf_md_d,df$vol[which(df[K]==c & df$lineage==l & df$number_in_lineage==j & df$cycle_stage==m)])
              print(j)
            } else {} ## vectors with vf for mother and daughter
            if (b %in% df$cycle_stage[which(df[K]==c & df$lineage==l & df$number_in_lineage==j)]
                && b %in% df$cycle_stage[which(df[K]==c & df$lineage==l & df$number_in_lineage==i)]){
              dt_md_d = c(dt_md_d,(df$Time[which(df[K]==c & df$lineage==l & df$number_in_lineage==j & df$cycle_stage == m)]-
                                     df$Time[which(df[K]==c & df$lineage==l & df$number_in_lineage==j & df$cycle_stage == b)])/60)
              dt_md_m = c(dt_md_m,(df$Time[which(df[K]==c & df$lineage==l & df$number_in_lineage==(i) & df$cycle_stage == m)]-
                                     df$Time[which(df[K]==c & df$lineage==l & df$number_in_lineage==(i) & df$cycle_stage == b)])/60)
              print(j)
              if (is.na(df$vol[which(df[K]==c & df$lineage==l & df$number_in_lineage==j & df$cycle_stage == m)])==F
                  && is.na(df$vol[which(df[K]==c & df$lineage==l & df$number_in_lineage==i & df$cycle_stage == m)])==F
                  && is.na(df$vol[which(df[K]==c & df$lineage==l & df$number_in_lineage==j & df$cycle_stage == b)])==F
                  && is.na(df$vol[which(df[K]==c & df$lineage==l & df$number_in_lineage==i & df$cycle_stage == b)])==F){
                dt_d=(df$Time[which(df[K]==c & df$lineage==l & df$number_in_lineage==j & df$cycle_stage == m)]-
                        df$Time[which(df[K]==c & df$lineage==l & df$number_in_lineage==j & df$cycle_stage == b)])/60
                dt_m=(df$Time[which(df[K]==c & df$lineage==l & df$number_in_lineage==i & df$cycle_stage == m)]-
                        df$Time[which(df[K]==c & df$lineage==l & df$number_in_lineage==i & df$cycle_stage == b)])/60
                dv_d =df$vol[which(df[K]==c & df$lineage==l & df$number_in_lineage==j & df$cycle_stage == m)]-
                  df$vol[which(df[K]==c & df$lineage==l & df$number_in_lineage==j & df$cycle_stage == b)]
                dv_m =df$vol[which(df[K]==c & df$lineage==l & df$number_in_lineage==i & df$cycle_stage == m)]-
                  df$vol[which(df[K]==c & df$lineage==l & df$number_in_lineage==i & df$cycle_stage == b)]
                dvdt_md_d = c(dvdt_md_d,dv_d/dt_d)
                dvdt_md_m = c(dvdt_md_m,dv_m/dt_m)
                print(j)
              } else {}
            } else {} ##vectors with dvdt for mother anddaughter and dt for mother and daughter
          } else {}
        }
        for (k in unique(df$number_in_lineage[which(df[K]==c & df$lineage==l & df$generation==(g+2))])){ ## mother/grand-daughter
          if (k %in% c((i*4):(i*4+3))
              && m %in% df$cycle_stage[which(df[K]==c & df$lineage==l & df$number_in_lineage==k)]){ ## grand-daughter exist and finished cell cycle
            if (is.na(df$vol[which(df[K]==c & df$lineage==l & df$number_in_lineage==k & df$cycle_stage == m)])==F
                && is.na(df$vol[which(df[K]==c & df$lineage==l & df$number_in_lineage==i & df$cycle_stage == m)])==F){
              vf_mgd_m = c(vf_mgd_m,df$vol[which(df[K]==c & df$lineage==l & df$number_in_lineage==i & df$cycle_stage==m)])
              vf_mgd_gd = c(vf_mgd_gd,df$vol[which(df[K]==c & df$lineage==l & df$number_in_lineage==k & df$cycle_stage==m)])
              print("test1")
              print(k)
            } else {} ## vectors with vf for mother and grand-daughter
            if (b %in% df$cycle_stage[which(df[K]==c & df$lineage==l & df$number_in_lineage==k)]
                && b %in% df$cycle_stage[which(df[K]==c & df$lineage==l & df$number_in_lineage==i)]){ ## vectors with dt for mother and grand-daughter
              dt_mgd_gd = c(dt_mgd_gd,(df$Time[which(df[K]==c & df$lineage==l & df$number_in_lineage==k & df$cycle_stage == m)]-
                                     df$Time[which(df[K]==c & df$lineage==l & df$number_in_lineage==k & df$cycle_stage == b)])/60)
              dt_mgd_m = c(dt_mgd_m,(df$Time[which(df[K]==c & df$lineage==l & df$number_in_lineage==(i) & df$cycle_stage == m)]-
                                     df$Time[which(df[K]==c & df$lineage==l & df$number_in_lineage==(i) & df$cycle_stage == b)])/60)
              print("test2")
              print(k)
              if (is.na(df$vol[which(df[K]==c & df$lineage==l & df$number_in_lineage==k & df$cycle_stage == m)])==F
                  && is.na(df$vol[which(df[K]==c & df$lineage==l & df$number_in_lineage==i & df$cycle_stage == m)])==F
                  && is.na(df$vol[which(df[K]==c & df$lineage==l & df$number_in_lineage==k & df$cycle_stage == b)])==F
                  && is.na(df$vol[which(df[K]==c & df$lineage==l & df$number_in_lineage==i & df$cycle_stage == b)])==F){ ## vectors with dvdt for mother and grand-daughter
                dt_d=(df$Time[which(df[K]==c & df$lineage==l & df$number_in_lineage==k & df$cycle_stage == m)]-
                        df$Time[which(df[K]==c & df$lineage==l & df$number_in_lineage==k & df$cycle_stage == b)])/60
                dt_m=(df$Time[which(df[K]==c & df$lineage==l & df$number_in_lineage==i & df$cycle_stage == m)]-
                        df$Time[which(df[K]==c & df$lineage==l & df$number_in_lineage==i & df$cycle_stage == b)])/60
                dv_d =df$vol[which(df[K]==c & df$lineage==l & df$number_in_lineage==k & df$cycle_stage == m)]-
                  df$vol[which(df[K]==c & df$lineage==l & df$number_in_lineage==k & df$cycle_stage == b)]
                dv_m =df$vol[which(df[K]==c & df$lineage==l & df$number_in_lineage==i & df$cycle_stage == m)]-
                  df$vol[which(df[K]==c & df$lineage==l & df$number_in_lineage==i & df$cycle_stage == b)]
                dvdt_mgd_gd = c(dvdt_md_d,dv_d/dt_d)
                dvdt_mgd_m = c(dvdt_md_m,dv_m/dt_m)
                print("test3")
                print(k)
              } else {} 
            } else {}
          } else {}
        }
        for (o in unique(df$number_in_lineage[which(df[K]==c & df$lineage==l & df$generation==(g+3))])){ ## mother/grand-grand-daughter
          if (o %in% c((i*8):(i*8+7))
              && m %in% df$cycle_stage[which(df[K]==c & df$lineage==l & df$number_in_lineage==o)]){ ## grand-grand daughter exists
            if (is.na(df$vol[which(df[K]==c & df$lineage==l & df$number_in_lineage==o & df$cycle_stage == m)])==F
                && is.na(df$vol[which(df[K]==c & df$lineage==l & df$number_in_lineage==i & df$cycle_stage == m)])==F){ ## vf exists for both
              vf_mgd_m = c(vf_mgd_m,df$vol[which(df[K]==c & df$lineage==l & df$number_in_lineage==i & df$cycle_stage==m)])
              vf_mgd_gd = c(vf_mgd_gd,df$vol[which(df[K]==c & df$lineage==l & df$number_in_lineage==o & df$cycle_stage==m)])
              print("test1")
              print(o)
            } else {}
            if (b %in% df$cycle_stage[which(df[K]==c & df$lineage==l & df$number_in_lineage==o)]
                && b %in% df$cycle_stage[which(df[K]==c & df$lineage==l & df$number_in_lineage==i)]){ ## vectors with dt for mother and grand-daughter
              dt_mggd_ggd = c(dt_mggd_ggd,(df$Time[which(df[K]==c & df$lineage==l & df$number_in_lineage==o & df$cycle_stage == m)]-
                                         df$Time[which(df[K]==c & df$lineage==l & df$number_in_lineage==o & df$cycle_stage == b)])/60)
              dt_mggd_m = c(dt_mggd_m,(df$Time[which(df[K]==c & df$lineage==l & df$number_in_lineage==(i) & df$cycle_stage == m)]-
                                       df$Time[which(df[K]==c & df$lineage==l & df$number_in_lineage==(i) & df$cycle_stage == b)])/60)
              print("test2")
              print(o)
              if (is.na(df$vol[which(df[K]==c & df$lineage==l & df$number_in_lineage==o & df$cycle_stage == m)])==F
                  && is.na(df$vol[which(df[K]==c & df$lineage==l & df$number_in_lineage==i & df$cycle_stage == m)])==F
                  && is.na(df$vol[which(df[K]==c & df$lineage==l & df$number_in_lineage==o & df$cycle_stage == b)])==F
                  && is.na(df$vol[which(df[K]==c & df$lineage==l & df$number_in_lineage==i & df$cycle_stage == b)])==F){ ## vectors with dvdt for mother and grand-grand-daughter
                dt_d=(df$Time[which(df[K]==c & df$lineage==l & df$number_in_lineage==o & df$cycle_stage == m)]-
                        df$Time[which(df[K]==c & df$lineage==l & df$number_in_lineage==o & df$cycle_stage == b)])/60
                dt_m=(df$Time[which(df[K]==c & df$lineage==l & df$number_in_lineage==i & df$cycle_stage == m)]-
                        df$Time[which(df[K]==c & df$lineage==l & df$number_in_lineage==i & df$cycle_stage == b)])/60
                dv_d =df$vol[which(df[K]==c & df$lineage==l & df$number_in_lineage==o & df$cycle_stage == m)]-
                  df$vol[which(df[K]==c & df$lineage==l & df$number_in_lineage==o & df$cycle_stage == b)]
                dv_m =df$vol[which(df[K]==c & df$lineage==l & df$number_in_lineage==i & df$cycle_stage == m)]-
                  df$vol[which(df[K]==c & df$lineage==l & df$number_in_lineage==i & df$cycle_stage == b)]
                dvdt_mggd_ggd = c(dvdt_md_d,dv_d/dt_d)
                dvdt_mggd_m = c(dvdt_md_m,dv_m/dt_m)
                print("test3")
                print(o)
              } else {} 
            } else {}
          }
        }
      }
    }
  }
}

## ENORMOUS loop to sort for sisters and cousins ####
# elements in the loop ####
K= which(colnames(df)=="cond_date")
cond=c("131122","131220")
s1_dt = NULL
s2_dt = NULL
s1_number_in_lineage_dt = NULL
s2_number_in_lineage_dt = NULL
sister_lineage_dt = NULL
c1_dt = NULL
c2_dt = NULL
c1_number_in_lineage_dt = NULL
c2_number_in_lineage_dt = NULL
cousin_generation_dt = NULL
cousin_lineage_dt = NULL

s1_vf = NULL
s2_vf = NULL
s1_number_in_lineage_vf = NULL
s2_number_in_lineage_vf = NULL
sister_lineage_vf = NULL
c1_vf = NULL
c2_vf = NULL
c1_number_in_lineage_vf = NULL
c2_number_in_lineage_vf = NULL
cousin_generation_vf = NULL
cousin_lineage_vf = NULL

s1_dvdt = NULL
s2_dvdt= NULL
s1_number_in_lineage_dvdt = NULL
s2_number_in_lineage_dvdt = NULL
sister_lineage_dvdt = NULL
c1_dvdt = NULL
c2_dvdt = NULL
c1_number_in_lineage_dvdt = NULL
c2_number_in_lineage_dvdt = NULL
cousin_generation_dvdt = NULL
cousin_lineage_dvdt = NULL

# loop ####
for (c in cond) {
  for (i in unique(df$lineage[which(df[K]==c)])){
    for (j in unique(df$generation[which(df[K]==c & df$lineage==i)])) {
      #print("generation")
      #print(j)
      for (k in unique(df$number_in_lineage[which(df[K]==c & df$lineage==i & df$generation==j)])){
        l=k+1
        while (max(unique(df$number_in_lineage[which(df[K]==c & df$lineage==i & df$generation==j)]))+1 > l & l > k) {
          if (m %in% df$cycle_stage[which(df[K]==c & df$lineage==i & df$number_in_lineage==k)]
              && m %in% df$cycle_stage[which(df[K]==c & df$lineage==i & df$number_in_lineage==l)]) { ## both divided : dt can be calculated
            if (k%%2==F && l==k+1){
              s1_dt <- c(s1_dt, df$Time[which(df[K]==c & df$lineage==i & df$number_in_lineage==k & df$cycle_stage==m)] - 
                           df$Time[which(df[K]==c & df$lineage==i & df$number_in_lineage==k & df$cycle_stage==b)])
              s2_dt <- c(s2_dt, df$Time[which(df[K]==c & df$lineage==i & df$number_in_lineage==l & df$cycle_stage==m)] - 
                           df$Time[which(df[K]==c & df$lineage==i & df$number_in_lineage==l & df$cycle_stage==b)])
              s1_number_in_lineage_dt <- c(s1_number_in_lineage_dt,k)
              s2_number_in_lineage_dt <- c(s2_number_in_lineage_dt,l)
              sister_lineage_dt <- c(sister_lineage_dt, i)
            } else {
              c1_dt <- c(c1_dt, df$Time[which(df[K]==c & df$lineage==i & df$number_in_lineage==k & df$cycle_stage==m)] - 
                           df$Time[which(df[K]==c & df$lineage==i & df$number_in_lineage==k & df$cycle_stage==b)])
              c2_dt <- c(c2_dt, df$Time[which(df[K]==c & df$lineage==i & df$number_in_lineage==l & df$cycle_stage==m)] - 
                           df$Time[which(df[K]==c & df$lineage==i & df$number_in_lineage==l & df$cycle_stage==b)])
              c1_number_in_lineage_dt <- c(c1_number_in_lineage_dt,k)
              c2_number_in_lineage_dt <- c(c2_number_in_lineage_dt,l)
              cousin_generation_dt <- c(cousin_generation_dt, j)
              cousin_lineage_dt <- c(cousin_lineage_dt, i)
              #print("cousin")
              #print(l)
            }
            if (is.na(df$vol[which(df[K]==c & df$lineage==i & df$number_in_lineage==k & df$cycle_stage==m)])==F
                && is.na(df$vol[which(df[K]==c & df$lineage==i & df$number_in_lineage==l & df$cycle_stage==m)])==F){ ## both divided and have vf: vf can be calculated
              if (k%%2==F && l==k+1){ ## for sisters :
                s1_vf <- c(s1_vf, df$vol[which(df[K]==c & df$lineage==i & df$number_in_lineage==k & df$cycle_stage==m)])
                s2_vf <- c(s2_vf, df$vol[which(df[K]==c & df$lineage==i & df$number_in_lineage==l & df$cycle_stage==m)])
                s1_number_in_lineage_vf <- c(s1_number_in_lineage_vf,k)
                s2_number_in_lineage_vf <- c(s2_number_in_lineage_vf,l)
                sister_lineage_vf <- c(sister_lineage_vf, i)
                #print("sister")
                #print(l)
                if (b %in% df$cycle_stage[which(df[K]==c & df$lineage==i & df$number_in_lineage==k)]
                    && b %in% df$cycle_stage[which(df[K]==c & df$lineage==i & df$number_in_lineage==l)]
                    && is.na(df$vol[which(df[K]==c & df$lineage==i & df$number_in_lineage==k & df$cycle_stage==b)])==F
                    && is.na(df$vol[which(df[K]==c & df$lineage==i & df$number_in_lineage==l & df$cycle_stage==b)])==F){
                  s1dv = (df$vol[which(df[K]==c & df$lineage==i & df$number_in_lineage==k & df$cycle_stage==m)]-
                            df$vol[which(df[K]==c & df$lineage==i & df$number_in_lineage==k & df$cycle_stage==b)])
                  s2dv = (df$vol[which(df[K]==c & df$lineage==i & df$number_in_lineage==l & df$cycle_stage==m)]-
                            df$vol[which(df[K]==c & df$lineage==i & df$number_in_lineage==l & df$cycle_stage==b)])
                  s1dt = (df$Time[which(df[K]==c & df$lineage==i & df$number_in_lineage==k & df$cycle_stage==m)]-
                            df$Time[which(df[K]==c & df$lineage==i & df$number_in_lineage==k & df$cycle_stage==b)])
                  s2dt = (df$Time[which(df[K]==c & df$lineage==i & df$number_in_lineage==l & df$cycle_stage==m)]-
                            df$Time[which(df[K]==c & df$lineage==i & df$number_in_lineage==l & df$cycle_stage==b)])
                  s1_dvdt <- c(s1_dvdt,s1dv/s1dt)
                  s2_dvdt <- c(s2_dvdt,s2dv/s2dt)
                  s1_number_in_lineage_dvdt <- c(s1_number_in_lineage_dvdt,k)
                  s2_number_in_lineage_dvdt <- c(s2_number_in_lineage_dvdt,l)
                  sister_lineage_dvdt <- c(sister_lineage_dvdt,i)
#                   rm(s1dv,s1dt,s2dt,s2dv)
                } else {}
              } else { ## for cousins : 
                c1_vf <- c(c1_vf, df$vol[which(df[K]==c & df$lineage==i & df$number_in_lineage==k & df$cycle_stage==m)])
                c2_vf <- c(c2_vf, df$vol[which(df[K]==c & df$lineage==i & df$number_in_lineage==l & df$cycle_stage==m)])
                c1_number_in_lineage_vf <- c(c1_number_in_lineage_vf,k)
                c2_number_in_lineage_vf <- c(c2_number_in_lineage_vf,l)
                cousin_generation_vf <- c(cousin_generation_vf, j)
                cousin_lineage_vf <- c(cousin_lineage_vf, i)
                #print("cousin")
                #print(l)
                if (is.na(df$vol[which(df[K]==c & df$lineage==i & df$number_in_lineage==k & df$cycle_stage==b)])==F
                    && is.na(df$vol[which(df[K]==c & df$lineage==i & df$number_in_lineage==l & df$cycle_stage==b)])==F){
                  c1dv = (df$vol[which(df[K]==c & df$lineage==i & df$number_in_lineage==k & df$cycle_stage==m)]-
                            df$vol[which(df[K]==c & df$lineage==i & df$number_in_lineage==k & df$cycle_stage==b)])
                  c2dv = (df$vol[which(df[K]==c & df$lineage==i & df$number_in_lineage==l & df$cycle_stage==m)]-
                            df$vol[which(df[K]==c & df$lineage==i & df$number_in_lineage==l & df$cycle_stage==b)])
                  c1dt = (df$Time[which(df[K]==c & df$lineage==i & df$number_in_lineage==k & df$cycle_stage==m)]-
                            df$Time[which(df[K]==c & df$lineage==i & df$number_in_lineage==k & df$cycle_stage==b)])
                  c2dt = (df$Time[which(df[K]==c & df$lineage==i & df$number_in_lineage==l & df$cycle_stage==m)]-
                            df$Time[which(df[K]==c & df$lineage==i & df$number_in_lineage==l & df$cycle_stage==b)])
                  c1_dvdt <- c(c1_dvdt,s1dv/s1dt)
                  c2_dvdt <- c(c2_dvdt,s2dv/s2dt)
                  rm(c1dv,c1dt,c2dt,c2dv)
                  c1_number_in_lineage_dvdt <- c(c1_number_in_lineage_dvdt,k)
                  c2_number_in_lineage_dvdt <- c(c2_number_in_lineage_dvdt,l)
                  cousin_lineage_dvdt <- c(cousin_lineage_dvdt,i)
                } else {}
              }
            }
          } else {}
          l=l+1
        }
      }
    }
  }
}





## dvdt across generations ####
# mother/daughter ####
x=dvdt_md_m
y=dvdt_md_d
test=lm(y~x)
testb=rlm(y~x)
pa <- ggplot(data.frame(x,y),aes(x=x,y=y))+#, color=factor(lineage_vf))) +
  geom_point() +
  labs(title= paste(title, "\n average_GR [µm3/hrs]"), y="daughter (g1)", 
       x=paste("mother (g0)",
               "\nlm:","R2=",signif(summary(test)$r.squared,a),", p=",signif(anova(test)$'Pr(>F)',a),
               ", n=",length(x),", rlm: a=",signif(coef(testb)[2],a), ", b=", signif(coef(testb)[1],a),
               "\n","a=",signif(coef(test)[2],a), ", b=", signif(coef(test)[1],a)))+ 
  #   scale_colour_discrete (name="lineage") +
  geom_abline(intercept = coef(testb)[1], slope = coef(testb)[2])+
  theme_bw()
sink(file=paste(dir,title,"_lm_corr_lineages_gr_g0_g1.txt",sep=""))
sink()
rm(test, testb)

# mother/grand-daughter ####
x=dvdt_mgd_m
y=dvdt_mgd_gd
test=lm(y~x)
testb=rlm(y~x)
pb <- ggplot(data.frame(x,y),aes(x=x,y=y))+#, color=factor(lineage_vf))) +
  geom_point() +
  labs(title= paste(title, "\n average_GR [µm3/hrs]"), y="grand-daughter (g2)", 
       x=paste("mother (g0)",
               "\nlm:","R2=",signif(summary(test)$r.squared,a),", p=",signif(anova(test)$'Pr(>F)',a),
               ", n=",length(x),", rlm: a=",signif(coef(testb)[2],a), ", b=", signif(coef(testb)[1],a),
               "\n","a=",signif(coef(test)[2],a), ", b=", signif(coef(test)[1],a)))+ 
  #   scale_colour_discrete (name="lineage") +
  geom_abline(intercept = coef(testb)[1], slope = coef(testb)[2])+
  theme_bw()
sink(file=paste(dir,title,"_lm_corr_lineages_gr_g0_g2.txt",sep=""))
sink()
rm(test, testb)

# mother/grand-grand-daughter ####
x=dvdt_mggd_m
y=dvdt_mggd_ggd
test=lm(y~x)
testb=rlm(y~x)
pc <- ggplot(data.frame(x,y),aes(x=x,y=y))+#, color=factor(lineage_vf))) +
  geom_point() +
  labs(title= paste(title, "\n average_GR [µm3/hrs]"), y="grand-grand-daughter (g3)", 
       x=paste("mother (g0)",
               "\nlm:","R2=",signif(summary(test)$r.squared,a),", p=",signif(anova(test)$'Pr(>F)',a),
               ", n=",length(x),", rlm: a=",signif(coef(testb)[2],a), ", b=", signif(coef(testb)[1],a),
               "\n","a=",signif(coef(test)[2],a), ", b=", signif(coef(test)[1],a)))+ 
  #   scale_colour_discrete (name="lineage") +
  geom_abline(intercept = coef(testb)[1], slope = coef(testb)[2])+
  theme_bw()
sink(file=paste(dir,title,"_lm_corr_lineages_gr_g0_g3.txt",sep=""))
sink()
rm(test, testb)

# daughter/daughter (vectors come from 1) just above) ####
x=dvdt_1
y=dvdt_2
test=lm(y~x)
testb=rlm(y~x)
pc <- ggplot(data.frame(x,y),aes(x=x,y=y))+#, color=factor(lineage_vf))) +
  geom_point() +
  labs(title= paste(title, "\n average_GR [µm3/hrs]"), y="grand-grand-daughter (g3)", 
       x=paste("mother (g0)",
               "\nlm:","R2=",signif(summary(test)$r.squared,a),", p=",signif(anova(test)$'Pr(>F)',a),
               ", n=",length(x),", rlm: a=",signif(coef(testb)[2],a), ", b=", signif(coef(testb)[1],a),
               "\n","a=",signif(coef(test)[2],a), ", b=", signif(coef(test)[1],a)))+ 
  #   scale_colour_discrete (name="lineage") +
  geom_abline(intercept = coef(testb)[1], slope = coef(testb)[2])+
  theme_bw()
sink(file=paste(dir,title,"_lm_corr_lineages_gr_g0_g3.txt",sep=""))
sink()
rm(test, testb)

# make figure for gr with also daughter/daughter
tiff(paste(dir,title,"_ratios_sisters_dv_dt_vi_vf.tiff",sep=""), width=1000,height=200)
grid.arrange(pa,pb,pc,pd,ncol=1)
dev.off()
rm(pa,pb,pc,pd)
