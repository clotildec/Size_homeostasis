#############   LINEAGE ANALYSIS VERSION 2 ##################################

#cette version permet d'avoir les données rangées dans un ordre aléatoire car elle n'utilise
# QUE les numéros d'identité des cellules

# export des données
data <- read.csv("C:/Users/Clotilde/Desktop/150612_hela-fucci_rosco_analysis_150630_clo_programme/150612_hela-fucci-rosco-vol_analysis_donnees_seuls.csv",sep=";",header=TRUE)
# data <- data[-which(data$vol>5000),]
d <- data
# data <- data[-which(data$lineage==3 & data$number_in_lineage==2),] # because ccl ==35hrs

### prepare a clean dataframe ####
#parameters to enter
dir <- "C:/Users/Clotilde/Desktop/150612_hela-fucci_rosco_analysis_150630_clo_programme/"
title <- "150612_hela_fucci_rosco_analysis_150630_ds0"
time_resolution <- 10 
data$cond_cell_type <- "HeLa-fucci"
data$cond_date <- "150612"
data$cond_vol <- "chamber_5.1_h2"
# data$cond_Glu <- "0"

# add generation column ####
generation=NULL
for (i in c(1:length(data$lineage))){
  if (data$number_in_lineage[i]==1){
    g=0
    } else if (data$number_in_lineage[i] %in% c(2,3)){
    g=1
    } else if (data$number_in_lineage[i] %in% c(4,5,6,7)){
      g=2
    } else if (data$number_in_lineage[i] %in% c(8,9,10,11,12,13,14,15)){
      g=3
    } else if (data$number_in_lineage[i] %in% c(16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31)){
      g=4
    } else if (data$number_in_lineage[i] %in% c(32,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63)) {
      g=5
    }
    generation=c(generation,g)
  }
data$generation <-generation
rm(generation)
# add volume different from NA column using volum_NA ####
vol=NULL
for (i in c(1:length(data$lineage))){
  if (data$volume_NA[i]==1){
    v=data$Vol[i]
    vol=c(vol,v)
    rm(v)
  }else{
    v=NA
    vol=c(vol,v)
    rm(v)
  }
}
data$vol <- vol 
rm(vol)

write.csv(data,paste(dir,title,"_df.csv",sep=""))

### call packages####
library(ggplot2)
library(grid)
library(gridExtra)
library(MASS)
library(xtable)

###################### ANALYSIS ############################################
# choix dataset to analyse
df <- data

########### controls for ccl and volume ########## 
#### 1) CCL #####
##### choose parameters for analysis ####
m = 1.6 # timepoint chosen for end of cell cycle
b = 240 # timepoint chosen for bgining of cell cycle
e = 3 # timepoitn where cell is lost alive (e for exit)
end = max(df$Time) # final frame of the movie
a=2
tri=df$lineage_complete
df$name_tri='lineages'

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

for (i in unique(data$lineage)){
  print("i")
  print(i)
  for (j in unique(c(data$number_in_lineage[which(data$lineage==i & data$number_in_lineage>1)]))){
    print("j")
    print(j)
    if (m %in% data$cycle_stage[which(data$lineage==i & data$number_in_lineage==j)]
        && b %in% data$cycle_stage[which(data$lineage==i & data$number_in_lineage==j)]){
      t1 = data$Time[which(data$lineage==i & data$number_in_lineage==j & data$cycle_stage==m)]-
        data$Time[which(data$lineage==i & data$number_in_lineage==j & data$cycle_stage==b)]
      ccl_complete = c(ccl_complete,t1)
      cell_number_complete = c(cell_number_complete,j)
      generation_complete = c(generation_complete, data$generation[which(data$lineage==i & data$number_in_lineage==j & data$cycle_stage==b)])
      lineage_complete = c(lineage_complete, i)
      ti_complete=c(ti_complete,data$Time[which(data$lineage==i & data$number_in_lineage==j & data$cycle_stage==b)])
      tf_complete = c(tf_complete,data$Time[which(data$lineage==i & data$number_in_lineage==j & data$cycle_stage==m)])
      rm (t1) 
      print(length(cell_number_complete))
      print(length(ccl_complete))
    } else if (b %in% data$cycle_stage[which(data$lineage==i & data$number_in_lineage==j)]){
      if (e %in% data$cycle_stage[which(data$lineage==i & data$number_in_lineage==j)]) {
        t2 = data$Time[which(data$lineage==i & data$number_in_lineage==j & data$cycle_stage==e)] -
          data$Time[which(data$lineage==i & data$number_in_lineage==j & data$cycle_stage==b)]
        ccl_uncomplete = c(ccl_uncomplete,t2)
        cell_number_uncomplete = c(cell_number_uncomplete,j)
        generation_uncomplete = c(generation_uncomplete, data$generation[which(data$lineage==i & data$number_in_lineage==j & data$cycle_stage==b)])
        lineage_uncomplete = c(lineage_uncomplete,i)
        rm (t2)
      } else {
        t2 = end -
          data$Time[which(data$lineage==i & data$number_in_lineage==j & data$cycle_stage==b)]
        ccl_uncomplete = c(ccl_uncomplete,t2)
        cell_number_uncomplete = c(cell_number_uncomplete,j)
        generation_uncomplete = c(generation_uncomplete, data$generation[which(data$lineage==i & data$number_in_lineage==j & data$cycle_stage==b)])
        lineage_uncomplete = c(lineage_uncomplete,i)
        rm (t2) 
      }
      print(length(cell_number_uncomplete))
      print(length(ccl_uncomplete))
    }
  }
}

##### plots ####
ccl <- c(ccl_complete/60,ccl_uncomplete/60)
type <- c(rep('complete',length(ccl_complete)),rep('uncomplete',length(ccl_uncomplete)))
cell_number = c(cell_number_complete,cell_number_uncomplete)
lineage <- c(lineage_complete, lineage_uncomplete)

###### hist  CCL complètes/non complètes ####
df <- data.frame(ccl,type,cell_number, lineage)
hist1 <- ggplot(df,aes(x=ccl,fill=as.factor(type)))+
  geom_histogram(position="dodge")+
  theme_bw()+
  labs(title=title,
       fill="ccl",
       x= paste('DT [hours]',"\n","complete: mean=",signif(mean(ccl_complete)/60,a), ", n=",length(ccl_complete), "\n uncomplete: n=", length(ccl_uncomplete),sep=""))
rm(df)

df <- data.frame(ccl_complete/60,generation_complete,ti_complete/60,tf_complete/60,lineage_complete,cell_number_complete)
hist2 <- ggplot(df,aes(ccl_complete/60, fill=as.factor(generation_complete)))+
  geom_histogram(position="dodge")+
  theme_bw()+
  labs(title=title,
       fill="generation",
       x= paste('DT [hours]',"\n","complete: mean=",mean(ccl_complete)/60, ", n=",length(ccl_complete), "\n uncomplete: n=", length(ccl_uncomplete),sep=""))
write.csv(df,paste(dir,title,"_ctrl_ccl_df.csv",sep=""))
rm(df)

###### plot  CCL time ####
df <- data.frame(ccl_complete/60,generation_complete,ti_complete/60,tf_complete/60,lineage_complete)
x=df$ti_complete.60
y=df$ccl_complete.60
test=lm(y~x)
testb=rlm(y~x)
pa <- ggplot(df,aes(x=ti_complete.60,y=ccl_complete.60))+#,color=factor(tri))) +
  geom_point() +
  labs(title= title, y="DT [hrs]", 
       x=paste("Time_i [hrs]","\n lm:","R2=",signif(summary(test)$r.squared,a),", p=",signif(anova(test)$'Pr(>F)',a),
               "a=",signif(coef(test)[2],a), ", b=", signif(coef(test)[1],2),
               "\n rlm: a=",signif(coef(testb)[2],a), ", b=", signif(coef(testb)[1],a)))+
  scale_colour_discrete (name="lineage") +
  geom_abline(intercept = coef(testb)[1], slope = coef(testb)[2])+
#   stat_smooth(method='lm',se=FALSE)+
  theme_bw()
sink(file=paste(dir,title,"ctrl_dt_ti_lm_rlm.txt",sep=""))
summary(test)
summary(testb)
sink()
rm(test,testb,y,x)

x=df$tf_complete.60
y=df$ccl_complete.60
test=lm(y~x)
testb=rlm(y~x)
pb <- ggplot(df,aes(x=tf_complete.60,y=ccl_complete.60))+#,color=factor(tri))) +
  geom_point() +
  labs(title= title, y="DT [hrs]", 
       x=paste("Time_f [hrs]","\n lm:","R2=",signif(summary(test)$r.squared,a),", p=",signif(anova(test)$'Pr(>F)',a),
               "a=",signif(coef(test)[2],a), ", b=", signif(coef(test)[1],2),
               "\n rlm: a=",signif(coef(testb)[2],a), ", b=", signif(coef(testb)[1],a)))+
  scale_colour_discrete (name="lineage") +
  #   geom_abline(intercept = coef(testb)[1], slope = coef(testb)[2])+
  #   stat_smooth(method='lm',se=FALSE)+
  theme_bw()
sink(file=paste(dir,title,"ctrl_dt_tf_lm_rlm.txt",sep=""))
summary(test)
summary(testb)
sink()
rm(test,testb,y,x)

##### make figure ####
tiff(paste(dir,title,"_hist_ccl_ctrl.tiff",sep=""), width=1000,height=500)
grid.arrange(hist1,pa,hist2,pb,ncol=2)
dev.off()
rm(df,tf_complete,ti_complete,ccl_complete,ccl_uncomplete,generation_complete,generation_uncomplete,lineage_complete,lineage_uncomplete,cell_number,ccl,cell_number_complete,cell_number_uncomplete)
rm(hist1,hist2,pa,pb)

#### 2) volume ####
# 1) all Vf and generation ####
m = 1.6 # timepoint chosen for end of cell cycle
vf = NULL
cell_number = NULL
generation = NULL
for (i in unique(data$lineage)){
  for (j in unique(c(data$number_in_lineage[which(data$lineage==i)]))){
    if (m %in% data$cycle_stage[which(data$lineage==i & data$number_in_lineage==j)]){
#       if (data$volume_NA[which(data$lineage==i & data$number_in_lineage==j & data$cycle_stage==m)]==1){
        v = data$vol[which(data$lineage==i & data$number_in_lineage==j & data$cycle_stage==m)]
        vf = c(vf,v)
        cell_number = c(cell_number,j)
        generation = c(generation, data$generation[which(data$lineage==i & data$number_in_lineage==j & data$cycle_stage==m)])
        rm (v)         
#       }else {}
    } else {}
  }
}
df <- data.frame (vf, cell_number, generation)
rm(vf, cell_number,generation)
hist3 <- ggplot(df,aes(x=vf,fill=as.factor(generation)))+
  geom_histogram(position="dodge",binwidth=50)+
  theme_bw()+
  labs(title=title,
       fill="generation",
       x= paste('Vf [µm3]','\n','mean=',signif(mean(df$vf),a), ", n=",length(df$vf),sep=""))
rm(df)

# 2) Vf, DV and time ####
###### tri ####
m = 1.6 # timepoint chosen for end of cell cycle
b = 240 # timepoint chosen for begining of cell cycle
vf = NULL
vi = NULL
tf = NULL
ti = NULL
number_in_lineage = NULL
lineage = NULL
generation = NULL
for (i in unique(data$lineage)){
  for (j in unique(c(data$number_in_lineage[which(data$lineage==i)]))){
      if(m %in% data$cycle_stage[which(data$lineage==i & data$number_in_lineage==j)]){
        vf_temp = data$vol[which(data$lineage==i & data$number_in_lineage==j & data$cycle_stage==m)]
        tf_temp = data$Time[which(data$lineage==i & data$number_in_lineage==j & data$cycle_stage==m)]
        g = data$generation[which(data$lineage==i & data$number_in_lineage==j & data$cycle_stage==m)]
          if(b %in% data$cycle_stage[which(data$lineage==i & data$number_in_lineage==j)]){
            vi_temp = data$vol[which(data$lineage==i & data$number_in_lineage==j & data$cycle_stage==b)]
            ti_temp = data$Time[which(data$lineage==i & data$number_in_lineage==j & data$cycle_stage==b)]
          } else {
            vi_temp = NA
            ti_temp = NA
          }
        vf = c(vf,vf_temp)
        vi = c(vi,vi_temp)
        tf = c(tf,tf_temp)
        ti = c(ti,ti_temp)
        number_in_lineage = c(number_in_lineage,j)
        lineage = c(lineage, i)
        generation = c(generation, g)
        rm(vf_temp, vi_temp, tf_temp, ti_temp)
      } else {}
  }
}
df <- data.frame (vf, vi, tf, ti, number_in_lineage,lineage, generation)
dv=NULL
for (i in c(1:length(df$vf))){
  if (is.na(df$vf[i])==F & is.na(df$vi[i])==F) {
    dv_temp = df$vf[i]- df$vi[i]
  }else{
    dv_temp = NA
  }
  dv = c(dv,dv_temp)
}

###### plots ####
df$dv = dv
df$dt = df$tf-df$ti
df$dvdt <- df$dv/df$dt
rm(vf, vi, tf, ti, number_in_lineage,lineage, generation,dv)

hist4 <- ggplot(df,aes(x=vi,fill=as.factor(generation)))+
  geom_histogram(position="dodge",binwidth=30)+
  theme_bw()+
  labs(title=title,
       fill="generation",
       x= paste('Vi [µm3]','\n', 'mean=',signif(mean(df$vi,na.rm=TRUE),a),", n=",length(df$vi[which(is.na(df$vi)==F)]),sep=""))

hist5 <- ggplot(df,aes(x=dv,fill=as.factor(generation)))+
  geom_histogram(position="dodge",binwidth=30)+
  theme_bw()+
  labs(title=title,
       fill="generation",
       x= paste('DV [µm3]','\n', 'mean=',signif(mean(df$dv,na.rm=TRUE),a),', n=',length(df$dv[which(is.na(df$dv)==F)]),sep=""))

test=lm(df$vf~df$ti)
testb=rlm(df$vf~df$ti)
p1 <- ggplot(df,aes(x=ti,y=vf))+#, color=factor(lineage))) +
  geom_point() +
  labs(title= title, y="Vf [µm3]", x=paste("Time_i [hours]","\n","R2=",signif(summary(test)$r.squared,a),", p=",signif(anova(test)$'Pr(>F)',a),
                                           ", n=",length(df$vf[which(is.na(df$vf)==F & is.na(df$ti)==F)]))) +
  scale_colour_discrete (name="lineage") +
#   geom_abline(intercept = coef(testb)[1], slope = coef(testb)[2])+
  theme_bw()
sink(file=paste(dir,title,"ctrl_lm_vf_ti.txt",sep=""))
summary(test)
summary(testb)
sink() 
rm(test)

test = lm(df$dv~df$ti)
testb = rlm(df$dv~df$ti)
p2 <- ggplot(df,aes(x=ti,y=dv))+#, color=factor(lineage))) +
  geom_point() +
  labs(title= title, y="DV [µm3]", x=paste("Time_i [hours]","\n","R2=",signif(summary(test)$r.squared,a),
                                           ", p=",signif(anova(test)$'Pr(>F)',a),
                                           ", n=",length(df$dv[which(is.na(df$dv)==F & is.na(df$ti)==F)]))) +
  scale_colour_discrete (name="lineage") +
#   geom_abline(intercept = coef(testb)[1], slope = coef(testb)[2])+
  theme_bw()
sink(file=paste(dir,title,"ctrl_lm_dv_ti.txt",sep=""))
summary(test)
summary(testb)
sink() 
rm(test)

test=lm(df$dvdt~df$ti)
p3 <- ggplot(df,aes(x=ti,y=dvdt))+#, color=factor(lineage))) +
  geom_point() +
  labs(title= title, y="DV/DT [µm3/hrs]", x=paste("Time_i [hours]","\n","R2=",signif(summary(test)$r.squared,a)
                                                  ,", p=",signif(anova(test)$'Pr(>F)',a)
                                                  ,", n=",length(df$dvdt[which(is.na(df$dv)==F)]))) +
  scale_colour_discrete (name="lineage") +
  theme_bw()
sink(file=paste(dir,title,"ctrl_lm_dvdt_ti.txt",sep=""))
summary(test)
sink() 
rm(test)

###### save ####
tiff(paste(dir,title,"_hist_plot_vol_ctrl.tiff",sep=""), width=1000,height=750)
grid.arrange(hist3,p1,hist4,p2,hist5,p3, ncol=2)
dev.off()
# 
# tiff(paste(dir,title,"_hist_vol_ctrl.tiff",sep=""), width=500,height=750)
# grid.arrange(hist3,hist4,hist5)
# dev.off()
# 
# tiff(paste(dir,title,"_plot_vol_ctrl.tiff",sep=""), width=500,height=750)
# grid.arrange(p1,p2,p3, ncol=1)
# dev.off()

write.csv(df,paste(dir,title,"_ctrl_vol_df.csv",sep=""))
rm(df, p1,p2,p3,hist3,hist4,hist5)

# 3) check what 's the best timepoint TO BE FINISHED #############
m = c(1.6,1.4,1) # timepoint chosen for end of cell cycle
b = c(200,240) # timepoint chosen for begining of cell cycle
t = c(b,m)
vf=NULL

for (i in c(1:(max(d$lineage)))){
  for (j in unique(c(data$number_in_lineage[which(data$lineage==i)]))){
    v=c(rep(1,length(t)))
    for (k in c(1:length(t))) {
      if (t[k] %in% data$cycle_stage[which(d$lineage==i & d$number_in_lineage==j)]){
        v[k] = data$vol[which(d$lineage==i & d$number_in_lineage==j & data$cycle_stage==t[k])]
      } else {
        v[k]=NA
      }
    }
    vf=rbind(vf,v)
  }
}
df <- data.frame(vf[,1])
df$v240 <- vf[,2]
df$v1.6 <- vf[,3]
df$v1.4 <- vf[,4]
df$v1 <- vf[,5]
df <- t(df)

########## absolute values dv, vf, vi, dt, dv/dt ########## 
#### tri ####
m = 1.6 # timepoint chosen for end of cell cycle
b = 240 # timepoint chosen for begining of cell cycle
vf = NULL
vi = NULL
tf = NULL
ti = NULL
number_in_lineage = NULL
lineage = NULL
generation = NULL
for (i in unique(data$lineage)){
  for (j in unique(c(data$number_in_lineage[which(data$lineage==i)]))){
    if(m %in% data$cycle_stage[which(data$lineage==i & data$number_in_lineage==j)]){
      vf_temp = data$vol[which(data$lineage==i & data$number_in_lineage==j & data$cycle_stage==m)]
      tf_temp = data$Time[which(data$lineage==i & data$number_in_lineage==j & data$cycle_stage==m)]
      g = data$generation[which(data$lineage==i & data$number_in_lineage==j & data$cycle_stage==m)]
      if(b %in% data$cycle_stage[which(data$lineage==i & data$number_in_lineage==j)]){
        vi_temp = data$vol[which(data$lineage==i & data$number_in_lineage==j & data$cycle_stage==b)]
        ti_temp = data$Time[which(data$lineage==i & data$number_in_lineage==j & data$cycle_stage==b)]
      } else {
        vi_temp = NA
        ti_temp = NA
      }
      vf = c(vf,vf_temp)
      vi = c(vi,vi_temp)
      tf = c(tf,tf_temp)
      ti = c(ti,ti_temp)
      number_in_lineage = c(number_in_lineage,j)
      lineage = c(lineage, i)
      generation = c(generation, g)
      rm(vf_temp, vi_temp, tf_temp, ti_temp)
    } else {}
  }
}
df <- data.frame (vf, vi, tf, ti, number_in_lineage,lineage, generation)
dv=NULL
for (i in c(1:length(df$vf))){
  if (is.na(df$vf[i])==F & is.na(df$vi[i])==F) {
    dv_temp = df$vf[i]- df$vi[i]
  }else{
    dv_temp = NA
  }
  dv = c(dv,dv_temp)
}

#### make df #### CAREFULL THINGS TO ENTER HERE ! ####
df$dv = dv
df$dt = (df$tf-df$ti)/60
df$dvdt = df$dv/df$dt
df$vfvi = df$vf/df$vi*100
df$dvvi = df$dv/df$vi*100
df$cond_cell_type <- "HT29"
df$cond_date <- "150605_sylvain"
df$cond_vol <- "chamber_5.2"
# df$cond_Glu <- "0"
rm(vf, vi, tf, ti, number_in_lineage,lineage, generation,dv,dv_temp)
write.csv(df, paste(dir,title,"_df_abs_val.csv",sep=""))

#### add parameter for sorting ####
##### for time_i ####
# enter parameters here + change line 5 in the loop , to add sorting-variable in df #####
n = 4
sort = df$ti[which(is.na(df$ti)==F)]
sort2 = df$ti
bin = (max(sort)-min(sort))/n
df$sort_ti = NA
bin_sort_ti=NULL
name_tri="time_i"

for (i in c(1:nrow(df))){
  for (j in c(1:n)){
    if (is.na(sort2[i])==F
        && (sort2[i] <= min(sort)+bin*j)==TRUE && (sort2[i] > min(sort)-1+bin*(j-1))==TRUE) {
        df$sort_ti[i] = j ## change this line !
    }else{}
  }
} 

# save a vector with bining values, make and save a table with informations about the sorting variable ####
bin_sort_ti=NULL
for (j in c(0:n)){
  bin_sort_ti=c(bin_sort_ti,(min(sort)+bin*j))
}
res=matrix(ncol=3,nrow=n)
for (i in c(1:n)){
  res[i,1]=min(sort)-1+bin*(i-1)
  res[i,2]=min(sort)+bin*i
  res[i,3]=length(df$sort_ti[which(df$sort_ti==i)])
}
res=as.data.frame(res)
colnames(res)=c('min','max','n')
table=tableGrob(res,equal.width = F)
grid.newpage()
h <- grobHeight(table)
w <- grobWidth(table)
title_temp <- textGrob(paste(title,"\nsorting by ",name_tri), y=unit(0.5,"npc") + 0.6*h, 
                       vjust=0, gp=gpar(fontsize=15))
gt <- gTree(children=gList(table, title_temp))
grid.draw(gt)
tiff(paste(dir,title,"_abs_val_table_sorting_",name_tri,".tiff",sep=""), width=400,height=300)
grid.draw(gt)
dev.off()

# write.csv(df, paste(dir,title,"_df_abs_values.csv",sep=""))

##### for average_gr ####
# enter parameters here + change line 5 in the loop , to add sorting-variable in df #####
n = 4
sort = df$dvdt[which(is.na(df$dvdt)==F)]
sort2 = df$dvdt
bin = (max(sort)-min(sort))/n
df$sort_dvdt = NA
bin_sort_dvdt=NULL
name_tri="average_GR"

for (i in c(1:nrow(df))){
  for (j in c(1:n)){
    if (is.na(sort2[i])==F
        && (sort2[i] <= min(sort)+bin*j)==TRUE && (sort2[i] > min(sort)-1+bin*(j-1))==TRUE) {
      df$sort_dvdt[i] = j ## change this line !
    }else{}
  }
} 

# save a vector with bining values, make and save a table with informations about the sorting variable ####
bin_sort_dvdt=NULL
for (j in c(0:n)){
  bin_sort_dvdt=c(bin_sort_ti,(min(sort)+bin*j))
}
res=matrix(ncol=3,nrow=n)
for (i in c(1:n)){
  res[i,1]=signif(min(sort)-1+bin*(i-1),a)
  res[i,2]=signif(min(sort)+bin*i,a)
  res[i,3]=length(df$sort_dvdt[which(df$sort_dvdt==i)])
}
res=as.data.frame(res)
colnames(res)=c('min','max','n')
table=tableGrob(res,equal.width = F)
grid.newpage()
h <- grobHeight(table)
w <- grobWidth(table)
title_temp <- textGrob(paste(title,"\nsorting by ",name_tri), y=unit(0.5,"npc") + 0.6*h, 
                       vjust=0, gp=gpar(fontsize=15))
gt <- gTree(children=gList(table, title_temp))
grid.draw(gt)
tiff(paste(dir,title,"_abs_val_table_sorting_",name_tri,".tiff",sep=""), width=400,height=300)
grid.draw(gt)
dev.off()

# write.csv(df, paste(dir,title,"_df_abs_values.csv",sep=""))

#### choose parameters for analysis ####
tri = df$sort_ti
name_tri = 'time_i'
a=2 # number of significant digits

#### plots ####
## p1 dv/dt vs vi ####
x=df$vi[which(is.na(df$vi)==FALSE & is.na(df$dvdt)==FALSE)]
y=df$dvdt[which(is.na(df$vi)==FALSE & is.na(df$dvdt)==FALSE)]
test = lm(y~x)
testb = rlm(y~x)
p1 <- ggplot(df,aes(x=vi,y=dvdt, color=factor(tri))) +
  geom_point() +
  labs(title= title, y="DV/DT [µm3/hrs]", 
       x=paste("Vol_i [µm3]","\n lm:","R2=",signif(summary(test)$r.squared,a),", p=",signif(anova(test)$'Pr(>F)',a),
               ", n=",length(df$dvdt[which(is.na(df$vi)==FALSE & is.na(df$dvdt)==FALSE)]), 
               "a=",signif(coef(test)[2],a), ", b=", signif(coef(test)[1],2),
               "\n rlm: a=",signif(coef(testb)[2],a), ", b=", signif(coef(testb)[1],a),
               "\n Pearson: r=", signif(cor(x,y),a), ", p=",signif(cor.test(x,y)$p.value,a)))+
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
x=df$vi[which(is.na(df$vi)==FALSE & is.na(df$vfvi)==FALSE)]
y=df$vfvi[which(is.na(df$vi)==FALSE & is.na(df$vfvi)==FALSE)]
test = lm(y~x)
testb = rlm(y~x)
p2 <- ggplot(df,aes(x=vi,y=vfvi, color=factor(tri))) +
  geom_point() +
  labs(title= title, y="vf/vi %", 
       x=paste("Vol_i [µm3]","\n lm: ","R2=",signif(summary(test)$r.squared,2),", p=",signif(anova(test)$'Pr(>F)',2),
               ", n=",length(df$vfvi[which(is.na(df$vi)==FALSE & is.na(df$vfvi)==FALSE)]),
               "a=",signif(coef(test)[2],a), ", b=", signif(coef(test)[1],a),
               "\n rlm: a=",signif(coef(testb)[2],a), ", b=", signif(coef(testb)[1],a),
               "\n Pearson: r=", signif(cor(x,y),a), "p=",signif(cor.test(x,y)$p.value,a)))+
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
x=df$vi[which(is.na(df$vi)==FALSE & is.na(df$dv)==FALSE)]
y=df$dv[which(is.na(df$vi)==FALSE & is.na(df$dv)==FALSE)]
test=lm(y~x)
testb=rlm(y~x)
p3 <- ggplot(df,aes(x=vi,y=dv, color=factor(tri))) +
  geom_point() +
  labs(title= title, y="DV [µm3]", 
       x=paste("Vol_i [µm3]","\n","lm: R2=",signif(summary(test)$r.squared,2),", p=",signif(anova(test)$'Pr(>F)',2),
               ", n=",length(df$dv[which(is.na(df$vi)==FALSE & is.na(df$dv)==FALSE)]),
               "a=",signif(coef(test)[2],2), ", b=", signif(coef(test)[1],2),
               "\n rlm: a=",signif(coef(testb)[2],2), ", b=", signif(coef(testb)[1],2),
               "\n Pearson: r=", signif(cor(x,y),a), "p=",signif(cor.test(x,y)$p.value,a)))+
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
x=df$vi[which(is.na(df$vi)==FALSE & is.na(df$dt)==FALSE)]
y=df$dt[which(is.na(df$vi)==FALSE & is.na(df$dt)==FALSE)]
test=lm(df$dt~df$vi)
testb=rlm(y~x)
p4 <- ggplot(df,aes(x=vi,y=dt, color=factor(tri))) +
  geom_point() +
  labs(title= title, y="DT [hrs]", 
       x=paste("Vol_i [µm3]","\n","lm: R2=",signif(summary(test)$r.squared,a),", p=",signif(anova(test)$'Pr(>F)',a),
               ", n=",length(df$dt[which(is.na(df$vi)==FALSE & is.na(df$dt)==FALSE)]),
               "a=",signif(coef(test)[2],a), ", b=", signif(coef(test)[1],a),
               "\n rlm: a=",signif(coef(testb)[2],2), ", b=", signif(coef(testb)[1],2),
               "\n Pearson: r=", signif(cor(x,y),a), "p=",signif(cor.test(x,y)$p.value,a)))+
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
x=df$vi[which(is.na(df$vi)==FALSE & is.na(df$vf)==FALSE)]
y=df$vf[which(is.na(df$vi)==FALSE & is.na(df$vf)==FALSE)]
test=lm(df$vf~df$vi)
testb=rlm(y~x)
p5 <- ggplot(df,aes(x=vi,y=vf, color=factor(tri))) +
  geom_point() +
  labs(title= title, y="Vf [µm3]", 
       x=paste("Vol_i [µm3]","\n","lm: R2=",signif(summary(test)$r.squared,a),", p=",signif(anova(test)$'Pr(>F)',a),
               ", n=",length(df$vf[which(is.na(df$vf)==F & is.na(df$vi)==F)]),
               "a=",signif(coef(test)[2],a), ", b=", signif(coef(test)[1],a),
               "\n rlm: a=",signif(coef(testb)[2],2), ", b=", signif(coef(testb)[1],2),
               "\n Pearson: r=", signif(cor(x,y),a), "p=",signif(cor.test(x,y)$p.value,a),
               "\n", "mean DV=",signif(mean(df$dv,na.rm=TRUE),a)))+
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
x=df$dt[which(is.na(df$dt)==FALSE & is.na(df$dv)==FALSE)]
y=df$dv[which(is.na(df$dt)==FALSE & is.na(df$dv)==FALSE)]
test=lm(y~x)
testb=rlm(y~x)
p6 <- ggplot(df,aes(x=dt,y=dv, color=factor(tri))) +
  geom_point() +
  labs(title= title, y="DV [µm3]", 
       x=paste("DT [hrs]","\n","R2=",signif(summary(test)$r.squared,a),", p=",signif(anova(test)$'Pr(>F)',a),
               ", n=",length(df$dt[which(is.na(df$dt)==FALSE & is.na(df$dv)==FALSE)]),
               "a=",signif(coef(test)[2],a), ", b=", signif(coef(test)[1],a),
               "\n rlm: a=",signif(coef(testb)[2],2), ", b=", signif(coef(testb)[1],2),
               "\n Pearson: r=", signif(cor(x,y),a), "p=",signif(cor.test(x,y)$p.value,a)))+
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
x=df$dt[which(is.na(df$dt)==FALSE & is.na(df$dvvi)==FALSE)]
y=df$dvvi[which(is.na(df$dt)==FALSE & is.na(df$dvvi)==FALSE)]
test=lm(y~x)
testb=rlm(y~x)
p7 <- ggplot(df,aes(x=dt,y=dvvi,color=factor(tri))) +
  geom_point() +
  labs(title= title, y="DV/Vi %", 
       x=paste("DT [hrs]","\n","R2=",signif(summary(test)$r.squared,a),", p=",signif(anova(test)$'Pr(>F)',a),
               ", n=",length(x),
               "a=",signif(coef(test)[2],a), ", b=", signif(coef(test)[1],a),
               "\n rlm: a=",signif(coef(testb)[2],2), ", b=", signif(coef(testb)[1],2),
               "\n Pearson: r=", signif(cor(x,y),a), "p=",signif(cor.test(x,y)$p.value,a)))+
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
x=df$vi[which(is.na(df$vi)==FALSE & is.na(df$dvvi)==FALSE)]
y=df$dvvi[which(is.na(df$vi)==FALSE & is.na(df$dvvi)==FALSE)]
test=lm(y~x)
testb=rlm(y~x)
p8 <- ggplot(df,aes(x=vi,y=dvvi, color=factor(tri))) +
  geom_point() +
  labs(title= title, y="DV/Vi %", 
       x=paste("vi [µm3]","\n","R2=",signif(summary(test)$r.squared,a),", p=",signif(anova(test)$'Pr(>F)',a),
               ", n=",length(x),
               "a=",signif(coef(test)[2],a), ", b=", signif(coef(test)[1],a),
               "\n rlm: a=",signif(coef(testb)[2],2), ", b=", signif(coef(testb)[1],2),
               "\n Pearson: r=", signif(cor(x,y),a), "p=",signif(cor.test(x,y)$p.value,a)))+
  scale_colour_discrete (guide=FALSE) +
#   geom_abline(intercept = coef(testb)[1], slope = coef(testb)[2])+
  stat_smooth(method="lm", se=FALSE)+
  theme_bw()
sink(file=paste(dir,title,"_lm_rlm_coeff-corr_dvvi_vi.txt",sep=""))
summary(test)
summary(testb)
cor.test(x,y)
sink()
rm(test,testb,y,x)

## make figure #####
tiff(paste(dir,title,"_abs_val_",name_tri,".tiff",sep=""), width=1000,height=1000)
grid.arrange(p2,p1,p3,p4,p5,p6,p7,p8,ncol=2)
dev.off()
# rm(df)
#rm(p1,p2,p3,p4,p5,p6,p7,p8)

####### comparison in lineages, correlations between sisters/mothers/cousins ######
#### 1) ratios between sisters for vf, vi, dv, dt ####
m = 1.6
b = 240
###### 1.1 vf/vi and dv/vi ####
### tri ####
vf_1 = NULL
vf_2 = NULL
vi_1 = NULL
vi_2 = NULL
dt_1 = NULL
dt_2 = NULL
lineage = NULL
number_in_lineage = NULL
for (i in c(1:(max(data$lineage)))){
  for (j in unique(data$number_in_lineage[which(data$lineage==i & data$number_in_lineage %%2==F)])){
    if (m %in% data$cycle_stage[which(data$lineage==i & data$number_in_lineage==j)]
        && m %in% data$cycle_stage[which(data$lineage==i & data$number_in_lineage==j+1)]
        && b %in% data$cycle_stage[which(data$lineage==i & data$number_in_lineage==j)]
        && b %in% data$cycle_stage[which(data$lineage==i & data$number_in_lineage==j+1)]
        && is.na(data$vol[which(data$lineage==i & data$number_in_lineage==j & data$cycle_stage==m)])==F
        && is.na(data$vol[which(data$lineage==i & data$number_in_lineage==j+1 & data$cycle_stage==m)])==F
        && is.na(data$vol[which(data$lineage==i & data$number_in_lineage==j & data$cycle_stage==b)])==F
        && is.na(data$vol[which(data$lineage==i & data$number_in_lineage==j+1 & data$cycle_stage==b)])==F){
      vf_1 = c(vf_1,data$vol[which(data$lineage==i & data$number_in_lineage==j & data$cycle_stage==m)])
      vf_2 = c(vf_2,data$vol[which(data$lineage==i & data$number_in_lineage==j+1 & data$cycle_stage==m)])
      vi_1 = c(vi_1,data$vol[which(data$lineage==i & data$number_in_lineage==j & data$cycle_stage==b)])
      vi_2 = c(vi_2,data$vol[which(data$lineage==i & data$number_in_lineage==j+1 & data$cycle_stage==b)])
      dt_1 = c(dt_1,data$Time[which(data$lineage==i & data$number_in_lineage==j & data$cycle_stage==m)]-
                 data$Time[which(data$lineage==i & data$number_in_lineage==j & data$cycle_stage==b)])
      dt_2 = c(dt_2,data$Time[which(data$lineage==i & data$number_in_lineage==j+1 & data$cycle_stage==m)]
               - data$Time[which(data$lineage==i & data$number_in_lineage==j+1 & data$cycle_stage==b)])
      lineage = c(lineage,i)
      number_in_lineage = c(number_in_lineage,j)
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

test=lm(df$vf_ratio~df$vi_ratio)
p1 <- ggplot(df,aes(x=vi_ratio,y=vf_ratio)) +
  geom_point() +
  labs(title= title, y="Vf ratio d1/d2 %", 
       x=paste("Vi ratio d1/d2 %","\n","R2=",summary(test)$r.squared,", p=",anova(test)$'Pr(>F)',", n=",length(df$vf_ratio[which(is.na(df$vf_ratio)==F)]),
               "\n","a=",coef(test)[2], ", b=", coef(test)[1]))+
  theme_bw()
sink(file=paste(dir,title,"_ratio_sisters_lm_vf_vi.txt",sep=""))
summary(test)
sink()
rm(test)

test=lm(df$dt_ratio~df$vi_ratio)
p2 <- ggplot(df,aes(x=vi_ratio,y=dt_ratio)) +
  geom_point() +
  labs(title= title, y="DT ratio d1/d2 %", 
       x=paste("Vi ratio d1/d2 %","\n","R2=",summary(test)$r.squared,", p=",anova(test)$'Pr(>F)',
               ", n=",length(df$dt_ratio[which(is.na(df$dt_ratio)==F)]),
               "\n","a=",coef(test)[2], ", b=", coef(test)[1]))+
  #scale_colour_discrete (name="lineage", labels=c("daughter 1","daughter 2")) +
  theme_bw()
sink(file=paste(dir,title,"_ratio_sisters_lm_dt_vi.txt",sep=""))
summary(test)
sink()
rm(test)

test=lm(df$dv_ratio~df$vi_ratio)
p3 <- ggplot(df,aes(x=vi_ratio,y=dv_ratio)) +
  geom_point() +
  labs(title= title, y="DV ratio d1/d2 %", 
       x=paste("Vi ratio d1/d2 %","\n","R2=",summary(test)$r.squared,", p=",anova(test)$'Pr(>F)',
               ", n=",length(df$dv_ratio[which(is.na(df$dv_ratio)==F)]),
               "\n","a=",coef(test)[2], ", b=", coef(test)[1]))+
  #scale_colour_discrete (name="lineage", labels=c("daughter 1","daughter 2")) +
  theme_bw()
sink(file=paste(dir,title,"_ratio_sisters_lm_dv_vi.txt",sep=""))
summary(test)
sink()
rm(test)

test=lm(df$dv_ratio~df$dt_ratio)
p4 <- ggplot(df,aes(x=dt_ratio,y=dv_ratio)) +
  geom_point() +
  labs(title= title, y="dv ratio d1/d2 %", 
       x=paste("dt ratio d1/d2 %","\n","R2=",summary(test)$r.squared,", p=",anova(test)$'Pr(>F)',
                ", n=",length(df$dv_ratio[which(is.na(df$dv_ratio)==F)]),
               "\n","a=",coef(test)[2], ", b=", coef(test)[1])) +
  #scale_colour_discrete (name="lineage", labels=c("daughter 1","daughter 2")) +
  theme_bw()
sink(file=paste(dir,title,"_ratio_sisters_lm_dv_dt.txt",sep=""))
summary(test)
sink()
rm(test)

###### 1.2 dt/vi (use the other one justs above) ####
dt_1 = NULL
dt_2 = NULL
vi_1 = NULL
vi_2 = NULL
lineage = NULL
number_in_lineage = NULL
for (i in c(1:(max(data$lineage)))){
  print(i)
  for (j in unique(data$number_in_lineage[which(data$lineage==i & data$number_in_lineage %%2==F)])){
    if (m %in% data$cycle_stage[which(data$lineage==i & data$number_in_lineage==j)]
        && m %in% data$cycle_stage[which(data$lineage==i & data$number_in_lineage==j+1)]
        && b %in% data$cycle_stage[which(data$lineage==i & data$number_in_lineage==j)]
        && b %in% data$cycle_stage[which(data$lineage==i & data$number_in_lineage==j+1)]
        && is.na(data$vol[which(data$lineage==i & data$number_in_lineage==j & data$cycle_stage==b)])==F
        && is.na(data$vol[which(data$lineage==i & data$number_in_lineage==j+1 & data$cycle_stage==b)])==F){
      dt_1 = c(dt_1,data$Time[which(data$lineage==i & data$number_in_lineage==j & data$cycle_stage==m)]-
                 data$Time[which(data$lineage==i & data$number_in_lineage==j & data$cycle_stage==b)])
      dt_2 = c(dt_2,data$Time[which(data$lineage==i & data$number_in_lineage==j+1 & data$cycle_stage==m)]
               - data$Time[which(data$lineage==i & data$number_in_lineage==j+1 & data$cycle_stage==b)])
      vi_1 = c(vi_1,data$vol[which(data$lineage==i & data$number_in_lineage==j & data$cycle_stage==b)])
      vi_2 = c(vi_2,data$vol[which(data$lineage==i & data$number_in_lineage==j+1 & data$cycle_stage==b)])
      lineage = c(lineage,i)
      number_in_lineage = c(number_in_lineage,j)
    }
  }
}
df <- data.frame(dt_1,dt_2,vi_1,vi_2, lineage, number_in_lineage)
df$dt_ratio <- df$dt_2/df$dt_1 *100
df$vi_ratio <- df$vi_2/df$vi_1 *100
rm(vf_1,vf_2,vi_1,vi_2, lineage, number_in_lineage)
p2 <- ggplot(df,aes(x=vi_ratio,y=dt_ratio, color=factor(lineage))) +
  geom_point() +
  labs(title= title, y="dt ratio d2/d1 %", x="Vi ratio d2/d1 %") +
  #scale_colour_discrete (name="lineage", labels=c("daughter 1","daughter 2")) +
  theme_bw()

###### save 1) ####
tiff(paste(dir,title,"_ratios_sisters_dv_dt_vi_vf.tiff",sep=""), width=750,height=500)
grid.arrange(p1,p2,p3,p4,ncol=2)
dev.off()
rm(p1,p2,p3,p4)

### 2) correlations sister/sister sister/mother sister/cousin ####
m = 1.6
b = 240

#### 2) lineages for vf and dt : comparison mother/daughter sister/sister cousin/cousin ####
m=1.6
b=240
###### 2.1 mother/daughter ####
######## tri ####
v1 = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15)
v2 = c(2,4,6,8,10,12,14,16,18,20,22,24,26,28,30)
v3 = c(3,5,7,9,11,13,15,17,19,21,23,25,27,29,31)
dff = data.frame(v1,v2,v3)
tf_m_ccl = NULL
tf_d_ccl = NULL
ti_m_ccl = NULL
ti_d_ccl = NULL
lineage_ccl = NULL
mother_number_ccl = NULL
daughter_number_ccl = NULL

vf_m = NULL
vf_d = NULL
lineage_vf = NULL
mother_number_vf = NULL
daughter_number_vf = NULL

for (i in unique(data$lineage)){
  for (j in c(1:length(dff[,1]))) {
    if (dff[j,1] %in% data$number_in_lineage[which(data$lineage==i)]
        && m %in% data$cycle_stage[which(data$lineage==i,data$number_in_lineage==dff[j,1])]){
      for (k in c(dff[j,2],dff[j,3])){
        if (k %in% data$number_in_lineage[which(data$lineage==i)]
            && m %in% data$cycle_stage[which(data$lineage==i & data$number_in_lineage==k)]){
          if (is.na(data$vol[which(data$lineage==i & data$number_in_lineage==dff[j,1] & data$cycle_stage==m)])==F
              && is.na(data$vol[which(data$lineage==i & data$number_in_lineage==k & data$cycle_stage==m)])==F){
            vf_m = c(vf_m, data$vol[which(data$lineage==i & data$number_in_lineage==dff[j,1] & data$cycle_stage==m)])
            vf_d = c(vf_d, data$vol[which(data$lineage==i & data$number_in_lineage==k & data$cycle_stage==m)])
            lineage_vf = c(lineage_vf,i)
            mother_number_vf = c(mother_number_vf,dff[j,1])
            daughter_number_vf = c(daughter_number_vf, k)
          } else {}
          if (b %in% data$cycle_stage[which(data$lineage==i & data$number_in_lineage==k)]
              && b %in% data$cycle_stage[which(data$lineage==i & data$number_in_lineage==dff[j,1])]){
            tf_m_ccl <- c(tf_m_ccl, data$Time[which(data$lineage==i & data$number_in_lineage==dff[j,1] & data$cycle_stage==m)])
            tf_d_ccl <- c(tf_d_ccl, data$Time[which(data$lineage==i & data$number_in_lineage==k & data$cycle_stage==m)])
            ti_m_ccl <- c(ti_m_ccl, data$Time[which(data$lineage==i & data$number_in_lineage==dff[j,1] & data$cycle_stage==b)])
            ti_d_ccl <- c(ti_d_ccl, data$Time[which(data$lineage==i & data$number_in_lineage==k & data$cycle_stage==b)])
            lineage_ccl <- c(lineage_ccl, i)
            mother_number_ccl <- c(mother_number_ccl,dff[j,1])
            daughter_number_ccl <- c(daughter_number_ccl, k)
          } else {}
        }
      }
    }
  }
}

######## plots ####
df_ccl <- data.frame ( tf_m_ccl, tf_d_ccl, ti_m_ccl, ti_d_ccl, lineage_ccl, mother_number_ccl, daughter_number_ccl)
df_ccl$dt_m <- (df_ccl$tf_m_ccl-df_ccl$ti_m_ccl)/60
df_ccl$dt_d <- (df_ccl$tf_d_ccl-df_ccl$ti_d_ccl)/60
rm(tf_m_ccl, tf_d_ccl, ti_m_ccl, ti_d_ccl, lineage_ccl, mother_number_ccl, daughter_number_ccl)
df_vf <- data.frame(vf_m,vf_d,lineage_vf,mother_number_vf,daughter_number_vf)
rm(vf_m,vf_d,lineage_vf,mother_number_vf,daughter_number_vf)

test = lm(df_ccl$dt_d~df_ccl$dt_m)
p1 <- ggplot(df_ccl,aes(x=dt_m,y=dt_d, color=factor(lineage_ccl))) +
  geom_point() +
  labs(title= title, y="DT daughter [hrs]", 
       x= paste("DT mother [hrs]",
                "\n","R2=",summary(test)$r.squared,", p=",anova(test)$'Pr(>F)',", n=",length(df_ccl$dt_m),
                "\n","a=",coef(test)[2], ", b=", coef(test)[1])) + 
  scale_colour_discrete (name="lineage") +
  theme_bw()
sink(file=paste(dir,title,"_lm_corr_lineages_dt_mother_daughter.txt",sep=""))
summary(test)
sink() 
rm(test)

test = lm(df_vf$vf_d~df_vf$vf_m)
pa <- ggplot(df_vf,aes(x=vf_m,y=vf_d, color=factor(lineage_vf))) +
  geom_point() +
  labs(title= title, y="Vf daughter [µm3]", 
       x=paste("Vf mother [µm3]",
               "\n","R2=",summary(test)$r.squared,", p=",anova(test)$'Pr(>F)',", n=",length(df_vf$vf_m),
               "\n","a=",coef(test)[2], ", b=", coef(test)[1]))+ 
  scale_colour_discrete (name="lineage") +
  theme_bw()
sink(file=paste(dir,title,"_lm_corr_lineages_vf_mother_daughter.txt",sep=""))
sink()
rm(test)

rm(dff,df_ccl,df_vf)

##### 2.2 sister/sister (works but use 2.3 for the plots)####
####### tri ####
vf_d1 = NULL
vf_d2 = NULL
dt_d1 = NULL
dt_d2 = NULL
lineage_ccl = NULL
lineage_vf = NULL
number_in_lineage_d1_ccl = NULL
number_in_lineage_d1_vf = NULL
for (i in unique(data$lineage)){
  for (j in unique(data$number_in_lineage[which(data$lineage==i & data$number_in_lineage %%2==F)])){
    if (m %in% data$cycle_stage[which(data$lineage==i & data$number_in_lineage==j)]
        && m %in% data$cycle_stage[which(data$lineage==i & data$number_in_lineage==j+1)]){
      if (is.na(data$vol[which(data$lineage==i & data$number_in_lineage==j & data$cycle_stage==m)])==F
          && is.na(data$vol[which(data$lineage==i & data$number_in_lineage==j+1 & data$cycle_stage==m)])==F){
        vf_d1 = c(vf_d1,data$vol[which(data$lineage==i & data$number_in_lineage==j & data$cycle_stage==m)])
        vf_d2 = c(vf_d2,data$vol[which(data$lineage==i & data$number_in_lineage==j+1 & data$cycle_stage==m)])
        lineage_vf = c(lineage_vf,i)
        number_in_lineage_d1_vf = c(number_in_lineage_d1_vf,j)
      } else {}
      if (b %in% data$cycle_stage[which(data$lineage==i & data$number_in_lineage==j)]
          && b %in% data$cycle_stage[which(data$lineage==i & data$number_in_lineage==j+1)]){
        dt_d1 = c(dt_d1,data$Time[which(data$lineage==i & data$number_in_lineage==j & data$cycle_stage==m)]-
                    data$Time[which(data$lineage==i & data$number_in_lineage==j & data$cycle_stage==b)])
        dt_d2 = c(dt_d2,data$Time[which(data$lineage==i & data$number_in_lineage==j+1 & data$cycle_stage==m)]
                  - data$Time[which(data$lineage==i & data$number_in_lineage==j+1 & data$cycle_stage==b)])
        lineage_ccl = c(lineage_ccl,i)
        number_in_lineage_d1_ccl = c(number_in_lineage_d1_ccl,j)
      }else {}
    } else {}
  }
}

####### plots ####
df_ccl <- data.frame(dt_d1,dt_d2, lineage_ccl, number_in_lineage_d1_ccl)
df_vf <- data.frame(vf_d1,vf_d2,lineage_vf,number_in_lineage_d1_vf)
rm(dt_d1,dt_d2, lineage_ccl, number_in_lineage_d1_ccl)
rm(vf_d1,vf_d2,lineage_vf,number_in_lineage_d1_vf)
p2 <- ggplot(df_ccl,aes(x=dt_d1,y=dt_d2, color=factor(lineage_ccl))) +
  geom_point() +
  labs(title= title, y="DT daughter [hrs]", x="DT mother [hrs]") +
  #scale_colour_discrete (name="lineage", labels=c("daughter 1","daughter 2")) +
  theme_bw()
sink(file=paste(dir,title,"_lm_corr_lineages_dt_daughter_daughter.txt",sep=""))
summary(lm(df_ccl$dt_d2~df_ccl$dt_d1))
sink() 

pb <- ggplot(df_vf,aes(x=vf_d1,y=vf_d2, color=factor(lineage_vf))) +
  geom_point() +
  labs(title= title, y="Vf sister 2 [µm3]", x="Vf sister 1 [µm3]") +
  theme_bw()
sink(file=paste(dir,title,"_lm_corr_lineages_vf_sister_sister.txt",sep=""))
summary(lm(df_vf$vf_d2~df_vf$vf_d1))
sink()

rm(df_ccl,df_vf)

##### 2.3 cousin/cousin and sister/sister ####
m=1.6
b=240
####### tri ####
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

for (i in unique(data$lineage)){
  for (j in unique(data$generation[which(data$lineage==i)])) {
    #print("generation")
    #print(j)
    for (k in unique(data$number_in_lineage[which(data$lineage==i & data$generation==j)])){
      l=k+1
      while (max(unique(data$number_in_lineage[which(data$lineage==i & data$generation==j)]))+1 > l & l > k) {
        if (m %in% data$cycle_stage[which(data$lineage==i & data$number_in_lineage==k)]
            && m %in% data$cycle_stage[which(data$lineage==i & data$number_in_lineage==l)]) {
          if (k%%2==F && l==k+1){
            s1_dt <- c(s1_dt, data$Time[which(data$lineage==i & data$number_in_lineage==k & data$cycle_stage==m)] - data$Time[which(data$lineage==i & data$number_in_lineage==k & data$cycle_stage==b)])
            s2_dt <- c(s2_dt, data$Time[which(data$lineage==i & data$number_in_lineage==l & data$cycle_stage==m)] - data$Time[which(data$lineage==i & data$number_in_lineage==l & data$cycle_stage==b)])
            s1_number_in_lineage_dt <- c(s1_number_in_lineage_dt,k)
            s2_number_in_lineage_dt <- c(s2_number_in_lineage_dt,l)
            sister_lineage_dt <- c(sister_lineage_dt, i)
          } else {
            c1_dt <- c(c1_dt, data$Time[which(data$lineage==i & data$number_in_lineage==k & data$cycle_stage==m)] - data$Time[which(data$lineage==i & data$number_in_lineage==k & data$cycle_stage==b)])
            c2_dt <- c(c2_dt, data$Time[which(data$lineage==i & data$number_in_lineage==l & data$cycle_stage==m)] - data$Time[which(data$lineage==i & data$number_in_lineage==l & data$cycle_stage==b)])
            c1_number_in_lineage_dt <- c(c1_number_in_lineage_dt,k)
            c2_number_in_lineage_dt <- c(c2_number_in_lineage_dt,l)
            cousin_generation_dt <- c(cousin_generation_dt, j)
            cousin_lineage_dt <- c(cousin_lineage_dt, i)
            #print("cousin")
            #print(l)
          }
          if (is.na(data$vol[which(data$lineage==i & data$number_in_lineage==k & data$cycle_stage==m)])==F
              && is.na(data$vol[which(data$lineage==i & data$number_in_lineage==l & data$cycle_stage==m)])==F){
            if (k%%2==F && l==k+1){
              s1_vf <- c(s1_vf, data$vol[which(data$lineage==i & data$number_in_lineage==k & data$cycle_stage==m)])
              s2_vf <- c(s2_vf, data$vol[which(data$lineage==i & data$number_in_lineage==l & data$cycle_stage==m)])
              s1_number_in_lineage_vf <- c(s1_number_in_lineage_vf,k)
              s2_number_in_lineage_vf <- c(s2_number_in_lineage_vf,l)
              sister_lineage_vf <- c(sister_lineage_vf, i)
              #print("sister")
              #print(l)
            } else {
              c1_vf <- c(c1_vf, data$vol[which(data$lineage==i & data$number_in_lineage==k & data$cycle_stage==m)])
              c2_vf <- c(c2_vf, data$vol[which(data$lineage==i & data$number_in_lineage==l & data$cycle_stage==m)])
              c1_number_in_lineage_vf <- c(c1_number_in_lineage_vf,k)
              c2_number_in_lineage_vf <- c(c2_number_in_lineage_vf,l)
              cousin_generation_vf <- c(cousin_generation_vf, j)
              cousin_lineage_vf <- c(cousin_lineage_vf, i)
              #print("cousin")
              #print(l)
            }
          }
        } else {}
        l=l+1
      }
    }
  }
}

####### plots ####
df_dt_cousin <- data.frame(c1_dt/60, c2_dt/60, c1_number_in_lineage_dt, c2_number_in_lineage_dt, cousin_generation_dt, cousin_lineage_dt)
df_dt_sister <- data.frame(s1_dt/60, s2_dt/60, s1_number_in_lineage_dt, s2_number_in_lineage_dt, sister_lineage_dt)
df_vf_cousin <- data.frame(c1_vf, c2_vf, c1_number_in_lineage_vf, c2_number_in_lineage_vf, cousin_generation_vf, cousin_lineage_vf)
df_vf_sister <- data.frame(s1_vf, s2_vf, s1_number_in_lineage_vf, s2_number_in_lineage_vf, sister_lineage_vf)
rm(c1_dt, c2_dt, c1_number_in_lineage_dt, c2_number_in_lineage_dt, cousin_generation_dt, cousin_lineage_dt,
     s1_dt, s2_dt, s1_number_in_lineage_dt, s2_number_in_lineage_dt, sister_lineage_dt,
     c1_vf, c2_vf, c1_number_in_lineage_vf, c2_number_in_lineage_vf, cousin_generation_vf, cousin_lineage_vf,s1_vf, s2_vf, s1_number_in_lineage_vf, s2_number_in_lineage_vf, sister_lineage_vf)

## vf ####
test = lm(df_vf_sister$s2_vf~df_vf_sister$s1_vf)
pb <- ggplot(df_vf_sister,aes(x=s1_vf,y=s2_vf, color=factor(sister_lineage_vf))) +
  geom_point() +
  labs(title= title, y="Vf sister 2 [µm3]", 
       x= paste("Vf sister 1 [µm3]","\n","R2=",summary(test)$r.squared,", p=",anova(test)$'Pr(>F)',
                ", n=",length(df_vf_sister$s1_vf),
                "\n","a=",coef(test)[2], ", b=", coef(test)[1])) +
  theme_bw()
sink(file=paste(dir,title,"_lm_corr_lineages_vf_sister_sister.txt",sep=""))
summary(test)
sink()
rm(test)

test = lm(df_vf_cousin$c2_vf~df_vf_cousin$c1_vf)
pc1 <- ggplot(df_vf_cousin,aes(x=c1_vf,y=c2_vf, color=factor(cousin_lineage_vf))) +
  geom_point() +
  labs(title= title, y="vf cousin 2 [µm3]", x= paste("Vf cousin 1 [µm3]","\n","R2=",summary(test)$r.squared,", p=",anova(test)$'Pr(>F)',", n=",length(df_vf_cousin$c1_vf),"\n","a=",coef(test)[2], ", b=", coef(test)[1])) +
  theme_bw()
pc2 <- ggplot(df_vf_cousin,aes(x=c1_vf,y=c2_vf, color=factor(cousin_generation_vf))) +
  geom_point() +
  labs(title= title, y="vf cousin 2 [µm3]", x= paste("Vf cousin 1 [µm3]","\n","R2=",summary(test)$r.squared,", p=",anova(test)$'Pr(>F)',", n=",length(df_vf_cousin$c1_vf),"\n","a=",coef(test)[2], ", b=", coef(test)[1])) +
  theme_bw()
pc3 <- ggplot(df_vf_cousin,aes(x=c1_vf,y=c2_vf, color=factor(c1_number_in_lineage_vf))) +
  geom_point() +
  labs(title= title, y="vf cousin 2 [µm3]", x= paste("Vf cousin 1 [µm3]","\n","R2=",summary(test)$r.squared,", p=",anova(test)$'Pr(>F)',", n=",length(df_vf_cousin$c1_vf),"\n","a=",coef(test)[2], ", b=", coef(test)[1])) +
  theme_bw()
sink(file=paste(dir,title,"_lm_corr_lineages_vf_cousin_cousin.txt",sep=""))
summary(test)
sink()
rm(test)
#tiff(paste(dir,title,"_vf_cousin_cousin_3different_labels.tiff",sep=""), width=500,height=750)
#grid.arrange(pc1,pc2,pc3,ncol=1)
#dev.off()

## dt ####
test= lm(df_dt_sister$s2_dt.60~df_dt_sister$s1_dt.60)
p2 <- ggplot(df_dt_sister,aes(x=df_dt_sister$s1_dt.60,y=df_dt_sister$s2_dt.60, color=factor(sister_lineage_dt))) +
  geom_point() +
  labs(title= title, y="DT sister 2 [hrs]", x=paste("DT sister1 [hrs]","\n","R2=",summary(test)$r.squared,", p=",anova(test)$'Pr(>F)',", n=",length(df_dt_sister$s1_dt.60),"\n","a=",coef(test)[2], ", b=", coef(test)[1])) +
  scale_colour_discrete (name="lineage") +
  theme_bw()
sink(file=paste(dir,title,"_lm_corr_lineages_dt_sister_sister.txt",sep=""))
summary(test)
sink() 
rm(test)

test = lm(df_dt_cousin$c2_dt.60~df_dt_cousin$c1_dt.60)
p3 <- ggplot(df_dt_cousin,aes(x=df_dt_cousin$c1_dt.60,y=df_dt_cousin$c2_dt.60, color=factor(cousin_lineage_dt))) +
  geom_point() +
  labs(title= title, y="DT cousin 2 [hrs]", x=paste("DT cousin 1 [hrs]","\n","R2=",summary(test)$r.squared,", p=",anova(test)$'Pr(>F)',", n=",length(df_dt_cousin$c1_dt.60),"\n","a=",coef(test)[2], ", b=", coef(test)[1])) +
  scale_colour_discrete (name="lineage") +
  theme_bw()
sink(file=paste(dir,title,"_lm_corr_lineages_dt_cousin_cousin.txt",sep=""))
summary(test)
sink() 
rm(test)

#### 2.4 figure with all plots for 2) ####
tiff(paste(dir,title,"_corr_lineages_vf_m_d_c.tiff",sep=""), width=500,height=750)
grid.arrange(pa,pb,pc1,ncol=1)
dev.off()

tiff(paste(dir,title,"_corr_lineages_dt_m_d_c.tiff",sep=""), width=500,height=750)
grid.arrange(p1,p2,p3,ncol=1)
dev.off()





