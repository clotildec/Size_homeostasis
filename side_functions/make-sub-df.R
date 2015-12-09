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
# by generation ####
df<-data[which(data$cond==0),]
df$dvdt_tot=df$dv_tot/df$dt_tot
df$dvdt_g1=df$dv_g1/df$dt_g1
df$dvdt_s_g2_m=df$dv_s_g2_m/df$dt_s_g2_m
t1 = unique(df$generation) ## enter the sequence of cycl_stage, in the order of appearance for one cell, from birth to division


t2= c(5,1.6,1.6)
names = c("birth","begining_g1_s","mitosis") ## enter all the timepoint you want to analyze for volume, time
names_phase = c("g1","s_g2_m","tot") ## enter all the phases you want to analyze for dv or dt, should be the name of the coumn in df1 -"dv" or "dt"
t_phase=c(240,5,240)
names_t_phase = c("birth","begining_g1_s","birth")


for (l in c(1:length(t1))) {
  d<-df[which(df$generation==t1[l]),]
  write.csv(d,paste(dir,title,"_sub_df_rosco20_generation",t1[l],".csv",sep=""))
}
