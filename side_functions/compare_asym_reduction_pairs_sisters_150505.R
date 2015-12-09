data <- read.csv("H:/150423_Raji_bilan/150505_combine Raji_131122-131220_red_asym/ 131122_131220_Raji_analysis_150505 _sisters_df.csv")
View(data)
dir="H:/150423_Raji_bilan/150505_combine Raji_131122-131220_red_asym/"
title="Raji_chambers_analysis_150505"

df <- data
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

df <- df[-which(df$vi_ratio>200),]
title=paste(title,"-dataset3",sep="")
write.csv(df,paste(dir,title,"df_pairs_sisters.csv"))

var_f <- abs(df$vf_1-df$vf_2)/((df$vf_1+df$vf_2)/2)
var_i <- abs(df$vi_1-df$vi_2)/((df$vi_1+df$vi_2)/2)
d <- data.frame(c(var_f,var_i),c(rep(2,length(var_f)),rep(1,length(var_i))))
colnames(d)=c("var","stage")
write.csv(d,paste(dir,title,"df_pairs_sisters_asym_red.csv"))
x=d$var[which(is.na(d$var)==F & d$stage==1)]
y=d$var[which(is.na(d$var)==F & d$stage==2)]
p2 <- ggplot(d[which(is.na(d$var)==F),],aes(factor(d$stage[which(is.na(d$var)==F)]),d$var[which(is.na(d$var)==F)]))+
  geom_boxplot()+
  labs(title=title,
       x=paste("t test, p=",signif(t.test(x,y)$p.value,a),", n=",length(x),
               "\nratio_i, m=", signif(mean(x),a),",sd=",signif(sd(x),a),
               "\nratio_f, m=",signif(mean(y),a),",sd=",signif(sd(y),a)),
       y="ratio volume in pair of sisters %")+
  scale_x_discrete (labels=c("initial","final"))+
  theme_bw()+
  theme(axis.text.x=element_text(size=12),
        text=element_text(size=12))

tiff(paste(dir,title,"_ratio_decrease_var_sister_boxplot.tiff",sep=""), width=250,height=250)
grid.arrange(p2)
dev.off()
