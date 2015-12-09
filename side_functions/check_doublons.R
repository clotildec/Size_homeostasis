m=0
for (k in cond){  
  end = max(df$Time[which(df$sort_cond==k)])
  for (i in unique(df$lineage[which(df$sort_cond==k)])){
    for (j in unique(c(df$number_in_lineage[which(df$lineage==i & df$sort_cond==k)]))){
      m=m+1
      for (l in c(1:length(t))){
        if (t[l] %in% df$cycle_stage[which(df$lineage==i & df$number_in_lineage==j & df$sort_cond==k)]){ ### for volume and time at entry into the phase
          temp2=df$Time[which(df$lineage==i & df$number_in_lineage==j & df$cycle_stage==t[(l)] & df$sort_cond==k)]
          if (length(temp2)>1){
            print("k")
            print(k)
            print("i")
            print(i)
            print("j")
            print(j)
            print("l")
            print(l) 
            print(temp2)
          } else {}
        }
      }
    }
  }
}
