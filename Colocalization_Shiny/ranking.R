library(plyr)
library(dplyr)

#'cytosol','mitochondria','ER','Nucleus','Vacuole','Plasma,membrane'

##Loads overlap of 10% brightest pixels
data = list()
data2 = list()
sorteddata = list()
names= c('1m3u','1pok','2cg4','2vyc','2iv1','1frw','1yac')
for (i in c(1:length(names))){
  data[[i]] = read.csv(paste0("/media/elusers/users/Tal/elmicro/2020_12_21_chaperone2/Colocalization-App2/Data/colocalization-data-",names[i],".csv"),stringsAsFactors = FALSE) 
  data2[[i]] = read.csv(paste0("/media/elusers/users/Tal/elmicro/2020_12_21_chaperone2/Colocalization-App2/Data/colocalization2-data-",names[i],".csv"),stringsAsFactors = FALSE)
  data2[[i]]=data2[[i]][c(1:90),]
  removecol = c('Peroxisome','golgi','lipid.droplet','SGA.hit')
  data2[[i]] = data2[[i]][,!(names(data2[[i]]) %in% removecol)]
  sorteddata[[i]] = rbind(data[[i]],data2[[i]])
  sorteddata[[i]] = sorteddata[[i]][order(sorteddata[[i]]$fiber, decreasing = TRUE),]
  sorteddata[[i]]$Localization = ''
  sorteddata[[i]]$Localization[which(sorteddata[[i]]$cytosol ==1)] = paste0(sorteddata[[i]]$Localization[which(sorteddata[[i]]$cytosol ==1)], 'Cytosol')
  sorteddata[[i]]$Localization[which(sorteddata[[i]]$Nucleus ==1)] = paste(sorteddata[[i]]$Localization[which(sorteddata[[i]]$Nucleus ==1)], 'Nucleus')
  sorteddata[[i]]$Localization[which(sorteddata[[i]]$mitochondria ==1)] = paste(sorteddata[[i]]$Localization[which(sorteddata[[i]]$mitochondria ==1)], 'Mitochondria')
  sorteddata[[i]]$Localization[which(sorteddata[[i]]$ER ==1)] = paste(sorteddata[[i]]$Localization[which(sorteddata[[i]]$ER ==1)], 'Endoplasmic-Reticulum')
  sorteddata[[i]]$Localization[which(sorteddata[[i]]$Vacuole ==1)] = paste(sorteddata[[i]]$Localization[which(sorteddata[[i]]$Vacuole ==1)], 'Vacuole')
  sorteddata[[i]]$Localization[which(sorteddata[[i]]$Plasma.membrane ==1)] = paste(sorteddata[[i]]$Localization[which(sorteddata[[i]]$Plasma.membrane ==1)], 'Plasma-membrane')
  
}
names(data)=names
names(data2)=names
names(sorteddata)=names

#Creates localization column


#Data frame of ranked proteins according to fibers only
ranked = data.frame(sorteddata$`1m3u`$names)
names(ranked)='1m3u'
for (i in c(2:length(sorteddata))){
  ranked[,names[i]] = sorteddata[[i]]$names
} 

#unique proteins names of top5 hits
prot_hits = c()
z = c() #Tracks the Z (colocalization/intensity) of each hit, sums if the hit repeats in different plates
for (lib in c(1:length(ranked[c(1:5),]))){
  for (prot in ranked[c(1:5),lib]) {
    if (!(prot %in% prot_hits)){
      prot_hits=c(prot_hits,list(prot,1,names(ranked)[lib]))
      z=c(z,sorteddata[[lib]][[which(sorteddata[[lib]]$names == prot)[[1]],'fiber']])
    }
    else{
      prot_hits[(which(prot_hits == prot)+1)] = as.numeric(prot_hits[(which(prot_hits == prot)+1)])+1
      prot_hits[(which(prot_hits == prot)+2)] = paste(prot_hits[(which(prot_hits == prot)+2)],names(ranked)[lib])
      z[(which(prot_hits == prot)+2)/3] = z[(which(prot_hits == prot)+2)/3] + 
        sorteddata[[lib]][[which(sorteddata[[lib]]$names == prot)[[1]],input$sCategory]]
    }
  }
}
top_hit_freq=data.frame(prot_hits[which(c(1:length(prot_hits)) %% 3 == 1)], stringsAsFactors=FALSE)
top_hit_freq = t(top_hit_freq)
top_hit_freq = as.data.frame(top_hit_freq, stringsAsFactors=FALSE)
names(top_hit_freq) = 'protein'
df = data.frame(prot_hits[which(c(1:length(prot_hits)) %% 3 == 2)], stringsAsFactors=FALSE)
df = t(df)
df = as.data.frame(df, stringsAsFactors=FALSE)
top_hit_freq$freq = df
top_hit_freq[,'fiber'] = z
top_hit_freq[,'fiberdiv'] = z/top_hit_freq$freq
dfscreen = data.frame(prot_hits[which(c(1:length(prot_hits)) %% 3 == 0)], stringsAsFactors=FALSE)
dfscreen = t(dfscreen)
dfscreen = as.data.frame(dfscreen, stringsAsFactors=FALSE)
top_hit_freq$screen = dfscreen


