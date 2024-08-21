#Loads csv results of colocalization plate. gathers the 4 strains into 1 column
prepare_plate2 = function (name,file,to_order='',onlypuncta,norm='No'){
  #change localization folder
  if ('Yes' %in% onlypuncta){
      myImgResources2 <- paste0("imgResources/",name,"-2-sel-puncta-R_8-bit-",seq_len(384),".png")
      myImgResources <- paste0("imgResources/",name,"-2-sel-puncta_8-bit-",seq_len(384),".png")
    }
  else{
      myImgResources2 <- paste0("imgResources/",name,"-2-sel-R_8-bit-", seq_len(384), ".png")
      myImgResources <- paste0("imgResources/",name,"-2-sel_8-bit-",seq_len(384),".png")
    }
  df = read.csv(file = file)
  if ('RFPnorm' %in% norm){
    df[,2] = df[,2] / df[,5]
    df[,3] = df[,3] / df[,5]
    df[,4] = df[,4] / df[,5]
    df[,5] = df[,5] / df[,5]
  }
  ncells = gather(df,key = 'strain', 'ncells',19:22)
  YFP = gather(df,key = 'strain', 'YFPmean',23:26)
  RFP = gather(df,key = 'strain', 'RFPmean',27:30)
  df = gather(df,key = 'strain', 'Colocalization',2:5)
  df["image_url"] = myImgResources
  df= df[,-(15:26)]
  df["image_url2"] = myImgResources2
  df["ncells"] = ncells$ncells
  df["YFPmean"] = YFP$YFPmean
  df["RFPmean"] = RFP$RFPmean
  df["text"] = paste('n:',df$ncells,'\nYFP:',df$YFPmean,'\nRFP:',df$RFPmean)
  if (to_order == 'colocalization'){
    order = c(4,19,31,83,18,57,21,32,20,50,41,51,40,33,37,30,1,15,80,17,84,47,52,35,49,96,95,94,93,91,92,90,44,71,72,6,70,86,88,2,8,66,
              25,24,64,81,54,69,27,23,53,87,26,77,29,36,85,63,67,76,79,61,13,39,28,74,55,65,68,62,73,82,5,7,10,75,3,9,34,46,22,43,42,56,
              11,16,78,12,59,38,45,14,48,58,60,89)}
  else if (to_order == 'correlation'){
    order = c(51,21,57,32,41,18,4,19,20,50,83,31,33,37,89,58,14,45,48,12,59,40,47,52,35,49,2,60,16,11,78,30,1,9,3,80,15,5,7,10,75,42,22,
              23,34,43,46,6,90,79,27,53,44,87,38,56,88,55,86,28,66,13,77,61,74,26,68,85,36,54,39,29,76,70,71,82,62,63,67,73,65,72,25,96,
              95,94,93,91,92,64,81,8,24,69,17,84)}
  else if (to_order == 'RFPmean'){
    order=c(1,80,43,48,42,32,24,58,59,12,17,14,22,34,96,95,94,93,91,92,25,36,49,35,26,47,69,74,63,73,72,13,67,29,65,71,52,55,39,53,64,61,
            62,79,86,44,88,77,27,76,40,56,54,23,8,21,37,38,41,85,70,28,57,51,81,60,33,46,66,68,7,11,30,84,90,45,10,89,5,31,15,2,3,16,9,50,
            20,4,18,82,6,78,87,19,75,83)}
  else if (to_order == 'N Wt'){
    order =order(df$ncells[289:384])
  }
  else if (to_order == 'Color Wt'){
    order=order(df$Colocalization[289:384])
  }
  else if (to_order == 'N misfolded'){
    order =order(df$ncells[97:192])
  }
  else if (to_order == 'Color misfolded'){
    order=order(df$Colocalization[97:192])
  }
  else if (to_order == 'N Fiber'){
    order =order(df$ncells[193:288])
  }
  else if (to_order == 'Color Fiber'){
    order=order(df$Colocalization[193:288])
  }
  else{order=c(1:96)}  
  order = c(order,order+96,order+192,order+288)
  df=df[order,]
  return(df)
}

# cbind.fill <- function(...){
#   nm <- list(...) 
#   nm <- lapply(nm, as.matrix)
#   n <- max(sapply(nm, nrow)) 
#   do.call(cbind, lapply(nm, function (x) 
#     rbind(x, matrix(, n-nrow(x), ncol(x))))) 
# }

require(plyr) # requires plyr for rbind.fill()
cbind.fill <- function(...) {                                                                                                                                                       
  transpoted <- lapply(list(...),t)                                                                                                                                                 
  transpoted_dataframe <- lapply(transpoted, as.data.frame)                                                                                                                         
  return (data.frame(t(rbind.fill(transpoted_dataframe))))                                                                                                                          
} 
