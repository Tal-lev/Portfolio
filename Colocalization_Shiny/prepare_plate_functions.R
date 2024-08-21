#Loads csv results of colocalization plate. gathers the 4 strains into 1 column
prepare_plate = function (name,file,to_order='',onlypuncta,norm='No'){
    Imgurl = "https://media.githubusercontent.com/media/Tal-lev/Portfolio/main/Colocalization_Shiny/myImages/"
  #change localization folder
  if ('Yes' %in% onlypuncta){
      myImgResources2 <- paste0(Imgurl,name,"-sel-puncta-R_8-bit-",seq_len(384),".png")
      myImgResources <- paste0(Imgurl,name,"-sel-puncta_8-bit-",seq_len(384),".png")
    }
  else{
      myImgResources2 <- paste0(Imgurl,name,"-sel-R_8-bit-", seq_len(384), ".png")
      myImgResources <- paste0(Imgurl,name,"-sel_8-bit-",seq_len(384),".png")
    }
  df = read.csv(file = file)
  if ('RFPnorm' %in% norm){
    df[,5] = df[,5] / df[,2]
    df[,4] = df[,4] / df[,2]
    df[,3] = df[,3] / df[,2]
    df[,2] = df[,2] / df[,2]
  }
  ncells = gather(df,key = 'strain', 'ncells',15:18)
  YFP = gather(df,key = 'strain', 'YFPmean',19:22)
  RFP = gather(df,key = 'strain', 'RFPmean',23:26)
  df = gather(df,key = 'strain', 'Colocalization',2:5)
  df["image_url"] = myImgResources
  df= df[,-(11:22)]
  df["image_url2"] = myImgResources2
  df["ncells"] = ncells$ncells
  df["YFPmean"] = YFP$YFPmean
  df["RFPmean"] = RFP$RFPmean
  df["text"] = paste('n:',df$ncells,'\nYFP:',df$YFPmean,'\nRFP:',df$RFPmean)
  if (to_order == 'colocalization'){
    order = c(74,40,39,62,42,1,15,95,55,56,34,65,30,4,5,63,16,93,14,25,2,26,3,52,13,32,48,81,44,92,35,94,58,70,11,83,50,
            87,66,43,91,60,89,10,59,88,79,90,51,64,49,69,7,57,75,17,67,76,22,27,46,20,33,68,72,82,84,86,31,85,19,24,47,
            8,6,54,9,45,96,37,53,28,38,77,23,78,36,12,80,41,18,21,29,71,61,73)}
  else if (to_order == 'correlation'){
    order=c(34,65,15,42,62,39,40,1,56,55,95,33,2,26,3,52,13,32,4,5,30,14,16,63,25,93,96,28,38,9,37,53,61,23,45,41,77,78,80,12,36,67,17,
            76,29,20,46,19,24,47,22,27,54,6,18,8,21,7,75,49,69,66,87,74,51,64,82,86,72,31,85,71,73,68,84,83,88,90,57,79,35,94,50,92,44,
            70,58,89,43,60,10,59,11,91,48,81)}
  else if (to_order == 'RFPmean'){
    order=c(40,48,87,13,66,37,34,89,53,49,59,88,15,90,29,64,69,79,46,67,84,3,21,31,25,85,86,76,9,54,36,96,8,47,80,12,24,27,78,77,45,6,22,
            33,68,18,19,72,41,38,63,71,16,23,61,14,20,73,32,2,5,58,70,26,82,10,17,93,51,42,28,57,83,7,52,4,11,94,35,92,60,62,75,65,91,50,
            30,1,74,55,44,56,39,95,43,81)}
  else if (to_order == 'N Wt'){
    order =order(df$ncells[1:96])
  }
  else if (to_order == 'Color Wt'){
    order=order(df$Colocalization[1:96])
  }
  else if (to_order == 'N Fiber'){
    order =order(df$ncells[97:192])
  }
  else if (to_order == 'Color Fiber'){
    order=order(df$Colocalization[97:192])
  }
  else if (to_order == 'N misfolded'){
    order =order(df$ncells[193:288])
  }
  else if (to_order == 'Color misfolded'){
    order=order(df$Colocalization[193:288])
  }
  else{order=c(1:96)}  
  order = c(order,order+96,order+192,order+288)
  df=df[order,]
  return(df)
}
