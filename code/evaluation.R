
library(caret)

setwd("~/extension_AD/data/result")


fn<-list.files(getwd())
fn
fn = fn[c(1:4, 12:21, 23:27)]
fn
fn = fn[-c(6,11,16)]
fn

for(j in 1:length(fn)){
  input = data.frame(read.csv(fn[j], header=T))
  input = input[order(input$timestamp),]
  start = min(which(input$leverage >=0))
  print(paste("data=", fn[j]))
  print(paste("nrow=", nrow(input)))
  start
  print(paste("start=", start))
  
  
  
  input$Fs_T1_c = rep(-1, nrow(input))
  input$Rec_T1_c = rep(-1, nrow(input))
  input$Pre_T1_c = rep(-1, nrow(input))
  input$Fs_T2_c = rep(-1, nrow(input))
  input$Rec_T2_c = rep(-1, nrow(input))
  input$Pre_T2_c = rep(-1, nrow(input))
  input$Fs_T3_c = rep(-1, nrow(input))
  input$Rec_T3_c = rep(-1, nrow(input))
  input$Pre_T3_c = rep(-1, nrow(input))
  input$Fs_T4_c = rep(-1, nrow(input))
  input$Rec_T4_c = rep(-1, nrow(input))
  input$Pre_T4_c = rep(-1, nrow(input))
  input$Fs_Tn_c = rep(-1, nrow(input))
  input$Rec_Tn_c = rep(-1, nrow(input))
  input$Pre_Tn_c = rep(-1, nrow(input))
  input$Fs_wT1_c = rep(-1, nrow(input))
  input$Rec_wT1_c = rep(-1, nrow(input))
  input$Pre_wT1_c = rep(-1, nrow(input))
  input$Fs_wT2_c = rep(-1, nrow(input))
  input$Rec_wT2_c = rep(-1, nrow(input))
  input$Pre_wT2_c = rep(-1, nrow(input))
  input$Fs_wT3_c = rep(-1, nrow(input))
  input$Rec_wT3_c = rep(-1, nrow(input))
  input$Pre_wT3_c = rep(-1, nrow(input))
  input$Fs_wT4_c = rep(-1, nrow(input))
  input$Rec_wT4_c = rep(-1, nrow(input))
  input$Pre_wT4_c = rep(-1, nrow(input))
  input$Fs_wTn_c = rep(-1, nrow(input))
  input$Rec_wTn_c = rep(-1, nrow(input))
  input$Pre_wTn_c = rep(-1, nrow(input))

  
  act = rep(0, nrow(input))
  act[which(input$anomaly_type!="normal")] = 1
  act= as.factor(act)
  
  pred_T1 = as.factor(as.numeric(input$leverage >= 0.05))
  pred_T2 = as.factor(as.numeric(input$leverage >= 0.1))
  pred_T3 = as.factor(as.numeric(input$leverage >= 0.15))
  pred_T4 = as.factor(as.numeric(input$leverage >= 0.2))
  pred_Tn = as.factor(input$tn)
  
  pred_wT1 = as.factor(as.numeric(input$w.leverage >= 0.05))
  pred_wT2 = as.factor(as.numeric(input$w.leverage >= 0.1))  
  pred_wT3 =  as.factor(as.numeric(input$w.leverage >= 0.15))
  pred_wT4 = as.factor(as.numeric(input$leverage >= 0.2))
  pred_wTn = as.factor(input$w.tn)
  MM= ceiling((nrow(input)-start)/500)
  for(i in 1:MM){
    print(paste(j,"th",i, "/", MM))

    #case level
    loc1 = input$Case[1:min(start+500*i-1, nrow(input))]
    loc2 = aggregate(1:length(loc1), by=list(loc1), FUN= max)[,2]
    loc2 = loc2[which(loc2>=start+500*(i-1))]
    cm1 <- confusionMatrix(pred_T1[loc2], act[loc2],positive = '1')
    cm1 <- as.vector(cm1[4])[[1]]
    cm2 <- confusionMatrix(pred_T2[loc2], act[loc2],positive = '1')
    cm2 <- as.vector(cm2[4])[[1]]
    cm3 <- confusionMatrix(pred_T3[loc2], act[loc2],positive = '1')
    cm3 <- as.vector(cm3[4])[[1]]
    cm3.2 <- confusionMatrix(pred_T4[loc2], act[loc2],positive = '1')
    cm3.2 <- as.vector(cm3.2[4])[[1]]
    cm4 <- confusionMatrix(pred_Tn[loc2], act[loc2],positive = '1')
    cm4 <- as.vector(cm4[4])[[1]]
    cm5 <- confusionMatrix(pred_wT1[loc2], act[loc2],positive = '1')
    cm5 <- as.vector(cm5[4])[[1]]
    cm6 <- confusionMatrix(pred_wT2[loc2], act[loc2],positive = '1')
    cm6 <- as.vector(cm6[4])[[1]]
    cm7 <- confusionMatrix(pred_wT3[loc2], act[loc2],positive = '1')
    cm7 <- as.vector(cm7[4])[[1]]
    cm7.2 <- confusionMatrix(pred_wT4[loc2], act[loc2],positive = '1')
    cm7.2 <- as.vector(cm7.2[4])[[1]]
    cm8 <- confusionMatrix(pred_wTn[loc2], act[loc2],positive = '1')
    cm8 <- as.vector(cm8[4])[[1]]
    
    input$Pre_T1_c[max(loc2)]  <- cm1[5]
    input$Rec_T1_c[max(loc2)]  <- cm1[6]
    input$Fs_T1_c[max(loc2)] <- cm1[7]
    input$Pre_T2_c[max(loc2)]  <- cm2[5]
    input$Rec_T2_c[max(loc2)]  <- cm2[6]
    input$Fs_T2_c[max(loc2)] <- cm2[7]
    input$Pre_T3_c[max(loc2)]  <- cm3[5]
    input$Rec_T3_c[max(loc2)]  <- cm3[6]
    input$Fs_T3_c[max(loc2)] <- cm3[7]
    input$Pre_T4_c[max(loc2)]  <- cm3.2[5]
    input$Rec_T4_c[max(loc2)]  <- cm3.2[6]
    input$Fs_T4_c[max(loc2)] <- cm3.2[7]
    input$Pre_Tn_c[max(loc2)]  <- cm4[5]
    input$Rec_Tn_c[max(loc2)]  <- cm4[6]
    input$Fs_Tn_c[max(loc2)] <- cm4[7]
    input$Pre_wT1_c[max(loc2)]  <- cm5[5]
    input$Rec_wT1_c[max(loc2)]  <- cm5[6]
    input$Fs_wT1_c[max(loc2)] <- cm5[7]
    input$Pre_wT2_c[max(loc2)]  <- cm6[5]
    input$Rec_wT2_c[max(loc2)]  <- cm6[6]
    input$Fs_wT2_c[max(loc2)] <- cm6[7]
    input$Pre_wT3_c[max(loc2)]  <- cm7[5]
    input$Rec_wT3_c[max(loc2)]  <- cm7[6]
    input$Fs_wT3_c[max(loc2)] <- cm7[7]
    input$Pre_wT4_c[max(loc2)]  <- cm7.2[5]
    input$Rec_wT4_c[max(loc2)]  <- cm7.2[6]
    input$Fs_wT4_c[max(loc2)] <- cm7.2[7]
    input$Pre_wTn_c[max(loc2)]  <- cm8[5]
    input$Rec_wTn_c[max(loc2)]  <- cm8[6]
    input$Fs_wTn_c[max(loc2)] <- cm8[7]
    
  }
  
  write.csv(input, paste(strsplit(fn[j],".csv")[[1]],"_eval.csv", sep=""), row.names = F)
}

