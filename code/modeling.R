
library(MASS)
library(caret)
library(fGarch)
library(fitdistrplus)
library(pracma)
library(BBmisc)

setwd("/media/sf_jonghyeon3/extension_AD/code")

# from = "/media/sf_jonghyeon3/extension_AD/data/stream/org"
from = "/media/sf_jonghyeon3/extension_AD/data/preprocessed_binet"
to = "/media/sf_jonghyeon3/extension_AD/data/stream/preprocess"

setwd(from)
fn<-list.files(getwd())
fn
input = data.frame(read.csv(fn[3], header=T))
normal= input[which(input$anomaly_type =="normal"),]
anomaly= input[which(input$anomaly_type !="normal"),]
normal_seq = aggregate(normal$Activity, by=list(normal$Case), FUN=paste0)
anomaly_seq = aggregate(anomaly$Activity, by=list(anomaly$Case), FUN=paste0)
delete_case= anomaly_seq[which(is.element(anomaly_seq$x , normal_seq$x)),'Group.1']
input = input[which(!is.element(input$Case, delete_case)),]
input$Event = 1:nrow(input)
input$Event = as.factor(input$Event)

####

fun_leverage = function(x){
  A<- ginv(t(x)%*%x)
  H_part1<- x%*%A  
  h_diag <- colSums(t(H_part1)*t(x))
  h_diag
}

####
streaming_score = function(input, Min = 100, Max = 0, until=0){
  preprocess_start <- Sys.time()
  
  pre<-input
  pre= pre[ with(pre, order(Case,timestamp)),]
  one= rep(1, nrow(pre))
  pre[,'start'] = ave(one, by= pre$Case, FUN= cumsum) -1
  pre[which(pre$start !=1),'start'] =0
  pre= pre[ with(pre, order(timestamp)),]
  pre[,'Event'] = as.factor(1:nrow(pre))
  
  pre[,'num_case'] = cumsum(pre$start)
  pre[,'leverage'] = rep(-1, nrow(pre))
  pre[,'w-leverage'] = rep(-1, nrow(pre))
  pre[,'t1'] = rep(0, nrow(pre))
  pre[,'w-t1'] = rep(0, nrow(pre))
  pre[,'t2'] = rep(0, nrow(pre))
  pre[,'w-t2'] = rep(0, nrow(pre))
  pre[,'t3'] = rep(0, nrow(pre))
  pre[,'w-t3'] = rep(0, nrow(pre))
  pre[,'tn']= rep(0, nrow(pre))
  pre[,'w-tn']= rep(0, nrow(pre))
  # pre[,'tg']= rep(0, nrow(pre))
  # pre[,'w-tg']= rep(0, nrow(pre))
  pre[,'time'] = rep(0, nrow(pre))
  event_num = nrow(pre)
  case_num= length(unique(pre$Case))
  # start_index  = which(pre$num_case == Min +1)[1]
  start_index  = which(pre$num_case == Min +1)[1]
  last_index = nrow(pre)
  print(paste("Start to calculate leverage score of ", start_index ,"-th event (total ",event_num," events)" ,sep=''))
  
  leverage_start <- Sys.time()
  
  pre2 = pre[1:start_index,]
  cur_len = length(unique(pre2$Case))
  if(Max==0){
    Max=case_num
  }else{
    if(Max < cur_len){
      delete_number = cur_len - Max
      pre2 = pre2[which(!is.element(pre2$Case,unique(pre2[which(pre2$num_case == delete_number),'Case']))),]
    }
  }
  
  data<- pre2[,c("Case","Activity")]  
  names(data)[1:2] <- c("ID", "ActivityID") 
  object_case = pre2$Case[nrow(pre2)]
  object_event = pre2$Event[nrow(pre2)]
  data$ID <- as.factor(data$ID)
  data$ActivityID <- as.factor(data$ActivityID)
  
  if(length(levels(data$ActivityID))>1){
    a<- model.matrix(~ActivityID, data = data)
    A<- as.numeric(data[,2])
    A[which(A!=1)] <- 0
    a<- cbind(ActivityID1 = A, a[,-1])
    onehot<- as.data.frame(a)
  }else{
    A<- as.numeric(data[,2])
    A[which(A!=1)] <- 0
    a<- cbind(ActivityID1 = A)
    onehot<- as.data.frame(a)
  }
  
  data1 <- onehot
  newdat <- cbind(data[,1], data1)
  newdat[,1] <- as.factor(newdat[,1])
  n<- length(levels((newdat[,1])))   # the number of cases
  m<- summary((newdat[,1]))[1]      # maximum trace length
  l<- unique((newdat[,1]))
  max<- m*(ncol(newdat)-1)
  c=unique(pre2[,c("Case","anomaly_type")]) #CHANGE
  label = as.character(c[,2])
  
  newdat2<- matrix(NA, nrow=n , ncol=max)
  for(j in 1:n){
    save2 <- as.vector(t(newdat[which(newdat[,1]==l[j]),-1]))  
    newdat2[j,1:length(save2)] <- save2
  }
  newdat2[which(is.na(newdat2))] <- 0 # zero-padding
  l_save=l
  newdat2_save= newdat2
  newdat3 = data.frame(cbind(Case=l, label= label, newdat2))
  
  
  #Caculate leverage
  x2= newdat3[,-(1:2)]
  x= as.matrix(sapply(x2, as.numeric))  
  h_diag <- fun_leverage(x)
  
  #Calculate weighted leverage
  length <- apply(x,1,sum)
  z_norm <- (length- mean(length))/sd(length)
  
  sigmoid_leng <- 1/(1+exp(-z_norm))
  if(-2.2822+max(length)^0.3422 <0 | length(unique(length))==1 ){
    h_diag2 = h_diag}else{h_diag2 <-h_diag*(1-sigmoid_leng)^(-2.2822+max(length)^0.3422) } #weighted leverage
  h_diag2 = h_diag2*sum(h_diag)/sum(h_diag2)
  
  pre[start_index, 'leverage'] = h_diag[which(l==object_case)]
  pre[start_index, 'w-leverage'] = h_diag2[which(l==object_case)]
  
  leverage_end <- Sys.time()
  
  print(paste("Anomaly score of", start_index ,"-th event = ", round( h_diag2[which(l==object_case)],5), " (CaseID=",object_case,")"  ,sep=''))
  pre[start_index, 'time'] =   (leverage_end-leverage_start)
  pre[start_index, 'tn'] = (h_diag[which(l==object_case)] > (mean(h_diag)+sd(h_diag)))
  pre[start_index, 'w-tn'] =  (h_diag2[which(l==object_case)] > (mean(h_diag2)+sd(h_diag2)))
  
  if(until==0){
    until = last_index
  }else{
    until= start_index+until
  }
  for(i in (start_index+1):until){ # last_index
    
    print(paste("Start to calculate leverage score of ", i ,"-th event (total ",event_num," events)" ,sep=''))
    leverage_start <- Sys.time()
    
    pre2 = rbind(pre2, pre[i,])
    cur_len = length(unique(pre2$Case))
    
    data<- pre2[,c("Case","Activity")]  
    names(data)[1:2] <- c("ID", "ActivityID") 
    object_case = pre2$Case[nrow(pre2)]
    object_event = pre2$Event[nrow(pre2)]
    data$ID <- as.factor(data$ID)
    data$ActivityID <- as.factor(data$ActivityID)
    
    if(length(levels(data$ActivityID))>1){
      a<- model.matrix(~ActivityID, data = data)
      A<- as.numeric(data[,2])
      A[which(A!=1)] <- 0
      a<- cbind(ActivityID1 = A, a[,-1])
      onehot<- as.data.frame(a)
    }else{
      A<- as.numeric(data[,2])
      A[which(A!=1)] <- 0
      a<- cbind(ActivityID1 = A)
      onehot<- as.data.frame(a)
    }
    
    data1 <- onehot
    newdat <- cbind(data[,1], data1)
    newdat[,1] <- as.factor(newdat[,1])
    n<- length(levels((newdat[,1])))   # the number of cases
    m<- summary((newdat[,1]))[1]      # maximum trace length

    
    l<-unique((newdat[,1]))
    max<- m*(ncol(newdat)-1)
    c=unique(pre2[,c("Case","anomaly_type")])
    label = as.character(c[,2])
    newdat2<- matrix(NA, nrow=n , ncol=max)
    newdat2[1:nrow(newdat2_save), 1:ncol(newdat2_save)] = newdat2_save
    if(is.element(object_case, l_save) ){
      j = which(l == object_case)
      save2 <- as.vector(t(newdat[which(newdat[,1]==l[j]),-1]))  
      newdat2[j,1:length(save2)] <- save2
      loc =j
    }else{
      j = which(l == object_case)
      save2 <- as.vector(t(newdat[which(newdat[,1]==l[j]),-1]))  
      newdat2[n,1:length(save2)] <- save2
      loc =n
      
      if(cur_len > Max ){
        newdat2 = newdat2[-(1:(cur_len - Max)),]
        loc = loc - (cur_len - Max)
        l= l[-(1:(cur_len - Max))]
        label= label[-(1:(cur_len - Max))]
      }   
    }
    
    newdat2[which(is.na(newdat2))] <- 0 # zero-padding
    
    l_save = l
    newdat2_save= newdat2
    newdat3 <-data.frame(cbind(Case=l, label= label, newdat2))
    
    #Caculate leverage
    x2= newdat3[,-(1:2)]
    x= as.matrix(sapply(x2, as.numeric))  
    h_diag <- fun_leverage(x)
    
    #Calculate weighted leverage
    length <- apply(x,1,sum)
    z_norm <- (length- mean(length))/sd(length)
    
    sigmoid_leng <- 1/(1+exp(-z_norm))
    if(-2.2822+max(length)^0.3422 <0 | length(unique(length))==1 ){
      h_diag2 = h_diag}else{h_diag2 <-h_diag*(1-sigmoid_leng)^(-2.2822+max(length)^0.3422) } #weighted leverage
    h_diag2 = h_diag2*sum(h_diag)/sum(h_diag2)
    
    pre[i, 'leverage'] = h_diag[loc]
    pre[i, 'w-leverage'] = h_diag2[loc]
    
    
    leverage_end <- Sys.time()
    

    print(paste("Anomaly score of", i ,"-th event = ", round( h_diag2[loc],5), " (CaseID=",object_case,")" ,sep=''))
    pre[i, 'time'] =   (leverage_end-leverage_start)
    pre[i, 'tn'] = (h_diag[loc] > (mean(h_diag)+sd(h_diag)))
    pre[i, 'w-tn'] =  (h_diag2[loc] > (mean(h_diag2)+sd(h_diag2)))
    
    # fit= try(tryCatch(fitdist(h_diag, distr="gamma", method='mle'),  error = 0))
    # if(is.error(fit)){
    #   fit.gamma <- fitdist(h_diag, distr="gamma", method='mme')
    # }else{fit.gamma <- fitdist(h_diag, distr="gamma", method='mle')}
    # 
    # fit2= try(tryCatch(fitdist(h_diag2, distr="gamma", method='mle'),  error = 0))
    # if(is.error(fit2)){
    #   fit2.gamma <- fitdist(h_diag2, distr="gamma", method='mme')
    # }else{fit2.gamma <- fitdist(h_diag2, distr="gamma", method='mle')}
    # 
    # pre[which(pre$Case == object_case & pre$Event == object_event), 'tg'] = (pgamma(h_diag[loc],
    #                                                                                 shape= fit.gamma$estimate[1],
    #                                                                                 rate= fit.gamma$estimate[2],
    #                                                                                 lower.tail = F) <0.10)
    # pre[which(pre$Case == object_case & pre$Event == object_event), 'w-tg'] = (pgamma(h_diag2[loc],
    #                                                                                 shape= fit2.gamma$estimate[1],
    #                                                                                 rate= fit2.gamma$estimate[2],
    #                                                                                 lower.tail = F) <0.10)
    # if(i == start_index + N3){
    #   write.csv(pre,"pre_0.3.csv", row.names = F)
    # }
    # if(i == start_index + N6){
    #   write.csv(pre,"pre_0.6.csv", row.names = F)
    # }
    
  }
  
  
  preprocess_end <- Sys.time()
  print(preprocess_end - preprocess_start)
  return(pre)
  # write.csv(pre,paste("t1",".csv",sep=''), row.names = F)
}



output = streaming_score(input,Min=1000,Max=3000, until =0)

setwd("/media/sf_jonghyeon3/extension_AD/data/result")
write.csv(output, "small_3000.csv",row.names = FALSE)

output = streaming_score(input,Min=1000,Max=10000, until =0)

setwd("/media/sf_jonghyeon3/extension_AD/data/result")
write.csv(output, "small_10000.csv",row.names = FALSE)
# 
# output = streaming_score(input, Min=100,Max=1000, until  = 0)
# write.csv(output, "helpdesk_1000.csv", row.names= FALSE)
# 
# output = streaming_score(input, Min=100, Max=2000, until = 0)
# write.csv(output, "helpdesk_2000.csv", row.names = FALSE)

}

write.csv(output, "bpic17-1_2000.csv",row.names = FALSE)




# evaluation
pre[, 't1'] = as.factor(as.numeric(pre[, 'leverage'] > 0.1))
pre[, 'w-t1'] = as.factor(as.numeric(pre[, 'w-leverage'] > 0.1))
pre[, 't2'] = as.factor(as.numeric(pre[, 'leverage'] > 0.15))
pre[, 'w-t2'] =as.factor(as.numeric(pre[, 'w-leverage'] > 0.15))
pre[, 't3'] = as.factor(as.numeric(pre[, 'leverage'] > 0.2))
pre[, 'w-t3'] = as.factor(as.numeric(pre[, 'w-leverage'] > 0.2))

part = pre[which(pre$leverage> -1), ]
p1 = part$time[ with(part, order(Timestamp))]
p1 = p1[1:(length(p1)-1)]
plot(p1, type='l', xlab='event index', ylab= 'time', main ='First 13,557 events calculated' )
axis(side = 2, at=seq(0,45,5))

summary(part$leverage)
summary(part$leverage[which(part$label==0)])
summary(part$leverage[which(part$label==1)])

summary(part$'w-leverage')
summary(part$'w-leverage'[which(part$label==0)])
summary(part$'w-leverage'[which(part$label==1)])


mean(part$leverage)
sd(part$leverage)
hist(part$leverage)
hist(part$leverage[which(part$leverage<0.2)])
summary(part$leverage[which(part$leverage<0.2)])
mean(part$leverage[which(part$leverage<0.2)])
sd(part$leverage[which(part$leverage<0.2)])
mean(part$leverage[which(part$leverage<0.2)]) + 2*sd(part$leverage[which(part$leverage<0.2)])

## Starting point : 이 방식은 앞에서 먼저 잘못 디텍팅 이벤트를 고려 못함.
act_start = part$label[which(part$start==1)]
act_start <- as.factor(act_start)

cm1 <- confusionMatrix(part$t1[which(part$start==1)], act_start,positive = '1')
cm1 <- as.vector(cm1[4])[[1]]
cm1[5]
cm1[6]
cm1[7]

## general evaluation
act = part$label
act <- as.factor(act)

cm1 <- confusionMatrix(part$t1, act,positive = '1')
cm1 <- as.vector(cm1[4])[[1]]
cm1[5]
cm1[6]
cm1[7]

cm2 <- confusionMatrix(part$t2, act,positive = '1')
cm2 <- as.vector(cm2[4])[[1]]
cm2[5]
cm2[6]
cm2[7]

cm3 <- confusionMatrix(part$t3, act,positive = '1')
cm3 <- as.vector(cm3[4])[[1]]
cm3[5]
cm3[6]
cm3[7]

cm4 <- confusionMatrix(part$'w-t1', act,positive = '1')
cm4 <- as.vector(cm4[4])[[1]]
cm4[5]
cm4[6]
cm4[7]

cm5 <- confusionMatrix(part$'w-t2', act,positive = '1')
cm5 <- as.vector(cm5[4])[[1]]
cm5[5]
cm5[6]
cm5[7]

cm6 <- confusionMatrix(part$'w-t3', act,positive = '1')
cm6 <- as.vector(cm6[4])[[1]]
cm6[5]
cm6[6]
cm6[7]




