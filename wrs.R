
fishersim = function(nm,nf,pm,pf){
  MP = rbinom(1,nm, pm)             #samples value for p
  FP = rbinom(1,nf,pf) 
  MQ = nm - MP                      #calculates rest of reads as q
  FQ = nf - FP
  mable = c(MP,MQ)                  # puts values into a column
  fable = c(FP,FQ)
  cable = cbind(mable,fable)        #makes a 2x2 table
  return (fisher.test(cable))
}



NfisherSim = function(N,nm,nf,pm,pf){
  mylist = c() 
  for (i in 1:N){
    mylist = append(mylist,c(fishersim(nm,nf,pm,pf)$p.value)) #collects N pvalues in a list
  } 
  return(mylist)    #returns list of N pvalues
}

tprGet = function(target,true,false){
  
  pcut = 1
  fpr = 1 
  tpr = 1 
  
  while(fpr >= target){
    fpr = (sum(false < pcut) / 10000)
    tpr = (sum(true < pcut) / 10000)
    pcut = pcut / 1.01
  }
  return(tpr)
}



FPRlines = function(deltap){
  
  depth = seq(0,500,5)
  list1 = c()
  list2 = c()
  list3 = c()
  
  for (dep in depth){
    true = NfisherSim(10000,dep,dep,.5,(.5+deltap))
    false = NfisherSim(10000,dep,dep,.5,.5)
    
    entry1 = tprGet(.1,true,false)
    entry2 = tprGet(.2,true,false)
    entry3 = tprGet(.5,true,false)
    
    list1 = append(list1,entry1)
    list2 = append(list2,entry2)
    list3 = append(list3,entry3)
    
  }
  
  plot(depth,list1,ylim=c(0,1), ylab='TPR', xlab='Depth', main = c('Delta P =', deltap), col='red', cex = .5)
  par(new=TRUE)
  plot(depth,list2,ylim=c(0,1), ylab='TPR', xlab='Depth', main = c('Delta P =', deltap), col='blue', cex = .5)
  par(new=TRUE)
  plot(depth,list3,ylim=c(0,1), ylab='TPR', xlab='Depth', main = c('Delta P =', deltap), col='green', cex = .5)
  legend(310,.4, cex=.65, c(".1 FPR",".2 FPR",".5 FPR"), lty=c(1,1), lwd=c(1,1),col=c("red","blue","green"))
  
  newlist = list(list1,list2,list3)
  return(newlist)
  
}

FPRonly = function(deltap){
  
  depth = seq(0,500,5)
  list1 = c()
  list2 = c()
  list3 = c()
  
  for (dep in depth){
    true = NfisherSim(10000,dep,dep,.5,(.5+deltap))
    false = NfisherSim(10000,dep,dep,.5,.5)
    
    entry1 = tprGet(.1,true,false)
    entry2 = tprGet(.2,true,false)
    entry3 = tprGet(.5,true,false)
    
    list1 = append(list1,entry1)
    list2 = append(list2,entry2)
    list3 = append(list3,entry3)
    
  }
  
  
  newlist = list(list1,list2,list3)
  return(newlist)
  
}

doit = function(){
  wun= FPRonly(.15)
   two = FPRonly(.1)
  tree = FPRonly(.05)
  fur = FPRonly(.2)
  
  list  = c()
  list = append(list,tree)
  list = append(list,two)
  list = append(list,wun)
  list  = append(list,fur)
  
  return(list)
}

doit()

graphs = function(){
par(mfrow=c(2,2))
FPRlines(.05)
FPRlines(.1)
FPRlines(.15)
FPRlines(.2) 

}

x= FDRlines(.2)

