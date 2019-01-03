
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


false = NfisherSim(100,100,100,.5,.5)
true = NfisherSim(100,100,100,.5,.7)


tprGet = function(target,true,false){
  
  pcut = 1
  fpr = 1 
  tpr = 1 
  
  while(fpr >= target){
    fpr = (sum(false < pcut) / 100)
    tpr = (sum(true < pcut) / 100)
    pcut = pcut / 1.01
  }
  return(tpr)
}


tprGet(.05,true,false)


FPRlines = function(deltap){
  
  depth = seq(0,500,2)
  list1 = c()
  list2 = c()
  list3 = c()
  
  for (dep in depth){
    true = NfisherSim(100,dep,dep,.5,(.5+deltap))
    false = NfisherSim(100,dep,dep,.5,.5)
    
    entry1 = tprGet(.1,true,false)
    entry2 = tprGet(.2,true,false)
    entry3 = tprGet(.5,true,false)
    
    list1 = append(list1,entry1)
    list2 = append(list2,entry2)
    list3 = append(list3,entry3)
    
  }
  
  plot(depth,list1,ylim=c(0,1), ylab='TPR', xlab='Depth', main = c('Delta P =', deltap), col='red', cex = .7)
  par(new=TRUE)
  plot(depth,list2,ylim=c(0,1), ylab='TPR', xlab='Depth', main = c('Delta P =', deltap), col='blue', cex = .7)
  par(new=TRUE)
  plot(depth,list3,ylim=c(0,1), ylab='TPR', xlab='Depth', main = c('Delta P =', deltap), col='green', cex = .7)
  legend(300,.4, cex=.8, c(".1 FDR",".2 FDR",".4 FDR"), lty=c(1,1), lwd=c(1,1),col=c("red","blue","green"))
  
  
}




graphs = function(){
  list  = c()
  par(mfrow=c(2,2))
  
  q = FPRlines(.05)
  w = FPRlines(.1)
  e = FPRlines(.15)
  r = FPRlines(.2)
  
  list = append(list,q)
  list = append(list,w)
  list = append(list,e)
  list = append(list,r)
  
  return(list)
  }

graphs()

