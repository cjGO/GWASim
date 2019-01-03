
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
  
  depth = seq(0,500,2)
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
  
  #plot(depth,list1,ylim=c(0,1), ylab='TPR', xlab='Depth', main = c('Delta P =', deltap), col='red', cex = .5)
  #par(new=TRUE)
  #plot(depth,list2,ylim=c(0,1), ylab='TPR', xlab='Depth', main = c('Delta P =', deltap), col='blue', cex = .5)
  #par(new=TRUE)
  #plot(depth,list3,ylim=c(0,1), ylab='TPR', xlab='Depth', main = c('Delta P =', deltap), col='green', cex = .5)
  #legend(310,.4, cex=.65, c(".1 FPR",".2 FPR",".5 FPR"), lty=c(1,1), lwd=c(1,1),col=c("red","blue","green"))
  listfin = c()
  listfin = append(listfin,list1)
  listfin = append(listfin,list2)
  listfin = append(listfin,list3)
  
  return(listfin) 
  
}
#q used for .05 alraedy
#ww .1
#ee
#ff

graphs = function(){
par(mfrow=c(2,2))
FPRlines(.05)
FPRlines(.1)
FPRlines(.15)
FPRlines(.2) 

}

graphs()


grapit = function(list1,list2,list3, deltap){

depth = seq(0,500,2)
  
plot(depth,list1,ylim=c(0,1), ylab='TPR', xlab='Depth', main = c('Delta P =', deltap), col='red', cex = .5)
par(new=TRUE)
plot(depth,list2,ylim=c(0,1), ylab='TPR', xlab='Depth', main = c('Delta P =', deltap), col='blue', cex = .5)
par(new=TRUE)
plot(depth,list3,ylim=c(0,1), ylab='TPR', xlab='Depth', main = c('Delta P =', deltap), col='green', cex = .5)
legend(260,.2,ncol=3, cex=.65, c(".1 FPR",".2 FPR",".5 FPR"),y.intersp = 0.5, lty=c(1,1), lwd=c(2,2),col=c("red","blue","green"))

}



par(mfrow=c(2,2))
grapit(qq1,qq2,qq3,.05)
grapit(ww1,ww2,ww3,.10)
grapit(ee1,ee2,ee3,.15)
grapit(ff1,ff2,ff3,.20)
mtext('Depth vs TPR for set DeltaP/FPR',outer=TRUE,cex=1.5,side = 3,line=-2)





