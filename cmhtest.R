
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




cmhsim = function(nm,nf,pm,pf){
  MP = rbinom(1,nm, pm)             #samples value for p
  FP = rbinom(1,nf,pf) 
  MQ = nm - MP                      #calculates rest of reads as q
  FQ = nf - FP
  mable = c(MP,MQ)                  # puts values into a column
  fable = c(FP,FQ)
  cable = c(mable,fable)
  
  return(cable)
}


cmharray = function(nm,nf,pm,pf){
  first = cmhsim(nm,nf,pm,pf)
  second = cmhsim(nm,nf,pm,pf)
  combo =array(c(first,second),dim = c(2, 2, 2))
  return(mantelhaen.test(combo)$p.value)  
}

Ncmhsim = function(N,nm,nf,pm,pf){
  mylist =  c()
  for( i in 1:N){
    mylist = append(mylist, c(cmharray(nm,nf,pm,pf)))
  }
  return(mylist)
}


Ncmhsim(100,65,65,.75,.94)


CMHcurve = function(deltaP){
  
  list50f = c()
  list50t = c()
  
  pcutoff=c()
  pcutoff= seq(0,1,.0005)
  
  for (depth in seq(25,100,25)){
    
    false = Ncmhsim(500,(depth),(depth),.5,.5)
    true = Ncmhsim(500,(depth),(depth),.5,(.5+ deltaP))
    
    for (pcut in pcutoff){
      
      fp = sum(false <= pcut)
      tp = sum(true <= pcut)
      
      fpr = (fp / 500 )
      tpr = (tp / 500 )
      
      list50f = append(list50f, fpr)
      list50t = append(list50t, tpr)
      
      }
    }
    
    sample50 = data.frame(list50f, list50t)
    
  plot(sample50, xlab = "FPR", ylab = "TPR", pch = 21, bg= 'red', ylim = c(0,1), xlim= c(0,1), cex=.4, main = c('Delta P = ', deltaP))
  
    listf = list()
    listf = append(listf,sample50)
    
    return(sample50)

  }


plot(sample50, xlab = "FPR", ylab = "TPR", pch = 21, bg= 'red', ylim = c(0,1), xlim= c(0,1), cex=.4, main = c('Delta P = ', deltaP))
par(new=TRUE)
plot(sample100, xlab = "FPR", ylab = "TPR", pch = 21, bg= 'blue',ylim = c(0,1), xlim= c(0,1), cex=.4) #type='n' removes points
par(new=TRUE)
plot(sample150, xlab = "FPR", ylab = "TPR", pch = 21, bg= 'green', ylim = c(0,1), xlim= c(0,1), cex=.4)
par(new=TRUE)
plot(sample200, xlab = "FPR", ylab = "TPR", pch = 21, bg= 'orange',ylim = c(0,1), xlim= c(0,1), cex=.4)

legend(.6,.4, cex=.8, c("25 depth","50 depth","75 depth","100 depth"), lty=c(1,1), lwd=c(1,1),col=c("red","blue","green","orange"))

#points(sample100$list100f, predict(fit2), type = 'l', col='blue',lwd=1)

}}



