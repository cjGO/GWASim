
TPRcurve = function(deltaP,pcut){
List = c()

for (depth in seq(5,1000,5)){
false = NfisherSim(1000,(depth),(depth),.5,.5)
true = NfisherSim(1000,(depth),(depth),.5,(.5+ deltaP))

fp = sum(false < pcut)
tp = sum(true < pcut)

fpr = (fp / 1000 )
tpr = (tp / 1000 )

List = append(List,tpr)
}

return(List)

}


x = TPRcurve(.1,.05)
x = TPRgraph(.1)


TPRgraph = function(deltaP){
line1 = TPRcurve(deltaP,.05)
line2 = TPRcurve(deltaP,(.05/1000))
xaxis = seq(5,1000,10)
plot(xaxis,line1, col = 'red', type='l', xlim=c(0,1000), ylim=c(0,1), main = c('Delta P =',deltaP), ylab='TPR', xlab='Depth')
lines(xaxis,line2, col = 'blue',lty=2)
legend(700,.2, cex =1, c(".05", ".05/1000"), lwd =c(1,1), col=c('red','blue'), lty=1:2, title='alpha level')
}



plot(sample50, xlab = "FPR", ylab = "TPR", pch = 21, bg= 'red', ylim = c(0,1), xlim= c(0,1), cex=.8, main = c('Delta P = ', deltaP))
par(new=TRUE)
plot(sample100, xlab = "FPR", ylab = "TPR", pch = 21, bg= 'blue',ylim = c(0,1), xlim= c(0,1), cex=.8) #type='n' removes points
par(new=TRUE)
plot(sample150, xlab = "FPR", ylab = "TPR", pch = 21, bg= 'green', ylim = c(0,1), xlim= c(0,1), cex=.8)
par(new=TRUE)
plot(sample200, xlab = "FPR", ylab = "TPR", pch = 21, bg= 'orange',ylim = c(0,1), xlim= c(0,1), cex=.8)

legend(.6,.4, cex=.8, c("25 depth","50 depth","75 depth","100 depth"), lty=1:2, lwd=c(1,1),col=c("red","blue","green","orange"))

#points(sample100$list100f, predict(fit2), type = 'l', col='blue',lwd=1)






TPRfin  = function(){
  attach(mtcars)
  par(mfrow=c(2,2))
  TPRgraph(.05)
  TPRgraph(.10)
  TPRgraph(.15)
  TPRgraph(.20)
}


TPRfin()










