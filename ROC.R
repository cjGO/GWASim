


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



dimSim = function(){ #compares the regression lines of five values of delta P across different depths
list1=c()
list2=c()
list3=c()
list4=c()
list5 = c()
for (delta in seq(.6,1,.1)){         # concantonates the median p-value for each depth at each delta p value
for (depth in seq(20,200,10)){
test = median(NfisherSim(30,(depth/2),(depth/2),.5,delta))
if (delta == .6){
list1 = append(list1,test)}
if (delta == .7){
list2 = append(list2,test)}
if (delta == .8){
list3 = append(list3,test)}
if (delta == .9){
list4 = append(list4,test)}
if (delta == 1){
list5 = append(list5,test)}
}
}

xaxis = seq(20,200,10)             # enters the lists into a data frame
sample1 = data.frame(xaxis, list1)
sample2 = data.frame(xaxis, list2)
sample3 = data.frame(xaxis, list3)
sample4 = data.frame(xaxis, list4)
sample5 = data.frame(xaxis, list5)

fit1 = lm(sample1$list1 ~ poly(sample1$xaxis,4,raw=TRUE)) # polynomial regression for each delta p
fit2 = lm(sample2$list2 ~ poly(sample2$xaxis,4,raw=TRUE))
fit3 = lm(sample3$list3 ~ poly(sample3$xaxis,4,raw=TRUE))
fit4 = lm(sample4$list4 ~ poly(sample4$xaxis,4,raw=TRUE))
fit5 = lm(sample5$list5 ~ poly(sample5$xaxis,4,raw=TRUE))

plot(xaxis,list1, xlab = "Depth", ylab = "median p", pch = 21, bg= 'red', ylim = c(0,.7), type='n', main = "Depth Simulator for Various delta P's")
par(new=TRUE)
plot(xaxis,list2, xlab = "Depth", ylab = "median p", pch = 21, bg= 'blue', ylim = c(0,.7), cex=.5, type='n')
par(new=TRUE)
plot(xaxis,list3, xlab = "Depth", ylab = "median p", pch = 21, bg= 'green', ylim = c(0,.7), cex=.5, type='n')
par(new=TRUE)
plot(xaxis,list4, xlab = "Depth", ylab = "median p", pch = 21, bg= 'orange', ylim = c(0,.7), cex=.5, type='n')
par(new=TRUE)
plot(xaxis,list5, xlab = "Depth", ylab = "median p", pch = 21, bg= 'black', ylim = c(0,.7), cex=.5, type='n')

points(sample1$xaxis, predict(fit1), type = 'l', col='red',lwd=1)
points(sample2$xaxis, predict(fit2), type = 'l', col='blue',lwd=1)
points(sample2$xaxis, predict(fit3), type = 'l', col='green',lwd=1)
points(sample2$xaxis, predict(fit4), type = 'l', col='orange',lwd=1)
points(sample2$xaxis, predict(fit5), type = 'l', col='black',lwd=1)

legend(150,.6, c(".1 delta p",".2 delta p",".3 delta p",".4 delta p",".5 delta p"), lty=c(1,1), lwd=c(1,1),col=c("red","blue","green","orange","black"))


return(list1)
}

dimSim()





ROCcurve = function(deltaP){

list50f = c()
list50t = c()
list100f = c()
list100t = c()
list150f = c()
list150t = c()
list200f = c()
list200t = c()

pcutoff =  10^-seq(0,7,.05)

for (depth in seq(50,200,50)){
false = NfisherSim(1000,(depth/2),(depth/2),.5,.5)
true = NfisherSim(1000,(depth/2),(depth/2),.5,(.5+ deltaP)) # delta  p  = .2
for (pcut in pcutoff){

fp = sum(false < pcut)
tp = sum(true < pcut)

fpr = (fp / 1000 )
tpr = (tp / 1000 )

if (depth == 50){
list50f = append(list50f, fpr)
list50t = append(list50t, tpr)
}
if (depth == 100){
list100f = append(list100f, fpr)
list100t = append(list100t, tpr)
}
if (depth == 150){
list150f = append(list150f, fpr)
list150t = append(list150t, tpr)
}
if (depth == 200){
list200f = append(list200f, fpr)
list200t = append(list200t, tpr)
}
}


sample50 = data.frame(list50f, list50t)
sample100 = data.frame(list100f, list100t)
sample150 = data.frame(list150f, list150t)
sample200 = data.frame(list200f, list200t)

}







fit1 = lm(sample50$list50t ~ poly(sample50$listf,4,raw=TRUE)) # polynomial regression for each delta p
fit2 = lm(sample100$list100t ~ poly(sample100$list100f,4,raw=TRUE))
fit3 = lm(sample150$list150t ~ poly(sample150$list150f,4,raw=TRUE))
fit4 = lm(sample200$list200t ~ poly(sample200$list200f,4,raw=TRUE))

plot(xaxis,list1, xlab = "FPR", ylab = "TPR", pch = 21, bg= 'red', ylim = c(0,.7), type='n', main = "Depth Simulator for Various delta P's")
par(new=TRUE)
plot(xaxis,list2, xlab = "FPR", ylab = "TPR", pch = 21, bg= 'blue', ylim = c(0,.7), cex=.5, type='n')
par(new=TRUE)
plot(xaxis,list3, xlab = "FPR", ylab = "TPR", pch = 21, bg= 'green', ylim = c(0,.7), cex=.5, type='n')
par(new=TRUE)
plot(xaxis,list4, xlab = "FPR", ylab = "TPR", pch = 21, bg= 'orange', ylim = c(0,.7), cex=.5, type='n')

}

plot(fpr,tpr)

}









ROCgraph = function{

ROCcurve(.1)
ROCcurve(.2)
ROCcurve(.3)
ROCcurve(.4)

}





