


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




ROCcurve = function(deltaP){

list50f = c()
list50t = c()
list100f = c()
list100t = c()
list150f = c()
list150t = c()
list200f = c()
list200t = c()

pcutoff= seq(0,1,.00005)
pcutoff = sort(pcutoff)

for (depth in seq(25,100,25)){
false = NfisherSim(10000,(depth),(depth),.5,.5)
true = NfisherSim(10000,(depth),(depth),.5,(.5+ deltaP))
for (pcut in pcutoff){

fp = sum(false < pcut)
tp = sum(true < pcut)

fpr = (fp / 10000 )
tpr = (tp / 10000 )

if (depth == 25){
list50f = append(list50f, fpr)
list50t = append(list50t, tpr)
}
if (depth == 50){
list100f = append(list100f, fpr)
list100t = append(list100t, tpr)
}
if (depth == 75){
list150f = append(list150f, fpr)
list150t = append(list150t, tpr)
}
if (depth == 100){
list200f = append(list200f, fpr)
list200t = append(list200t, tpr)
}
}

sample50 = data.frame(list50f, list50t)
sample100 = data.frame(list100f, list100t)
sample150 = data.frame(list150f, list150t)
sample200 = data.frame(list200f, list200t)

}

fit1 = lm(sample150$list150t ~ poly(sample150$list150f,4,raw=TRUE)) # polynomial regression for each delta p : not used

plot(sample50, xlab = "FPR", ylab = "TPR", pch = 21, bg= 'red', ylim = c(0,1), xlim= c(0,1), cex=.8, main = c('Delta P = ', deltaP))
par(new=TRUE)
plot(sample100, xlab = "FPR", ylab = "TPR", pch = 21, bg= 'blue',ylim = c(0,1), xlim= c(0,1), cex=.8) #type='n' removes points
par(new=TRUE)
plot(sample150, xlab = "FPR", ylab = "TPR", pch = 21, bg= 'green', ylim = c(0,1), xlim= c(0,1), cex=.8)
par(new=TRUE)
plot(sample200, xlab = "FPR", ylab = "TPR", pch = 21, bg= 'orange',ylim = c(0,1), xlim= c(0,1), cex=.8)

legend(.6,.4, cex=.8, c("25 depth","50 depth","75 depth","100 depth"), lty=c(1,1), lwd=c(1,1),col=c("red","blue","green","orange"))

#points(sample100$list100f, predict(fit2), type = 'l', col='blue',lwd=1)

}






ROCgraph = function(){

attach(mtcars)
par(mfrow=c(2,2))

(ROCcurve(.05))
(ROCcurve(.1))
(ROCcurve(.15))
(ROCcurve(.2))

}

ROCgraph()



