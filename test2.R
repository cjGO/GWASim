test = function(){

list50f = c()
list50t = c()
list100f = c()
list100t = c()

deltaP = .2

depth=100

false = NfisherSim(100,(depth/2),(depth/2),.5,.5)
true = NfisherSim(100,(depth/2),(depth/2),.5,(.5+ deltaP)) # delta  p  = .2

for (pcut in seq(.01,.05,.0025)){

fp = sum(false < pcut)
tp = sum(true < pcut)

fpr = (fp / 100 )
tpr = (tp / 100 )

if (depth == 50){
list50f = append(list50f, fpr)
list50t = append(list50t, tpr)
}
if (depth == 100){
list50f = append(list100f, fpr)
list50t = append(list100t, tpr)
}
if (depth == 150){
list50f = append(list150f, fpr)
list50t = append(list150t, tpr)
}
if (depth == 200){
list50f = append(list200f, fpr)
list50t = append(list200t, tpr)
}}}
