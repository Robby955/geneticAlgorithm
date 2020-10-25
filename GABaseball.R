
baseball=read.table('baseball.txt',header=TRUE)

log_salary=log(baseball.dat$salary)
baseball.pred=baseball.dat[,-1]


N= nrow(baseball) # Will be 337
M= ncol(baseball.pred) # Will be 27
P= 20 # Number of individuals per generation.
numGenerations= 100 # We will use 100 generation
mu=0.015 # We will use this to increase diversity
genAic= matrix(0,nrow=P,ncol=1) # Used to store each gen AIC 
phi= matrix(0,nrow=P,ncol=1) # Used to store ranks of each cand

chromosome= matrix(0,nrow=P,ncol=M) # Save the predicotrs for each generation.

chromosome.next= matrix(0,nrow=P,ncol=M) # Save the predicotrs for each generation.


best_chromosome=NULL

chromosomeAIC=matrix(0,nrow=P,ncol=1) # Will be the AIC from each chromsome in that generation.
cacheAIC=matrix(0,nrow=P,ncol=numGenerations)
bestAIC=0 # Best global AIC
bestAIC_EachGen=rep(0,numGenerations) # Best AIC from this generation of individuals






set.seed(20)
for(i in 1:P){
  chromosome[i,] = rbinom(M,1,.5)
  chromosome.variables = baseball.pred[,chromosome[i,]==1]
  getModel= lm(log_salary~.,chromosome.variables)
  chromosomeAIC[i] = extractAIC(getModel)[2]
  cacheAIC[i,1] = chromosomeAIC[i]
  if(chromosomeAIC[i] < bestAIC){
    chrom = chromosome[i,]
    bestAIC = chromosomeAIC[i]
  }
}


r = rank(-chromosomeAIC)
phi = 2*r/(P*(P+1))
bestAIC_EachGen[1]=bestAIC







## MAIN
for(j in 1:numGenerations-1){
  
  for(i in 1:(P/2)){
    parent_1 = chromosome[sample(1:P,1,prob=phi),]
    parent_2 = chromosome[sample(1:P,1),]
    pos = sample(1:(M-1),1)
    mutate = rbinom(M,1,mu)
    chromosome.next[i,] = c(parent_1[1:pos],parent_2[(pos+1):M])
    chromosome.next[i,] = (chromosome.next[i,]+mutate)%%2
    mutate = rbinom(M,1,mu)
    chromosome.next[P+1-i,] = c(parent_2[1:pos],parent_1[(pos+1):M])
    chromosome.next[P+1-i,] = (chromosome.next[P+1-i,]+mutate)%%2
  }
  
  chromosome = chromosome.next
  
  for(i in 1:P){
    chromosome.variables = baseball.pred[,chromosome[i,]==1]
    getModel = lm(log_salary~.,chromosome.variables)
    chromosomeAIC[i] = extractAIC(getModel)[2]
    cacheAIC[i,j+1] = chromosomeAIC[i]
    if(chromosomeAIC[i] < bestAIC){
      best_chromosome = chromosome[i,]
      bestAIC = chromosomeAIC[i]
    }
  }
  bestAIC_EachGen[j+1]=bestAIC
  r = rank(-chromosomeAIC)
  phi = 2*r/(P*(P+1))
}

best_chromosome
bestAIC



#png(file='newbest.png')
plot(-cacheAIC,xlim=c(0,numGenerations),ylim=c(50,425),type="n",ylab="-AIC",
     xlab="Generation",main=" Negative AIC Value per generation")
for(i in 1:numGenerations){points(rep(i,P),-cacheAIC[,i],pch=18,col='red')}
#dev.off()

