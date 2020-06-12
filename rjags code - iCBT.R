
# load libraries ----------------------------------------------------------
library(netmeta)
library(rjags)
library(MCMCvis)

# load data ---------------------------------------------------------------
rm(list=ls()) 
pairwise.all=readRDS("treatment.level")
for.cNMA=readRDS("for.cNMA")


# pairwise MA -------------------------------------------------------------

#### pairwise comparison: CBT vs WL
dev.off()
pw1=pairwise.all[pairwise.all$comp=="CBT vs WL",]
pm1=metagen(data=pw1, TE=TE,seTE=seTE,studlab=studlab, byvar=type)

forest(pm1, prediction=T, label.right = "Favours WL",
       label.left = "Favours CBT", colgap.forest.left ="3.5cm", 
       leftcols=c("studlab"), xlim = c(-10,4), resid.hetstat=F, digits.tau2=2)
#### pairwise comparison: CBT vs TAU
pw2=pairwise.all[pairwise.all$comp=="CBT vs TAU",]
pm2=metagen(data=pw2, TE=TE,seTE=seTE,studlab=studlab, byvar=type)
forest(pm2, prediction=T, label.right = "Favours TAU",
       label.left = "Favours CBT", colgap.forest.left ="3.5cm", 
       leftcols=c("studlab"), xlim = c(-10,4), resid.hetstat=F,digits.tau2=2)


# NMA treatment level -----------------------------------------------------
#### NMA aggregated level, no components
n2=netmeta(data=pairwise.all, studlab = studlab,TE=TE,seTE=seTE,treat1 = treat1,treat2=treat2)
dev.off()
netgraph(n2, seq = "optimal", col = "black", plastic = FALSE,
         points = TRUE, pch = 21, cex.points = 3, col.points = "black",
         bg.points = "gray", thickness = "se.random",
         multiarm = FALSE, number.of.studies = TRUE)
nl2=netleague(n2, digits=1)
nma_primary=as.matrix(nl2$random)
netrank(n2)
netsplit(n2)

# cNMA --------------------------------------------------------------------

#### this analysis uses only patients with complete outcomes, and gives slighly different results than the 
#### analysis that uses the multiply imputed datasets from studies with IPD
model1.string <-  "
model {
####   AD 
for(i in 1:Ns) { 
w[i,1]<- 0
delta[i,1]<- 0                                             

for (k in 1:na[i])  { 
se[i,k]<-sd[i,k]/sqrt(n[i,k])
prec[i,k]<- (1/(se[i,k]*se[i,k]))
y[i,k]~dnorm(phi[i,k],prec[i,k])
phi[i,k]<-u[i]+delta[i,k]
}

for (k in 2:na[i]) {
delta[i,k] ~ dnorm(md[i,k],precd[i,k])          
md[i,k] <- mean[i,k]  + sw[i,k]                  
precd[i,k] <- 2*prect*(k-1)/k                                   
w[i,k] <- (delta[i,k] - mean[i,k] )         
sw[i,k] <-sum(w[i,1:(k-1)])/(k-1)        

##consistency equations
mean[i,k] <-A1[i,k]-B1[i]
A1[i,k]<-
d[1]*(1-equals(c1[i,k],0)) + 
d[3]*(1-equals(c3[i,k],0)) + d[4]*(1-equals(c4[i,k],0)) +
d[5]*(1-equals(c5[i,k],0)) + d[6]*(1-equals(c6[i,k],0)) + 
d[7]*(1-equals(c7[i,k],0)) + d[8]*(1-equals(c8[i,k],0)) + 
d[9]*(1-equals(c9[i,k],0)) + d[10]*(1-equals(c10[i,k],0))+
d[11]*(1-equals(c11[i,k],0)) + d[12]*(1-equals(c12[i,k],0))+
d[13]*(1-equals(c13[i,k],0)) + d[14]*(1-equals(c14[i,k],0))+ 
d[15]*(1-equals(c15[i,k],0)) + d[16]*(1-equals(c16[i,k],0))+
d[17]*(1-equals(c17[i,k],0))
}
B1[i]<-
d[1]*(1-equals(c1[i,1],0))  + 
d[3]*(1-equals(c3[i,1],0))  + d[4]*(1-equals(c4[i,1],0)) +
d[5]*(1-equals(c5[i,1],0))  + d[6]*(1-equals(c6[i,1],0)) + 
d[7]*(1-equals(c7[i,1],0))  + d[8]*(1-equals(c8[i,1],0)) + 
d[9]*(1-equals(c9[i,1],0))  + d[10]*(1-equals(c10[i,1],0)) + 
d[11]*(1-equals(c11[i,1],0))+ d[12]*(1-equals(c12[i,1],0))+
d[13]*(1-equals(c13[i,1],0))+ d[14]*(1-equals(c14[i,1],0))+
d[15]*(1-equals(c15[i,1],0)) +d[16]*(1-equals(c16[i,1],0))+
d[17]*(1-equals(c17[i,1],0))
}
##prior distribution for log-odds in baseline arm of study i
for (i in 1:Ns) {	u[i] ~ dnorm(0,.001)	}
##prior distribution for heterogeneity	
tau ~ dnorm(0,0.1)I(0,)                                      
prect<- 1/tau.sq
tau.sq<- pow(tau,2)

##prior distribution for basic parameters		
for(k in 1:Nc) {d[k] ~ dnorm(0,.001)}
}
"
model1.spec<-textConnection(model1.string) 
data <- list(y=for.cNMA$y,n=for.cNMA$n, sd=for.cNMA$sd, Nc=for.cNMA$Nc, c1=for.cNMA$c1, c3=for.cNMA$c3, 
             c4=for.cNMA$c4, c5=for.cNMA$c5, c6=for.cNMA$c6, c7=for.cNMA$c7, c8=for.cNMA$c8, 
             c9=for.cNMA$c9, c10=for.cNMA$c10, c11=for.cNMA$c11, c12=for.cNMA$c12, c13=for.cNMA$c13, 
             c14=for.cNMA$c14, Ns=for.cNMA$Ns, c15=for.cNMA$c15, c16=for.cNMA$c16, c17=for.cNMA$c17, na=for.cNMA$na)
jags.m=0
jags.m <- jags.model(model1.spec, data = data, n.chains =2, n.adapt = 5000)
params <- c("tau", "d"   ) 
closeAllConnections()
samps<- coda.samples(jags.m, params, n.iter =25000)
MCMCtrace(samps,pdf = FALSE, params = "d") 
A1= MCMCsummary(samps)
rownames(A1)[1:17]=components
rownames(A1)[1:17]=c("wl", "dt", "pl", "pe", "cr", "ba", "is", "ps", "re", "w3", 
                     "bi", "rp", "hw", "ff", "ae", "he", "tg")
round(A1[,c(4,3,5)], digits=3)

