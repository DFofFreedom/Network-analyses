##Load R packages that are required

library(asnipe)
library(ape)
library(igraph)
library(sna)
library(ergm)
library(ergm.count)
library(tnet)
library(assortnet)


##Read in test data
##(Paths will need changing)

#true network
setwd("C:/Users/matth/Desktop/Dave_Julian simulations folder/Network-analyses-jevansbranch/test/8_0.4_0.4_-0.6_0.25_0.2_1000/truenet")
truenet<-read.csv("1.csv",check.names=FALSE)

#observed network
setwd("C:/Users/matth/Desktop/Dave_Julian simulations folder/Network-analyses-jevansbranch/test/8_0.4_0.4_-0.6_0.25_0.2_1000/obsnet")
obsnet<-read.csv("1.csv",check.names=FALSE)

#association network
setwd("C:/Users/matth/Desktop/Dave_Julian simulations folder/Network-analyses-jevansbranch/test/8_0.4_0.4_-0.6_0.25_0.2_1000/obsgbimat")
gbimat<-read.csv("1.csv",check.names=FALSE)

#data on individuals
setwd("C:/Users/matth/Desktop/Dave_Julian simulations folder/Network-analyses-jevansbranch/test/8_0.4_0.4_-0.6_0.25_0.2_1000/popdat")
indiv_dat<-read.csv("1.csv",check.names=FALSE)

#location of groups
setwd("C:/Users/matth/Desktop/Dave_Julian simulations folder/Network-analyses-jevansbranch/test/8_0.4_0.4_-0.6_0.25_0.2_1000/obsgbigroups")
grouplocs1<-read.csv("1.csv",check.names=FALSE)
names(grouplocs1)<-"Gr"

##--------------------------------------------------------------------------

##Create association networks

#create association network from gbi matrix
gbinet<-get_network(gbimat, data_format = "GBI",identities=colnames(gbimat))

#generate a second "empty" association network
#This will be used to create a count (rather than association index) version of the association network
gbinet2<-matrix(0,nr=nrow(indiv_dat),nc=nrow(indiv_dat))
colnames(gbinet2)<-rownames(gbinet2)<-names(gbimat)

#fill with association data
for(i in 1:nrow(gbinet2)){
for(j in 1:ncol(gbinet2)){
gbinet2[i,j]<-sum(gbimat[,i]==1&gbimat[,j]==1)
}
}

#calculate numberof times each individual is observed in the association network (to use as an explanatory variable)
nres<-diag(gbinet2)

#then set diagonal of the matrix to zero
diag(gbinet2)<-0

##--------------------------------------------------------------------------------

#this works out home range centroids for each individual and distances between them

grouplocs1$x<-rep(NA,nrow(grouplocs1))
grouplocs1$y<-rep(NA,nrow(grouplocs1))

for(i in 1:nrow(grouplocs1)){
grouplocs1$x[i]<-indiv_dat$x[min(which(indiv_dat$groups==grouplocs1$Gr[i]))]
grouplocs1$y[i]<-indiv_dat$y[min(which(indiv_dat$groups==grouplocs1$Gr[i]))]
}

indiv.groups<-list()
for(i in 1:nrow(indiv_dat)){
indiv.groups[[i]]<-which(gbimat[,i]==1)
}

indiv.centroids<-data.frame(indiv_dat$indivs)
names(indiv.centroids)<-"id"

indiv.centroids$x<-rep(NA,nrow(indiv.centroids))
indiv.centroids$y<-rep(NA,nrow(indiv.centroids))

for(i in 1:nrow(indiv.centroids)){
indiv.centroids$x[i]<-mean(grouplocs1$x[indiv.groups[[i]]])
indiv.centroids$y[i]<-mean(grouplocs1$y[indiv.groups[[i]]])
}

dist.centroids<-as.matrix(dist(indiv.centroids[,2:3]))

##-------------------------------------------------------------------------------------------------------

##Here we run some "diagnostics" to make sure networks are similar to each other in the expected manner

#First set is for weighted matrices

A<-netlm(list(gbinet),list(truenet),mode="graph",nullhyp=c("qapspp"))
B<-netlm(list(obsnet),list(truenet),mode="graph",nullhyp=c("qapspp"))
C<-netlm(list(gbinet),list(obsnet),mode="graph",nullhyp=c("qapspp"))

#---------------------------------

#Second set is for binary matrices

A_bin<-netlm(list(sign(gbinet)),list(sign(truenet)),mode="graph",nullhyp=c("qapspp"))
B_bin<-netlm(list(sign(obsnet)),list(sign(truenet)),mode="graph",nullhyp=c("qapspp"))
C_bin<-netlm(list(sign(gbinet)),list(sign(obsnet)),mode="graph",nullhyp=c("qapspp"))

#-----------------------------------------------------------------------------------------------------------

##remove any individuals that were not observed in any of the networks
##(going back through this I can't remember why this was necessary but it obviously was....)
##Seems to make sense as a thing to do anyway

miss.true<-which(colSums(truenet)==0)
miss.obs<-which(colSums(obsnet)==0)
miss.gbi<-which(colSums(gbinet2)==0)


if(length(miss.true)>0){
truenet<-truenet[!miss.true,!miss.true]
}

if(length(miss.obs)>0){
obsnet<-obsnet[!miss.obs,!miss.obs]
}

if(length(miss.gbi)>0){
gbinet<-gbinet[-miss.gbi,-miss.gbi]
gbinet2<-gbinet2[-miss.gbi,-miss.gbi]
gbimat<-gbimat[,-miss.gbi]
dist.centroids2<-dist.centroids[-miss.gbi,-miss.gbi]
}

if(length(miss.true)==0){
miss.true<-nrow(indiv_dat)+1
}

if(length(miss.obs)==0){
miss.obs<-nrow(indiv_dat)+1
}

##---------------------------------------------------------------------------------------------------------------------

#do some plotting and some calculation/coomparison of metrics

#create igraph network objects
gbi.NET<-graph.adjacency(gbinet,mode="undirected",weighted=TRUE,diag=FALSE)
obs.NET<-graph.adjacency(as.matrix(obsnet),mode="undirected",weighted=TRUE,diag=FALSE)
tru.NET<-graph.adjacency(as.matrix(truenet),mode="undirected",weighted=TRUE,diag=FALSE)

#plot association network using igraph
gbi.NET<-igraph::set.vertex.attribute(gbi.NET, "group", index=V(gbi.NET), indiv_dat$groups[-miss.gbi])
gbi.NET<-igraph::set.vertex.attribute(gbi.NET, "sex", index=V(gbi.NET), indiv_dat$sex[-miss.gbi])	   		
V(gbi.NET)$color=V(gbi.NET)$group #assign the "Group" attribute as the vertex color
igraph::plot.igraph(gbi.NET,edge.width=(E(gbi.NET)$weight)^0.25,vertex.shape=c("circle","square")[factor(V(gbi.NET)$sex)],layout=layout.fruchterman.reingold(gbi.NET),vertex.size=8,vertex.label=NA,margin=c(0,0,0,0))


#----------------------------------

#plot observed network using igraph
obs.NET<-igraph::set.vertex.attribute(obs.NET, "group", index=V(obs.NET), indiv_dat$groups[-miss.obs])
obs.NET<-igraph::set.vertex.attribute(obs.NET, "sex", index=V(obs.NET), indiv_dat$sex[-miss.obs])	   	
V(obs.NET)$color=V(obs.NET)$group #assign the "Group" attribute as the vertex color
igraph::plot.igraph(obs.NET,edge.width=(E(obs.NET)$weight)^0.25,vertex.shape=c("circle","square")[factor(V(obs.NET)$sex)],layout=layout.fruchterman.reingold(obs.NET),vertex.size=8,vertex.label=NA,margin=c(0,0,0,0))


#----------------------------------

#plot true network using igraph
tru.NET<-igraph::set.vertex.attribute(tru.NET, "group", index=V(tru.NET), indiv_dat$groups[-miss.true])
tru.NET<-igraph::set.vertex.attribute(tru.NET, "sex", index=V(tru.NET), indiv_dat$sex[-miss.true])	   
V(tru.NET)$color=V(tru.NET)$group #assign the "Group" attribute as the vertex color
igraph::plot.igraph(tru.NET,edge.width=(E(tru.NET)$weight)^0.25,vertex.shape=c("circle","square")[factor(V(tru.NET)$sex)],layout=layout.fruchterman.reingold(tru.NET),vertex.size=8,vertex.label=NA,margin=c(0,0,0,0))


##-------------------------------------------------------------------------------------------------------------------------
##-------------------------------------------------------------------------------------------------------------------------
##-------------------------------------------------------------------------------------------------------------------------

##fit some ERGMs to the data

#set up true network as a network object (that's right, another sna package!)
tru.NET2.edgelist<-as.tnet(as.matrix(truenet))
tru.NET2<-network(tru.NET2.edgelist[,1:2],directed=FALSE)
set.edge.attribute(tru.NET2,"weight",as.vector(tru.NET2.edgelist[,3]))

#set up attributes
set.vertex.attribute(tru.NET2,"group",as.vector(indiv_dat$groups[-miss.true]))
set.vertex.attribute(tru.NET2,"sex",as.vector(indiv_dat$sex[-miss.true]))

#create a shared group matrix (binary) to use as a dyadic covariate in the model
sh.gr<-array(NA,dim=dim(truenet))
for(i in 1:ncol(sh.gr)){
 for(j in 1:nrow(sh.gr)){
  ifelse(indiv_dat$groups[indiv_dat$indivs==colnames(truenet)[i]]==indiv_dat$groups[indiv_dat$indivs==colnames(truenet)[j]],sh.gr[i,j]<-1,sh.gr[i,j]<-0)
 }
}
diag(sh.gr)<-0

#create a matrix of distances between groups to use as a covariate in the model
dist.gr<-array(NA,dim=dim(truenet))
for(i in 1:ncol(dist.gr)){
 for(j in 1:nrow(dist.gr)){
  dist.gr[i,j]<-dist(rbind(c(indiv_dat$x[indiv_dat$indivs==colnames(truenet)[i]],indiv_dat$y[indiv_dat$indivs==colnames(truenet)[i]]),c(indiv_dat$x[indiv_dat$indivs==colnames(truenet)[j]],indiv_dat$y[indiv_dat$indivs==colnames(truenet)[j]])))
 }
}

#Run the model (shared group effect only)
modA<-ergm(tru.NET2~sum+nonzero+nodefactor("sex")+nodematch("sex")+edgecov(sh.gr),reference=~Poisson,response="weight")

#Run an alternative model (shared group effect and distance between group effect)
modA_2<-ergm(tru.NET2~sum+nonzero+nodefactor("sex")+nodematch("sex")+edgecov(sh.gr)+edgecov(dist.gr),reference=~Poisson,response="weight")

#--------------------------------------------------------------

#set up observed network as a network object
obs.NET2.edgelist<-as.tnet(as.matrix(obsnet))
obs.NET2<-network(obs.NET2.edgelist[,1:2],directed=FALSE)
set.edge.attribute(obs.NET2,"weight",as.vector(obs.NET2.edgelist[,3]))

#add attributes
set.vertex.attribute(obs.NET2,"group",as.vector(indiv_dat$groups[-miss.obs]))
set.vertex.attribute(obs.NET2,"sex",as.vector(indiv_dat$sex[-miss.obs]))

#shared group matrix
sh.gr<-array(NA,dim=dim(obsnet))
for(i in 1:ncol(sh.gr)){
 for(j in 1:nrow(sh.gr)){
  ifelse(indiv_dat$groups[indiv_dat$indivs==colnames(obsnet)[i]]==indiv_dat$groups[indiv_dat$indivs==colnames(obsnet)[j]],sh.gr[i,j]<-1,sh.gr[i,j]<-0)
 }
}
diag(sh.gr)<-0

#distance between groups matrix
dist.gr<-array(NA,dim=dim(obsnet))
for(i in 1:ncol(dist.gr)){
 for(j in 1:nrow(dist.gr)){
  dist.gr[i,j]<-dist(rbind(c(indiv_dat$x[indiv_dat$indivs==colnames(obsnet)[i]],indiv_dat$y[indiv_dat$indivs==colnames(obsnet)[i]]),c(indiv_dat$x[indiv_dat$indivs==colnames(obsnet)[j]],indiv_dat$y[indiv_dat$indivs==colnames(obsnet)[j]])))
 }
}

#Run first model (shared group only)
modB<-ergm(obs.NET2~sum+nonzero+nodefactor("sex")+nodematch("sex")+edgecov(sh.gr),reference=~Poisson,response="weight")

#Run alternative model (shared group + distance)
modB_2<-ergm(obs.NET2~sum+nonzero+nodefactor("sex")+nodematch("sex")+edgecov(sh.gr)+edgecov(dist.gr),reference=~Poisson,response="weight")

#--------------------------------------------------------------

#association network as network package object
gbi.NET2.edgelist<-as.tnet(as.matrix(gbinet2))
gbi.NET2<-network(gbi.NET2.edgelist[,1:2],directed=FALSE)
set.edge.attribute(gbi.NET2,"weight",as.vector(gbi.NET2.edgelist[,3]))

#add attributes
set.vertex.attribute(gbi.NET2,"group",as.vector(indiv_dat$groups[-miss.gbi]))
set.vertex.attribute(gbi.NET2,"sex",as.vector(indiv_dat$sex[-miss.gbi]))
set.vertex.attribute(gbi.NET2,"nres",nres[-miss.gbi])

#shared group matrix
sh.gr<-array(NA,dim=dim(gbinet2))
for(i in 1:ncol(sh.gr)){
 for(j in 1:nrow(sh.gr)){
  ifelse(indiv_dat$groups[indiv_dat$indivs==colnames(gbinet2)[i]]==indiv_dat$groups[indiv_dat$indivs==colnames(gbinet2)[j]],sh.gr[i,j]<-1,sh.gr[i,j]<-0)
 }
}
diag(sh.gr)<-0

##distance matrix pre-calculated above using home range info

#Run two models with/without home range info
#also have an effect for the number of resightings of each individual
modC3<-ergm(gbi.NET2~sum+nonzero+nodefactor("sex")+nodecov("nres")+nodematch("sex"),reference=~Poisson,response="weight")
modC4<-ergm(gbi.NET2~sum+nonzero+nodefactor("sex")+nodecov("nres")+nodematch("sex")+edgecov(dist.centroids2),reference=~Poisson,response="weight")

##------------------------------------------------------------------------------------------------------------------------------
##------------------------------------------------------------------------------------------------------------------------------
##------------------------------------------------------------------------------------------------------------------------------

##Now for the randomisations approach

#calculate observed homophily for each network
assort.gog<-assortment.discrete(gbinet2,indiv_dat$sex[-miss.gbi])$r
assort.tru<-assortment.discrete(truenet,indiv_dat$sex[-miss.true])$r
assort.obs<-assortment.discrete(obsnet,indiv_dat$sex[-miss.obs])$r


#-----------------------------

#Model sex differences in degree
#I'm going to use the GLM/GLMM approach from appendix 2 of Farine and Whitehad 2015. 
#No random effect needed yet, however could add one for group/location

#Calculate weighted degree
str.tru<-colSums(truenet)
str.obs<-colSums(obsnet)
str.gog<-colSums(gbinet2)

#Plot weighted degree
plot(str.obs~str.tru)

#Look at histogram of weighted degree (and examine mean/variance)
hist(log(str.tru))
var(str.gog)
mean(str.gog)

#model weighted degree in each network using Poisson GLMMs
tru.mod<-glm(str.tru~indiv_dat$sex[-miss.true],family=poisson)
obs.mod<-glm(str.obs~indiv_dat$sex[-miss.obs],family=poisson)
gog.mod<-glm(str.gog~indiv_dat$sex[-miss.gbi],family=poisson)

##---------------------------------------------------------------------------------------------------
##---------------------------------------------------------------------------------------------------

##And now generate randomised networks using datastream permutations for association-based networks

#create empty vectors to store results
r.assorts<-rep(NA,10000)
r.effs<-rep(NA,10000)

#set up an object that's useful for the randomisations
perm<-list(gbimat,gbinet)

#start loop through datastream permutations
for(ii in 1:10000){
 
 #generate permuted data
 random_networks_gbi<-network_swap(association_data=perm[[1]],association_matrix=perm[[2]], swaps=10) 
 newgbi<-random_networks_gbi[[2]]
 r.gbinet<-random_networks_gbi[[1]]

 #recreate the object used for re-permuting data
 perm<-list(newgbi,r.gbinet)

 #calculate count matrix for associations
 r.gbinet2<-matrix(0,nr=nrow(indiv_dat[-miss.gbi,]),nc=nrow(indiv_dat[-miss.gbi,]))
 colnames(r.gbinet2)<-rownames(r.gbinet2)<-names(newgbi)
 for(i in 1:nrow(r.gbinet2)){
  for(j in 1:ncol(r.gbinet2)){
   r.gbinet2[i,j]<-sum(newgbi[,i]==1&newgbi[,j]==1)
  }
 }

 #calculate assortativity of randomised network
 tmp.assort<-assortment.discrete(r.gbinet2,indiv_dat$sex[-miss.gbi])
 
 #model differences in degree in randomised network
 str.tmp<-colSums(r.gbinet2)
 tmp.mod<-glm(str.tmp~indiv_dat$sex[-miss.gbi],family=poisson)

 #store m vs f estimate from model
 r.effs[ii]<-coef(tmp.mod)[2]
 #store assortativity coefficient
 r.assorts[ii]<-tmp.assort$r
 print(ii)
}

#------------------------------------------------------------------------------

#Use node swaps (in sna) to do randomisations for true and observed networks

#empty vectors for results
r.assorts.tr<-rep(NA,10000)
r.effs.tr<-rep(NA,10000)

trnet<-truenet

for(ii in 1:10000){
 #permute network
 random_networks_tr<-rmperm(trnet) 

 trnet<-random_networks_tr
 
 #calculate assortativity and degree effect
 tmp.assort<-assortment.discrete(trnet,indiv_dat$sex[-miss.true])
 str.tmp<-colSums(trnet)
 tmp.mod<-glm(str.tmp~indiv_dat$sex[-miss.true],family=poisson)

 #store temporary results
 r.effs.tr[ii]<-coef(tmp.mod)[2]
 r.assorts.tr[ii]<-tmp.assort$r
 
 print(ii)
}

#empty vectors for results
r.assorts.obs<-rep(NA,10000)
r.effs.obs<-rep(NA,10000)

onet<-obsnet

for(ii in 1:10000){

 #permute network
 random_networks_o<-rmperm(onet) 
 onet<-random_networks_o

 #calculate assortativity and degree effect
 tmp.assort<-assortment.discrete(onet,indiv_dat$sex[-miss.obs])
 str.tmp<-colSums(onet)
 tmp.mod<-glm(str.tmp~indiv_dat$sex[-miss.obs],family=poisson)

 #store results
 r.effs.obs[ii]<-coef(tmp.mod)[2]
 r.assorts.obs[ii]<-tmp.assort$r

 print(ii)
}

##-------------------------------------------------------------------------------------------------------------------------------------------------
##-------------------------------------------------------------------------------------------------------------------------------------------------
##-------------------------------------------------------------------------------------------------------------------------------------------------

##output needs to be ERGM estimates, errors and p values, glm estimates errors and p values, randomisation 2.5 and 97.5 confidence ints and p val
##for all three networks

##ERGMoutputs

ergm_res<-list()

ergm_res[[1]]<-summary(modA_2)$coefs
ergm_res[[2]]<-summary(modB_2)$coefs
ergm_res[[3]]<-summary(modC4)$coefs

#-----------------------------------------------

##randomisations outputs
#assortativity

ran_res<-list()

tru.assort<-c("T",assort.tru,quantile(r.assorts.tr,c(0.005,0.025,0.975,0.995)),sum(r.assorts.tr<assort.tru)/10001)
obs.assort<-c("O",assort.obs,quantile(r.assorts.obs,c(0.005,0.025,0.975,0.995)),sum(r.assorts.obs<assort.obs)/10001)
gog.assort<-c("G",assort.gog,quantile(r.assorts,c(0.005,0.025,0.975,0.995)),sum(r.assorts<assort.gog)/10001)
names(tru.assort)<-names(obs.assort)<-names(gog.assort)<-c("Network","Value","q0.5","q2.5","q97.5","q99.5","P")

assort.out<-data.frame(rbind(tru.assort,obs.assort,gog.assort))

ran_res[[1]]<-assort.out

#-----------------------------------------------

#sex effect

tru.eff<-c("T",coef(tru.mod)[2],quantile(r.effs.tr,c(0.005,0.025,0.975,0.995)),sum(r.effs.tr<coef(tru.mod)[2])/10001)
obs.eff<-c("O",coef(obs.mod)[2],quantile(r.effs.obs,c(0.005,0.025,0.975,0.995)),sum(r.effs.obs<coef(obs.mod)[2])/10001)
gog.eff<-c("G",coef(gog.mod)[2],quantile(r.effs,c(0.005,0.025,0.975,0.995)),sum(r.effs<coef(gog.mod)[2])/10001)
names(tru.eff)<-names(obs.assort)<-names(gog.assort)<-c("Network","Value","q0.5","q2.5","q97.5","q99.5","P")

effs.out<-data.frame(rbind(tru.eff,obs.eff,gog.eff))

ran_res[[2]]<-effs.out



##-----------------------------------------------------------
##-----------------------------------------------------------
##-----------------------------------------------------------
##-----------------------------------------------------------
##-----------------------------------------------------------
##-----------------------------------------------------------


###Some abandoned code from when we were going to model edge values as well

#now GLMMs of edge variables to test homophily and sex effects together
#currently random effects of each individual in the dyad. Additional random effects could be added

edges<-gbinet2[upper.tri(gbinet2)]
indiv1<-numeric()

for(i in 2:ncol(gbinet2)){
tmp<-seq(1,i-1,1)
indiv1<-c(indiv1,tmp)
}

indiv2<-numeric()
for(i in 2:ncol(gbinet2)){
tmp<-rep(i,i-1)
indiv2<-c(indiv2,tmp)
}

sex1<-factor(levels=c("F","M"))
sex2<-factor(levels=c("F","M"))
for(i in 1:length(indiv1)){
sex1[i]<-indiv_dat$sex[indiv_dat$indivs==indiv1[i]]
sex2[i]<-indiv_dat$sex[indiv_dat$indivs==indiv2[i]]
}

samesex<-as.numeric(sex1==sex2)

gog.edgelist<-data.frame(indiv1,indiv2,sex1,sex2,samesex,edges)

##discover that edges are horribly zero-inflated. Back off for now and use randomisations instead




