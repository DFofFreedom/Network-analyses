network.generator<-function(groups,mean.group.size,max.group.size,d.eff,o.dens,i.dens,plot=F){

#and use this to calculate the overall size of the population
n.indivs<-mean.group.size*groups

#create the individuals
indivs<-seq(1,n.indivs,1)

##and sample individuals into groups
#size of possible groups is equal to max group size, then sample the overall population size from this
poss.groups<-rep(seq(1,groups,1),each=max.group.size)
indiv.groups<-sample(poss.groups,n.indivs,replace=F)

#arrange in data frame and order by group
inds<-data.frame(indivs,indiv.groups)

###SPATIAL LOCATION OF GROUPS/CLUSTERS###

#create a dataframe of group locations
#points end up on a grid in space
group.id<-seq(1,groups,1)
group.x<-rnorm(length(group.id),3,1)
group.y<-rnorm(length(group.id),3,1)
xs<-rep(seq(1,groups,1),each=groups)
ys<-rep(seq(1,groups,1),groups)
rows<-sample(1:length(xs),groups,replace=F)
for(i in 1:groups){
group.x[i]<-xs[rows[i]]
group.y[i]<-ys[rows[i]]
}

#create dataframe of spatial information
spat<-data.frame(group.id,group.x,group.y)

#plot to check spatial distribution
#dev.new()
#plot(group.x,group.y)

#create distance matrix for groups
dists<-as.matrix(dist(spat[,2:3]))
rownames(dists)<-colnames(dists)<-group.id

#standardise and invert (so 1 is closest and 0.001 is furthest away)
   #edited this to a neater version and re-named
close<-1/(1+dists)
diag(close)<-1

#-----------------------------------------------------------------------------------------------------------------

#####NETWORK STUFF#####
#create empty network in association matrix form
net.d<-array(NA,dim=rep(nrow(inds),2))
colnames(net.d)<-rownames(net.d)<-inds[,1]

#create network info
#current negbinoms are:
#within group - size=i.dens and prob=0.3 
#out of group - size=o.dens x distance (remember higher is closer) then to the power of d.eff. 

for (i in 1:(nrow(inds)-1)){
   for(j in (i+1):nrow(inds)){
      if(inds[,2][i]==inds[,2][j]){
         net.d[i,j]<-round(rnbinom(1,size=i.dens,prob=0.3))
      }
      else{
         net.d[i,j]<-round(rnbinom(1,size=o.dens*(close[which(colnames(close)==inds[,2][i]),
                                                        which(colnames(close)==inds[,2][j])])^d.eff,
                                   prob=0.3))
      }
   }
}    

#this sets the diagonal to zero i.e. no self loops in network
net3.d<-t(net.d)
diag(net3.d)<-0

#this symmetrises the association matrix
for (i in 1:nrow(inds)){
for (j in 1:nrow(inds)){
if(j<i){
net3.d[j,i]<-net.d[j,i]
}
}
}


#plots the graph of each network made if plot=TRUE
if(plot==T){
   dev.new()
   par(mfrow=c(1,1))
   net2.d<-graph.adjacency(net.d,mode="undirected",weighted=TRUE,diag=FALSE)
   net2.d<-set.vertex.attribute(net2.d, "group", index=V(net2.d), inds[,2])
   V(net2.d)$group
   #V(net2.d)$color=V(net2.d)$group #assign the "Group" attribute as the vertex color
   V(net2.d)$color="blue"
   plot(net2.d,edge.width=(E(net2.d)$weight)^0.25,layout=layout.fruchterman.reingold(net2.d),vertex.size=8,vertex.label=NA,margin=c(0,0,0,0))
}

#fill in lists with network data
net4.d<-net3.d
for(i in 1:nrow(net4.d)){
for(j in 1:ncol(net4.d)){
if(is.na(net4.d[i,j])==TRUE){
net4.d[i,j]<-0
}
}
}

pop.dat<-list(inds,net4.d)

return(pop.dat)

#end function
}
