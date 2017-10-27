network.generator<-function(groups,mean.group.size,max.group.size,d.eff,o.dens,i.dens,plot=F){
	#i.dens = density of within group associations?
	#o.dens = density of outside group associations?
	#d.eff = 
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
	inds<-inds[order(indiv.groups),]
	
	#Generate sex vector. This sounds rude.
	inds$sex=sample(c("M","F"),nrow(inds),replace=T)

	###SPATIAL LOCATION OF GROUPS/CLUSTERS###

	#create a dataframe of group locations
	#points end up on a grid in space
	group.id<-seq(1,groups,1)

	group.x=rep(seq(1,groups,1),each=groups)[sample(1:groups^2,groups,replace=F)]
	group.y=rep(seq(1,groups,1),each=groups)[sample(1:groups^2,groups,replace=F)]

	#create dataframe of spatial information
	spat<-data.frame(group.id,group.x,group.y)

	#plot to check spatial distribution
	if(plot==T){
		dev.new()
		plot(group.x,group.y,pch=16,col=1:groups)
	}

	#create distance matrix for groups
	dists<-as.matrix(dist(spat[,2:3]))
	rownames(dists)<-colnames(dists)<-group.id

	#standardise and invert (so 1 is closest and 0.001 is furthest away)
	dists2<-dists/max(dists)
	dists3<-1.001-dists2
	diag(dists3)<-1

	inds$x=group.x[inds[,2]]
	inds$y=group.x[inds[,2]]
	#-----------------------------------------------------------------------------------------------------------------

	#####NETWORK STUFF#####
	#create empty network in association matrix form
	net.d<-array(NA,dim=rep(nrow(inds),2))
	colnames(net.d)<-rownames(net.d)<-inds[,1]

	#create network info
	#current negbinoms are:
	#within group - size=i.dens and prob=0.3 and then multiply output by 10
	#out of group - size=o.dens x distance (remember higher is closer) then to the power of d.eff. 
	#For out of group interactions the edge weight for each indiv is calculated and added together. This is because of the old code this has been adapted from and could be removed

	dyads=which(upper.tri(net.d),arr.ind=T)
	dsex=sapply(1:2,function (x) inds[dyads[,x],"sex"])
	dsites=sapply(1:2,function (x) inds[dyads[,x],2])

	distsv=sapply(1:nrow(dsites),function (x) dists3[dsites[x,1],dsites[x,2]])

	samesex=dsex[,1]==dsex[,2]
	samesite=dsites[,1]==dsites[,2]

	
	#What is the effect of samesex? What is the effect of sex? 
	
	net.d[rbind(dyads[samesite,],dyads[samesite,c(2,1)])]=sapply(which(samesite),function (x) round(rnbinom(1,size=i.dens*(distsv[x])^d.eff,prob=0.3))+round(rnbinom(1,size=i.dens*(distsv[x])^d.eff,prob=0.3)))
	net.d[rbind(dyads[!samesite,],dyads[!samesite,c(2,1)])]=sapply(which(!samesite),function (x) round(rnbinom(1,size=o.dens*(distsv[x])^d.eff,prob=0.3))+round(rnbinom(1,size=o.dens*(distsv[x])^d.eff,prob=0.3)))
	diag(net.d)=0
		

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

	pop.dat<-list(ind_data=inds,network=net.d,distmat=dists3)

	return(pop.dat)

	#end function
}