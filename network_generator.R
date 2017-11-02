network.generator<-function(groups,mean.group.size,max.group.size,d.eff,o.dens,i.dens,sex.eff=NA,m.i.eff=NA,m.o.eff=NA,plot=F){
	#####
	#Network generator: a function that does STUFF.
	#i.dens = density of within group associations?
	#o.dens = density of outside group associations?
	#d.eff = d effect description
	#m.i.eff = effect of being a male on within group assoc
	#m.o.eff = effect of being a male on outside group assoc
	#sex.eff = same sex effect description
	#and so forth
	#####

	if(is.na(m.i.eff)){
		m.i.eff=0
	}

	if(is.na(m.o.eff)){
		m.o.eff=m.i.eff # if no outside sex effect included, same as inside
	}

	if(is.na(sex.eff)){
		sex.eff=1
	}

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

	##could try this instead? Neater
	#dists3= 1/(1+dists) # neater, but does not give exactly the same results
	

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
	ismale=dsex[,1]=="M"	
	
	#What is the effect of samesex? What is the effect of sex? 
	#added this utility function to make code more readable
	#move elsewhere eventually
	symmat = function (bool,dyads){
		return(rbind(dyads[bool,],dyads[bool,c(2,1)]))
	}
	
	#Wasn't sure of your rationale for change the equations - have retained.
	#o.dens and i.dens standing for the effects for females
	net.d[symmat(samesite,dyads)]=sapply(which(samesite),function (x) round(rnbinom(1,size=i.dens*(distsv[x])^d.eff,prob=0.3))+
		round(rnbinom(1,size=i.dens*(distsv[x])^d.eff,prob=0.3)))
		
	net.d[symmat(!samesite,dyads)]=sapply(which(!samesite),function (x) round(rnbinom(1,size=o.dens*(distsv[x])^d.eff,prob=0.3))+
		round(rnbinom(1,size=o.dens*(distsv[x])^d.eff,prob=0.3)))
	
	#and then you an m.eff which adds/substracts a particular value to these for males
	net.d[symmat(ismale&samesite,dyads)]=net.d[symmat(ismale&samesite,dyads)]-m.i.eff

	net.d[symmat(ismale&!samesite,dyads)]=net.d[symmat(ismale&!samesite,dyads)]-m.o.eff
	
	#between versus same-sex. Baseline = FM, then have separate effects for FF and MM
	#set the same to keep things simple for now
	net.d[symmat(samesex,dyads)]= net.d[symmat(samesex,dyads)]*sex.eff

	diag(net.d)=0
	net.d[net.d<0]=0 #Possibly not needed - but for testing purposes works fine

	
		

	#plots the graph of each network made if plot=TRUE
	if(plot==T){
		dev.new()
		par(mfrow=c(1,1))
		net2.d<-graph.adjacency(net.d,mode="undirected",weighted=TRUE,diag=FALSE)
		net2.d<-set.vertex.attribute(net2.d, "group", index=V(net2.d), inds$indiv.groups)
		net2.d<-set.vertex.attribute(net2.d, "sex", index=V(net2.d), inds$sex)	   
		
		#V(net2.d)$color=V(net2.d)$group #assign the "Group" attribute as the vertex color
		#V(net2.d)$color="blue"
		#V(net2.d)$color=factor(V(net2.d)$sex) #assign the "sex" attribute as the vertex color
		plot(net2.d,edge.width=(E(net2.d)$weight)^0.25,layout=layout.fruchterman.reingold(net2.d),vertex.size=8,vertex.label=NA,margin=c(0,0,0,0))
	}

	pop.dat<-list(ind_data=inds,network=net.d,distmat=dists3)

	return(pop.dat)
}