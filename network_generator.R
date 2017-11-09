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

	require(igraph)

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
	inds<-data.frame(indivs,groups=indiv.groups)
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
	
	#and then you an m.eff which adds/substracts a particular value to these for males (depending on if effect is +/-)
	net.d[symmat(ismale&samesite,dyads)]=sapply(which(ismale&samesite),function (x) net.d[dyads[x,1],dyads[x,2]]+
		round((m.i.eff/abs(m.i.eff))*rnbinom(1,size=abs(m.i.eff),prob=0.3)))

	net.d[symmat(ismale&!samesite,dyads)]=sapply(which(ismale&!samesite),function (x) net.d[dyads[x,1],dyads[x,2]]+
		round((m.o.eff/abs(m.i.eff))*rnbinom(1,size=abs(m.o.eff),prob=0.3)))

	#between versus same-sex. Baseline = FM, then have separate effects for FF and MM
	#set the same to keep things simple for now
	# have kept as a simple scaling factor
	net.d[symmat(samesex,dyads)]= net.d[symmat(samesex,dyads)]*sex.eff

	diag(net.d)=0
	net.d[net.d<0]=0 #Possibly not needed - but for testing purposes works fine

	#plots the graph of each network made if plot=TRUE
	if(plot==T){
		dev.new()
		par(mfrow=c(1,1))
		net2.d<-graph.adjacency(net.d,mode="undirected",weighted=TRUE,diag=FALSE)
		net2.d<-set.vertex.attribute(net2.d, "group", index=V(net2.d), inds$groups)
		net2.d<-set.vertex.attribute(net2.d, "sex", index=V(net2.d), inds$sex)	   
		
		V(net2.d)$color=V(net2.d)$group #assign the "Group" attribute as the vertex color
		plot(net2.d,edge.width=(E(net2.d)$weight)^0.25,vertex.shape=c("circle","square")[factor(V(net2.d)$sex)],layout=layout.fruchterman.reingold(net2.d),vertex.size=8,vertex.label=NA,margin=c(0,0,0,0))
	}

	pop.dat<-list(ind_data=inds,network=net.d,distmat=dists3)

	return(pop.dat)
}

networkobs<-function(pop.dat,timesteps){
	#add observer effort parameter at some point

	#get objects
	inds=pop.dat$ind_data
	network=pop.dat$network
	distmat=pop.dat$distmat
	
	#data frame of dyads
	dyads=which(upper.tri(net.d),arr.ind=T)
	assocs=network[dyads]
	names1=as.numeric(row.names(network)[dyads[,1]])
	names2=as.numeric(colnames(network)[dyads[,2]])
	dyads=data.frame(names1,names2,dyads,assocs=assocs)
	
	#generate GBI mat. Select a random individual - check chance of being seen with other
	#indiv based on their associations. Repeat this for any individual added.
	gbimat=t(sapply(1:timesteps,function (x) makeevent(inds,dyads)))
	colnames(gbimat)=inds$indivs

	#what is the effect of observer effort? Randomly remove individuals from a grouping event?
	# Could calculate an individual's degree and base removal on combination of this and 
	# observer effort - less central individuals less likely to be seen when observer effort
	# is low?

	# Or should the probability of addition of an individual to a grouping event be dependant
	# on observer effort?
	
	#get observed interactions (for interaction based network)
	#
	
	#TODO: INTEGRATE INTERACTIONS FUNCTION HERE ONCE IT'S DONE

	return(list(gbimat=gbimat,interactions=interactions)
}

makeevent<-function(inds,dyads){
	#randomly pick an individual to start
	gbirow=matrix(0,1,length(inds$indivs))
	colnames(gbirow)=inds$indivs
	seed=sample(inds$indivs,1)
	gbirow[colnames(gbirow)==seed]=1	
	fassoc=focalassoc(seed,dyads,T)
	todo=fassoc
	checked=seed
	gbirow[colnames(gbirow)%in%fassoc]=1
	
	#check the associations of the individuals we are adding until we are no longer adding members
	while(length(todo)>0){
		fassoc=unlist(sapply(todo,function (x) focalassoc(x,dyads)))
		gbirow[colnames(gbirow)%in%fassoc]=1		
		checked=c(checked,todo)
		todo=unique(fassoc[!fassoc%in%checked])
	}
	return(gbirow)	
}

focalassoc<-function(focal,dyads,forcemulti=F){
	####
	# focal = ID of individual whose associations we are using
	# dyads = all possible dyads
	# forcemulti = Force function to return at least one interaction?
	####
	potentials=dyads[dyads$names1==focal|dyads$names2==focal,]
	#swap the names around for readability
	potentials$names2[potentials$names1!=focal]=potentials$names1[potentials$names1!=focal]
	potentials$names1[potentials$names1!=focal]=focal
	
	#Might need a better way of coming up with probability. Currently normalise these associations
	#by the maximum assoc strength of a network + 1. Depending on the effects being added to 
	#associations this can lead to some very small probablities

	if(forcemulti==T){
		gbirow=rep(0,nrow(potentials))
		while(sum(gbirow)==0){
			gbirow=rbinom(potentials$assoc, 1, potentials$assoc/(max(dyads$assoc)+1))
		}
	} else {
		gbirow=rbinom(potentials$assoc, 1, potentials$assoc/(max(dyads$assoc)+1))

	}
	return(potentials$names2[gbirow==1])
}

getinteractions(inds,dyads,gbirow,t,intfreq,obseff){
	#inds dataframe
	#dyads dataframe
	#gbirow = a row of a group by individual matrix (so we can apply this to one event at a time)
	#t = current timestep (the row number of the GBI matrix)
	#intfreq = what is the frequency of interactions within a grouping event
	#obseff = chance of seeing these interactions

	#get indviduals in current grouping event
	gindivs=inds$indivs[gbirow==1]
	pdyads=data.frame(t(combn(gindivs,2)))
	
	#get the association of these dyads
	pdyads$assoc=sapply(1:nrow(pdyads),function (x) dyads$assoc[(dyads$names1==pdyads[x,1]&dyads$names2==pdyads[x,2])|(dyads$names2==pdyads[x,1]&dyads$names1==pdyads[x,2])])

	#NEXT: a number of binomial(?) trials depending on freq of interaction, assoc strength and observer effort?
}