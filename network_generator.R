#SYMMAT UTILITY FUNC
symmat = function (bool,dyads){
	return(rbind(dyads[bool,],dyads[bool,c(2,1)]))
}


###NETWORK GENERATION FUNCTION####
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
	###SETUP BASE POPULATION####
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

	###SPATIAL LOCATION OF GROUPS/CLUSTERS####

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
	ismale=dsex[,1]=="M"|dsex[,2]=="M"
	

	symmat = function (bool,dyads){
		return(rbind(dyads[bool,],dyads[bool,c(2,1)]))
	}
	
	#WITHIN GROUP EDGES####

	#FF
	net.d[symmat(samesite&!ismale,dyads)]=sapply(which(samesite&!ismale,dyads),function (x) round(rnbinom(1,size=i.dens*(distsv[x])^d.eff,prob=0.3))+
		round(rnbinom(1,size=i.dens*(distsv[x])^d.eff,prob=0.3)))
	#MF
	net.d[symmat(samesite&ismale&!samesex,dyads)]=sapply(which(samesite&ismale&!samesex),function (x) round(rnbinom(1,size=i.dens*(distsv[x])^d.eff,prob=0.3))+
		round(rnbinom(1,size=m.i.eff+(i.dens*(distsv[x])^d.eff),prob=0.3)))
	#MM
	net.d[symmat(samesite&ismale&samesex,dyads)]=sapply(which(samesite&ismale&samesex),function (x) round(rnbinom(1,size=m.i.eff+(i.dens*(distsv[x])^d.eff),prob=0.3))+
		round(rnbinom(1,size=m.i.eff+(i.dens*(distsv[x])^d.eff),prob=0.3)))
	
	#OUTSIDE GROUP EDGES####

	#FF
	net.d[symmat(!samesite&!ismale,dyads)]=sapply(which(!samesite&!ismale,dyads),function (x) round(rnbinom(1,size=o.dens*(distsv[x])^d.eff,prob=0.3))+
		round(rnbinom(1,size=o.dens*(distsv[x])^d.eff,prob=0.3)))
	#MF
	net.d[symmat(!samesite&ismale&!samesex,dyads)]=sapply(which(!samesite&ismale&!samesex),function (x) round(rnbinom(1,size=o.dens*(distsv[x])^d.eff,prob=0.3))+
		round(rnbinom(1,size=m.o.eff+(o.dens*(distsv[x])^d.eff),prob=0.3)))
	#MM
	net.d[symmat(!samesite&ismale&samesex,dyads)]=sapply(which(!samesite&ismale&samesex),function (x) round(rnbinom(1,size=m.o.eff+(o.dens*(distsv[x])^d.eff),prob=0.3))+
		round(rnbinom(1,size=m.o.eff+(o.dens*(distsv[x])^d.eff),prob=0.3)))	
	
	#SEX HOMOPHILY EFFECT####
	net.d[symmat(samesex,dyads)]= round(net.d[symmat(samesex,dyads)]*sex.eff)

	diag(net.d)=0
	
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

###NETWORK SAMPLING FUNCTION####		      
		      
networkobs<-function(pop.dat,timesteps,obseff,intfreq,probnorm=NA){
	####
	#pop.dat - pop.dat object produced by network_generator
	#timesteps - number of timesteps
	#obseff - number from 0 to 1 representing how much is seen
	#intfreq - how many interactions are seen in a timestep
	#probnorm - a value to normalise associations
	# 		if not given is obtained from the network being considered
	####

	#get objects
	inds=pop.dat$ind_data
	network=pop.dat$network
	distmat=pop.dat$distmat
	
	#data frame of dyads
	dyads=which(upper.tri(network),arr.ind=T)
	assocs=network[dyads]
	names1=as.numeric(row.names(network)[dyads[,1]])
	names2=as.numeric(colnames(network)[dyads[,2]])
	dyads=data.frame(names1,names2,dyads,assocs=assocs)

	#observed contact network - adds noise based on observation effort. 
	#network
	obsnetwork=network
	obsnetwork[symmat(T,as.matrix(dyads[,3:4]))]=sapply(dyads$assoc,function (x) assocnoise(x,obseff))
	
	#generate GBI mat. Select a random individual - check chance of being seen with other
	#indiv based on their associations. Repeat this for any individual added.
	gbimat=t(sapply(1:timesteps,function (x) makeevent(inds,dyads,probnorm)))
	colnames(gbimat)=inds$indivs

	#randomly see only certain groups based on observation effort (could make it so larger groups are more likely to be seen - but, simple for now)
	obsgbimat=gbimat[sort(sample(1:timesteps,round(obseff*timesteps))),]

	#similarly we could throw away individuals from the GBI here based on their sociality - but keeping it simple for now

	
	#get observed interactions (for interaction based network)
	interactions=lapply(1:nrow(obsgbimat),function (x) data.frame(getinteractions(inds,dyads,gbimat[x,],intfreq,obseff),timestep=x))
	interactions=do.call(rbind,interactions)
	interactions$sex1=inds$sex[match(interactions$name1,inds$indivs)]	
	interactions$sex2=inds$sex[match(interactions$name2,inds$indivs)]	
	return(list(truegbimat=gbimat,obsgbi=obsgbimat,obsnetwork=obsnetwork,interactions=interactions))
}

assocnoise<-function(x,obseff){
	if(x==0){
		return(0)
	} else {
		xvec=seq(0,1,length.out=x+1)
		probs=(-(obseff-1)*(obseff*(xvec-1)+1))^(2-(2*xvec))
		obsassoc=sample(0:x,1,1,prob=probs)
		return(obsassoc)
	}
}

makeevent<-function(inds,dyads,probnorm=NA,verbose=F){
	if(is.na(probnorm)){
		probnorm=max(dyads$assoc)
	}
	#randomly pick an individual to start
	gbirow=matrix(0,1,length(inds$indivs))
	colnames(gbirow)=inds$indivs
	seed=sample(inds$indivs,1)
	seeddegree=abs(round(rnorm(1,sum(focalpotentials(seed,dyads)$assocs>0),2)))
	gbirow[colnames(gbirow)==seed]=1
	#scalev=(seeddegree-(sum(gbirow)-1))/probnorm
	scalev=1
	#not usins scalev as it stands. Far too stingy with group sizes!
	if(verbose==T){
		print(paste("seed is",seed))
		print(paste("seed degree is",seeddegree))
		print(paste("groupsize probability scaling factor is",scalev))
	}	
	fassoc=focalassoc(seed,dyads,gbirow,scalev,probnorm,T)
	todo=fassoc
	checked=seed
	if(verbose==T){
		print(paste("todo:",paste(todo,collapse=" ")))
		print(paste("checked:",paste(checked,collapse=" ")))
	}	
	gbirow[colnames(gbirow)%in%fassoc]=1
	
	#check the associations of the individuals we are adding until we are no longer adding members
	while(length(todo)>0){
		#scalev=(seeddegree-(sum(gbirow)-1))/probnorm
		scalev=1
		#randomise order of todo so we're not biased by the order they were added to matrix
		todo=todo[sample(1:length(todo))]
		#ISSUE: I love sapply, but this means that several members can be added to the group at once
		#without considering the size of the group or the members being added at the same time.
		#fassoc=unlist(sapply(todo,function (x) focalassoc(x,dyads,gbirow,scalev,probnorm)))
		if(verbose==T){
			print(paste("groupsize:",sum(gbirow)))
			print(paste("groupsize probability scaling factor is",scalev))
			print(paste("focal individual:",todo[1]))
		}	
		
		fassoc=focalassoc(todo[1],dyads,gbirow,scalev,probnorm)
		gbirow[colnames(gbirow)%in%fassoc]=1		
		checked=c(checked,todo[1])
		todo=c(todo[!todo%in%checked],unique(fassoc[!fassoc%in%checked&!fassoc%in%todo]))
		if(verbose==T){
			print(paste("todo:",paste(todo,collapse=" ")))
			print(paste("checked:",paste(checked,collapse=" ")))
		}	

	}
	return(gbirow)	
}

focalpotentials=function(focal,dyads){
	potentials=dyads[dyads$names1==focal|dyads$names2==focal,]
	#swap the names around for readability
	potentials$names2[potentials$names1!=focal]=potentials$names1[potentials$names1!=focal]
	potentials$names1[potentials$names1!=focal]=focal
	return(potentials)
}


focalassoc<-function(focal,dyads,currevent,scalev,probnorm,forcemulti=F){
	####
	# focal = ID of individual whose associations we are using
	# dyads = all possible dyads
	# forcemulti = Force function to return at least one interaction?
	# probnorm = value used for the normalisation of associations probabilities - if not given is set to 
	# 	
	####
	

	
	if(is.na(probnorm)){
		probnorm=max(dyads$assoc)
	}
	
	potentials=focalpotentials(focal,dyads)
	
	eventmembers=colnames(currevent)[currevent==1&(colnames(currevent)!=focal)]
	probs=potentials$assoc/(probnorm+1)
	potids=potentials$names2[probs>0]
	#consider these potentials vs existing group members - if they have a zero association with 
	#an individual already in the group reduce the probability of them being in the GBI
	probs[probs>0][sapply(potids, function (i) {potassocs=focalpotentials(i,dyads);sum(potassocs[potassocs$assocs==0,"names2"]%in%eventmembers)>0})]=0.001
	#probabilities are scaled depending on the sociality of the seed
	probs=probs*scalev
	probs=probs^2
	if(forcemulti==T){
		gbirow=rep(0,nrow(potentials))
		while(sum(gbirow)==0){

			gbirow=rbinom(potentials$assoc, 1, probs)
		}
	} else {
		gbirow=rbinom(potentials$assoc, 1,probs )

	}
	return(potentials$names2[gbirow==1])
}

getinteractions<-function(inds,dyads,gbirow,intfreq,obseff){
	#inds dataframe
	#dyads dataframe
	#gbirow = a row of a group by individual matrix (so we can apply this to one event at a time)
	#intfreq = what is the frequency of interactions within a grouping event
	#obseff = chance of seeing these interactions

	#get indviduals in current grouping event
	gindivs=inds$indivs[gbirow==1]
	pdyads=data.frame(t(combn(gindivs,2)))

	#get the association of these dyads
	pdyads$assoc=sapply(1:nrow(pdyads),function (x) dyads$assoc[(dyads$names1==pdyads[x,1]&dyads$names2==pdyads[x,2])|(dyads$names2==pdyads[x,1]&dyads$names1==pdyads[x,2])])
	#convert these into probabilities (so they add up to 1)
	pdyads$prob=pdyads$assoc/sum(pdyads$assoc)
	#get all interactions that have occured in this grouping event
	#individuals with higher association have higher probability of being observed interacting
	allinters=sample(1:nrow(pdyads),intfreq,1,prob=pdyads$prob)
	
	#randomly sample this based on the observation effort
	seenindex=sample(allinters,intfreq*obseff)
	seeninters=pdyads[seenindex,1:2]
	names(seeninters)[1:2]=c("name1","name2")	
	return(seeninters)
}
