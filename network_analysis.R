##Load R packages that are required

library(asnipe)
library(ape)
library(igraph)
library(sna)
library(ergm)
library(ergm.count)
library(tnet)
library(assortnet)

do_analysis<-function(exportdir1,exportdir2,filename,nreps,startrep){
	#Prep export folders
	if(!dir.exists(file.path(exportdir2,filename))){
		dir.create(file.path(exportdir2,filename))
	}
	resultfolders=c("ergmresult","randconresult","randeffresult")
	for(folder in resultfolders){
		fulldir=file.path(exportdir2,filename,folder,sep="/")
		if(!dir.exists(fulldir)){
			dir.create(fulldir)
		}
	}
	
	for(currrep in (startrep+1):(startrep+nreps)){
		gbimat=get_gbimat(exportdir1,filename,currrep)
		indiv_dat=get_indiv_dat(exportdir1,filename,currrep)
		dist.centroids=get_dist_centroids(exportdir1,filename,currrep)
		obsnet=get_obs_net(exportdir1,filename,currrep)
		#just to be sure
		row.names(obsnet)=colnames(obsnet)
		row.names(dist.centroids)=colnames(dist.centroids)=row.names(obsnet)

		#remove any individuals that were not observed in any of the networks
		miss.obs=get_missing(obsnet)
		miss.gbi=get_missing(gbimat)

		obsnet=obsnet[!miss.obs,!miss.obs]
		gbimat=gbimat[,!miss.gbi]

		gbinet2=gbimat_to_net(gbimat)
		dist.centroids2=dist.centroids[!miss.gbi,!miss.gbi]
		nres=colSums(gbimat)

		obs.ergm=do_ergm(obsnet,indiv_dat,miss.obs,dist.centroids2)
		obs.ergm=data.frame(effect=row.names(obs.ergm),obs.ergm,type="OBS",rep=currrep)
		
		gbi.ergm=do_ergm(gbinet2,indiv_dat,miss.gbi,dist.centroids2,nres)
		gbi.ergm=data.frame(effect=row.names(gbi.ergm),gbi.ergm,type="GBI",rep=currrep)
		
		obs.rand=do_randomisations(obsnet,indiv_dat,miss.obs,10000)
		gbi.rand=do_randomisations(gbinet2,indiv_dat,miss.gbi,10000,cgbimat=gbimat)

		obs.rand[[1]]=data.frame(obs.rand[[1]],type="OBS",rep=currrep)
		obs.rand[[2]]=data.frame(obs.rand[[2]],type="OBS",rep=currrep)

		gbi.rand[[1]]=data.frame(gbi.rand[[1]],type="GBI",rep=currrep)
		gbi.rand[[2]]=data.frame(gbi.rand[[2]],type="GBI",rep=currrep)

		ergm.result=rbind(obs.ergm,gbi.ergm)
		rand.con.result=rbind(obs.rand[[1]],gbi.rand[[1]])
		rand.eff.result=rbind(obs.rand[[2]],gbi.rand[[2]])

		#Do exports
		write.csv(ergm.result,paste(exportdir2,"/",filename,"/ergmresult/",currrep,".csv",sep=""),row.names=F)
		write.csv(rand.con.result,paste(exportdir2,"/",filename,"/randconresult/",currrep,".csv",sep=""),row.names=F)
		write.csv(rand.eff.result,paste(exportdir2,"/",filename,"/randeffnresult/",currrep,".csv",sep=""),row.names=F)
	}
}

do_randomisations<-function(net,indiv_dat2,miss,nrand=1000,swapn=1,cgbimat=NA){
	#remove missing individual from indiv_dat
	indiv_dat2=indiv_dat2[!miss,]
	
	assort<-assortment.discrete(net,indiv_dat2$sex)$r
	wdeg<-colSums(net)
	net.mod<-glm(wdeg~indiv_dat2$sex,family=poisson)
	r.assorts<-rep(NA,nrand)
	r.effs<-rep(NA,nrand)

	if(!is.na(sum(cgbimat))){
		perm<-list(cgbimat,net)
		#start loop through datastream permutations
		for(ii in 1:nrand){
			random_networks_gbi<-network_swap(association_data=perm[[1]],association_matrix=perm[[2]], swaps=swapn) 
			newgbi<-random_networks_gbi[[2]]
			r.gbinet<-gbimat_to_net(random_networks_gbi[[2]])
			perm<-list(newgbi,r.gbinet)
			tmp.assort<-assortment.discrete(r.gbinet,indiv_dat2$sex)
			tmp.wdeg<-colSums(r.gbinet)
			tmp.mod<-glm(tmp.wdeg~indiv_dat2$sex,family=poisson)
			r.assorts[ii]<-tmp.assort$r
			r.effs[ii]<-coef(tmp.mod)[2]
		}
	} else {
		rnet=net
		for(ii in 1:nrand){
 			#node swaps
 			rnet<-rmperm(rnet) 
			tmp.assort<-assortment.discrete(rnet,indiv_dat2$sex)
 			tmp.wdeg<-colSums(rnet)
			tmp.mod<-glm(tmp.wdeg~indiv_dat2$sex,family=poisson)
 			r.effs[ii]<-coef(tmp.mod)[2]
 			r.assorts[ii]<-tmp.assort$r
 		}
	}
	return(list(data.frame(realassort=assort,
		data.frame(t(quantile(r.assorts,c(0.005,0.025,0.975,0.995),na.rm=T))),
		pval=sum(r.assorts<assort)/(nrand+1)),
		data.frame(realeff=coef(net.mod)[2],
		data.frame(t(quantile(r.effs,c(0.005,0.025,0.975,0.995),na.rm=T))),
		pval=sum(r.effs<coef(net.mod)[2])/(nrand+1)))
	)
}

do_ergm<-function(net,indiv_dat2,miss,dist.centroids2,nres2=NA){
	#remove missing individual from indiv_dat
	indiv_dat2=indiv_dat2[!miss,]
	
	#set up observed network as a network object
	NET2.edgelist<-as.tnet(as.matrix(net))
	NET2<-network(NET2.edgelist, directed=F, ignore.eval=F, names.eval="weight")

	#add attributes
	set.vertex.attribute(NET2,"group",as.vector(indiv_dat2$groups))
	set.vertex.attribute(NET2,"sex",as.vector(indiv_dat2$sex))
	if(!is.na(sum(nres2))){
		set.vertex.attribute(NET2,"nres",as.vector(nres2))
	}

	#shared group matrix - Julian likes vectorisation version
	mat1=matrix(rep(indiv_dat2$groups,each=ncol(net)),nrow=nrow(net))
	mat2=matrix(rep(indiv_dat2$groups,each=ncol(net)),nrow=nrow(net),byrow=T)
	sh.gr=mat1==mat2
	sh.gr[sh.gr]=1

	#sh.gr<-array(NA,dim=dim(net))
	#for(i in 1:ncol(sh.gr)){
 	#	for(j in 1:nrow(sh.gr)){
  	#		ifelse(indiv_dat$groups[indiv_dat$indivs==colnames(net)[i]]==indiv_dat$groups[indiv_dat$indivs==colnames(net)[j]],sh.gr[i,j]<-1,sh.gr[i,j]<-0)
 	#	}
	#}
	
	diag(sh.gr)<-0
	
	#distance between groups matrix - Julian likes vectorisation still
	mat1x=matrix(rep(indiv_dat2$x,each=ncol(net)),nrow=nrow(net))
	mat1y=matrix(rep(indiv_dat2$y,each=ncol(net)),nrow=nrow(net))
	mat2x=matrix(rep(indiv_dat2$x,each=ncol(net)),nrow=nrow(net),byrow=T)
	mat2y=matrix(rep(indiv_dat2$y,each=ncol(net)),nrow=nrow(net),byrow=T)

	dist.gr=((mat1x-mat2x)^2+(mat1y-mat2y)^2)^0.5
	
	#dist.gr<-array(NA,dim=dim(obsnet))
	#for(i in 1:ncol(dist.gr)){
 	#	for(j in 1:nrow(dist.gr)){
  	#	dist.gr[i,j]<-dist(rbind(c(indiv_dat$x[indiv_dat$indivs==colnames(obsnet)[i]],indiv_dat$y[indiv_dat$indivs==colnames(obsnet)[i]]),c(indiv_dat$x[indiv_dat$indivs==colnames(obsnet)[j]],indiv_dat$y[indiv_dat$indivs==colnames(obsnet)[j]])))
 	#	}
	#}

	if(is.na(sum(nres2))){
		modB_2<-ergm(NET2~sum+nonzero+nodefactor("sex")+nodematch("sex")+edgecov(sh.gr)+edgecov(dist.gr),reference=~Poisson,response="weight")
	} else {
		modB_2<-ergm(NET2~sum+nonzero+nodefactor("sex")+nodecov("nres")+nodematch("sex")+edgecov(dist.centroids2),reference=~Poisson,response="weight")
	}
	resultdf=summary(modB_2)$coefs
	return(resultdf)

}



get_missing=function(net){
	miss=(colSums(net)==0)

	return(miss)
}



get_gbimat<-function(exportdir1,filename,currrep){
	currfp=file.path(exportdir1,filename,"obsgbimat",paste(currrep,".csv",sep=""))
	gbimat<-read.csv(currfp,check.names=FALSE)
	return(gbimat)
}

get_indiv_dat<-function(exportdir1,filename,currrep){
	currfp=file.path(exportdir1,filename,"popdat",paste(currrep,".csv",sep=""))
	indiv_dat<-read.csv(currfp,check.names=FALSE)
	return(indiv_dat)
}

gbimat_to_net<-function(gbimat){
	gbinetT= t(as.matrix(gbimat)) %*% as.matrix(gbimat)
	diag(gbinetT)=0
	return(gbinetT)
}

get_obs_net<-function(exportdir1,filename,currrep){
	currfp=file.path(exportdir1,filename,"obsnet",paste(currrep,".csv",sep=""))
	obsnet<-read.csv(currfp,check.names=FALSE)
	return(obsnet)
}

get_dist_centroids<-function(exportdir1,filename,currrep){
	#this works out home range centroids for each individual and distances between them
	currfp=file.path(exportdir1,filename,"obsgbigroups",paste(currrep,".csv",sep=""))
	grouplocs1<-read.csv(currfp,check.names=FALSE)
	names(grouplocs1)<-"Gr"
	gbimat=get_gbimat(exportdir1,filename,currrep)

	currfp=file.path(exportdir1,filename,"popdat",paste(currrep,".csv",sep=""))
	indiv_dat<-read.csv(currfp,check.names=FALSE)


	grouplocs1$x<-rep(NA,nrow(grouplocs1))
	grouplocs1$y<-rep(NA,nrow(grouplocs1))
	
	grouplocs1$x=sapply(1:nrow(grouplocs1),function (i){
		indiv_dat$x[min(which(indiv_dat$groups==grouplocs1$Gr[i]))]
	})
	grouplocs1$y=sapply(1:nrow(grouplocs1),function (i){
		indiv_dat$y[min(which(indiv_dat$groups==grouplocs1$Gr[i]))]
	})
	
	indiv.groups=lapply(1:nrow(indiv_dat),function (i){
		which(gbimat[,i]==1)
	})

	indiv.centroids<-data.frame(indiv_dat$indivs)
	names(indiv.centroids)<-"id"

	indiv.centroids$x=sapply(1:nrow(indiv.centroids),function (i){
		mean(grouplocs1$x[indiv.groups[[i]]])
	})

	indiv.centroids$y=sapply(1:nrow(indiv.centroids),function (i){
		mean(grouplocs1$y[indiv.groups[[i]]])
	})

	dist.centroids<-as.matrix(dist(indiv.centroids[,2:3]))
	return(dist.centroids)
}

