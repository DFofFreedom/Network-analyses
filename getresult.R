
exportdir="C:\\Users\\Julian\\Dropbox\\POSTDOC\\Network_generator\\2results"
toread=dir(exportdir)
rlist3=lapply(1:3,function(x){data.frame()})
for(parset in toread){
	resulttypes=dir(file.path(exportdir,parset))
	pars=strsplit(parset,"_")[[1]]
	pars=(t(data.frame(as.numeric(pars))))
	colnames(pars)=c("d.eff","i.dens","o.dens","m.i.eff","sex.eff","obs.eff.c","timesteps")
	
	for(i in 1:3){
		currpath=(file.path(exportdir,parset,resulttypes[i]))
		rlist1=lapply(file.path(currpath,dir(currpath)),read.csv)
		rlist2=lapply(rlist1,function(x){
			data.frame(parset=parset,pars,x)
		})
		rlist3[[i]]=rbind(rlist3[[i]],do.call(rbind,rlist2))
	}
	
}

#sort pval effect direction - does this need to happen for randomised assortivity too?
rlist3[[3]]$pval=ifelse(rlist3[[3]]$realeff>0,1-rlist3[[3]]$pval,rlist3[[3]]$pval)



for(i in 1:3){
	write.csv(rlist3[[i]],paste("all",resulttypes[i],".csv",sep=""))
}