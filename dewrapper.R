setwd("C:/Users/Julian/Dropbox/POSTDOC/Network_generator")
exportdir='test'
ptoread=dir(exportdir)
plist=strsplit(ptoread,"_")
ddf=do.call(rbind,lapply(plist,FUN=as.numeric))
colnames(ddf)=c("d.eff","i.dens","o.dens","m.i.eff","sex.eff","obs.eff","timesteps")
ddf=data.frame(filename=ptoread,ddf)

#check if these combos have data
ddf$DATA=F
for(filename in ptoread){
	typestoread=dir(file.path(exportdir,filename))
	typelist=as.list(rep(NA,6))
	names(typelist)=typestoread
	for(type in typestoread){
		repstoread=dir(file.path(exportdir,filename,type))
		if(length(repstoread)>0){
			ddf$DATA[ddf$filename==filename]=T
		}
	}
	

}

#example of selecting and reading in data
ddf2=ddf[ddf$d.eff==4&ddf$DATA==T,]

#Read data into this dataframe.
ddf2$data=I(as.list(rep(NA,nrow(ddf2))))
for(filename in ddf2$filename){
	typestoread=dir(file.path(exportdir,filename))
	typelist=as.list(rep(NA,6))
	names(typelist)=typestoread
	for(type in typestoread){
		repstoread=dir(file.path(exportdir,filename,type))
		repslist=as.list(rep(NA,length(repstoread)))
		for(csvfile in repstoread){
			repslist[[which(repstoread==csvfile)]]=read.csv(file.path(exportdir,filename,type,csvfile),check.names=F)
		}
		typelist[[which(typestoread==type)]]=repslist
	}
	
	ddf2$data[[which(ddf2$filename==filename)]]=typelist
}






