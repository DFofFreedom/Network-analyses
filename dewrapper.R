setwd("C:/Users/jevan8/Dropbox/POSTDOC/Network_generator")
exportdir='test'
ptoread=dir(exportdir)
plist=strsplit(ptoread,"_")
ddf=do.call(rbind,lapply(plist,FUN=as.numeric))
colnames(ddf)=c("d.eff","i.dens","o.dens","m.i.eff","sex.eff","obs.eff","timesteps")
ddf=data.frame(filename=ptoread,ddf)


datalist=as.list(rep(NA,nrow(ddf)))

for(filename in ptoread){
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
	
	datalist[[which(ddf$filename==filename)]]=typelist
}






