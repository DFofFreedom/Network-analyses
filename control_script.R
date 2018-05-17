source("network_generator.R")
source("network_analysis.R")
checkdata=function(exportdir,currpar2){
	if(dir.exists(as.character(exportdir))){
		#directory does not exist, code definitely not run	
		return(data.frame(DATA=T,donereps=0))
	}
	
	#Directory exists - are there the correct number of folders
	ptoread=dir(exportdir)
	plist=strsplit(ptoread,"_")
	ddf=do.call(rbind,lapply(plist,FUN=as.numeric))
	colnames(ddf)=c("d.eff","i.dens","o.dens","m.i.eff","sex.eff","obs.eff","timesteps")
	ddf=data.frame(filename=ptoread,ddf)

	if(nrow(ddf)!=nrow(currpar2)){
		#wrong number of subfolders, code needs to be rerun
		return(data.frame(DATA=T,donereps=0))
	}

	#check if these combos have data
	ddf$DATA=F
	ddf$donereps=0
	typestoread=dir(file.path(exportdir,ptoread[1]))
	ddf2=data.frame(type=rep(typestoread,each=nrow(ddf)),ddf)
	for(filename in ptoread){

		for(type in typestoread){
			repstoread=dir(file.path(exportdir,filename,type))
			#if we have not done all the reps, need to run this set
			if(length(repstoread)==currpar2$nreps[1]){
				
				ddf2$DATA[ddf$filename==filename&ddf$type==type]=T
			}
			ddf2$donereps=length(repstoread)
		}
		
	}
	return(data.frame(DATA=any(ddf2$DATA==F),donereps=min(ddf2$donereps)))
}

parameters=read.csv("parametersets.csv")
job<- as.numeric(commandArgs(trailingOnly = TRUE)[1])
exportdir1=paste(job,"nets")
exportdir2=paste(job,"results")
currpar=parameters[job,]

#get a vector of obs eff
obs.effvec=as.numeric(strsplit(as.character(currpar$obs.eff)," ")[[1]])

#version of currpar with all obs.effs for checking for completion
currpar2=do.call(rbind,(lapply(obs.effvec,function (i) {
	cpar1=currpar
	cpar1$obs.eff=i
	return(cpar1)
})))

#First, check we have not already run this parameter set
#Useful if we need to restart the jobs
outputcheck=checkdata(exportdir1,currpar2)
#If we have not already generated this parameter set or not finished it
if(outputcheck$DATA){
	#generate networks
	currpar3=currpar
	currpar3$nreps=currpar$nreps-outputcheck$donereps
	currpar3$startrep=outputcheck$donereps
	currpar3$exportdir=exportdir1
	currpar3=as.list(currpar3)
	currpar3$obs.eff=obs.effvec
	do.call(do_networksim,currpar3)
}

#if networks have been simulated, check if the analysis files exist
outputcheck2=checkdata(exportdir2,currpar2)
if(outputcheck2$DATA){
	if(!dir.exists(as.character(exportdir2))){
		dir.create(as.character(exportdir2))
	}
	#analysis not yet done, we need to do it
	toanalyse=dir(exportdir1)
	for(filename in toanalyse){
		nreps=currpar$nreps-outputcheck$donereps
		startrep=outputcheck$donereps
		do.call(do_analysis,list(exportdir1,exportdir2,filename,nreps,startrep))
	}
}

#Done
