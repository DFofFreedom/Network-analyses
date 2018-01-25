library(asnipe)
library(ape)

setwd("C:/Users/matth/Desktop/Dave_Julian simulations folder/Network-analyses-jevansbranch")
source("network_generator.R")
#all combos of variables below

groups = 10
mean.group.size = 10
max.group.size=10
d.eff=c(0,4,8)
i.dens=c(0.4,0.8,1.2)
o.dens=c(0.4,0.2,0.1)
m.i.eff=c(-0.6,-0.4,-0.2,0,0.2,0.4,0.6)
m.o.eff=NA
sex.eff=c(0.25,0.5,0.75,1,0.25,1.5,1.75,2)
obs.eff=c(0.2,0.4,0.6,0.8,1)
timesteps=c(100,1000)
intfreq=10
floaterprob=0.01 #?
probnorm=NA


number_of_cores=2
number_of_cores=detectCores()#probably not reccomended for all as this uses ALL available cores.
exportdir='test'
nreps=2

iterdf1=expand.grid(groups=groups,mean.group.size=mean.group.size,max.group.size=max.group.size,
	d.eff=d.eff,i.dens=i.dens,o.dens=o.dens,m.i.eff=m.i.eff,m.o.eff=m.o.eff,sex.eff=sex.eff,
	obs.eff=obs.eff,timesteps=timesteps,intfreq=intfreq,floaterprob=floaterprob,probnorm=probnorm,
	nreps=nreps,exportdir=exportdir)

#combine combos into a data frame and make it an iterator
library(iterators)
iterdf=iter(iterdf1[1:10,],by = "row")

require(doParallel)
require(parallel)
registerDoParallel(cores=number_of_cores)
require(foreach)	

foreach(i = iterdf)%dopar%{
do.call(do_networksim,i)
}
#just a note when the iterator is done, you must call it again to run the loop again

#do.call(do_networksim,nextElem(iterdf))
do.call(do_networksim,as.list(iterdf1[3,]))






