setwd("C:/Users/Julian/Dropbox/POSTDOC/Network_generator")
source("network_generator.R")

simulated.networks=network.generator(5, 10, 20,
                  0.9, 0.2, 0.8,
                  sex.eff=0.1, m.i.eff=NA, m.o.eff=NA,
                  plot=T)

obs.sim.networks=networkobs(simulated.networks,
           timesteps = 1000, obseff = 1, intfreq = 10,
           floaterprob=0.01,probnorm=NA)

#difference between observed network and true network
sum(abs(obs.sim.networks$obsnetwork-simulated.networks$network))

#plot this for fun (might want to force they layout algorithm so as to have comparable plots)
dev.new()
par(mfrow=c(1,1))

net2.d<-graph.adjacency(obs.sim.networks$obsnetwork,mode="undirected",weighted=TRUE,diag=FALSE)
net2.d<-set.vertex.attribute(net2.d, "group", index=V(net2.d), simulated.networks$ind_data$groups)
net2.d<-set.vertex.attribute(net2.d, "sex", index=V(net2.d), simulated.networks$ind_data$sex)	   

l=layout_(net2.d, nicely())	
V(net2.d)$color=V(net2.d)$group #assign the "Group" attribute as the vertex color
plot(net2.d,edge.width=(E(net2.d)$weight)^0.25,vertex.shape=c("circle","square")[factor(V(net2.d)$sex)],layout=l,vertex.size=8,vertex.label=NA,margin=c(0,0,0,0))


#Generate network from the gbi
library(asnipe)
gbinet=get_network(obs.sim.networks$truegbi, data_format = "GBI",identities=colnames(obs.sim.networks$truegbi))
#need to have sim network be comparable to make this comparions

library(ape)
mantel.test(simulated.networks$network,gbinet,graph=T)

net2.d<-graph.adjacency(gbinet,mode="undirected",weighted=TRUE,diag=FALSE)
net2.d<-set.vertex.attribute(net2.d, "group", index=V(net2.d), simulated.networks$ind_data$groups)
net2.d<-set.vertex.attribute(net2.d, "sex", index=V(net2.d), simulated.networks$ind_data$sex)	   
dev.new()
V(net2.d)$color=V(net2.d)$group #assign the "Group" attribute as the vertex color
plot(net2.d,edge.width=(E(net2.d)$weight)^2.5,vertex.shape=c("circle","square")[factor(V(net2.d)$sex)],layout=l,vertex.size=8,vertex.label=NA,margin=c(0,0,0,0))


#TESTING INTERACTION NETWORK
g=graph.data.frame(obs.sim.networks$interactions[,1:2],directed=FALSE)
g2=get.adjacency(g,sparse=FALSE)
#note that I've never used this kind of data, so there may be a better
#way of doing this?

net2.d=graph.adjacency(g2,mode="undirected",weighted=TRUE)
net2.d<-set.vertex.attribute(net2.d, "group", index=V(net2.d), simulated.networks$ind_data$groups[match(V(net2.d),simulated.networks$ind_data$indivs)])
net2.d<-set.vertex.attribute(net2.d, "sex", index=V(net2.d), simulated.networks$ind_data$sex[match(V(net2.d),simulated.networks$ind_data$indivs)])	   
dev.new()
V(net2.d)$color=V(net2.d)$group #assign the "Group" attribute as the vertex color
plot(net2.d,edge.width=(E(net2.d)$weight)^0.1,vertex.shape=c("circle","square")[factor(V(net2.d)$sex)],layout=l[match(V(net2.d),colnames(obs.sim.networks$obsnetwork)),],vertex.size=8,vertex.label=NA,margin=c(0,0,0,0))


pop.dat=simulated.networks
timesteps=100
obseff=1
intfreq=20