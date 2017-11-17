setwd("C:/Users/Julian/Dropbox/POSTDOC/Network_generator")
source("network_generator.R")

simulated.networks=network.generator(5, 10, 20,
                  0.9, 0.2, 0.8,
                  sex.eff=0.1, m.i.eff=NA, m.o.eff=NA,
                  plot=T)

obs.sim.networks=networkobs(simulated.networks,
           timesteps = 100, obseff = 1, intfreq = 10,
           probnorm=NA)

#difference between observed network and true network
sum(abs(obs.sim.networks$obsnetwork-simulated.networks$network))

#plot this for fun (might want to force they layout algorithm so as to have comparable plots)
dev.new()
par(mfrow=c(1,1))
net2.d<-graph.adjacency(obs.sim.networks$obsnetwork,mode="undirected",weighted=TRUE,diag=FALSE)
net2.d<-set.vertex.attribute(net2.d, "group", index=V(net2.d), simulated.networks$ind_data$groups)
net2.d<-set.vertex.attribute(net2.d, "sex", index=V(net2.d), simulated.networks$ind_data$sex)	   
		
V(net2.d)$color=V(net2.d)$group #assign the "Group" attribute as the vertex color
plot(net2.d,edge.width=(E(net2.d)$weight)^0.25,vertex.shape=c("circle","square")[factor(V(net2.d)$sex)],layout=layout.fruchterman.reingold(net2.d),vertex.size=8,vertex.label=NA,margin=c(0,0,0,0))


#Generate network from the gbi
library(asnipe)
gbinet=get_network(obs.sim.networks$truegbi, data_format = "GBI",identities=colnames(obs.sim.networks$truegbi))
#need to have sim network be comparable to make this comparions


sum(abs(gbinet-normnet))


dev.new()
par(mfrow=c(1,1))
net2.d<-graph.adjacency(gbinet,mode="undirected",weighted=TRUE,diag=FALSE)
net2.d<-set.vertex.attribute(net2.d, "group", index=V(net2.d), simulated.networks$ind_data$groups)
net2.d<-set.vertex.attribute(net2.d, "sex", index=V(net2.d), simulated.networks$ind_data$sex)	   
		
V(net2.d)$color=V(net2.d)$group #assign the "Group" attribute as the vertex color
plot(net2.d,edge.width=(E(net2.d)$weight)^0.25,vertex.shape=c("circle","square")[factor(V(net2.d)$sex)],layout=layout.fruchterman.reingold(net2.d),vertex.size=8,vertex.label=NA,margin=c(0,0,0,0))

