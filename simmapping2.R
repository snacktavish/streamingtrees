library(phytools)

to.matrix<-function(x,seq) sapply(seq,"==",x)*1

trees<-read.nexus("all_noEPI_1.best.tre")
seqs<-as.data.frame(read.nexus.data("rbcL_noEPI_1.nex"))
seqs<- data.frame(lapply(seqs, as.character), stringsAsFactors=FALSE)
trees <- root(trees, "Pinus")


for (x in 1:10){
    ancseq=matrix('',nrow=trees$Nnode,ncol=dim(seqs)[1])
  for (i in 1:dim(seqs)[1]) {	
   xy = as.character(seqs[i,])
   names(xy)=colnames(seqs)
   PP<-to.matrix(xy,c('a','g','t','c'))
   for (ii in 1:dim(PP)[1]){
     if (sum(PP[ii,])==0){
     	PP[ii,]=0.25
     }
     }
   simtrees=make.simmap(trees,PP,nsim=1,model="ER")
   xyz=describe.simmap(simtrees)
   for (ii in 1:length(xyz$states)){
      ancseq[ii,i]=xyz$states[ii]	
   }  
   }
    rownames(ancseq)=names(xyz$states)
    totdat=cbind(seqs,t(ancseq))
   write.nexus.data(as.data.frame(totdat),paste("ancestorT",x,".nex",sep=''), format = "dna", datablock = TRUE,interleaved = T)
}

plotSimmap(simtrees,node.numbers=T)
fulltree=read.nexus("all_noEPI.best.tre")
fulltree <- root(fulltree, "Pinus")
plot.phylo(fulltree)

trees2$node.label<-names(xyz$states)
plot.phylo(trees,show.node.label=T)
write.tree(trees,"rep_1.tre")
