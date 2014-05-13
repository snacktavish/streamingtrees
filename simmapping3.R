library(phytools)
library(seqinr)

treefile<-"fig.best.tre"
seqfile<-"fig.fst"
prefix<-"figAnc"

to.matrix<-function(x,seq) sapply(seq,"==",x)*1

trees<-read.nexus(treefile)
alig<-read.alignment(seqfile,"fasta")
seqs=matrix(nrow=length(alig$seq),ncol=length(s2c(alig$seq[[x]][1])))
rownames(seqs)<- alig$nam
for (x in 1:length(alig$seq)){
  seqs[x,]=s2c(alig$seq[[x]][1])
}

#trees <- root(trees, "Pinus")
for (x in 1:2){
  ancseq=matrix('',nrow=trees$Nnode,ncol=dim(seqs)[2])
  for (site in 1:dim(seqs)[2]) {	
   xy = seqs[,site]
   PP<-to.matrix(xy,c('a','g','t','c'))
   for (ii in 1:dim(PP)[1]){
     if (sum(PP[ii,])==0){
     	PP[ii,]=0.25
     }
     }
   simtrees=make.simmap(trees,PP,nsim=1,Q="mcmc")
   xyz=describe.simmap(simtrees)
   for (node in 1:length(xyz$states)){
      ancseq[node,site]=xyz$states[node]	
   }  
   
  }
    rownames(ancseq)=names(xyz$states)
    totdat=rbind(seqs,ancseq)
  write.fasta(as.data.frame(t(totdat)),names=rownames(totdat),file.out="ancseq.fas")
}


plotSimmap(simtrees,node.numbers=T,)
fulltree=read.nexus("all_noEPI.best.tre")
fulltree <- root(fulltree, "Pinus")
plot.phylo(fulltree)

trees2$node.label<-names(xyz$states)
plot.phylo(trees,show.node.label=T)
write.tree(trees,"rep_1.tre")


getDescendants<-function(tree,node,curr=NULL){
  if(is.null(curr)) curr<-vector()
  daughters<-tree$edge[which(tree$edge[,1]==node),2]
  curr<-c(curr,daughters)
  w<-which(daughters>=length(tree$tip))
  if(length(w)>0) for(i in 1:length(w))
    curr<-getDescendants(tree,daughters[w[i]],curr)
  return(curr)
}