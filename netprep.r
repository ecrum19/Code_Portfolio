#!/usr/bin/env Rscript
library("optparse")
library("igraph")
library("Matrix")
 
option_list = list(
	make_option(c("-d", "--panpath"), type="character", default=NULL, 
              help="path to pangenome directory", metavar="character"),
	make_option(c("-T", "--threads"), type="double", default=2,
			  help="number of theads to use when possible [default= %default]", metavar="double")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

netprep<-function(panpath,threads)
{
	setwd(panpath)
	filelist=list.files()
	panpos=grep("PAN.db",filelist)
	if(length(panpos)==0)
	{
		print("No pangenome output recognized. Please ensure that your pangenome database file ends with PAN.db")
		break
	}
	if(length(panpos)!=0)
	{
		panfile=filelist[panpos[1]]
		print(paste("pangenome database file found. ",panfile," will be used in downstream analysis. If this is not the correct file, please cancel the script and ensure only one PAN.db file is present in the directory.",sep=''))
			
		print("Exporting anvio information table")
		system(paste('anvi-export-table ',panfile,' --table self -o self.txt',sep=''))
		info=readLines('self.txt')
		infosplit=strsplit(info,'\t')
		external=TRUE
		if(length(infosplit[[3]])==2){external=FALSE}

		print("Relabeling sequence headers and building network edgelist")
		anvirelabel(panpath,external)
		mcl2edgelist('mcl-clusters-relabeled.txt')
	}
}

anvirelabel<-function(panpath,external=TRUE)
{
	setwd(panpath)
	#get the original genomes in same order as hashes
	dbinfo=readLines('self.txt')				#assumes the output of anvi-export-table is in the panpath directory and named "self.txt"
	if(external==TRUE){genomelist=strsplit(dbinfo[4],'\t')[[1]][2]}		#pulls out the comma separated list of genomes from the anvio db list.
	if(external==FALSE){genomelist=strsplit(dbinfo[3],'\t')[[1]][2]}	#uses the appropriate row if from internal genome db and not external.
	genomelist=unlist(strsplit(genomelist,','))
	
	#Extract hash names
	allseqs=readLines("combined-aas.fa")	#read in every sequence (with hash labels)
	allaccpos=grep(">",allseqs)
	allaccs=allseqs[allaccpos]
	allaccs=unique(unlist(strsplit(allaccs,'>')))[-1]
	allseqsvec=allseqs[-allaccpos]
	names(allseqsvec)=allaccs				#now we have the sequences named by their accessions, instead of in fasta format. This will create some indexing efficiency later on and avoid a nested loop.
	hashlist=character()
	accsplit=strsplit(allaccs,'_')
	for(i in 1:length(accsplit)){hashlist[i]=paste(accsplit[[i]][1:(length(accsplit[[i]])-1)],collapse="_")}
	unihashlist=unique(hashlist)
		
	#make a final table matching hashes to genome names
	hashmat=cbind(unihashlist,genomelist)
	colnames(hashmat)=c("hash","genome")
	rownames(hashmat)=hashmat[,1]
	names(genomelist)=unihashlist
	names(unihashlist)=genomelist		#now, we can index the hashes by the genomes and vice versa. This is faster to access than using the whole hash matrix
	#Relabel the sequence file
	newaccs=character()
	for(i in 1:length(accsplit)){newaccs[i]=paste(hashmat[accsplit[[i]][1],2],accsplit[[i]][2],sep="_")}
	#allseqs[accpos]=paste(">",newaccs,sep='')
	newseqvec=allseqsvec
	names(newseqvec)=newaccs
	newseqs=allseqs
	newseqs[allaccpos]=paste(">",newaccs,sep='')
	write.table(newseqs,'combined-aas-relabeled.fa',row.names=F,col.names=F,quote=F)
	
	#Write out MCL clusters as fasta files using the genome info, instead of hashes
	clustlist=strsplit(readLines("mcl-clusters.txt"),"\t")
	newmcl=character()
	for(i in 1:length(clustlist))			#this is slow partly because it's a nested for loop, but really because the allseqs object is just a big character vector. Could speed this up (and avoid that second loop) if restructure allseqs so that the deflines are index names for sequences.
	{
		tempn=2*length(clustlist[[i]])		#total number of positions for a new character vector with the cluster
		tempaccpos=match(clustlist[[i]],allaccs)		#gives index of the accessions in the accpos vector for each cluster gene. Takes advantage of original accpos being in hashes
		tempmcl=names(newseqvec[tempaccpos])
		newmcl=c(newmcl,paste(tempmcl,collapse='\t'))
	}
	write.table(newmcl,'mcl-clusters-relabeled.txt',row.names=F,col.names=F,quote=F)
}

acc2genome<-function(accs)
{
	accsplit=strsplit(accs,"_")
	genomes=character()
	for(i in 1:length(accsplit))
	{
		genomes[i]=paste(accsplit[[i]][1:length(accsplit[[i]])-1],collapse='_')
	}
	return(genomes)
}

mcl2edgelist<-function(mcl)
{
	clust=strsplit(readLines(mcl),'\t')
	gcvec=paste('gc_',c(0:(length(clust)-1)),sep='')
	gclist=list()
	for(i in 1:length(clust))
	{
		gclist[[i]]=rep(gcvec[i],length(clust[[i]]))
	}
	clustgenomes=lapply(clust,acc2genome)
	presedge=cbind(unlist(clustgenomes),unlist(gclist))
	genomelist=unique(presedge[,1])
	
	pres=Matrix(0,nrow=length(genomelist),ncol=length(gclist),sparse=TRUE)
	colnames(pres)=gcvec
	rownames(pres)=genomelist
	for(i in 1:length(genomelist))
	{
		temp=which(presedge[,1]==genomelist[i])
		pres[i,presedge[temp,2]]=1
	}
	
	padj=pres%*%t(pres)   # No weight
	gadj=t(pres)%*%pres
	
	# For trimming number of edges
	psize=as.matrix(apply(pres,1,sum))      # sums over every row in pres to find the number of genes in each genome
	psizeprod=psize%*%t(psize)              # makes a matrix with same dimensions of padj where every entry is the product of genome sizes between genomes i and j
	padjw=padj/sqrt(psizeprod)              # makes weighted padj by dividing by sqrt of psizeprod
	
	# Trimming edges happens here
	padj2=padjw
	minw <- 0.3
	#padj2n <- sapply(padj, as.numeric)
	#hist(padj2n)
	padj2[padj2<minw]=0
	
	
	# Original
	#pnet=graph.adjacency(sign(padj2),diag=FALSE,mode='upper')

	# for adding edge weights to the output file
	pnet=graph.adjacency(padj2,diag=FALSE,mode='upper', weighted=TRUE)
	
	gnet=graph.adjacency(sign(gadj),diag=FALSE,mode='upper')
	
	pedge=get.data.frame(pnet)
	gedge=as_edgelist(gnet)
	
	write.table(pedge,'genome_edgelist3_weight.txt',row.names=F,col.names=F,quote=F,sep='\t')
	write.table(gedge,'genecluster_edgelist.txt',row.names=F,col.names=F,quote=F,sep='\t')
}

netprep(opt$panpath, opt$threads)