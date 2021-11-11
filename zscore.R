# average all z-scored bins
args<-commandArgs(TRUE)
dir = args[1]
filenames = list.files(dir,".*.bin.cov$")

sample_fn=args[2]
samples=read.table(sample_fn, stringsAsFactor=F)
# matching filenames
tmp = c()
for(s in samples[,1]){
	idx = grep(s, filenames, fixed = TRUE)
	if(length(idx)>0){
		print( paste("sample=",s) )
		print( idx )
		print( paste("file=",filenames[idx[1]]) )
		tmp=c(tmp, filenames[idx[1]])
	}
}
filenames=tmp
print(filenames)

ids = c()
for(f in filenames){
   if(dir.exists(paste(dir, f, sep="/"))){next}
   print( paste("reading", f) )
   x=read.table(paste(dir, f, sep="/"), header=T, stringsAsFactors=F, sep="\t")
   id = paste(x[,1],x[,2],sep="_")
   print( length(id) )
   if( length(id)<10000 ) next;
   if(length(ids)>10000) { ids = unique(sort(c(ids,id))) }
   else {ids=id}
}


#
print("making a data matrix")

dat=matrix(nrow=length(ids), ncol=length(filenames)+2)
rownames(dat) = ids;
colnames(dat) = c("chr","pos",filenames)
for(f in filenames){
   if(dir.exists(paste(dir, f, sep="/"))){next}
   x=read.table(paste(dir, f, sep="/"), header=T, stringsAsFactors=F, sep="\t")
   id = paste(x[,1],x[,2],sep="_")
   rownames(x) = id
   print( paste(f, sd(x[ids,3], na.rm=T)) )
   if(sd(x[ids,3], na.rm=T)<300) {dat[ids,f] = x[ids,3]}
   else { print(paste("skip",f)); dat[ids,f]=rep(NA,length(ids)) }
}

avg=matrix(ncol=3, nrow=length(ids))
for(i in 1:nrow(dat)){
   vs = strsplit(ids[i], "_")
   avg[i,1] = vs[[1]][1]
   avg[i,2] = as.integer(vs[[1]][2])
   #dat[i,1] = avg[i,1]
   #dat[i,2] = avg[i,2]
   avg[i,3] = mean(as.numeric(dat[i,3:ncol(dat)]), na.rm=T)
}

write.table(file=paste(sample_fn,"dat.txt",sep=""), dat, quote=FALSE, col.names = FALSE, row.names = FALSE, sep="\t" );
write.table(file=paste(sample_fn,"avg.txt",sep=""), avg, quote=FALSE, col.names = FALSE, row.names = FALSE, sep="\t" );

ziliang@ziliangdeMacBook-Pro manuscript % cat ./Z-segments/gene.R

# read in segments
segments = read.table("segments.txt", sep="\t", header=T, row.names=1, check.names=FALSE, stringsAsFactors=F)

# calculate values
args<-commandArgs(TRUE)
filenames=c()
for(dir in args){
	filenames = c(filenames, list.files(dir,".*.bin.cov$", full.names = TRUE))
}

dat = matrix(nrow=0, ncol=nrow(segments))
for(f in filenames){
   print(paste("processing",f))
   x=read.table(f, header=T, stringsAsFactors=F, sep="\t")
   m = median(x[x[,1]!="chrX"&x[,1]!="chrY",3], na.rm=T)
   print( paste(f,m) )
   m = 0

   v = c()
   for(i in 1:nrow(segments)){
      idx = (gsub("chr0","chr", x[,1])==gsub("chr0", "chr",segments[i,1])) & x[,2]>=segments[i,2] & x[,2]<=segments[i,3]
      v = c(v, median(x[idx,3], na.rm=T)-m)
   }
   dat = rbind(dat, v)
}

rownames(dat) = basename(filenames)
colnames(dat) = rownames(segments)

write.table(file="dat.txt", dat, sep="\t")

###
