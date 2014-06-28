printTable<-function(x, filter_range=2, cols=100 , fileName="attractorTable.csv", printCSV==TRUE){

	filtered<-rep(NA,nrow(x))
	for(i in 1:nrow(x))
		filtered[i]<-sort(x[i,],decreasing=T)[filter_range]
	eligible<- which(filtered>0.5)


	filtered<-rep(NA,nrow(x))
	for(i in 1:nrow(x))
		filtered[i]<-sort(x[i,],decreasing=T)[1]-sort(x[i,],decreasing=T)[2]
	eligible<- intersect(which(filtered<0.5), eligible)

	attractor_table<-matrix(" ", 1, cols)
	geneSet<-colnames(x)

	for(i in 1:length(eligible)  ){
		attractor_table<-rbind( attractor_table , geneSet[ order(x[eligible[i],], decreasing=T)[1:cols]])
		attractor_table<-rbind( attractor_table ,x[ eligible[i] , order(x[eligible[i],], decreasing=T)[1:cols] ])
		attractor_table<-rbind( attractor_table, matrix(" ", 1, cols) )
	}

	if(printCSV==TRUE)
		write.table(attractor_table, file=fileName,quote=F,sep=",")
		
	return(attractor_table)
}