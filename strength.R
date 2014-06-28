

getGeneSymbols = function(innames){
	outnames = sapply(innames, function(x){

		if(regexpr("\\?", x) > 0){
			o = strsplit(x, "\\|")[[1]][2]
		}else{
			o = strsplit(x, "\\|")[[1]][1]
		}
		return (o)

	}
	)
}

strength=rep(0,1,100)

for(j in 2:length(strength) ){

memberList<-names(sort(attractor,decreasing=T))[1:j]

similarityMatrix<-matrix(NA,length(memberList),length(memberList) )
rownames(similarityMatrix)=memberList
colnames(similarityMatrix)=memberList


for(i in 1:length(memberList)){
  	similarityMatrix[i,]<-getAllMIWz(ge[memberList,], ge[memberList[i],])
	}


strength[j]<- (sum(similarityMatrix^5)-nrow(similarityMatrix)) / (length(similarityMatrix)-nrow(similarityMatrix) )

}

plot(1:length(strength), strength)

