# R script for merging the results generated from Amazon EC2 cluster
# Wei-Yi Cheng
# 2013.03.07
# Revised by Tai-Hsien Ou Yang for deterministic results
# 2014.04.17
# arguments for script:
#	args[1] : directory of the results (contains exclusively the .rda files)
# Example
# > Rscript mergeResultsFull.R resultDir/
#
#

args = commandArgs(TRUE)

filePath = "./"
if(length(args) > 0){
	filePath = args[1]
}

fileNames = list.files(pattern="*.attractors.*.rda", path = filePath)
nf = length(fileNames)

x = NULL
for(i in 1:nf){
	cat("Processing", fileNames[i], "...\n");flush.console()
	load(file.path(filePath,fileNames[i]))
	if(is.null(as)) next
	if(is.null(x)){
		x = as
		next
	}

	x<-rbind(x, as)

} #End of for
gc()


task <- 1:nrow(x)
as_merged <- NULL

while(length(task) > 0){
if(  (length(x)/ncol(as)) > 1){
    out= x[1 ,]
    as_merged=rbind(as_merged, out)
    rownames(as_merged)[nrow(as_merged)]= rownames(x)[task[1]]
  
  #find the rows with distance to row-in-focus larger than the threshold 
  un <- apply(x, 1, function(xx){ 
          max(abs(xx - out)) > 1E-4
        })

  #if there are rows with distance less than the threshold, remove them
  if(prod(un) == 0){ 
    task <- setdiff(1:nrow(x), which(un == 0) )
    x<-x[ task ,]
    cat("Removed", which(un == 0) ,"seeds, done!\n");flush.console()
    next
  }
  }else{
    task=NULL
    as_merged=rbind(as_merged, x)
    rownames(as_merged)[nrow(as_merged)]= names(x)[which.max(x)]
  } #End of if
} #End of while

strength = apply(as_merged, 1, function(xx){sort(xx, decreasing=T)[10]})
as_merged = as_merged[order(strength, decreasing=T),]
cat(nrow(as_merged), "attractors in total.\n\n")
x<-as_merged

save(x, file="attractorMatrix.rda")
