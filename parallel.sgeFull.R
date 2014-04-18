# R script for Amazon EC2 cluster
# Wei-Yi Cheng
# 2013.03.07
# Revised by Tai-Hsien Ou Yang for deterministic results
# 2014.04.17
# arguments for script:
#	args[1] : input gene expression data matrix
#	args[2] : worker id
#	args[3] : total workers
#	args[4] : job codename (optional)
# For example, to run parFindAttractor on worker 1 in a 10-worker cluster:
# > Rscript parallel.sge.R data.rda 1 10 JOBID
#
#

#require(devtools)
#require(caret)
#require(RCurl)
require(cafr2)

args <- commandArgs(TRUE)



#Create Folds from caret
"createFolds" <- 
  function(y, k = 10, list = TRUE, returnTrain = FALSE)
{

  if(is.numeric(y))
    {
      ## Group the numeric data based on their magnitudes
      ## and sample within those groups.

      ## When the number of samples is low, we may have
      ## issues further slicing the numeric data into
      ## groups. The number of groups will depend on the
      ## ratio of the number of folds to the sample size.
      ## At most, we will use quantiles. If the sample
      ## is too small, we just do regular unstratified
      ## CV
      cuts <- floor(length(y)/k)
      if(cuts < 2) cuts <- 2
      if(cuts > 5) cuts <- 5
      y <- cut(
               y, 
               unique(
                      quantile(y,
                               probs =
                               seq(0, 1, length = cuts))), 
               include.lowest = TRUE)
    }


  if(k < length(y))
    {
      ## reset levels so that the possible levels and 
      ## the levels in the vector are the same
      y <- factor(as.character(y))
      numInClass <- table(y)
      foldVector <- vector(mode = "integer", length(y))
      
      ## For each class, balance the fold allocation as far 
      ## as possible, then resample the remainder.
      ## The final assignment of folds is also randomized. 
      for(i in 1:length(numInClass))
        {
          ## create a vector of integers from 1:k as many times as possible without 
          ## going over the number of samples in the class. Note that if the number 
          ## of samples in a class is less than k, nothing is producd here.
          seqVector <- rep(1:k, numInClass[i] %/% k)
          ## add enough random integers to get  length(seqVector) == numInClass[i]
          if(numInClass[i] %% k > 0) seqVector <- c(seqVector, sample(1:k, numInClass[i] %% k))
          ## shuffle the integers for fold assignment and assign to this classes's data
          foldVector[which(y == dimnames(numInClass)$y[i])] <- sample(seqVector)
        }
    } else foldVector <- seq(along = y)

  if(list)
    {
      out <- split(seq(along = y), foldVector)
      names(out) <- paste("Fold", gsub(" ", "0", format(seq(along = out))), sep = "")
      if(returnTrain) out <- lapply(out, function(data, y) y[-data], y = seq(along = y))
    } else out <- foldVector
  out
}

createMultiFolds <- function(y, k = 10, times = 5)
  {
    prettyNums <- paste("Rep", gsub(" ", "0", format(1:times)), sep = "")
    for(i in 1:times)
      {
        tmp <- createFolds(y, k = k, list = TRUE, returnTrain = TRUE)
        names(tmp) <- paste("Fold",
                            gsub(" ", "0", format(seq(along = tmp))),
                            ".",
                            prettyNums[i],
                            sep = "")
        out <- if(i == 1) tmp else c(out, tmp)
          
      }
    out
  }

source("parAttractorScanningFull.R")

args <- commandArgs(TRUE)

#================================
# Load data
#===============================

env <- new.env()
nm <- load(args[1], env)[1]
ge <- env[[nm]]

codename = args[1]
wid = 1
numWorkers=1
dirname = "output"
if(length(args)>1){
	wid = as.numeric(args[2])
	numWorkers = as.numeric(args[3])
}
if(length(args)>3){
	codename = args[4]
	dirname=codename
}

cat("=======================================\n")
cat("JOB NAME:\t", codename, "\n", sep="")
cat("FILE NAME:\t", args[1], "\n", sep="")
cat("WID:\t", wid, "\n",sep="")
cat("TOTAL WORKERS:\t",  numWorkers, "\n", sep="")
cat("=======================================\n")
flush.console()

m = nrow(ge)
set.seed(913503)
taskList = createFolds(1:m, k=numWorkers)

cat("Finding attractors...\n");flush.console()

as = parAttractorScanningFull(ge, maxIter=500 , taskList=taskList, wid=wid)
dir.create(dirname)
save(as, file=paste(dirname,"/",codename,".attractors.", sprintf("%04d", wid), ".rda", sep=""))

cat("Done.\n")
cat("=======================================\n");flush.console()

