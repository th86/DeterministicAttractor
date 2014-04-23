library("cafr")
CAFrun <- function(data, vec, a=5, maxIter = 100, epsilon=1E-14, bin = 6, so = 3,rankBased = FALSE,  negateMI = TRUE, verbose=TRUE, sorting=TRUE){
  m <- nrow(data)
  n <- ncol(data)
  
  if(rankBased){
    vec <- rank(vec)
    dataIn <- t( apply(data, 1, rank) )
  }
  miOut <- .Call("cafR2C", data, vec, a, maxIter, epsilon, m, n, bin, so, as.integer(negateMI), as.integer(verbose))
  
  if(miOut[1] == -999) return (NULL)
  
  names(miOut) <- rownames(data)
  if(sorting){
    return (sort(miOut, decreasing=T))
  }else{
    return (miOut)
  }
}


attractorScanningFull <- function(data, a=5, maxIter=100, epsilon=1E-14, bin=6, so=3, rankBased=FALSE, negateMI=TRUE){

  m <- nrow(data)
  n <- ncol(data)
  genes <- rownames(data)
  task <- 1:m
  c <- 1
  as <- NULL
while(length(task) > 0){
    i <- task[1]
    cat(genes[i], " ( ", c, " / ", length(task), ") ... ", sep="");flush.console()
   out <- CAFrun(data, data[i,], a=a, maxIter=maxIter, epsilon=epsilon, bin=bin, so=so, rankBased=rankBased, negateMI=negateMI, verbose=FALSE, sorting=FALSE)

    ##remove non-converged seeds
    if(is.null(out)){
      cat("not converged.\n");flush.console()
      task <- task[-1] #will affect i as well
      next
    }
    
    ##remove seeds heuristically
    #killIdx <- which(out >= out[i]) 
    #task <- setdiff(task, killIdx) 

    #Remove the seed just searched
    task <- task[-1]

    d <- out[order(out)[m-9]] #top - second > 0.5 is dominant, not to be included
    if(d > 0.4){

    #  if(!is.null(as)){
    #    un <- apply(as, 1, function(x){
    #      max(abs(x - out)) > 1E-14 #Distance between other rows and the attractor just found
    #    })
    #    if(prod(un) == 0){ #product of all values =0, means the attractor found is identical to one in the attractor matrix 
    #      killIdx <- which(out >= out[i])
    #      task <- setdiff(task, killIdx)
    #      rownames(as)[which(un==0)] <- genes[i]
    #      #if(verbose) {cat("done!\n");flush.console()}
    #      next
    #    }
    #  }

      as <- rbind(as, out)        #Bind the attractor found to the attractor matrix
      rownames(as)[c] <- genes[i] #Assign the rowname
      c <- c + 1
      cat("done!\n");flush.console()
    }else{
      cat("dominant.\n");flush.console()
    }

  } #End of while



#Merge attractors
cat("Merging Tasks\n");flush.console()


task <- 1:nrow(as)
as_merged <- NULL

while(length(task) > 0){
  #cat( task, rownames(as), "\n")
 if(  (length((as))/nrow(data)) > 1){
 out= as[ 1 ,]
 as_merged=rbind(as_merged, out)
 rownames(as_merged)[nrow(as_merged)]= rownames(as)[1]

  #find the rows with distance to row-in-focus larger than the threshold 
  un <- apply(as, 1, function(x){ 
          max(abs(x - out)) > 1E-4
        })

  #if there are rows with distance less than the threshold, remove them
  if(prod(un) == 0){ 
    task <- setdiff(1:nrow(as), which(un == 0) )
    #cat(length(task),"\n")
    #cat(nrow(as),"\n")
    as<-as[ task ,]

    cat("Removed seeds", which(un == 0) ,", done!\n");flush.console()
    next
  }
  }else{
    task=NULL
    as_merged=rbind(as_merged, as)
    rownames(as_merged)[nrow(as_merged)]= names(as)[which.max(as)]
  } #End of if
} #End of while

#if( length(as)==ncol(as_merged)  ){
#    as_merged=rbind(as_merged, as)
#    rownames(as_merged)[nrow(as_merged)]=names(as)[1]
#  }


  if(!is.null(as_merged)){
    colnames(as_merged) <- rownames(data)  #May cause NA here if column number changes because of identical values
  }

strength = apply(as_merged, 1, function(xx){sort(xx, decreasing=T)[10]})
as_merged = as_merged[order(strength, decreasing=T),]
cat(nrow(as_merged), "attractors in total.\n\n")
x<-as_merged

  return (as_merged)
}
