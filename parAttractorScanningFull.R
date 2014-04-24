#Revised and commented by Tai-Hsien Ou Yang

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



parAttractorScanningFull <- function(data, taskList=list(1:nrow(data)), wid=1,  a=5, maxIter=100, epsilon=1E-14, bin=6, so=3, rankBased=FALSE, negateMI=TRUE){
  cat("Worker", wid, "initialized.\n");flush.console()
  m <- nrow(data)
  n <- ncol(data)
  genes <- rownames(data)
  task <- taskList[[wid]]
  c <- 1
  as <- NULL
  while(length(task) > 0){
    i <- task[1]
    cat("Worker ", wid, " : ", genes[i], " ( ", c, " / ", length(task), ")\n", sep="");flush.console()
    out <- CAFrun(data, data[i,], a=a, maxIter=maxIter, epsilon=epsilon, bin=bin, so=so, rankBased=rankBased, negateMI=negateMI, verbose=FALSE, sorting=FALSE)

    ##remove non-converged seeds
    if(is.null(out)){
      #if(verbose) {cat("not converged.\n");flush.console()}
      task <- task[-1] #will affect i as well
      next
    }
    
    ##remove seeds heuristically
    #killIdx <- which(out >= out[i]) 
    #task <- setdiff(task, killIdx) 

    #Remove the seed just searched
    task <- task[-1]

    d <- out[order(out)[m]] - out[order(out)[m-1]] #top - second > 0.5 is dominant, not to be included
    if(d <= 0.5){

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
      #if(verbose) {cat("done!\n");flush.console()}
    }else{
      #if(verbose) {cat("dominant.\n");flush.console()}
    }

  } #End of while


  if(!is.null(as)){
    colnames(as) <- rownames(data)  #May cause NA here if column number changes because of identical values
  }


  cat("Worker", wid, "finished.\n");flush.console()
  return (as)
}
