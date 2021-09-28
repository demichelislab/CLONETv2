## define imports
#' @import parallel sets ggplot2 ggrepel arules dbscan
#' @export error_table
#' @export seg_tb_toy
#' @export pileup_tumor_toy
#' @export pileup_normal_toy
#' @export snv_reads_toy
#' @export bt_toy
#' @export pl_table_toy
#' @export adm_table_toy
#' @export scna_clonality_table_toy
#' @export allele_specific_cna_table_toy
#' @export snv_clonality_table_toy
NULL

fromListToDF <- function(inputList){
  if (is.null(inputList)){return(NULL)}

  #check if some is null and remove
  nullPositions <- which(sapply(inputList,is.null))
  if (length(nullPositions) > 0){
    inputList <- inputList[-nullPositions]
  }


  #inputList <- globalAdm.list
  firstEl <- inputList[[1]][1,]


  #inputList <- lapply(inputList, function(x){ cat(x$gene,i,"\n"); i<-i+1;matrix(unlist(x), ncol=ncol(x))} )
#   for(i in c(1:length(inputList))){
#     #cat(inputList[[i]]$Gene.id[1],i,"\n")
#     matrix(unlist(inputList[[i]]), ncol=ncol(inputList[[i]]))
#   }
#
  inputList <- lapply(inputList, function(x){ matrix(unlist(x), ncol=ncol(x))} )

  #outDF <- as.data.frame(do.call(rbind, inputList),stringsAsFactors=F)
  outDF <- as.data.frame(do.call(rbind, inputList),stringsAsFactors=F)
  colnames(outDF) <-names(firstEl)

   for(idx in c(1:ncol(outDF))){
     if (class(firstEl[[idx]]) == "logical"){
       if (is.na(firstEl[[idx]])){
         class(outDF[[idx]]) <- "numeric"
       }else if (outDF[1,idx] == "TRUE" ||outDF[1,idx] == "FALSE"  ){
        outDF[,idx] <- as.logical(outDF[,idx]) * 1
       }
       class(outDF[[idx]]) <- "numeric"
     }else if (class(firstEl[[idx]]) == "factor"){
     	 class(firstEl[[idx]]) == "character"
     }else {
       class(outDF[[idx]]) <- class(firstEl[[idx]])
     }
   }
#    tt<-lapply(c(1:ncol(outDF)),function(idx){
#      class(outDF[[idx]]) <- class(firstEl[[idx]])
#    })


  #sapply(firstEl, class)
  #sapply(outDF, class)
  #f <- which(sapply(firstEl, class)=="factor")
#   for(i in f) {
#        lev <- levels(firstEl[[i]])
#        outDF[[i]] <- factor(as.integer(outDF[[i]]), levels=seq_along(lev), labels=lev)
#    }

  return(outDF)
}

setBetween0and1 <- function(x){
  return(max(0,min(1,x)))
}



## Function to extend beta table with cnA and cnB
## It is required that column log2.plCorr is defined
extendBetaTableWithCopyNumbers <- function(BetaTable, G, ncores = 1){
  #i<-1
  extendBetaLine <- function(i, BetaTable, G ){
    cat(i,"\n");
    BetaLine <- BetaTable[i,,drop=F]

    cnAB.beta <- from_BetaLogR_to_cnAB(G = G, LogR = BetaLine$log2.plCorr, Beta = BetaLine$beta)
    cnAB.beta.min <- from_BetaLogR_to_cnAB(G = G, LogR = BetaLine$log2.plCorr, Beta = BetaLine$beta.min)
    cnAB.beta.max <- from_BetaLogR_to_cnAB(G = G, LogR = BetaLine$log2.plCorr, Beta = BetaLine$beta.max)

    BetaLine$cnA <- cnAB.beta$cnA
    BetaLine$cnB <- cnAB.beta$cnB
    BetaLine$cnA.betaMin <- cnAB.beta.min$cnA
    BetaLine$cnB.betaMin <- cnAB.beta.min$cnB
    BetaLine$cnA.betaMax <- cnAB.beta.max$cnA
    BetaLine$cnB.betaMax <- cnAB.beta.max$cnB
    return(BetaLine)

  }

  BetaTableExt.list <- mclapply(seq(1,nrow(BetaTable),1),extendBetaLine, BetaTable=BetaTable, G=G,mc.preschedule = T,mc.cores = ncores)
  BetaTableExt <- fromListToDF(BetaTableExt.list)
  return(BetaTableExt)
}


# for each line of the beta table compute min and max beta  accordingly to error table
addErrorToBetaTable <- function(BetaTable, errorTable,ncores=1){


  betaLine <- 1
  addErrorToBetaLine <- function(betaLine, BetaTable, errorTable ){
    betas <-  BetaTable$beta[betaLine]
    nsnps <- BetaTable$nsnp[betaLine]
    cov <- BetaTable$cov[betaLine]

    idxError <- which(errorTable$n.info.snps <= nsnps & errorTable$mean.cov <= cov)
    if ( length(idxError) == 0 ){ idxError <- 1 }
    local.errors <- errorTable$adm.estimation.error[max(idxError)]

    BetaTable.line <- BetaTable[betaLine,,drop=F]
    BetaTable.line$beta.error <- local.errors
    BetaTable.line$beta.min <- max(0,betas-local.errors)
    BetaTable.line$beta.max <- min(1,betas+local.errors)
    return(BetaTable.line)

  }

  extBetaTable.list <- mclapply(seq(1,nrow(BetaTable),1), addErrorToBetaLine,  BetaTable=BetaTable, errorTable=errorTable, mc.preschedule = T, mc.cores = ncores )
  extBetaTable <- do.call(rbind,extBetaTable.list)

  return(extBetaTable)
  #   plot(extBetaTable$beta,extBetaTable$log2)
  #   plot(extBetaTable$beta.min,extBetaTable$log2,pch=20,col="red4")
  #   points(extBetaTable$beta.max,extBetaTable$log2,pch=20,col="blue4",new=T)

}


getSegmentsPos <- function(chr, initialPos, finalPos, segments,flank = 0){
  #select segments that intersect initialPos and finalPs
  suppressWarnings(goodSegs <- which(
      ( as.numeric(gsub("chr","",segments[,1])) == as.numeric(gsub("chr","",chr)) |
        gsub("chr","",segments[,1]) == gsub("chr","",chr)  ) & (
      (as.numeric(segments[,2]) >= (as.numeric(initialPos) - flank) & as.numeric(segments[,2]) <= (as.numeric(finalPos) + flank) ) |
      (as.numeric(segments[,3]) >= (as.numeric(initialPos) - flank) & as.numeric(segments[,3]) <= (as.numeric(finalPos) + flank) ) |
      (as.numeric(segments[,2]) <= (as.numeric(initialPos) - flank) & as.numeric(segments[,3]) >= (as.numeric(finalPos) + flank) ))))

  return(goodSegs)
}
