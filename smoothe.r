#!/usr/bin/Rscript
options(stringsAsFactors = FALSE, scipen = 10)
source("libtrack-io.R")
library(parallel)
cpuThreads <- detectCores()
frameSize = 200
rdsDir <- function(x = "") paste("rds", x, sep = "/")
smooth <- function(file, frameSize = 200, chromoRows = rownames(chromoData) ){
    dt <- readRDS(rdsDir(file))
    ## colnames(dt)[ which(colnames(dt) == "value") ] <- file
    attach(dt)
    ## frameSteps:   0|111|22|3333|000
    ## frameBorders: 1, 4, 6, 10; start: 1, 2, ...; end: 1, 4, ...
    frameStep <- as.integer(start / frameSize)
    frameBorder <- which(diff(frameStep) != 0)
    frameStart <- c(0, frameBorder) + 1
    frameEnd <- c( frameBorder, length(start))
    sValue <- mapply(function(s, e) median(value[s:e]), frameStart, frameEnd)
    sSeqid <- seqid[frameStart]
    sStart <- frameStep[frameStart] * frameSize
    if (TRUE){
        names(sValue) <- paste(sSeqid, sStart, sep = ".")
        detach(dt)
        return(sValue[chromoRows])

    } else {
        sEnd <- sStart + frameSize - 1
        sDt <- data.frame(seqid = sSeqid, start = sStart, end = sEnd, value = sValue)
        return(sDt)
    }
}
files <- list.files(rdsDir())
chromoArm <- list(  # From flybase.org: dmel-all-r5.50.gff.gz: chromosome_arm
                  "2L" = 23011544,
                  "2R" = 21146708,
                  "3L" = 24543557,
                  "3R" = 27905053,
                  "X"  = 22422827)
chromoData <- do.call(
    rbind,
    mapply(
        function(seqid, tail) {
            start <- seq(frameSize, tail - frameSize, frameSize)
            data.frame(seqid, start, end = start + frameSize - 1, row.names = start)
        }, names(chromoArm), chromoArm,
        SIMPLIFY = FALSE))
cl <- makeCluster(cpuThreads)
clusterExport(cl, ls())
sDt <- na.omit( 
    cbind(
        chromoData,
        parSapply(cl, files, smooth, simplify = "data.frame")))
stopCluster(cl)
saveRDS(sDt, "Sm200V082.rds")
sDtPrc <- prcomp(sDt[, -(1:3)])
saveRDS(cbind(sDt[,1:3],sDtPrc$x), "Sm200V082Prc.rds")
