#!/usr/bin/Rscript
## (C) 2013 Fyodor Goncharov <fedor@mcb.nsc.ru>
## modEncode File Downloader
## Download wig files from modENCODE and serialize them.
options(stringsAsFactors=FALSE);
download<- function(...) while(try(download.file(...)) != 0) { cat("Error while file downloaded from modENCODE. I'm try restart process.", "\n")}
source("libtrack-io.R")
xDir <- function(x) function(name = "") paste(x, "/", name, sep="");
wigDir <- xDir("wig");
rdsDir <- xDir("rds");
files <- read.csv("modencodeFiles.csv");
dir.create(wigDir())
with(files, mapply(function(url, name){download(url, wigDir(name))}, URL, Title))
wigFiles <- list.files(wigDir())
dir.create(rdsDir())
library(parallel)                       #  IO bottleneck :'/
cl <- makeCluster(detectCores())
clusterExport(cl, ls())
parSapply(cl, wigFiles, function(file) saveRDS(read.wig(paste(wigDir(file))), rdsDir(file)))
stopCluster(cl)
