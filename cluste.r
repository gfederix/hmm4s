#!/usr/bin/Rscript --max-ppsize=500000
source("libtrack-io.R")
library(RHmm)
options(stringsAsFactors = FALSE)
COLORS <- c("cyan", "blue", "darkolivegreen" ,"magenta", "coral")
STEP <- 200
pc <- readRDS("Sm200V082Prc.rds")
x <- sample(nrow(pc), 1000)
mask <- which(diff(pc$start) != STEP )
nstart <- c(1, mask + 1)
nend <-   c(mask, length(pc$start))
obs.list <-  mapply(function(s,e) list(as.matrix(pc[s:e, c("PC1","PC2")])),
                    nstart, nend)

par(mfrow=c(1,1))
STATE <- 4
hmm <- HMMFit(obs = obs.list, nStates = STATE)
v <- unlist(sapply( obs.list, function (x) {viterbi(hmm, x)$states}, simplify=FALSE ))
hmm.track <- pc[, c("seqid", "start", "end")]
hmm.track$color <- COLORS[v]
hmm.track$value <- v
hmm.track$name  <- v
write.bed(hmm.track, "hmm")

