read.dummy <- function(file, header=FALSE, sep="", quote = NULL,
                       col.names=c("seqid", "start", "end","value"), skip=0,
                       colClasses = NA, na.strings = c("NA"), ...){
  res <- read.table(file, header = header, sep = sep, quote = quote, skip=skip,
                    na.strings = na.strings, colClasses = colClasses,
                    stringsAsFactor=FALSE, ...)
  if(!header) names(res) <- col.names
   if(substr(as.character(res[1,"seqid"]), 1, 3) == "chr"){
    seqid     <- as.character(res$seqid)
    res$seqid <- substr(seqid, nchar("chr") + 1, nchar(seqid))
  }
  return(res)
}
read.tab <- function(file,UCSC=TRUE){
  read.dummy(file, sep='\t', header=TRUE)
}
read.wig <- function(file,UCSC=TRUE){
  read.dummy(file, col.names=c("seqid", "start", "end","value"), skip=1)
}
read.bed <- function(file,UCSC=TRUE){
  read.dummy(file, col.names=c("seqid", "start", "end","name", "score","strand","thickStart","thickEnd","itemRgb"), skip=1)
}

read.gff <- function(file){
  gff <- read.dummy(file, sep="\t", na.strings='.',
                     comment.char="#",quote=NULL,
                     col.names=c("seqid","source","type","start","end","score","strand","phase","attributes"),
                     colClasses=c("factor","factor","factor","integer","integer",NA, "factor", "factor","character")
                     )
  gff$name <- sub(".*Name=(.+?);.*", "\\1",gff$attr)
  rownames(gff) <- sub(".*ID=(.+?);.*", "\\1",gff$attr)
  ## gff$seqid <- paste("chr", gff$seqid, sep = "") 
  return(gff)
}

write.dummy <- function(data, file='', header=NULL, sep="\t", zip="gz", UCSC=TRUE,
                        quote = FALSE, row.names = FALSE, col.names = FALSE,
                        append = FALSE){
  file <- gsub(" ","_",file)
  options(scipen=100)
  is.ucsc <- substr(as.character(data[1,"seqid"]), 1, 3) == "chr"
  if (UCSC & ! is.ucsc){
    data$seqid <- paste("chr",as.character(data$seqid),sep="")
  } else {
    if (!UCSC & is.ucsc){
      seqid     <- as.character(data$seqid)
      res$seqid <- substr(seqid, nchar("chr") + 1, nchar(seqid))
    }
  }
  if(is.character(file)) {
    if(zip == "gz"|(is.logical(zip) & zip == TRUE)) {
      file <- gzfile(file <- paste(file,"gz",sep='.'),"wt")
    }
  }
  if (is.character(header)) {
    write(header, sep = '', file = file)
  }
  write.table(data, sep=sep,
              quote = quote, row.names = row.names, col.names = col.names,
              append = append,
              file = file)
  close.connection(file)
}
write.delim <- function(x,file='') write.table(x,file,sep="\t",quote=FALSE, row.names=FALSE, col.names=TRUE)
write.wig <- function(data, file,
                      name = sub("^(.*)\\.wig$", "\\1", file),
                      description = NA, value = "value",
                      type = 'wiggle_0', UCSC = TRUE){
  write.dummy(data[, c("seqid","start","end",value)],
              header = paste('track type=', type,' name="',name,'" ',
              ifelse(is.na(description), "", paste('description="',description,'"')),
              sep = ""),
              sep="\t",
              quote = FALSE, row.names = FALSE, col.names = FALSE,
              append = TRUE,
              file = file)
}
write.bed <- function(data, file,
                      name = sub("^(.*)\\.bed$", "\\1", file),
                      description = NA,
                      value = "value",
                      value2name = NA,
                      value2rgb = NA,
                      UCSC=TRUE){
  options(scipen=100)
  ## if(UCSC){
  ##   data$seqid <- paste("chr", data$seqid, sep="")
  ## }
  for(col in c("name","source", "strand", "feature")){
    if (! col %in% colnames(data))
      data[, name]  <- "."
    }
  ## data$name    <- NA
  header <- paste('track name="',name,'"',
                  ifelse(is.na(description),"",
                     paste(' description="',description,'"')),
                  sep="")
  if(!is.na(value2name[1])){
    data$name <- sapply(data[,value],function(val){value2name[val]})
  }

  if("color" %in% names(data) ){
    data$RGB <- sapply(data$color, function(col){paste(col2rgb(col),collapse=",")})}
  if(!is.na(value2rgb[1])){
    rgb <- sapply(value2rgb, function(col){paste(col2rgb(col),collapse=",")})
    data$RGB <- sapply(data[,value], function(val){rgb[val]})
  }
  if("RGB" %in% names(data) ){header <- paste(header,'itemRgb="On"')}

  for (x in c("name","strand")){
    if (! x %in% colnames(data)) data[,x] <- '.'
  }
  values <- c("seqid", "start","end",
                      "name",value,"strand",
                      "start","end")
  if("RGB" %in% colnames(data)) values <- c(values, "RGB")
  write.dummy(data[, values],UCSC=UCSC,
              header = header,
              quote = FALSE, row.names = FALSE, col.names = FALSE,
              append = TRUE, sep = "\t",
              file = paste(file, ".bed", sep=""))
}

write.gff <- function(name,track,data){
  output.kmeans <- paste("kmeans-", name, "-banding.chr"
                         ,chromosome,".gff", sep="")
  data$source  <- "FedorKmeans"
  data$feature <- "4stateBanding"
  data$strand  <- "."
  data$frame   <- "."
  data$group   <- NA
  data[data$kmeans == 1,"group"] <- "hole"
  data[data$kmeans == 2,"group"] <- "band"
  data[data$kmeans == 3,"group"] <- "interband"
  data[data$kmeans == 4,"group"] <- "superinterband"
    
  write(paste('track name="kmeans.',name,'" ',
              'description="',paste(track,collapse=","),'"',sep=""),
        file=output.kmeans)
  write.table(data[,c("seqid","source","feature",
                      "start","end","kmeans",
                      "strand","frame","group")],
              quote=FALSE,row.names=FALSE, col.names=FALSE,
              append=TRUE, sep="\t",
              file=output.kmeans)
}

matchRow <- function(data, chr, point, end=NULL){
  data <- split(data, data$seqid)
  ## data <- lapply(data, function(x) split(x, round(x$start, -3)))
  if (is.null(end)) {                     #
    .match <- function(chr, point)
      {
        r <- data[[chr]] # [[round(point, -3)]]
        if(is.null(r)) return (NA)
        ## while( (l <- nrow(r)) > 10*4){
        ## for(i in 1:2) {
          l <- nrow(r)
          median <- as.integer(l/2)
          if (point < r[ median , "start"])
            r <- r[ 1:(median-1), ]
          else
            r <- r[median:nrow(r),]
         ## }
        name <- row.names(r[r$start <= point & r$end > point, , drop=FALSE ])
        if(length(name) == 1) return(name)
        return(NA)
      }
    mapply(.match, chr, point)
  }
}

write.fasta <- function(data, name, postfix = ".fasta")
  {
    head <- names(data)
    write(paste(paste(">", head, sep="")
                , data
                , collapse="\n", sep="\n")
          , paste(name, postfix, sep=""))

  }
getDNA <- function(dna, chr, start , end, names = NULL )
{
  rez <- mapply(function(chr, start, end)
         {
           as.character( subseq( dna[[chr]], start, end))
         }
         , chr, start, end)
  if (! is.null(names)) names(rez) <- names
  return(rez)
}
