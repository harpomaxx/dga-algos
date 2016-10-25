# Based on DGA from https://github.com/jayjacobs/dga
#

library(data.table)
library(lattice)
datadir="/home/harpo/Dropbox/ongoing-work/dga/data"
resultdir="/home/harpo/Dropbox/ongoing-work/dga/results/domainname"
samples=c()
samples[1]=paste(datadir,"/datadriven/legit-dga_domains.csv",sep="")
samples[2]=paste(datadir,"/datadriven/legit-dga_domains-pablo.csv",sep="")
samples[3]=paste(datadir,"/datadriven/legit-dga_domains-pablo-short.csv",sep="")

# Find ngrams
calculate_goodngrams<-function(id=0,ngramid=0,recalculate=FALSE){
  datasample=as.data.frame(fread(samples[id],sep=",",header=T))
  if (recalculate==T){
    legitname <- datasample$domain[datasample$class=="legit"]
    onegood=ngram(legitname,1)
    twogood=ngram(legitname,2)
    threegood=ngram(legitname,3)
    fourgood=ngram(legitname,4)
    fivegood=ngram(legitname,5)
    good345=ngram(legitname,c(3,4,5))
    save(onegood,file=paste("data/",basename(samples[id]),"_onegood.rda",sep=""))
    save(twogood,file=paste("data/",basename(samples[id]),"_twogood.rda",sep=""))
    save(threegood,file=paste("data/",basename(samples[id]),"_threegood.rda",sep=""))
    save(fourgood,file=paste("data/",basename(samples[id]),"_fourgood.rda",sep=""))
    save(fivegood,file=paste("data/",basename(samples[id]),"_fivegood.rda",sep=""))
    save(good345,file=paste("data/",basename(samples[id]),"_good345.rda",sep=""))
    
  }else{
    
    if( ngramid==0){
      load(paste("data/",basename(samples[id]),"_onegood.rda",sep=""))
      load(paste("data/",basename(samples[id]),"_twogood.rda",sep=""))
      load(paste("data/",basename(samples[id]),"_threegood.rda",sep=""))
      load(paste("data/",basename(samples[id]),"_fourgood.rda",sep=""))
      load(paste("data/",basename(samples[id]),"_fivegood.rda",sep=""))
      load(paste("data/",basename(samples[id]),"_good345.rda",sep=""))
    }else{
      load(paste("data/",basename(samples[ngramid]),"_onegood.rda",sep=""))
      load(paste("data/",basename(samples[ngramid]),"_twogood.rda",sep=""))
      load(paste("data/",basename(samples[ngramid]),"_threegood.rda",sep=""))
      load(paste("data/",basename(samples[ngramid]),"_fourgood.rda",sep=""))
      load(paste("data/",basename(samples[ngramid]),"_fivegood.rda",sep=""))
      load(paste("data/",basename(samples[ngramid]),"_good345.rda",sep=""))
    
    }
  }
  # calc distances for each domain
  print("Calculating getngram")
  datasample$onegram <- getngram(onegood, datasample$domain)
  print(".")
  datasample$twogram <- getngram(twogood, datasample$domain)
  print(".")
  datasample$threegram <- getngram(threegood, datasample$domain)
  print(".")
  datasample$fourgram <- getngram(fourgood, datasample$domain)
  print(".")
  datasample$fivegram <- getngram(fivegood, datasample$domain)
  print(".")
  datasample$gram345 <- getngram(good345, datasample$domain)
  print(".")
  
  return(datasample)
}


# GIven a datasetID create features
#
preprocess_dataset <-function(datasetid,recalculate=T,datasetngramid=0){
  print(samples[datasetid])
  if ( datasetngramid==0){
    datasample=calculate_goodngrams(datasetid,recalculate = recalculate)
    filename=paste(samples[datasetid],".rda",sep="")
  }else{  
    datasample=calculate_goodngrams(datasetid,ngramid= datasetngramid,recalculate=F  )
    filename=paste(samples[datasetid],"_",samples[datasetngramid],".rda",sep="")
  }
  
  # Calculate Entropy and length
  datasample$entropy=entropy(datasample$domain)
  datasample$length=nchar(datasample$domain)

    # calculate dictionary frequency
  print("Calculating dict freq")
 
  datasample$dict <- wmatch(datasample$domain)
  print(paste("saving file",filename))
  save(datasample,file=filename)
  print("done")
  
  return (datasample)
  
  
} 

create_matrix_plot<-function(ds,size=1000){
  samp=sample(nrow(ds),size)
  ds=ds[samp,c(3,5:13)]
  #head(ds)
  splom(~ds[,2:10],data=ds,
        groups=ds$class,
        diag.panel = function(x, ...){
          yrng <- current.panel.limits()$ylim
          d <- density(x, na.rm=TRUE)
          d$y <- with(d, yrng[1] + 0.95 * diff(yrng) * y / max(y) )
          panel.lines(d)
          diag.panel.splom(x,...)
        },cex=0.5,xlab="",ylab="",auto.key = TRUE,cex.labels=0.1,pscales = 0,alpha=0.1,varname.cex=0.5
  )
}



#' Calculate the proportion of a string that matches dictionary words
#' 
#' Given a word (or vector of character/words), this will do substring
#' matches of the string against a dictionary of known words and 
#' return the percentage of letters in the word that are explained
#' by the dictionary.
#' 
#' @param text character object (or vector) of words to match
#' @param cnum logical whether to count numbers as a match or not, if not counted and the word contains numbers it will never be 100\% match.
#' @import stringr
#' @export
library(stringr)
wmatch <- function(text, cnum = T) {
  load(file="data/dwords.rda")
  matched <- lapply(text, function(txt) rep(F, nchar(txt)))
  if(cnum) {
    outs <- str_locate_all(text, "[0-9]")
    for(i in seq_along(outs)) {
      matched[[i]][outs[[i]][, 1]] <- TRUE
    }
  }
  top <- min(max(nchar(text)), length(dwords))
  for(tlen in seq(top, 3)) {
    ngs <- ngram.name(text, n=tlen)
    for(i in seq_along(ngs)) {
      pos <- which(ngs[[i]] %in% dwords[[tlen]])
      for(x in pos) {
        matched[[i]][x:(x+tlen-1)] <- T
      }
    }
  }
  sapply(matched, mean)
}

#' Tranform a vector of DNS entries to a clean data frame
#' 
#' given a vector of DNS entries, this will call urltools::suffix_extract and do
#' the following clean up tasks:
#' * remove any invalid entries (no top level domain matched)
#' * remove any duplciated entries
#' * strip any domain less than \code{strip.len} characters long.
#' 
#' @param dnames the input vector of names
#' @param strip.len numeric on which to remove variables
#' @export
#' @import urltools
cleandns <- function(dnames, strip.len=6) {
  tld <- suffix_extract(dnames)
  # pull invalid names
  tld <- tld[!is.na(tld$tld), ]
  # remove duplicates
  tld <- tld[!duplicated(tld), ]
  # yank domains < 6 characters
  tld <- tld[sapply(tld$domain, nchar) > strip.len, ]
  tld
}

#' Get a sparse matrix of ngrams from a word, given an existing ngram matrix
#' 
#' Given a vector of words and the length of grams to slice, this will
#' return a matrix of counts each ngram appears in the words.
#' 
#' @param fit existing ngram counts (vector) to match
#' @param newtxt the new word to match ngrams and count
#' @export
getngram <- function(fit, newtxt) {
  n <- unique(sapply(colnames(fit), nchar))
  sapply(newtxt, function(domain) {
    cnm.list <- table(unlist(ngram.name(domain, n)))
    matched <- names(cnm.list) %in% colnames(fit)
    colmatch <- colnames(fit) %in% names(cnm.list)
    out <- matrix(0, ncol=ncol(fit), dimnames=list(NULL, colnames(fit)))
    out[1, colmatch] <- cnm.list[matched]
    as.vector(fit %*% t(out))    
    
  })
}

# Calculate the entropy for a given string
entropy <- function(instr) {
  if (mode(instr)!="character")
    stop("Expected character array, got ", mode(instr))
  sapply(instr, function(x) {
    cts <- table(unlist(strsplit(x, NULL)))
    lns <- nchar(x)
    -sum((cts/lns) * log2(cts/lns))    
  })
}



ngram <- function(instr, n, minct=0.0005) {
  cnm.list <- ngram.name(instr, n)
  # get a count of appearances in words
  cnm.count <- table(unlist(lapply(cnm.list, unique)))
  # get a unique list of names matching at least <minct>*<words> amount
  cnm <- unique(names(cnm.count)[cnm.count>=(minct*length(instr))])
  cnm.counts <- table(unlist(cnm.list))
  outs <- cnm.counts[names(cnm.counts) %in% cnm]
  #matrix(scale(outs, center=F), nrow=1, dimnames=list(NULL, cnm))
  matrix(log10(outs), nrow=1, dimnames=list(NULL, cnm))
  
}


ngram.name <- function(instr, n) {
  lapply(instr, function(x) {
    first <- unlist(strsplit(x, NULL))
    lns <- nchar(x)
    unlist(lapply(n[n<=lns], function(i) {
      sapply(seq(i, lns), function(p) paste0(first[(p-(i-1)):p], collapse="") )  
    }))
  })
}

#
# GIven two dataset with legit and dga samples it will build new one with classes and NA fixed
#
create_dataset<-function(legit_data,dga_data,result_data,sample=0){
  dga=read.csv(file=dga_data,sep=" ",header=FALSE,stringsAsFactors = F)
  legit=read.csv(legit_data,sep=" ",header=FALSE,stringsAsFactors = F)
  
  legitdomain=suffix_extract(legit$V1)[,c(1,3,4)]
  dgadomain=suffix_extract(dga$V1)[,c(1,3,4)]

  legitdomain$class=rep("legit",nrow(legitdomain))
  dgadomain$class=rep("dga",nrow(dgadomain))
  
  legitdomain$subclass=rep("legit",nrow(legitdomain))
  dgadomain$subclass=rep("dga",nrow(dgadomain))
  
  a=rbind(legitdomain,dgadomain)
  # replace NA with domain name
  a$domain[is.na(a$domain)]=a$host[is.na(a$domain)]
  
  names(a)=c("host","domain","tld","class","subclass")
  if (sample!=0){
    a=a[sample(nrow(a),sample),]
    write.csv(a,paste(result_data,".sample",sep="",sep=""),col.names = T,quote = F,row.names = F)
  }else
    write.csv(a,result_data,col.names = T,quote = F,row.names = F)
  
  return (a)
}