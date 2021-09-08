library(ggplot2)
library(reshape2)
library(plyr)

# CAUTION: for PTC there can be several rows for the same year and concept(s) due to removal of PTC type
loadRawData <- function(dir='data/21-extract-discoveries/ND.min100',indivOrJoint='indiv',sources=c('KD','PTC'), removePTCTypes=TRUE, minYear=1950,maxYear=2020) {
  ldply(sources, function(source) {
    f <- paste(dir,indivOrJoint,source,sep='/')
    print(f)
    d<-read.table(f,sep='\t',quote='',comment.char = '')
    d$source <- rep(source,nrow(d))
    d <- if (indivOrJoint == 'indiv') {
      colnames(d)[1:3] <- c('year', 'concept', 'freq')
      d[,c('year', 'concept', 'freq','source')]
    } else {
      colnames(d)[1:4] <- c('year', 'c1', 'c2', 'freq')
      d
    }
    if (source == 'PTC' & removePTCTypes) {
      if (indivOrJoint == 'indiv') {
        d<-removeTypePTC(d, 'concept')
      } else {
        if (indivOrJoint == 'joint') {
          d<-removeTypePTC(d)
        }
      }
    }
    d[d$year>=minYear & d$year<=maxYear,]
  })
}


# it's long (40 seconds for 8m joint rows) but at least it doesn't crash
removeTypePTC <- function(d,cols=c('c1','c2')) {
  res<-llply(cols, function(col) {
    l<-strsplit(as.character(d[,col]),'@', fixed=TRUE)
    res <- unlist(llply(l, function(row) {
      ll <- length(row)
      if (ll>2) {
        # remove last element then join the others
        paste(row[-c(ll)],collapse='@')
      } else {
        row[[1]]
      }
    }))
  })
  d[,cols] <- res
  d
}


loadTotalFiles <- function(dataPath='data/21-extract-discoveries', by='by-doc', minYear=1950,maxYear=2020) {
  f <- paste(dataPath,'totals',by,'KD',sep='/')
  kd<-read.table(f,sep='\t')
  f <- paste(dataPath,'totals',by,'PTC',sep='/')
  ptc<-read.table(f,sep='\t')
  kd$source <- rep('KD',nrow(kd))
  ptc$source <- rep('PTC',nrow(ptc))
  d<-rbind(kd,ptc)
  colnames(d) <- c('year', 'unique_concepts','nb','concepts_mentions','source')
  d[d$year>=minYear & d$year<=maxYear,]
}


computeProbsAndRanks <- function(df, totals, conceptCols=c('concept')) {
#  print(nrow(df))
#  print(colnames(df))
  ddply(df, c('source', 'year'), function(dataYear) {
    print(paste(dataYear[1,'source'],dataYear[1,'year']))
    totalYear <- totals[totals$source==dataYear[1,'source'] & totals$year==dataYear[1,'year'],'nb']
    resYear <- ddply(dataYear, conceptCols, function(dataYearConcepts) {
      f <- sum(dataYearConcepts$freq) # in case several rows
      data.frame(freq=f,prob=f/totalYear)
    })
    resYear$rank <- rank(resYear$freq)
    resYear$relRank <- resYear$rank / nrow(resYear)
    resYear
  })
} 
