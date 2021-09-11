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


loadTotalFiles <- function(dataPath='data/21-extract-discoveries', by='by-doc', minYear=1950,maxYear=2020, filterCols=c('source','year','nb')) {
  f <- paste(dataPath,'totals',by,'KD',sep='/')
  kd<-read.table(f,sep='\t')
  f <- paste(dataPath,'totals',by,'PTC',sep='/')
  ptc<-read.table(f,sep='\t')
  kd$source <- rep('KD',nrow(kd))
  ptc$source <- rep('PTC',nrow(ptc))
  d<-rbind(kd,ptc)
  colnames(d) <- c('year', 'unique_concepts','nb','concepts_mentions','source')
  d[d$year>=minYear & d$year<=maxYear,filterCols]
}


# OBSOLETE
# for old method (with rank added)
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


# OBSOLETE
# old method, same as first python script
applyWindow <- function(df, halfWindowSize=5,conceptCols=c('concept'), targetCols=c('prob','relRank')) {
  years <- sort(unique(df$year))
  minYear <- min(years)
  maxYear <- max(years)
  ddply(df, conceptCols, function(conceptDF) {
#    print(head(conceptDF))
    ldply(seq(minYear+halfWindowSize,maxYear-halfWindowSize+1), function(currentYear) {
      # note that pastDF can be empty
#      print(currentYear)
      pastDF <- conceptDF[conceptDF$year<currentYear & conceptDF$year>=currentYear-halfWindowSize,]
#      print(head(pastDF))
      futureDF <- conceptDF[conceptDF$year>=currentYear & conceptDF$year<currentYear+halfWindowSize,]
      ldply(targetCols, function(targetCol) {
        pastMean <- if (nrow(pastDF) > 0) { mean(pastDF[,targetCol]) } else 0
        futureMean <- if (nrow(pastDF) > 0) { mean(futureDF[,targetCol]) } else 0
        r<-data.frame(year=currentYear, target_col=targetCol, prev.mean=pastMean, next.mean=futureMean)
#        print("RES")
#        print(r)
        r
      })
    })
  })
}


# OBSOLETE
# adds zeros when then is no data for a year
completeFreqDataWithZeros <- function(df, conceptCols=c('concept'), minYear=1950, maxYear=2020) {

  cols <- c('source', conceptCols)
  ddply(df, cols, function(subDF) {
    ldply(seq(minYear,maxYear), function(year) {
      d <- subDF[subDF$year == year,]
      f <- if (nrow(d)==0) 0 else sum(d[,'freq'])   # in case several rows due to PTC types
      data.frame(year=year,freq=f)
    })
  })
    
}

# moving average on a vector 'x' with window size 'window' (centered).
# The first/last few elements in the resulting vector are filled with NA, so that the size is the same as the input vector.
# from https://stackoverflow.com/questions/743812/calculating-moving-average
# this solution gives an error when used insde ddply, I wasn't able to solve it.
#ma <- function(x, n = 5){
#  stats::filter(x, rep(1 / n, n), sides = 2)
#}

# moving average on a vector 'x' with window size 'window' (centered).
ma <- function(x,n=5,padWithNA=FALSE) {
  cx <- c(0,cumsum(x))
  r <- (cx[(n+1):length(cx)] - cx[1:(length(cx) - n)]) / n
  if (padWithNA) {
    nas <- rep(NA, (n-1)/2)
    c(nas,r,nas)
  } else {
    r
  }
}


# OBSOLETE
# adds a new column with moving average for 'col' named '<col>.ma<window>'
# df <- completeFreqDataWithZeros(...)
addMovingAverageCol <- function(df, col='freq', window=5, idCols=c('source','concept'),sanityCheckNbYears=71) {
  colName <- paste0(col,'.ma',window)
  ddply(df, idCols, function(subDF) {
#    print(subDF[1,idCols])
    if (!is.na(sanityCheckNbYears) & !is.null(sanityCheckNbYears)) {
      if (length(unique(subDF$year)) != sanityCheckNbYears) {
        print(head(subDF))
        stop(paste('Error: sanity check failed, the number of years should be',sanityCheckNbYears,'but this subset has',length(unique(subDF$year))))
      }
    }
    resDF <- subDF[order(subDF$year),] # making sure the years are consecutive
    resDF[,colName] <- ma(resDF[,col],window,padWithNA = TRUE)
#    print(head(resDF))
    resDF
  })
}


# input = data frame with year, 'targetCol'
# output = data frame with year, 'targetCol', prev.mean, next.mean
# df is supposed to contain a single time series, complete (single freq for every year)
computePrevVsNext <- function(df, halfWindowSize=2, targetCol='freq', minYear=1950, maxYear=2020, padWithNA=FALSE) {
  r<-ldply(seq(minYear+halfWindowSize,maxYear-halfWindowSize+1), function(currentYear) {
#    print(currentYear)
    prevDF <- df[df$year<currentYear & df$year>=currentYear-halfWindowSize,]
    nextDF <- df[df$year>=currentYear & df$year<currentYear+halfWindowSize,]
#    print(prevDF)
#    print(nextDF)
    prevMean <- if (nrow(prevDF) > 0) { mean(prevDF[,targetCol]) } else 0
    nextMean <- if (nrow(nextDF) > 0) { mean(nextDF[,targetCol]) } else 0
    data.frame(year=currentYear, target_col=targetCol, prev.mean=prevMean, next.mean=nextMean)
  })
  r<-if (padWithNA) {
    nas <- rep(NA, halfWindowSize-1)
    data.frame(year=df$year, target_col=df[,targetCol], prev.mean=c(NA,nas,r$prev.mean,nas), next.mean=c(NA,nas,r$next.mean,nas))
  } else {
    r
  }
  colnames(r)[2] <- targetCol
  r
}


# new method which prepares tha aggregated data for all the methods depending on options
# the data frame doesn't need to be completed with zeros first, i.e. df <- loadRawData(...)
#
# aggregMethods: ma = moving average, pvsn = previous vs. next
# half_window: 
#  - for ma uses window=half_window*2+1 (odd number, window centered on the year)
#  - for pvn uses half_window for previous years and current year + (half_window-1) next years
#
applyAggregMethods <- function(freqDF, totalsDF, idCols=c('source','concept'), half_windows=c(2,3,4),aggregMethods=c('ma','pvn'), minYear=1950, maxYear=2020) {
  fullYearsDF <- data.frame(year=seq(minYear,maxYear))
  totals <- list()
  # prepare totals with moving average or pvn (in a list)
  if ('ma' %in% aggregMethods) {
    totals$ma <- ddply(totalsDF, 'source', function(seqDF) {
      ldply(half_windows, function(half_window) {
        data.frame(half_window=rep(half_window,nrow(seqDF)),year=seqDF$year,nb=ma(seqDF[,'nb'],half_window*2+1,padWithNA = TRUE))
      })
    })
  }
  if ('pvn' %in% aggregMethods) {
    totals$pvn <- ddply(totalsDF, 'source', function(seqDF) {
      ldply(half_windows, function(half_window) {
        r <- computePrevVsNext(seqDF,half_window,'nb', minYear, maxYear,padWithNA = TRUE)
        r$half_window <- rep(half_window,nrow(r))
        r
      })
    })
  }
  ddply(freqDF, idCols, function(subDF) {
    #    print(subDF[1,idCols])
    # 1. aggregate in case multiple rows for same year due to PTC types
    aggDF <- ddply(subDF,'year', function(singleYearDF) { data.frame(freq=sum(singleYearDF$freq)) })
    # 2. complete with 0 for missing years
    stdDF <- merge(aggDF,fullYearsDF,all=TRUE)
    stdDF[is.na(stdDF$freq),'freq'] <- 0
    # 3. making sure the years are consecutive (probably not needed)
    stdDF <- stdDF[order(stdDF$year),]  
    ldply(half_windows, function(half_window) {
      resDF <- stdDF
      if ('ma' %in% aggregMethods) {
        resDF$ma <- ma(resDF[,'freq'],half_window*2+1,padWithNA = TRUE)
        selectedTotalSeries <- totals$ma[totals$ma$half_window == half_window & totals$ma$source == subDF[1,'source'],]
        resDF$ma.total <- selectedTotalSeries$nb
      }
      if ('pvn' %in% aggregMethods) {
        pvnDF <- computePrevVsNext(stdDF,half_window,'freq', minYear, maxYear,padWithNA = TRUE)
        resDF$prev.mean <- pvnDF$prev.mean
        resDF$next.mean <- pvnDF$next.mean
        selectedTotalSeries <- totals$pvn[totals$pvn$half_window == half_window & totals$pvn$source == subDF[1,'source'],]
        resDF$prev.mean.total <- selectedTotalSeries$prev.mean
        resDF$next.mean.total <- selectedTotalSeries$next.mean
      }
      resDF
    })
  })
}
