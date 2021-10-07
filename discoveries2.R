library(ggplot2)
library(reshape2)
library(plyr)

# CAUTION: for PTC there can be several rows for the same year and concept(s) due to removal of PTC type
loadRawData <- function(dir='data/21-extract-discoveries/',indivOrJoint='indiv',suffix='ND.min100',sources=c('KD','PTC','MED'), removePTCTypes=TRUE, addMeshPrefixMED=TRUE, minYear=1950,maxYear=2020, debugPath=FALSE) {
  ldply(sources, function(source) {
    f <- paste(dir,source, paste(indivOrJoint,suffix,sep='.'),sep='/')
    if (debugPath) {
      f
    } else {
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
      if (source=='MED' & addMeshPrefixMED) {
        if (indivOrJoint == 'indiv') {
          d$concept <- paste0('MESH:',d$concept)
        } else {
          if (indivOrJoint == 'joint') {
            d$c1 <- paste0('MESH:',d$c1)
            d$c2 <- paste0('MESH:',d$c2)
          }
        }
      }
      d[d$year>=minYear & d$year<=maxYear,]
    }
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


loadTotalFiles <- function(dataPath='data/21-extract-discoveries', sources=c('KD','MED','PTC'), minYear=1950,maxYear=2020, filterCols=c('source','year','nb')) {
  ldply(sources, function(source) {
    f <- paste(dataPath, source, 'indiv.full.total', sep='/')
    df <- read.table(f,sep='\t')
    colnames(df) <- c('year', 'unique_concepts','nb','concepts_mentions')
    df$source <- rep(source,nrow(df))
    df[df$year>=minYear & df$year<=maxYear,filterCols]
  })
}

loadTotalFilesOLD <- function(dataPath='data/21-extract-discoveries', by='by-doc', minYear=1950,maxYear=2020, filterCols=c('source','year','nb')) {
  f <- paste(dataPath,'totals',by,'KD',sep='/')
  kd<-read.table(f,sep='\t')
  f <- paste(dataPath,'totals',by,'PTC',sep='/')
  ptc<-read.table(f,sep='\t')
  f <- paste(dataPath,'totals',by,'MED',sep='/')
  med<-read.table(f,sep='\t')
  kd$source <- rep('KD',nrow(kd))
  ptc$source <- rep('PTC',nrow(ptc))
  med$source <- rep('MED',nrow(med))
  d<-rbind(kd,ptc,med)
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
      resDF$half_window <- rep(half_window, nrow(resDF))
      resDF
    })
  },.progress = "text")
}

calculateIndicators <- function(aggregDF, probRate, probDiff, freqDiff, indicators=c('prob.rate','prob.diff','product.log','product'),replaceNonFiniteWithNA=TRUE) {
  ldply(indicators, function(indicator) {
    resDF <- aggregDF
    if (indicator=='prob.rate') {
      resDF$trend <- probRate
    } else {
      if (indicator=='prob.diff') {
        resDF$trend <- probDiff
      } else {
        if (indicator=='product.log') {
          resDF$trend <- rep(NA, nrow(resDF))
          # note: freqDiff must be higher than 1 in order to prevent negative value
          regularCase <- !is.na(freqDiff) & freqDiff > 1
          resDF[regularCase,]$trend <- probRate[regularCase] * log2(freqDiff[regularCase])
        } else {
          if (indicator=='product') {
            resDF$trend <- probDiff * probRate
            negResult <- !is.na(probDiff) & !is.na(probRate) & probDiff<0 & probRate <0
            resDF$trend[negResult] <- -1 * resDF$trend[negResult] 
          } else {
            stop(paste('Error: invalid indicator id',indicator))
          }
        }
      }
    }
    if (replaceNonFiniteWithNA) {
      resDF[!is.finite(resDF$trend),]$trend <- rep(NA, length(resDF$trend[!is.finite(resDF$trend)]))
    }
    resDF$indicator <- rep(indicator,nrow(resDF))
    resDF
  })
}

# aggregDF <- applyAggregMethods(...)
computeSurgeIndicators <- function(aggregDF,indicators=c('prob.rate','prob.diff','product.log','product')) {
  minYear <- min(aggregDF$year)
#  resDF <- data.frame()
  maDF <- if ('ma' %in% colnames(aggregDF)) {
    maAsProb <- aggregDF$ma/aggregDF$ma.total
    maPrevValues <- head(maAsProb,-1) # remove last
    maNextValues <- tail(maAsProb,-1) # remove first
    probRate <- (maNextValues-maPrevValues) / maPrevValues
    probRate <- c(NA, probRate)
    probDiff <- maNextValues-maPrevValues
    probDiff <- c(NA, probDiff)
    maPrevFreq <- head(aggregDF$ma,-1) # remove last
    maNextFreq <- tail(aggregDF$ma,-1) # remove first
    freqDiff <- c(NA, maNextFreq-maPrevFreq)
    maDF <- calculateIndicators(aggregDF, probRate, probDiff, freqDiff, indicators)
    maDF$aggreg.method <- rep('ma',nrow(maDF))
    maDF
  }
  pvnDF <- if ('prev.mean' %in% colnames(aggregDF)) {
    prevAsProb <- aggregDF$prev.mean / aggregDF$prev.mean.total
    nextAsProb <- aggregDF$next.mean / aggregDF$next.mean.total
    probRate <- (nextAsProb-prevAsProb) / prevAsProb
    probDiff <- nextAsProb-prevAsProb
    freqDiff <- aggregDF$next.mean-aggregDF$prev.mean
    pvnDF <- calculateIndicators(aggregDF, probRate, probDiff, freqDiff, indicators)
    pvnDF$aggreg.method <- rep('pvn',nrow(maDF))
    pvnDF
  }
  rbind(maDF,pvnDF)
}

# from raw data:
# df <- loadRawData(indivOrJoint = 'joint')
filterConceptFromJoint <- function(df, cui='C0002736', mesh='MESH:D000690') {
  rbind(df[df$source == 'KD' & (df$c1 %in% cui | df$c2 %in% cui),],
        df[df$source == 'PTC' & (df$c1 %in% mesh | df$c2 %in% mesh),])
}


pickRandom <- function(df, idCols=c('source','c1','c2'), n=1, fastNotReallyRandom=FALSE) {
  if (fastNotReallyRandom) {
    selected <- df[sample(nrow(df),n),idCols]
    
  } else {
    values <- unique(df[,idCols])
    print(paste("Unique id cols = ",nrow(values)))
    n <- min(nrow(values),n)
    selected <- values[sample(nrow(values),n),idCols]
  }
#  print(selected)
  merge(selected,df)
}

# plot functions below for a single relation, e.g. picked by pickRandom

plotTrendByYear <- function(df) {
 ggplot(df, aes(year, trend)) + geom_col() + facet_grid(aggreg.method+indicator ~ half_window ,scales = "free") 
}

plotRawByYear <- function(df) {
  raw <- unique(df[,c('year','freq')])
  ggplot(raw,aes(year,freq))+geom_col()
}

plotMAByYear <- function(df) {
  maDF <- unique(df[,c('year','ma','half_window')])
  ggplot(maDF,aes(year,ma))+geom_col() + facet_grid(. ~ half_window ,scales = "free") 
}

plotPVNByYear <- function(df) {
  pvn <- df[df$aggreg.method=='pvn' & df$indicator=='prob.rate',c('year','prev.mean','next.mean','half_window')]
  d <- melt(pvn,id.vars=c('year','half_window'),variable.name='side',value.name='mean')
  ggplot(d,aes(year,mean,fill=side))+geom_col(position='identity',alpha=.5) + facet_grid(. ~ half_window ,scales = "free") 
}


calculateThresholdTopOutliers <- function(v0) {
  v <- v0[!is.na(v0)]
  if (length(v)>0) {
    quartiles <- quantile(v,c(.25,.75))
    iqr <- quartiles[2] - quartiles[1] # Q3-Q1
    quartiles[2] + 1.5*iqr
  } else {
    NA
  }
}

computeDiscoveryYear <- function(df, idCols=c('source','concept','half_window','indicator','aggreg.method'),discMethods=c('max.trend','earliest.outlier')) {
  minYear <- min(df$year)
  maxYear <- max(df$year)
  ddply(df, idCols, function(singleCaseDF) {
    if (nrow(singleCaseDF) != maxYear-minYear+1) {
      stop(paste('Error: sanity check failed, expected',maxYear-minYear+1,'years for single case but found',nrow(singleCaseDF)))
    }
    ldply(discMethods, function(discMethod) {
      discTrend <- NA
      discYear <- NA
      if (nrow(singleCaseDF[!is.na(singleCaseDF$trend),]) > 0) {
#        print(singleCaseDF[!is.na(singleCaseDF$trend),])
        if (discMethod == 'max.trend') {
          discTrend <- max(singleCaseDF$trend,na.rm = TRUE)
          discYear <- singleCaseDF[!is.na(singleCaseDF$trend) & singleCaseDF$trend==discTrend,'year']
        } else {
          if (discMethod == 'earliest.outlier') {
            t <- calculateThresholdTopOutliers(singleCaseDF$trend)
            if (!is.na(t)) {
              topOutliers <-singleCaseDF[!is.na(singleCaseDF$trend) & singleCaseDF$trend>=t,]
              if (nrow(topOutliers)>0) {
                discYear <- min(topOutliers$year) 
                discTrend <- topOutliers[topOutliers[,'year']==discYear,'trend']
              }
            }
          } else {
            stop(paste('Error: invalid discovery method id:',discMethod))
          }
        }
      }
      data.frame(disc.year=discYear, disc.trend=discTrend, disc.method=discMethod)
    })   
  })  
#  }, .progress=TRUE)  # CAUSES BUG!!!
}



# this one gives an error
computeDiscoveryYearBUG8 <- function(df, idCols=c('source','concept','half_window','indicator','aggreg.method'),discMethods=c('max.trend','earliest.outlier')) {
  ddply(df, idCols, function(singleCaseDF) {
    x <- singleCaseDF$trend
    length(x)
#  })  
  }, .progress=TRUE)  
}

# this one doesn't give an error
computeDiscoveryYearBUG9 <- function(df, idCols=c('source','concept','half_window','indicator','aggreg.method'),discMethods=c('max.trend','earliest.outlier')) {
  ddply(df, idCols, function(singleCaseDF) {
    x <- singleCaseDF$trend
    length(x)
#  })  
  })  
}


groupOrder <- function(disc_res, groupCols=c('source','half_window','indicator','aggreg.method','disc.method'), writeToFilesPrefix=NULL) {
  ddply(disc_res, groupCols, function(groupDF) {
    res<-groupDF[order(-groupDF$disc.trend),]
    if (!is.null(writeToFilesPrefix)) { 
      id = paste(groupDF[1,groupCols], collapse='.')
      name = paste(writeToFilesPrefix,id,'tsv',sep='.')
      write.table(res, name, quote=FALSE,sep='\t', row.names=FALSE)
    }
    res
  })
}


# temporary first study with annotated data
# aggJointDF <- read.table('21-extract-discoveries/recompute-with-ND-group/ND.agg-joint.top-pmi-npmi.tsv',sep='\t',quote='',header=TRUE)
reformatAggregatedJointTopPMI <- function(aggJointDF,aggIndivDF,aggTotalDF) {
  med <- aggJointDF[aggJointDF$source=='MED',]
  colnames(med)[colnames(med)=='freq.x'] <- 'all.years.freq.joint'
  colnames(med)[colnames(med)=='freq.y'] <- 'freq.joint'
  med
}