library(ggplot2)
library(reshape2)
library(plyr)
library(rpart)
library(rpart.plot)
library(cowplot)
library(lsa)

mesh.ALS <- 'MESH:D000690'
mesh.FTD <- 'MESH:D057180'
mesh.adult <- 'MESH:D000328'
mesh.elderly <- 'MESH:D000368'
mesh.brain <- 'MESH:D001921'

subsetForC9ORF72Study <- c(mesh.ALS,mesh.FTD,mesh.adult,mesh.elderly,mesh.brain)

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

computeDiscoveryYear <- function(df, idCols=c('source','concept','half_window','indicator','aggreg.method'),discMethods=c('max.trend','earliest.outlier'), discardNonFiniteTrendValues=TRUE) {
  minYear <- min(df$year)
  maxYear <- max(df$year)
  ddply(df, idCols, function(singleCaseDF) {
    if (nrow(singleCaseDF) != maxYear-minYear+1) {
      stop(paste('Error: sanity check failed, expected',maxYear-minYear+1,'years for single case but found',nrow(singleCaseDF)))
    }
    if (discardNonFiniteTrendValues) {
      singleCaseDF <- singleCaseDF[is.finite(singleCaseDF$trend),]
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
reformatAggregatedJointTopPMI <- function(aggJointDF,aggIndivDF,discDF=NULL,joint_cols=c('c1','c2','all.years.freq.joint','freq.joint','status','year','ma.joint'),indiv_cols=c('concept','year','freq','ma','ma.total'),withNames=FALSE) {
  topJointMed <- aggJointDF[aggJointDF$source=='MED',]
  colnames(topJointMed)[colnames(topJointMed)=='freq.x'] <- 'all.years.freq.joint'
  colnames(topJointMed)[colnames(topJointMed)=='freq.y'] <- 'freq.joint'
  colnames(topJointMed)[colnames(topJointMed)=='ma'] <- 'ma.joint'
#  colnames(topJointMed)[colnames(topJointMed)=='ma.total'] <- 'ma.total.joint'
  topJointMed[,'ma.total'] <- NULL
  topJointMed[,'NA.'] <- NULL
  if (withNames) {
    joint_cols <- c(joint_cols, 'term1','group1','term2', 'group2')
  }
  topJointMed <- topJointMed[,joint_cols]
  fullIndivMed <- aggIndivDF[aggIndivDF$source=='MED',]
  fullIndivMed <- fullIndivMed[,indiv_cols]
  joint1 <- merge(topJointMed,fullIndivMed,by.x=c('year','c1'),by.y=c('year','concept'))
  colnames(joint1)[colnames(joint1)=='freq'] <- 'freq.c1'
  colnames(joint1)[colnames(joint1)=='ma'] <- 'ma.c1'
  joint1[,'ma.total'] <- NULL # removing because ma.total is the same for c1 and c2
  res <- merge(joint1,fullIndivMed,by.x=c('year','c2'),by.y=c('year','concept'))
  colnames(res)[colnames(res)=='freq'] <- 'freq.c2'
  colnames(res)[colnames(res)=='ma'] <- 'ma.c2'
  if (!is.null(discDF)) {
    disc0 <- discDF[discDF$source=='MED',c('c1','c2','indicator','disc.year','disc.trend')]
    disc0$disc.prev.year <- disc0$disc.year-1
    res1 <- merge(res, disc0,by.x=c('c1','c2','year'),by.y=c('c1','c2','disc.prev.year'))
    res1$position <- rep(1,nrow(res1))
    res1$disc.year <- NULL
    res2 <- merge(res, disc0,by.x=c('c1','c2','year'),by.y=c('c1','c2','disc.year'))
    res2$position <- rep(2,nrow(res2))
    res2$disc.prev.year <- NULL
    res <- rbind(res1,res2)
    res$ma.c1Gc2 <- res$ma.joint / res$ma.c2
    res$ma.c2Gc1 <- res$ma.joint / res$ma.c1
    res$diffCondi <- abs(res$ma.c1Gc2 - res$ma.c2Gc1)
  }
  res  
}

study1 <- function(x, bins=4) {
  y <- x[,c('c1','c2',"status",'position','indicator','ma.c1Gc2',"ma.c2Gc1",'ma.c1','ma.c2')]
  y <- ddply(y, c('c1','c2','indicator','status'), function(s) {
    if (nrow(s)!=2) {
      print(s)
      stop("bug")
    }
    s1 <- s[s$position==1,]
    s2 <- s[s$position==2,]
    data.frame(prev.c1Gc2=s1$ma.c1Gc2, prev.c2Gc1=s1$ma.c2Gc1,next.c1Gc2=s2$ma.c1Gc2, next.c2Gc1=s2$ma.c2Gc1,
               prev.c1=s1$ma.c1,prev.c2=s1$ma.c2,next.c1=s2$ma.c1,next.c2=s2$ma.c2)
  })
  if (bins>0) {
    b <- seq(0,1,1/bins)
    y$bin.prev.c1Gc2 <- cut(y$prev.c1Gc2, b)
    y$bin.prev.c2Gc1 <- cut(y$prev.c2Gc1, b)
    y$bin.next.c1Gc2 <- cut(y$next.c1Gc2, b)
    y$bin.next.c2Gc1 <- cut(y$next.c2Gc1, b)
    y$concat <- paste(y$bin.prev.c1Gc2, y$bin.prev.c2Gc1, y$bin.next.c1Gc2, y$bin.next.c2Gc1)
  }
  y
}

# not really CV
study2 <- function(study1DF,returnDF=FALSE, includeCondi=TRUE, includeStatsTime=FALSE, includeStatsPairConcepts=FALSE,cv=FALSE) {
  x <- study1DF[study1DF$status!='unclear',c('prev.c1Gc2', 'prev.c2Gc1','next.c1Gc2', 'next.c2Gc1','status')]
  cols <- c('status')
  x$prevMin <- pmin(x$prev.c1Gc2, x$prev.c2Gc1)
  x$prevMax <- pmax(x$prev.c1Gc2, x$prev.c2Gc1)
  x$nextMin <- pmin(x$next.c1Gc2, x$next.c2Gc1)
  x$nextMax <- pmax(x$next.c1Gc2, x$next.c2Gc1)
  if (includeCondi) {
      cols <- c(cols, 'prevMin', 'prevMax', 'nextMin','nextMax')
  }
  if (includeStatsTime) {
    x$diffMinPrevNext <- x$nextMin - x$prevMin
    x$ratioMinPrevNext <- x$nextMin / x$prevMin
    x$diffMaxPrevNext <- x$nextMax - x$prevMax
    x$ratioMaxPrevNext <- x$nextMax / x$prevMax
    cols <- c(cols, 'diffMinPrevNext', 'ratioMinPrevNext','diffMaxPrevNext','ratioMaxPrevNext')
  }
  if (includeStatsPairConcepts) {
    x$diffPrevMinMax <- x$prevMax - x$prevMin
    x$ratioPrevMinMax <- x$prevMax / x$prevMin
    x$diffNextMinMax <- x$nextMax - x$nextMin
    x$ratioNextMinMax <- x$nextMax / x$nextMin
    cols <- c(cols, 'diffPrevMinMax', 'ratioPrevMinMax', 'diffNextMinMax', 'ratioNextMinMax')
  }
  d <- x[,cols]
  #  d2 <- x
#  colnames(d2)[colnames(d2)=='prev.c1Gc2'] <- 'prev.c2Gc1.tmp'
#  colnames(d2)[colnames(d2)=='prev.c2Gc1'] <- 'prev.c1Gc2'
#  colnames(d2)[colnames(d2)=='prev.c2Gc1.tmp'] <- 'prev.c2Gc1'
#  colnames(d2)[colnames(d2)=='next.c1Gc2'] <- 'next.c2Gc1.tmp'
#  colnames(d2)[colnames(d2)=='next.c2Gc1'] <- 'next.c1Gc2'
#  colnames(d2)[colnames(d2)=='next.c2Gc1.tmp'] <- 'next.c2Gc1'
#  d <- rbind(d,d2)
  if (returnDF) {
    d
  } else {
    if (cv) {
      r<-ldply(1:100, function(fold) {
        res <- trainTestDF(d)
        data.frame(accu=res$accuracy)        
      })
      print(paste("average accu:",mean(r$accu)))
    } else {
      res <- trainTestDF(d)
      print(res$conf_mat)
      print(paste('accuracy=',res$accuracy))
      rpart.plot(res$tree)
    }
  }
}


# returns a list 
trainTestDF <- function(d) {
  trainIdx <- sample(nrow(d),nrow(d)*.6)
  v<-1:nrow(d)
  trainBool <-  v %in% trainIdx
  testIdx <- v[!trainBool]
  tree <- rpart(status ~ ., data=d[trainIdx,])
  predicted = predict(tree, newdata=d[testIdx,], type="class")
  conf_mat <- table(predicted, d[testIdx,]$status)
  accu <- (conf_mat['discovery','discovery']+conf_mat['trivial','trivial']) / nrow(d[testIdx,])
  list(accuracy=accu, tree=tree, conf_mat=conf_mat)
}

pickpair <- function(df, categories=c('discovery','trivial')) {
  df0 <- df
  if (!is.null(categories)) {
    df0 <- df[df$status %in% categories,]
  }
  pairs <- unique(df0[,c('c1','c2')])
  n <- sample(nrow(pairs),1)
  pickedpair <- pairs[n,]
#  print("PAIR")
#  print(pickedpair)
  df0[df0$c1==pickedpair$c1 & df0$c2==pickedpair$c2,]
  
}

# df <- reformatAggregatedJointTopPMI(top_joint,aggregated_indiv,top_disc,withNames = TRUE)
#
multiplots_onecase <- function(df, top_discDF,top_disc_indivDF, sanitycheck_nbyears=71) {
  pairs <- unique(df[,c('c1','c2')])
  if (nrow(pairs)>1) {
    stop('BUG more than one pair of concepts')
  }
  df0 <- df[order(df$year),]
  print(cosine(replaceNAValues(df0$ma.c1),replaceNAValues(df0$ma.c2)))
  if (nrow(df) != sanitycheck_nbyears) {
    stop(paste('problem: expecting',sanitycheck_nbyears,'years'))
  }
  pair_disc <- top_discDF[top_discDF$source=='MED' & top_discDF$c1==pairs$c1 & top_discDF$c2==pairs$c2,]
#  print(pair_disc)
  n1 <- paste0(pair_disc$term1,' (',pair_disc$group1,")")
  n2 <- paste0(pair_disc$term2,' (',pair_disc$group2,")")
  name <- paste0(n1,' - ',n2,': ',pair_disc$status)

  indivDF <- melt(df,measure.vars = c('ma.c1','ma.c2'),value.name = 'ma.freq')
#  print(head(indivDF))
  firstYearDF <- ddply(indivDF,'variable',function(s) {
    nonzero <- s[!is.na(s$ma.freq) & s$ma.freq>0,]
    firstyear <- min(nonzero$year)
    nonzero[nonzero$year==firstyear,]
  })
  top_disc_indivDF <- top_disc_indivDF[top_disc_indivDF$aggreg.method=='ma',] 
  mapped_concepts <- ddply(firstYearDF,'variable',function(s) {
    if (s$variable == 'ma.c1') {
      data.frame(concept=s$c1)
    } else {
      data.frame(concept=s$c2)
    }
  })
  top_disc_indivDF<-merge(mapped_concepts, top_disc_indivDF,by='concept')
#  print(top_disc_indivDF)
  
  indivDF$ma.prob <- indivDF$ma.freq / indivDF$ma.total 
  indivDF$ma.condi <-indivDF$ma.joint / indivDF$ma.freq
  g_indiv_both <- ggplot(indivDF, aes(year,ma.freq,fill=variable))+geom_col(position='identity',alpha=.5) + geom_vline(data=pair_disc,aes(xintercept=disc.year,linetype=indicator)) + geom_vline(data=firstYearDF,aes(xintercept=year,colour=variable),linetype="dotted") +geom_vline(data=top_disc_indivDF,aes(xintercept=disc.year,colour=variable,linetype=indicator))+ theme(legend.position = "none")+xlab("")
#  g_indiv1 <- ggplot(df, aes(year,ma.c1))+geom_col()+geom_vline(data=pair_disc,aes(xintercept=disc.year,linetype=indicator))+ theme(legend.position = "none")
#  g_indiv2 <- ggplot(df, aes(year,ma.c2))+geom_col()+geom_vline(data=pair_disc,aes(xintercept=disc.year,linetype=indicator))+ theme(legend.position = "none")
  g_condi_both <- ggplot(indivDF, aes(year,ma.condi,colour=variable))+geom_line(alpha=.5,size=2)+geom_vline(data=pair_disc,aes(xintercept=disc.year,linetype=indicator))+ geom_vline(data=firstYearDF,aes(xintercept=year,colour=variable),linetype="dashed")+geom_vline(data=top_disc_indivDF,aes(xintercept=disc.year,colour=variable,linetype=indicator))+ theme(legend.position = "none")+xlab("")
  g_joint <- ggplot(df, aes(year, ma.joint))+geom_col()+xlab("")
  
  g_title <- ggdraw() + 
    draw_label(
      name,
      size=8
#      fontface = 'bold',
#      x = 0,
#      hjust = 0
    )
  plot_grid(g_title,
            g_indiv_both,
            g_condi_both,
            g_joint,
            labels = c('','Indiv','Condi','Joint'),
            label_x = 0.2,
            nrow = 4,
            rel_heights = c(0.15, 1,1,1))
}


study3 <- function(study1DF) {
  d <- study1DF
  d$prev.diff <- abs(d$prev.c1 - d$prev.c2)
  ddply(d, c('c1','c2','indicator','status'),function(s) {
    if (nrow(s)!=1) {
      stop('bug')
    }
    if (s$prev.c1 > s$prev.c2) {
      diff_cond_time_largest <- s$next.c1Gc2 - s$prev.c1Gc2
    } else {
      diff_cond_time_largest <- s$next.c2Gc1 - s$prev.c2Gc1
    }
    data.frame(diff_cond_time_largest=diff_cond_time_largest,prev.diff=s$prev.diff)
  })
}

study3prop <- function(study3DF) {
  d <- study3DF
  d$posi <- d$diff_cond_time_largest>0
  ddply(d,c('indicator','posi'),function(s) { 
    disc <- nrow(s[s$status=='discovery',])
    tot <- nrow(s)
    data.frame(disc=disc,tot=tot, prop_disc=disc/tot) 
  })
}


get_top_indiv <- function(aggJointDF, fullAggIndivDF) {
  fullIndivMed <- fullAggIndivDF[fullAggIndivDF$source=='MED',]
  targets <- unique(c(aggJointDF$c1, aggJointDF$c2))
  fullIndivMed[fullIndivMed$concept %in% targets,]
}


replaceNAValues <- function(v, replaceWith=0) {
  naValues <- is.na(v)
  v[naValues] <- rep(0,length(naValues[naValues]))
  v
}


similarityYears <- function(df, idCols=c('c1','c2'), valcol1='ma.c1',valcol2='ma.c2') {
  ddply(df, idCols, function(s) {
    if (nrow(s) != 71) {
      stop('bug')
    }
    val1 <- replaceNAValues(s[order(s$year),valcol1])
    val2 <- replaceNAValues(s[order(s$year),valcol2])
    data.frame(mae=mean(abs(val1-val2)),cos=cosine(val1,val2))
#    data.frame(cos=cosine(val1,val2))
  })
}

simThresholdDF <- function(discDFWithSim, returnF1DF=FALSE) {
  d <- discDFWithSim[discDFWithSim$status!='unclear',]
  values <- sort(unique(d$cos))
  res<-ldply(values, function(t) {
    c <- d[d$cos<=t,]
    nbDisc <- nrow(c[c$status=='discovery',])
    nbTriv <- nrow(c[c$status=='trivial',])
    data.frame(threshold=c(t,t),total=c(nrow(c),nrow(c)),status=c("discovery",'trivial'),nb=c(nbDisc, nbTriv),prop=c(nbDisc/nrow(c), nbTriv/nrow(c)))
  })
  optimsep<-ldply(values,function(t) {
    tp <- nrow(d[d$cos<=t & d$status=='discovery',])
    fp <- nrow(d[d$cos<=t & d$status!='discovery',])
    tn <- nrow(d[d$cos>t & d$status!='discovery',])
    fn <- nrow(d[d$cos>t & d$status=='discovery',])
    prec <- tp/(tp+fp)
    rec <- tp/(tp+fn)
    f1 <- 2 * prec * rec / (prec+rec)
    data.frame(threshold=c(t,t,t), score=c('prec','rec','f1'),perf=c(prec,rec,f1))
  })
  if (returnF1DF) {
    optimsep
  } else {
    res
  }
}

bestThreshold <- function(d,valCol='diff_cond_time_largest', positiveLowerThanThreshold=TRUE,positiveClass='discovery', labelCol='status') {
  values <- sort(unique(d[,valCol]))
  optimsep<-ldply(values,function(t) {
    tp <- nrow(d[d[,valCol]<=t & d[,labelCol]==positiveClass,])
    fp <- nrow(d[d[,valCol]<=t & d[,labelCol]!=positiveClass,])
    fn <- nrow(d[d[,valCol]>t & d[,labelCol]==positiveClass,])
    tn <- nrow(d[d[,valCol]>t & d[,labelCol]!=positiveClass,])
    if (!positiveLowerThanThreshold) {
      tmp <- tp
      tp <- fn
      fn <- tmp
      tmp <- tn
      tn <- fp
      fp <- tmp
    }
    prec <- tp/(tp+fp)
    rec <- tp/(tp+fn)
    f1 <- 2 * prec * rec / (prec+rec)
    data.frame(threshold=c(t,t,t), score=c('prec','rec','f1'),perf=c(prec,rec,f1))
  })
  maxf1<-max(optimsep[optimsep$score=='f1','perf'],na.rm = TRUE)
#  print(maxf1)
  maxt <- optimsep[!is.na(optimsep$perf) & optimsep$score=='f1' & optimsep$perf==maxf1 ,'threshold']
#  print(maxt)
  print(optimsep[optimsep$threshold==maxt,])
  optimsep
}

merge2columns <- function(mainDF, dictDF, mainDFcols=c('c1','c2'),dictDFcol='concept', extraIdCols=c(),suffixMain='.joint',suffixesPair=c('.c1','.c2')) {
  tmp <- merge(mainDF,dictDF,by.x=c(extraIdCols,mainDFcols[1]),by.y=c(extraIdCols,dictDFcol),suffixes=c(suffixMain,''))
  merge(tmp, dictDF,by.x=c(extraIdCols,mainDFcols[2]),by.y=c(extraIdCols,dictDFcol),suffixes=suffixesPair)
}

mergeAggJointAggIndiv <- function(agg_joint, agg_indiv,groupCols=c('source','half_window','year')) {
  # removing the ma.total col from indiv in order to avoid having it 3 times (it's the same in iindiv and joint)
  agg_indiv_no_total <- agg_indiv[,colnames(agg_indiv)[colnames(agg_indiv) !='ma.total'],]
  merge2columns(agg_joint, agg_indiv_no_total,extraIdCols = groupCols)
}


# df may contain several cases, default is to pick one at random
# df must contain columns all the ma columns:
# joint_rich_format<-mergeAggJointAggIndiv(aggregated_joint,aggregated_indiv)
# disc_joint and disc_indiv are optional, if used the disc years of the indiv/joint concepts are plotted
multiplot_joint <- function(joint_rich_format, sanitycheck_nbyears=71,pickRandomPair=TRUE,idCols=c('c1','c2'), termCols=c('term.c1','term.c2'), disc_joint=NULL,disc_indiv=NULL, assocGraphCols=NULL) {
  pairs <- unique(joint_rich_format[,idCols])
  if (nrow(pairs)>1) {
    if (pickRandomPair) {
      pairs <- unique(joint_rich_format[,idCols])
      n <- sample(nrow(pairs),1)
      pickedpair <- pairs[n,]
      joint_rich_format <- joint_rich_format[joint_rich_format[,idCols[1]]==pickedpair[,idCols[1]] & joint_rich_format[,idCols[2]]==pickedpair[,idCols[2]],]
    } else {
      stop('BUG more than one pair of concepts and pickRandomPair is FALSE')
    }
  }
  df0 <- joint_rich_format[order(joint_rich_format$year),]
  if (nrow(df0) != sanitycheck_nbyears) {
    stop(paste('problem: expecting',sanitycheck_nbyears,'years'))
  }
  if (!is.null(termCols) & termCols[1] %in% colnames(df0)) {
    name1 <- paste0(df0[1,idCols[1]],' [',df0[1,termCols[1]], ']')
    name2 <- paste0(df0[1,idCols[2]],' [',df0[1,termCols[2]], ']')
  } else {
    name1 <- paste0(df0[1,idCols[1]])
    name2 <- paste0(df0[1,idCols[2]])
  }
  name <- paste0(name1, ' - ', name2)
  
  indivDF <- melt(df0,measure.vars = c('ma.c1','ma.c2'),value.name = 'ma.freq')
  #  print(head(indivDF))
  firstYearDF <- ddply(indivDF,'variable',function(s) {
    nonzero <- s[!is.na(s$ma.freq) & s$ma.freq>0,]
    firstyear <- min(nonzero$year)
    nonzero[nonzero$year==firstyear,]
  })
  mapped_concepts <- ddply(firstYearDF,'variable',function(s) {
    if (s$variable == 'ma.c1') {
      data.frame(concept=s$c1)
    } else {
      data.frame(concept=s$c2)
    }
  })
  if (!is.null(disc_joint)){
    disc_joint <- disc_joint[disc_joint$source=='MED' & disc_joint$c1==pairs$c1 & disc_joint$c2==pairs$c2,]
  }
  if (!is.null(disc_indiv)) {
    disc_indiv <- disc_indiv[disc_indiv$aggreg.method=='ma',] 
    disc_indiv <- merge(mapped_concepts, disc_indiv,by='concept')
  }
  #  print(top_disc_indivDF)
  
  indivDF$ma.prob <- indivDF$ma.freq / indivDF$ma.total 
  indivDF$ma.condi <-indivDF$ma.joint / indivDF$ma.freq
  g_indiv_both <- ggplot(indivDF, aes(year,ma.freq,fill=variable))+geom_col(position='identity',alpha=.5)  + geom_vline(data=firstYearDF,aes(xintercept=year,colour=variable),linetype="dotted")
  g_condi_both <- ggplot(indivDF, aes(year,ma.condi,colour=variable))+geom_line(alpha=.5,size=2)+ geom_vline(data=firstYearDF,aes(xintercept=year,colour=variable),linetype="dotted")
  g_joint <- ggplot(df0, aes(year, ma.joint))+geom_col()+xlab("")
  if (!is.null(disc_joint)) {
    g_indiv_both <- g_indiv_both + geom_vline(data=disc_joint,aes(xintercept=disc.year,linetype=indicator))
    g_condi_both <- g_condi_both + geom_vline(data=disc_joint,aes(xintercept=disc.year,linetype=indicator))
  }
  if (!is.null(disc_indiv)) {
    g_indiv_both <- g_indiv_both  +geom_vline(data=disc_indiv,aes(xintercept=disc.year,colour=variable,linetype=indicator))
    g_condi_both <- g_condi_both  +geom_vline(data=disc_indiv,aes(xintercept=disc.year,colour=variable,linetype=indicator))
  }
  g_indiv_both <- g_indiv_both  + theme(legend.position = "none")+xlab("")
  g_condi_both <- g_condi_both  + theme(legend.position = "none")+xlab("")
  g_title <- ggdraw() + 
    draw_label(
      name,
      size=8
      #      fontface = 'bold',
      #      x = 0,
      #      hjust = 0
    )
  
  if (!is.null(assocGraphCols)) {
      assocDF <- melt(df0[,c('year',assocGraphCols)], id.vars='year',variable.name='assoc.measure')
    g_assoc <- ggplot(assocDF,aes(year,value,colour=assoc.measure))+geom_line()+ theme(legend.position = "none")+xlab("")
    plot_grid(g_title,
              g_indiv_both,
              g_condi_both,
              g_joint,
              g_assoc,
              labels = c('','Indiv','Condi','Joint','Assoc'),
              label_x = 0.2,
              nrow = 5,
              rel_heights = c(0.15, 1,1,1,1))
  } else {
  plot_grid(g_title,
            g_indiv_both,
            g_condi_both,
            g_joint,
            labels = c('','Indiv','Condi','Joint'),
            label_x = 0.2,
            nrow = 4,
            rel_heights = c(0.15, 1,1,1))
  }
}


#
# if normalized is FALSE returns the binary MI.
# if normalized is TRUE, returns a list (mi, nmi): the first is regular (binary) MI and the second is normalized MI (NMI)
#
binaryMI <- function(pA, pB, pA_B, normalized=FALSE) {
  pNA_B <- pB - pA_B
  pA_NB <- pA - pA_B
  pNA_NB <- 1 - ( pA_B + pNA_B + pA_NB)
  mi <- rep(0, length(pA_B))
  mi[pNA_NB>0] <- mi[pNA_NB>0] + pNA_NB[pNA_NB>0] * log2( pNA_NB[pNA_NB>0] / ((1-pA[pNA_NB>0]) * (1-pB[pNA_NB>0])) )
  mi[pNA_B>0] <- mi[pNA_B>0] + pNA_B[pNA_B>0] * log2( pNA_B[pNA_B>0] / ((1-pA[pNA_B>0]) * pB[pNA_B>0]) )
  mi[pA_NB>0] <- mi[pA_NB>0] + pA_NB[pA_NB>0] * log2( pA_NB[pA_NB>0] / (pA[pA_NB>0] * (1-pB[pA_NB>0])) )
  mi[pA_B>0] <- mi[pA_B>0] + pA_B[pA_B>0] * log2( pA_B[pA_B>0] / (pA[pA_B>0] * pB[pA_B>0]) )
  if (normalized) {
    norm <- rep(0, length(pA_B))
    norm[pNA_NB>0] <- norm[pNA_NB>0] + pNA_NB[pNA_NB>0] * log2( pNA_NB[pNA_NB>0]  )
    norm[pNA_B>0] <- norm[pNA_B>0] + pNA_B[pNA_B>0] * log2( pNA_B[pNA_B>0]  )
    norm[pA_NB>0] <- norm[pA_NB>0] + pA_NB[pA_NB>0] * log2( pA_NB[pA_NB>0]  )
    norm[pA_B>0] <- norm[pA_B>0] + pA_B[pA_B>0] * log2( pA_B[pA_B>0] )
    list(mi=mi, nmi=mi/-norm)
  } else {
    mi
  }
}

default_measures = c('scp', 'pmi', 'npmi', 'mi', 'nmi', 'pmi2', 'pmi3')

# - df must contain columns all the ma columns:
#   joint_rich_format <- mergeAggJointAggIndiv(aggregated_joint,aggregated_indiv)
# - there must not be any NA:
#   integratedRelationsDF <- joint_rich_format[!is.na(joint_rich_format$ma.total),]
#
calculateAssociation <- function(integratedRelationsDF, filterMeasures=default_measures,docFreq1Col='ma.c1',docFreq2Col='ma.c2', jointFreqCol='ma.joint',totalCol='ma.total') { 
  jointProb <- integratedRelationsDF[,jointFreqCol] / integratedRelationsDF[,totalCol]
  pA <- integratedRelationsDF[,docFreq1Col] / integratedRelationsDF[,totalCol]
  pB <- integratedRelationsDF[,docFreq2Col] / integratedRelationsDF[,totalCol]
  pmi <- log2( jointProb / (pA * pB) )
  if ('scp' %in% filterMeasures) {
    integratedRelationsDF$scp <- jointProb ^ 2 / (pA * pB)
  }
  if ('pmi' %in% filterMeasures) {
    integratedRelationsDF$pmi <- pmi
  }
  if ('npmi' %in% filterMeasures) {
    integratedRelationsDF$npmi <- - pmi / log2(jointProb)
  }
  if ('nmi' %in% filterMeasures) {
    l <- binaryMI(pA, pB, jointProb, normalized=TRUE)
    integratedRelationsDF$mi <- l$mi
    integratedRelationsDF$nmi <- l$nmi
  } else {
    if ('mi' %in% filterMeasures) {
      integratedRelationsDF$mi <- binaryMI(pA, pB, jointProb)
    }
  }
  if ('pmi2' %in% filterMeasures) {
    integratedRelationsDF$pmi2 <- log2( (jointProb^2) / (pA * pB) )
  }
  if ('pmi3' %in% filterMeasures) {
    integratedRelationsDF$pmi3 <- log2( (jointProb^3) / (pA * pB) )
  }
  integratedRelationsDF
}



computeSurgeIndicators2 <- function(aggregDF, valueCols='ma.joint', discardNonFiniteValues=FALSE) {
  minYear <- min(aggregDF$year)
  #  resDF <- data.frame()
  ldply(valueCols,function(valueCol) {
    resDF <- rbind(aggregDF,aggregDF)
    if (! valueCol %in% colnames(aggregDF)) {
      stop(paste0('Error: no column "',valueCol,'" (value column) in dataframe'))
    }
    values <- aggregDF[,valueCol]
    prevValues <- head(values,-1) # remove last
    nextValues <- tail(values,-1) # remove first
    rate <- (nextValues-prevValues) / prevValues
    rate <- c(NA, rate)
    diff <- nextValues-prevValues
    diff <- c(NA, diff)
    resDF$trend <- c(rate,diff)
    resDF$indicator <- c(rep(paste(valueCol,'rate',sep='.'),length(rate)),rep(paste(valueCol,'diff',sep='.'),length(rate)))
    if (discardNonFiniteValues) {
      resDF[is.finite(resDF$trend),]
    } else {
      resDF
    }
  })
}
  

showTopN <- function(df, N, idCols=c('source','half_window','indicator','disc.method'),orderBy='disc.trend') {
  ddply(df, idCols, function(s) {
    tail(s[order(s[,orderBy]),], N)
  })
}


studyDiscYears <- function(jointDF, discDF,mergeBy=c('c1','c2','half_window','source','indicator')) {
  discDF$disc.prevyear <- discDF$disc.year -1
  before <- merge(jointDF, discDF, by.x=c(mergeBy,'year'), by.y=c(mergeBy,'disc.prevyear'))
  before$disc.year <- NULL
  after <- merge(jointDF, discDF, by.x=c(mergeBy,'year'), by.y=c(mergeBy,'disc.year'))
  after$disc.prevyear <- NULL
  before$when <- rep('before', nrow(before))
  after$when <- rep('after', nrow(after))
  rbind(before,after)
}
