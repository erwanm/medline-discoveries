
library(data.table)
library(ggplot2)


# EM November 2021
#
# New version of the "discoveries" code:
#  - using data.table library for efficiency
#  - meant to be usable in a Shiny app for convenient visualization
#  - 3 main parts: detect surges across years; filter out trivial cases with high conditional prob; order by [N]PMI
#  - dropping support for sources other than Medline MeSH descriptors (KD, PTC)
#  - dropping "previous vs next" option (and possibly some indicator options as well)
#  - adding option to use other values than raw frequency for detecting surges (but preliminary experiments with PMI don't give good results) 
#

# Terminology: 'dynamic' means data by year ('year' column present), 'static' means data across all years (no 'year' column)
#


# returns the indiv or joint data as a data table, with result DT key set as the concept or pair of concepts
loadDynamicData <- function(dir='data/21-extract-discoveries/recompute-with-ND-group/MED',indivOrJoint='indiv',suffix='ND.min100', minYear=1950, maxYear=2020) {
  f <- paste(dir,paste(indivOrJoint,suffix,sep='.'),sep='/')
  if (indivOrJoint == 'indiv') {
    d<-fread(f,col.names=c('year','concept','freq'),drop=4)
    setkey(d, concept)
  } else {
    d<-fread(f,col.names=c('year','c1','c2','freq'))
    setkey(d, c1, c2)
  }
  d[d$year>=minYear & d$year<=maxYear,]
}


loadDynamicTotalFile <- function(dir='data/21-extract-discoveries/recompute-with-ND-group/MED', filename='indiv.full.total', minYear=1950, maxYear=2020) {
  f <- paste(dir, filename, sep='/')
  d<-fread(f, col.name=c('year','total'),drop=c(2,4))
  d[d$year>=minYear & d$year<=maxYear,]
}


# loads both the joint and indiv static data as well as the total file
#
# - By default (merged=TRUE) a single datatable is returned after merging the joint and indiv data, so that every pair
#   has the full information.
# - If merged=FALSE then the joint and indiv data are returned separately as a list.
#
loadStaticData <- function(dir='data/21-extract-discoveries/recompute-with-ND-group/MED/across-all-years/',jointFile='joint.min100.ND',indivFile='indiv.min100.ND.terms', totalFile='indiv.full.total',merged=TRUE) {
  joint <- fread(paste(dir, jointFile,sep='/'),col.names=c('c1','c2','freq'))
  setkey(joint,c1,c2)
  indiv <- fread(paste(dir, indivFile,sep='/'),col.names=c('concept','freq','term','group'),drop=3)
  setkey(indiv,concept)
  totaldocs <- read.table(paste(dir, totalFile,sep='/'))[1,1]
  if (merged) {
    mergeStaticJointWithIndivData(joint, indiv, totaldocs)
  } else {
    list(joint=joint,indiv=indiv,total=totaldocs)
  }
}


#
mergeStaticJointWithIndivData <- function(jointDT, indivDT, addTotalCol=NULL) {
  d0 <- merge(jointDT, indivDT, by.x='c2',by.y='concept', suffixes=c('.joint',''))
  d0 <- merge(d0, indivDT, by.x='c1', by.y='concept', suffixes=c('.c2','.c1'))
  if (!is.null(addTotalCol)) {
    d0[,total := addTotalCol,]
  }
  d0
}

pickRandomDynamic <- function(d, n=1, minFreqYear=0) {
  if (minFreqYear>0) {
    dmin <- d[freq>=minFreqYear,]
    u <- unique(dmin,by=key(dmin))
  } else {
    u <- unique(d,by=key(d))
  }
  print(paste('selecting among n pairs = ',nrow(u)))
  r <- pickOneRandomDynamic(u,d)
  if (n>1) {
    for (i in seq(2,n)) {
      r<-rbind(r,pickOneRandomDynamic(u,d))
    }
  }
  setkeyv(r,key(d))
  r
}

pickOneRandomDynamic <- function(uniq, data) {
  n <- sample(nrow(uniq),1)
  picked <- uniq[n,]
  data[c1==picked$c1 & c2==picked$c2,]
}

filterConcepts <- function(d, concepts) {
  d[c1 %in% concepts | c2 %in% concepts,]
}


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


# complete the sequence of years for every concept or pair: if a year between the start year and end year    
# is not present, it is added with frequency zero.
# additionally an option can be used to 'pad' before the first year with zeros (for moving average).
#
# from https://stackoverflow.com/questions/69815130/data-table-is-it-possible-to-merge-sd-and-return-a-new-sub-data-table-by-gro/
#
fillIncompleteYears <- function(dt,idCols=c('concept'),padBeforeStartYear=0,filterMinYear=NA) {
  # creates a sequence of years from startYear - padBeforeStartYear to endYear
  res<- dt[ dt[, .(year = seq(min(year)-padBeforeStartYear, max(year))), by = idCols],
           on = c(idCols,'year'),
        ][is.na(freq), freq := 0][]
  setkeyv(res,key(dt))
  if (!is.na(filterMinYear)) {
    res[year>=filterMinYear,]
  } else {
    res
  }
}


selectTotalYears <- function(totalsDT, minY,maxY) {
  totalsDT[year>=minY & year<=maxY,ma.total,]
}


# moving average for one of the dynamic tables, either joint or indiv (as mainDT). 
#
# - adds 2 columns: 'ma' (moving average) and 'ma.total' (total for the year)
# - window: odd number, window centered on the year
# - by default window=1, which means that 'ma' is the same as the raw frequency
# - The key of mainDT must have been set.
#
computeMovingAverage <- function(mainDT, totalsDT, window=1) {
  mainDT <- fillIncompleteYears(mainDT,idCols=key(mainDT),padBeforeStartYear=ceiling((window-1)/2)+1,filterMinYear=min(mainDT$year))
  totalsDT[,ma.total := ma(total, window, padWithNA = TRUE),]
  mainDT[,c('ma','ma.total') := list(ma(freq, window, padWithNA = TRUE), selectTotalYears(totalsDT,min(year),max(year))),by=key(mainDT)]
  mainDT[,ma := ma(freq, window, padWithNA = TRUE),by=key(mainDT)]
  mainDT[,ma.total := selectTotalYears(totalsDT,min(year),max(year)),by=key(mainDT)]
  mainDT
}


# returns x * log(y) if y>1, NA otherwise (in order to prevent negative log result)
calculateProdLog <- function(x,y) {
  res <- rep(NA, length(x))
  regularCase <- is.finite(y) & y > 1
  res[regularCase] <- x[regularCase] * log2(y[regularCase])
  res
}


# - input as data.table, modified by reference
# - possible indicators: 'prob.rate', 'prob.diff, 'prod.log'
#
# d <- computeMovingAverage(...)
#
computeTrend <- function(d, indicator='prob.rate') {
  d[,prob:=ma/ma.total,]
  d[, c('ma.prev', 'prob.prev') := list( c(NA,head(ma,-1)), c(NA,head(prob,-1)) ), by=key(d)]
  if ('prob.rate' == indicator) {
    d[, trend := (prob-prob.prev)/prob.prev,]
  }
  if ('prob.diff' == indicator) {
    d[, trend := prob-prob.prev,]
  }
  if ('prod.log' == indicator) {
    freqDiff <- d$ma - d$ma.prev
    probRate <- (d$prob - d$prob.prev) / d$prob.prev
    d[, trend := calculateProdLog(probRate, freqDiff),]
  }
  d
}


# returns the threshold for upper outliers for the values in v0: Q3 + 1.5 IQR
#
calculateThresholdTopOutliers <- function(v0, discardNegativeValues=FALSE) {
  if (discardNegativeValues) {
    v <- v0[!is.na(v0) & v0>0]
  } else {
    v <- v0[!is.na(v0)]
  }
  if (length(v)>0) {
    quartiles <- quantile(v,c(.25,.75))
    iqr <- quartiles[2] - quartiles[1] # Q3-Q1
    quartiles[2] + 1.5*iqr
  } else {
    NA
  }
}


#
# - outliers can be taken either globally (across all the 'trend' values) or locally (among the 'trend' values by key)
# trendDT <- computeTrend(..)
#
detectSurges <- function(trendDT, globalOutliers=FALSE, discardNegativeTrend=FALSE) {
  d <- trendDT
  if (globalOutliers) {
    d[, surge := is.finite(trend) & trend>=calculateThresholdTopOutliers(trend, discardNegativeValues=discardNegativeTrend),]
  } else {
    d[, surge := is.finite(trend) & trend>=calculateThresholdTopOutliers(trend, discardNegativeValues=discardNegativeTrend), by=key(d)]
  }
  d
}


computeAndSaveSurgesData <- function(dir='data/21-extract-discoveries/recompute-with-ND-group/MED', outputFilePrefix='./', indivOrJoint='joint', ma_windows=c(1,3,5),indicators=c('prob.rate','prod.log'), outlier_methods=c('local','global')) {
  dynamic_joint <- loadDynamicData(dir,indivOrJoint)
  dynamic_total <- loadDynamicTotalFile(dir)
  for (w in ma_windows) {
    for (i in indicators) {
      for (o in outlier_methods) {
        f <- paste0(outputFilePrefix,paste(w,i,o,'tsv',sep='.'))
        print(paste('processing and saving to', f))
        relations <- computeMovingAverage(dynamic_joint,dynamic_total, window=w)
        computeTrend(relations, indicator=i)
        if (o == 'local') {
          surges <- detectSurges(relations, globalOutliers=FALSE)
        } else {
          if (o == 'global') {
            surges <- detectSurges(relations, globalOutliers=TRUE)
          } else {
            stop('Error: invalid value for method_outliers, must be "global" or "local".')
          }
        }
        fwrite(surges[surge==TRUE,],f,sep='\t')
      }
    }
  }
}


#
# - 'oneSurgeByKey': 'no' for keeping all the surges, 'first' for picking the earliest, 'max' for picking the max trend
#
# surgesDT <- detectSurges(..)
#
filterSurges <- function(surgesDT, oneSurgeByKey='no',minTrend=-Inf) {
  d <- surgesDT[surge==TRUE & trend>=minTrend,]
  if (oneSurgeByKey != 'no') {
    if (oneSurgeByKey == 'max') {
      d<-d[,.SD[trend==max(trend),],by=key(d)]
    }
    if (oneSurgeByKey == 'first') {
      d<-d[,.SD[year==min(year),],by=key(d)]
    }
  }
  d
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

#
calculateAssociation <- function(dt, filterMeasures=default_measures,col1='freq.c1',col2='freq.c2',colJoint='freq.joint',colTotal='total') { 
  v1 <- as.name(col1)
  v2 <- as.name(col2)
  joint<-as.name(colJoint)
  total <- as.name(colTotal)
  dt[,prob.joint := eval(joint) / eval(total),]
  dt[,prob.c1 := eval(v1) / eval(total),]
  dt[,prob.c2 := eval(v2) / eval(total),]
  dt[,prob.C1GivenC2 := eval(joint) / eval(v2),]
  dt[,prob.C2GivenC1 := eval(joint) / eval(v1),]
  dt[,pmi := log2( prob.joint / (prob.c1 * prob.c2) ),]
  if ('scp' %in% filterMeasures) {
    dt[, scp := prob.joint^2 / (prob.c1*prob.c2),]
  }
  if ('npmi' %in% filterMeasures) {
    dt[,npmi := pmi / log2(prob.joint),]
  }
  if ('nmi' %in% filterMeasures) {
    l <- binaryMI(dt[,prob.c1], dt[,prob.c2], dt[,prob.joint], normalized=TRUE)
    dt[,mi := l$mi,]
    dt[,nmi := l$nmi,]
  } else {
    if ('mi' %in% filterMeasures) {
      dt[,mi := l$mi,]
      integratedRelationsDF$mi <- binaryMI(dt[,prob.c1], dt[,prob.c2], dt[,prob.joint])
    }
  }
  if ('pmi2' %in% filterMeasures) {
    dt[, pmi2 := log2( (prob.joint^2) / (prob.c1 * prob.c2) ), ]
  }
  if ('pmi3' %in% filterMeasures) {
    dt[, pmi3 := log2( (prob.joint^3) / (prob.c1 * prob.c2) ), ]
  }
}


# Merge dynamic datatable with static data and computes association measures.
#
# This function can be used with any dt as 'relations', but it is advisable to use it only with the
# smaller table obtained by filtering computeSurges:
#
# surges <- detectSurges(..)
# surges <- surges[surge==TRUE,]
# addStaticDataToRelations(surges, ...)
#
# The static data should be in the default 'merged' format obtained from loadStaticData;
# staticData <- loadStaticData(...)
#
addStaticAssociationToRelations <- function(relationsDT, staticData, filterMeasures = c('pmi','npmi')) {
  d <- merge(relationsDT, staticData, by=c('c1','c2'),suffixes=c('.dynamic','.static'))
  calculateAssociation(d,filterMeasures=filterMeasures)
  d
}



#
# relationsDT <- computeMovingAverage(dynamic_joint,dynamic_total, window=5)
# indivDT <- computeMovingAverage(dynamic_indiv,dynamic_total, window=5)
# 
addDynamicAssociationToRelations <- function(relationsDT, indivDT, filterMeasures = c('pmi','npmi')) {
  d <- merge(relationsDT, indivDT, by.x=c('year','ma.total','c2'),by.y=c('year','ma.total','concept'), suffixes=c('.joint',''))
  d <- merge(d, indivDT, by.x=c('year','ma.total','c1'), by.y=c('year','ma.total','concept'), suffixes=c('.c2','.c1'))
  calculateAssociation(d,filterMeasures=filterMeasures,col1='ma.c1',col2='ma.c2',colJoint = 'ma.joint',colTotal = 'ma.total')
  d
}


# pairData <- pickRandomDynamic(....)
# staticData <- loadStaticData(....)
displaySeveralPairsData <- function(pairsData, staticData, totalsDT,window=5,yearRange=c(1988,2018),excludeConceptsFromName=NULL,ncol=3) {
  maDT <- computeMovingAverage(pairsData,totalsDT,window=window)
  d <- merge(maDT, staticData,by=c('c1','c2'))
  d[,n1 := paste0(term.c1,' (',c1,")"),]
  d[,n2 := paste0(term.c2,' (',c2,")"),]
  d[,fulldescr := paste0(n1,' - ',n2),]
  print(unique(d[,fulldescr]))
  if (is.null(excludeConceptsFromName)) {
    d[,relation := paste0(term.c1,' - ',term.c2),]
  } else {
    d[,tmp1 := if (c1 %in% excludeConceptsFromName) "" else term.c1,by=c1]
    d[,tmp2 := if (c2 %in% excludeConceptsFromName) "" else term.c2,by=c2]
    d[,relation := paste0(tmp1,tmp2),]
  }
  ggplot(d,aes(year,freq))+geom_col(alpha=.4)+geom_line(aes(year,ma))+xlim(yearRange)+facet_wrap(.~relation,scales='free_y',ncol=ncol)
}

displaySurges <- function(pairsData, staticData, totalsDT,windows=c(1,3,5), indicators=c('prob.rate', 'prob.diff','prod.log'),withTitle=TRUE) {
  l <- list()
  for (window in windows) {
    for (indicator in indicators) {
      d <- computeMovingAverage(pairsData,totalsDT,window=window)
      computeTrend(d, indicator)  
      surges <- detectSurges(d, globalOutliers=FALSE)
      surges[,indicator:=indicator,]
      surges[,window:=window,]
      l[[length(l)+1]]<-surges
    }
  }
  r<-rbindlist(l)
  g <- ggplot(r,aes(year,prob,fill=surge))+geom_col()+facet_grid(window~indicator)+scale_fill_manual(values = c("#AAAAAA","#FF766D"))+ theme(legend.position="none")
  if (withTitle) {
    d <- merge(pairsData, staticData,by=c('c1','c2'))
    d[,n1 := paste0(term.c1,' (',c1,")"),]
    d[,n2 := paste0(term.c2,' (',c2,")"),]
    d[,fulldescr := paste0(n1,' - ',n2),]
    name <- unique(d[,fulldescr])
    print(name)
    g <- g +ggtitle(name)
  } else {
    g
  }
}


