
library(data.table)
library(ggplot2)
library(plyr)
library(cowplot)
library(scales)


# EM November 2021
#
# New version of the "discoveries" code:
#  - using data.table library for efficiency
#  - meant to be usable in a Shiny app for convenient visualization
#  - dropping support for sources other than Medline MeSH descriptors (KD, PTC)
#  - dropping "previous vs next" option (and possibly some indicator options as well)
#  - adding option to use other values than raw frequency for detecting surges (but preliminary experiments with PMI don't give good results) 
#

# Terminology: 'dynamic' means data by year ('year' column present), 'static' means data across all years (no 'year' column)
#

default_measures = c('prob.joint','scp', 'pmi', 'npmi', 'mi', 'nmi', 'pmi2', 'pmi3')


##### LOAD/SAVE FUNCTIONS

# returns the indiv or joint data as a data table, with result DT key set as the concept or pair of concepts
loadDynamicData <- function(dir='data/input',indivOrJoint='indiv',suffix='.min100.ND', minYear=1950, maxYear=2020) {
  f <- paste(dir,paste0(indivOrJoint,suffix),sep='/')
  if (indivOrJoint == 'indiv') {
    d<-fread(f,col.names=c('year','concept','freq'),drop=4)
    setkey(d, concept)
  } else {
    d<-fread(f,col.names=c('year','c1','c2','freq'))
    setkey(d, c1, c2)
  }
  d[d$year>=minYear & d$year<=maxYear,]
}


loadDynamicTotalFile <- function(dir='data/input', filename='indiv.total', minYear=1950, maxYear=2020) {
  f <- paste(dir, filename, sep='/')
  d<-fread(f, col.name=c('year','total'),drop=c(2,4))
  d[d$year>=minYear & d$year<=maxYear,]
}


# loads both the joint and indiv static data as well as the total file
#
# - By default (merged=TRUE) a single datatable is returned after merging the joint and indiv data, so that every pair
#   has the full information.
# - If merged=FALSE then the joint and indiv data are returned separately as a list.
# - By default the filenames are obtained using the suffix as follows:
#     - indiv<suffix>.terms
#     - joint<suffix>
#     - indiv.total
#   Otherwise the three filenames can be provided manually with filenames=c(indiv, joint, total). 'suffix' is ignored in this case.
#
loadStaticData <- function(dir='data/input/static/',suffix='.min100.ND',merged=TRUE, filenames=NULL) {
  if (is.null(filenames)) {
    indivFile <- paste0('indiv',suffix,'.terms')
    jointFile <- paste0('joint',suffix)
    totalFile <- 'indiv.total'
  } else {
    indivFile <- flienames[1]
    jointFile <- flienames[2]
    totalFile <- flienames[3]
  }    
  joint <- fread(paste(dir, jointFile,sep='/'),col.names=c('c1','c2','freq'))
  setkey(joint,c1,c2)
  indiv <- fread(paste(dir, indivFile,sep='/'),col.names=c('concept','freq','term','group'),drop=3)
  setkey(indiv,concept)
  totaldocs <- read.table(paste(dir, totalFile,sep='/'))[1,1]
  if (merged) {
    r<-mergeStaticJointWithIndivData(joint, indiv, totaldocs)
    setkey(r,c1,c2)
    r
  } else {
    list(joint=joint,indiv=indiv,total=totaldocs)
  }
}



computeAndSaveSurgesData <- function(dir='data/input', outputDir='data/output/', suffix='.min100.ND',ma_windows=c(1,3,5),measures=default_measures, indicators=c('rate','diff')) {
  if (!dir.exists(outputDir)) {
    print(paste('creating directory', outputDir))
    dir.create(outputDir,recursive = TRUE)
  }
  dynamic_joint <- loadDynamicData(dir,suffix,indivOrJoint = 'joint')
  # for debugging:
  # dynamic_joint<-pickRandomDynamic(dynamic_joint,n=1000)
  dynamic_indiv <- loadDynamicData(dir,suffix,indivOrJoint = 'indiv')
  dynamic_total <- loadDynamicTotalFile(dir)
  for (w in ma_windows) {
    joint.ma <- computeMovingAverage(dynamic_joint,dynamic_total, window=w)
    indiv.ma <- computeMovingAverage(dynamic_indiv,dynamic_total, window=w)
    for (m in measures) {
      for (i in indicators) {
        f <- paste(outputDir,paste0(paste(m,i,w,sep='.'),suffix,'.tsv'),sep='/')
        print(paste('processing and saving to', f))
        relations<-addDynamicAssociationToRelations(joint.ma,indiv.ma,measures = m)
        computeTrend(relations, indicator=i,measure=m)
        threshold <- calculateThresholdInflectionPoint(relations$trend)
        detectSurges(relations, globalThreshold=threshold)
        addNextNonZeroYear(relations)
        fwrite(relations[surge==TRUE,],f,sep='\t')
      }
    }
  }
}



loadSurgesData <- function(dir='data/output', suffix='.min100.ND', ma_window=1,measure='prob.joint', indicator='diff',dropMeasuresCols=FALSE) {
  f <- paste(dir,paste0(paste(measure,indicator,ma_window,sep='.'),suffix,'.tsv'),sep='/')
  print(f)
  if (dropMeasuresCols) {
    suppressWarnings(d<-fread(f, drop=default_measures))
  } else {
    d<-fread(f)
    
  }
  setkey(d,c1,c2)
  d
}


loadGoldDiscoveries <- function(dir='data',file='ND-discoveries-year.tsv') {
  fread(paste(dir,file,sep='/'),drop=c('source','notes'))
}


readMultipleSurgesFiles <- function(dir='data/output',suffix='.min100.ND',ma_windows=c(1,3,5),measures=default_measures, indicators=c('diff','rate'), firstYear=TRUE) {
  l <- list()
  for (w in ma_windows) {
    for (m in measures) {
      for (i in indicators) {
        #        print(paste('w=',w,'m=',m,'i=',i))
        d<-loadSurgesData(dir,suffix,w,m,i,dropMeasuresCols=TRUE)
        setkey(d,c1,c2)
        if (firstYear) {
          d <- d[,.SD[year==min(year)],by=key(d)]
        }
        d[,ma.window := w]
        d[,measure := m]
        d[,indicator:=i]
        l[[length(l)+1]]<-d
      }
    }
  }
  dt<-rbindlist(l)
  setkey(dt, c1, c2)
  dt
}




##### MAIN PROCESS

# moving average on a vector 'x' with window size 'window' (centered).
ma <- function(x,n=5,padWithNA=FALSE) {
  if (length(x)<n) { # not enough values to calculate 
    if (padWithNA) { 
      rep(NA, length(x))
    } else {
      stop('Error in ma: not enough values in vector ')
    }
  } else {
    cx <- c(0,cumsum(x))
    r <- (cx[(n+1):length(cx)] - cx[1:(length(cx) - n)]) / n
    if (padWithNA) {
      nas <- rep(NA, (n-1)/2)
      c(nas,r,nas)
    } else {
      r
    }
  }
}


# complete the sequence of years for every concept or pair: if a year between the start year and end year    
# is not present, it is added with frequency zero.
# additionally an option can be used to 'pad' before the first year with zeros (for moving average).
#
# from https://stackoverflow.com/questions/69815130/data-table-is-it-possible-to-merge-sd-and-return-a-new-sub-data-table-by-gro/
#
fillIncompleteYears <- function(dt,idCols=c('concept'),padBeforeStartYear=0,filterMinYear=NA) {
#  print(paste('DEBUG fillIncompleteYears:',min(dt$year),max(dt$year),padBeforeStartYear,filterMinYear))
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
  mainDT <- fillIncompleteYears(mainDT,idCols=key(mainDT),padBeforeStartYear=ceiling((window-1)/2)*2+1,filterMinYear=min(totalsDT$year))
  totalsDT[,ma.total := ma(total, window, padWithNA = TRUE),]
#  browser()
#  print(mainDT)
#  print(totalsDT)
#  mainDT[,c('ma','ma.total') := list(ma(freq, window, padWithNA = TRUE), selectTotalYears(totalsDT,min(year),max(year))),by=key(mainDT)]
  mainDT[,ma := ma(freq, window, padWithNA = TRUE),by=key(mainDT)]
  mainDT[,ma.total := selectTotalYears(totalsDT,min(year),max(year)),by=key(mainDT)]
  mainDT
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
    dt[,npmi := pmi / -log2(prob.joint),]
  }
  if ('nmi' %in% filterMeasures | 'mi' %in% filterMeasures) {
    nas <- is.na(dt[,prob.c1]) | is.na(dt[,prob.c2]) |  is.na(dt[,prob.joint])
    l <- binaryMI(dt[!nas,prob.c1], dt[!nas,prob.c2], dt[!nas,prob.joint], normalized=TRUE)
    if ('nmi' %in% filterMeasures) {
      res <- rep(NA,nrow(dt))
      res[!nas] <- l$nmi
      dt[,nmi := res,]
    }
    if ('mi' %in% filterMeasures) {
      res <- rep(NA,nrow(dt))
      res[!nas] <- l$mi
      dt[,mi := res,]
    }
  }
  if ('pmi2' %in% filterMeasures) {
    dt[, pmi2 := log2( (prob.joint^2) / (prob.c1 * prob.c2) ), ]
  }
  if ('pmi3' %in% filterMeasures) {
    dt[, pmi3 := log2( (prob.joint^3) / (prob.c1 * prob.c2) ), ]
  }
}



# - input as data.table, modified by reference
# - possible indicators: 'rate' or 'diff'
# if using ma freq as value, must pre-calculate prob:   d[,prob:=ma/ma.total,]
#
computeTrend <- function(d, indicator='rate',measure='prob.joint') {
  if (! indicator %in% c('rate','diff')) {
    stop(paste('Error: invalid indicator "',indicator,'"'))
  }
  val<-as.name(measure)
  d[, prev.value := c(NA,head(eval(val),-1)), by=key(d)]
  if ('rate' == indicator) {
    d[, trend := (eval(val)-prev.value)/abs(prev.value),]
  }
  if ('diff' == indicator) {
    d[, trend := eval(val)-prev.value,]
  }
#  if ('product' == indicator) {
#    d[, trend := (eval(val)-prev.value)^2/prev.value,]
#  }
}


# returns the threshold for upper outliers for the values in v0: Q3 + k IQR (k=3 by default, 1.5 if farOut=FALSE
#
calculateThresholdTopOutliers <- function(v0, discardNegativeValues=FALSE,farOut=TRUE) {
  if (discardNegativeValues) {
    v <- v0[!is.na(v0) & v0>0]
  } else {
    v <- v0[!is.na(v0)]
  }
  if (farOut) {
    k <- 3
  } else {
    k <- 1.5
  }
  if (length(v)>0) {
    quartiles <- quantile(v,c(.25,.75))
    iqr <- quartiles[2] - quartiles[1] # Q3-Q1
    quartiles[2] + k*iqr
  } else {
    NA
  }
}


minMaxScale <- function(values) {
  noNA <- values[!is.na(values)]
  maxtrend<-max(noNA)
  mintrend<-min(noNA)
  (noNA-mintrend)/(maxtrend-mintrend)
}

# Inflection point which maximizes the top left rectangle in the curve made of the 'flat' distribution.
# See analysis Rmd.
#
calculateThresholdInflectionPoint <- function(values) {
  fValues <- values[is.finite(values)]
  norma <- minMaxScale(fValues)
  relrank <- rank(norma, ties.method = 'random')/length(norma)
  areaQuantileTrend <- relrank*(1-norma)
  fValues[which(areaQuantileTrend==max(areaQuantileTrend))]
}

#
# - By default calculates the outlier threshold locally (i.e. for every concept/pair)
# - Recommended to provide a global threshold instead, typically by calculating it like this:
#
#   computeTrend(trendDT, ..)
#   detectSurges(trendDT, globalThreshold=calculateThresholdTopOutliers(trendDT$trend))
#
# - discardNegativeTrend is used only for the local outlier threshold, ignored otherwise
# 
#
detectSurges <- function(trendDT, globalThreshold=NA, discardNegativeTrend=FALSE) {
  d <- trendDT
  if (is.na(globalThreshold)) {
    d[, surge := (is.finite(trend) & (trend>=calculateThresholdTopOutliers(trend, discardNegativeValues=discardNegativeTrend))), by=key(d)]
  } else {
    d[, surge := (is.finite(trend) & (trend>=globalThreshold)),]
  }
  d
}


# OBSOLETE
# for a single concept/relation
# returns only the 'surge' column
adjustZeroFreqSurgesSingleRelation.OBSOLETE <- function(dt,window) {
  half_window <- ceiling((window-1)/2)
  zeroFreqSurgeYears <- dt[surge==TRUE & freq.joint==0,year]
  firstNonZeroYears <- unique(
    laply(zeroFreqSurgeYears, function(start_year) {
      nonZeroYearsInWindow <- dt[year>=start_year & year <= start_year+half_window & freq.joint>0,year]
      if (length(nonZeroYearsInWindow)>0) {
        min(nonZeroYearsInWindow)
      } else {
        NA
      }
    })
    ,na.rm=TRUE)
  newSurge <- dt$surge
  zeroFreqSurgeYearsIndexes <- dt[surge==TRUE & freq.joint==0, which=TRUE]
  newSurge[zeroFreqSurgeYearsIndexes] <- rep(FALSE,length(zeroFreqSurgeYearsIndexes))
  firstNonZeroIndexes <- dt[year %in% firstNonZeroYears,which=TRUE]
  newSurge[firstNonZeroIndexes] <- rep(TRUE,length(firstNonZeroIndexes))
  newSurge
}


# OBSOLETE TOO (too slow)
# for a single concept/relation
getNextNonZeroYear.DEPRECATED <- function(dt) {
  nz <- dt[freq.joint>0,year]
  laply(dt[,year], function(y0) { 
    l <- nz[nz>=y0]
    if (length(l)>0) {
      min(l) 
    } else {
      NA
    }
  })
} 
# receives a dt with colymns 'year' and 'freq.joint'
addNextNonZeroYear.DEPRECATED <- function(dt) {
  dt[,next.nonzero.year := getNextNonZeroYear(.SD),by=key(dt)]
}

# only for surge years!
addNextNonZeroYear <- function(d, maxiter=10) {
  # non zero data
  dnz <- d[freq.joint>0,.(year,c1,c2,freq.joint)]
  # zero frequency 
  dzf <- d[freq.joint==0 & surge==TRUE,.(year,c1,c2)]
  dzf[,year.adjusted:=year]
  print(paste('before while:',nrow(dzf)))
  i<-1
  while (i<maxiter & nrow(dzf)>0) {
    dzf[,year.adjusted:=year.adjusted+1]
    newdzf <- merge(dzf,dnz,by.x=c('c1','c2','year.adjusted'),by.y=c('c1','c2','year'))
    d <- merge(d, newdzf[freq.joint>0,.(year,year.adjusted,c1,c2)],by=c('c1','c2','year'),all.x=TRUE)
    dzf <- newdzf[freq.joint==0, .(year,year.adjusted,c1,c2)]
    print(paste('end iteration:',i,nrow(dzf)))
    i<-i+1
  }
  d
}



##### UTILITIES FUNCTIONS

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
  indexes <- sample(nrow(u),n)
  key <- key(d)
  picked <- u[indexes,..key]
  merge(d,picked)
}

filterConcepts <- function(d, concepts) {
  d[c1 %in% concepts | c2 %in% concepts,]
}


#
# relationsDT <- computeMovingAverage(dynamic_joint,dynamic_total, window=5)
# indivDT <- computeMovingAverage(dynamic_indiv,dynamic_total, window=5)
# 
addDynamicAssociationToRelations <- function(relationsDT, indivDT, measures = 'prob.joint') {
  d <- merge(relationsDT, indivDT, by.x=c('year','ma.total','c2'),by.y=c('year','ma.total','concept'), suffixes=c('.joint',''))
  d <- merge(d, indivDT, by.x=c('year','ma.total','c1'), by.y=c('year','ma.total','concept'), suffixes=c('.c2','.c1'))
  setkey(d,c1,c2)
  calculateAssociation(d,filterMeasures=measures,col1='ma.c1',col2='ma.c2',colJoint = 'ma.joint',colTotal = 'ma.total')
  d
}


# mostly for visualization purposes
#
addRelationName <- function(relationsDT, staticData,excludeConceptsFromName=NULL) {
  d <- merge(relationsDT, staticData,by=c('c1','c2'))
  d[,n1 := paste0(term.c1,' (',c1,")"),]
  d[,n2 := paste0(term.c2,' (',c2,")"),]
  d[,fulldescr := paste0(n1,' - ',n2),]
  #  print(unique(d[,fulldescr]))
  if (is.null(excludeConceptsFromName)) {
    d[,relation := paste0(term.c1,' - ',term.c2),]
  } else {
    d[,tmp1 := if (c1 %in% excludeConceptsFromName) "" else term.c1,by=c1]
    d[,tmp2 := if (c2 %in% excludeConceptsFromName) "" else term.c2,by=c2]
    d[,relation := paste0(tmp1,tmp2),]
  }
  d
  #  d[,c(colnames(relationsDT),'relation'), with = FALSE]
}


# returns a tidy version of the df after separating the elements in column 'colname'
# The rows which have multiple elements get duplicated, each with one of the elements
#
# source: https://stackoverflow.com/questions/13773770/split-comma-separated-strings-in-a-column-into-separate-rows
#
tidyUpMultiStringCol <- function(dt, colname, sep=' ') {
  otherCols <- colnames(dt)[colnames(dt) != colname]
  setDT(dt)[, lapply(.SD, function(x) unlist(tstrsplit(x, sep, fixed=TRUE))), by = otherCols]
  
}


# table must have group.c1 and group.c2 as single group: use tidyUpMultiStringCol() on both columns
#
selectRelationsGroups <- function(relationsDT, groups1=c('DISO','CHEM','GENE','ANAT'), groups2=c('DISO','CHEM','GENE','ANAT')) {
  relationsDT[group.c1 %in% groups1 & group.c2 %in% groups2,] 
}


filterCondProb <- function(surgesDT, conditional.threshold=.7, onlyFirstSurge=FALSE,yearMin=NA,yearMax=NA) {
  d <- surgesDT
  if (onlyFirstSurge) {
    d[,first.surge:=(year==min(year)),by=key(d)]
    d<-d[first.surge==TRUE,]
  }
  d<-d[prob.C1GivenC2<=conditional.threshold & prob.C2GivenC1<=conditional.threshold,]
  if (!is.na(yearMin)) {
    d<-d[year>=yearMin,]
  }
  if (!is.na(yearMax)) {
    d<-d[year<=yearMax,]
  }
  d
}

roundIfNotZero <- function(x,asPercentage) {
  res <- rep('',length(x))
  if (asPercentage) {
    res[x>0] <- round(x[x>0],digits=2)*100
  } else {
    res[x>0] <- round(x[x>0],digits=2)
  }
  res
}



##### VISUALIZATION


# pairData <- pickRandomDynamic(....)
# staticData <- loadStaticData(....)
displaySeveralPairsData <- function(pairsData, staticData, totalsDT,window=5,yearRange=c(1988,2018),excludeConceptsFromName=NULL,ncol=3) {
  maDT <- computeMovingAverage(pairsData,totalsDT,window=window)
  d <- addRelationName(maDT, staticData, excludeConceptsFromName)
  ggplot(d,aes(year,freq))+geom_col(alpha=.4)+geom_line(aes(year,ma))+xlim(yearRange)+facet_wrap(.~relation,scales='free_y',ncol=ncol)
}

displayMultiMeasure <- function(pairsData,indivData, totals, staticData,windows=c(1,3,5),measures=c('prob.joint','pmi'),yearRange=c(1988,2018),excludeConceptsFromName=NULL) {
  l <- list()
  for (window in windows) {
    ma.joint <- computeMovingAverage(pairsData,totals, window=window)
    ma.indiv <- computeMovingAverage(indivData,totals, window=window)
    d<-addDynamicAssociationToRelations(ma.joint,ma.indiv,measures = measures)
    d <- addRelationName(d, staticData, excludeConceptsFromName)
    d[,window:=window]
    l[[length(l)+1]]<-d
  }
  r<-rbindlist(l)
  cols <- c('year','relation','window',measures)
  r <- melt(r[,..cols],id.vars = c('year','relation','window'),variable.name='measure')
  ggplot(r,aes(year,value,alpha=as.factor(window)))+geom_line()+facet_grid(measure ~ relation,scales='free_y')+xlim(yearRange)+ scale_alpha_ordinal(range = c(0.3, 1))
}


displayMultiTrend <- function(pairsData,indivData, totals, staticData, indicator='rate',window=3,measure='prob.joint',yearRange=c(1988,2018),excludeConceptsFromName=NULL) {
  ma.joint <- computeMovingAverage(pairsData,totals, window=window)
  ma.indiv <- computeMovingAverage(indivData,totals, window=window)
  d<-addDynamicAssociationToRelations(ma.joint,ma.indiv,measures = measure)
  d <- addRelationName(d, staticData, excludeConceptsFromName)
  computeTrend(d,indicator,measure)
  raw <- d[,c('year','relation')]
  raw[,var:=measure]
  raw[,value:=d[,eval(as.name(measure))]]
  trend <- d[,c('year','relation')]
  trend[,var:='trend']
  trend[,value:=d[,trend]]
  r<-rbind(raw,trend)
  ggplot(r,aes(year,value))+geom_col()+facet_grid(var~relation,scales='free_y')+xlim(yearRange)
}


# receives a data table with ma  already calculated
# returns a lis of 4 graphs with/without log for x/y 
displayTrendDistribution1 <- function(ma.joint, ma.indiv, indicators=c('diff','rate'),measures=c('prob.joint','pmi','npmi','mi','nmi')) {
    l <- list()
    thresholds <- list()
    for (measure in measures) {
      for (indicator in indicators) {
        d<-addDynamicAssociationToRelations(ma.joint,ma.indiv,measures = measure)
        computeTrend(d,indicator,measure)
        l[[length(l)+1]] <- data.table(trend=d$trend, measure=measure,indicator=indicator)
        t1 <- calculateThresholdTopOutliers(d$trend,farOut = FALSE)
        t2 <- calculateThresholdTopOutliers(d$trend,farOut = TRUE)
        thresholds[[length(thresholds)+1]] <- data.table(threshold=c(t1,t2),outlier.type=c('k1.5', 'k3.0'),measure=measure,indicator=indicator)
      }
    }
    r<-rbindlist(l)
    t <- rbindlist(thresholds)
    graphlist <- list()
    g <- ggplot(r,aes(trend))+geom_histogram()+facet_wrap(measure~indicator,scales='free')+geom_vline(data=t,aes(xintercept=threshold,colour=outlier.type))
    graphlist[['x0.y0']] <- g
    graphlist[['x0.y1']] <- g+scale_y_log10()
    graphlist[['x1.y0']] <- g+scale_x_log10()
    graphlist[['x1.y1']] <- g+scale_x_log10()+scale_y_log10()
    graphlist
}


# receives a data table with ma  already calculated
displayTrendDistribution2 <- function(ma.joint, ma.indiv, indicators=c('diff','rate'),measures=c('prob.joint','pmi','npmi','mi','nmi'),withRect=TRUE, normaTrend=TRUE,withThresholds=TRUE) {
  l <- list()
  thresholds <- list()
  rects <- list()
  for (measure in measures) {
#    print(measure)
    for (indicator in indicators) {
#      print(indicator)
      d<-addDynamicAssociationToRelations(ma.joint,ma.indiv,measures = measure)
      computeTrend(d,indicator,measure)
      d<-d[is.finite(trend),]
      d[,trend.norma:=minMaxScale(trend)]
      d[,relrank:=rank(trend.norma, ties.method = 'random')/nrow(d)]
      d[,areaQuantileTrend := relrank*(1-trend.norma),]
      t1 <- calculateThresholdTopOutliers(d$trend,farOut = FALSE)
      rr1 <- d[abs(d$trend-t1)==min(abs(d$trend-t1)),relrank]
      t2 <- calculateThresholdTopOutliers(d$trend,farOut = TRUE)
      rr2 <- d[abs(d$trend-t2)==min(abs(d$trend-t2)),relrank]
      t3 <- calculateThresholdInflectionPoint(d$trend)
      rr3 <- d[abs(d$trend-t3)==min(abs(d$trend-t3)),relrank]
      thresholds[[length(thresholds)+1]] <- data.table(threshold=c(rr1,rr2,rr3),outlier.type=c('k1.5', 'k3.0','infl.pt.'),measure=measure,indicator=indicator)
      l[[length(l)+1]] <- data.table(trend=d$trend,trend.norma=d$trend.norma,relrank=d$relrank, measure=measure,indicator=indicator)
      if (withRect) {
        p<-d[areaQuantileTrend==max(areaQuantileTrend),]
        rects[[length(rects)+1]] <- data.table(xmin=0,xmax=p$relrank,ymin=p$trend.norma,ymax=1,measure=measure,indicator=indicator)
      }
    }
  }
  r<-rbindlist(l)
  t <- rbindlist(thresholds)
  rect<-rbindlist(rects)
  if (normaTrend) {
    yVar <- 'trend.norma'
  } else {
    yVar <- 'trend'
  }
  g <- ggplot(r,aes_string('relrank',yVar))+geom_point()+facet_grid(measure~indicator,scales = 'free_y')
  if (withThresholds) {
    g <- g+geom_vline(data=t,aes(xintercept=threshold,colour=outlier.type),size=1.5,linetype = "longdash")
  }
  if (withRect) {
    g<-g+geom_rect(data=rect,aes(NULL,NULL,xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,alpha=.1))+ guides(alpha = 'none')
  }
  g
}


displayTwoPlotsDistribTrend <- function(ma.joint, ma.indiv, indicator='diff',measure='nmi',withQuantileTest=FALSE,fontsize=16) {
  d<-addDynamicAssociationToRelations(ma.joint,ma.indiv,measures = measure)
  computeTrend(d,indicator,measure)
  d<-d[is.finite(trend),]
  d[,quantile:=rank(trend, ties.method = 'random')/nrow(d)]
  g1 <- ggplot(d,aes(trend))+geom_histogram(bins=60)+scale_y_log10(labels = scales::comma_format(accuracy=1,big.mark = ",", decimal.mark = "."))+theme(text=element_text(size=fontsize))
  g2 <- ggplot(d,aes(quantile,trend))+geom_point(size=.75)+theme(text=element_text(size=fontsize))
  if (withQuantileTest) {
    qdt <- data.table(q=seq(0,1,.05),v=quantile(d$trend,probs=seq(0,1,.05)))
    g2 <- g2 + +geom_point(data=qdt, aes(q,v),colour='red')
  }
  plot_grid(g1,
            g2,
            labels = NULL,
            label_x = 0.2,
            ncol = 2)
}


#
# single relation
#
displaySurges <- function(pairsData, indivData, totals, staticData,valueCol='prob.joint',windows=c(1,3,5), indicators=c('diff', 'rate'),withTitle=TRUE,excludeConceptsFromName=NULL) {
  if (nrow(unique(pairsData,by=key(pairsData)))>1) {
    stop('Error: "displaySuges" cannot handle more than one relation')
  }
  l <- list()
  for (window in windows) {
    for (indicator in indicators) {
      ma.joint <- computeMovingAverage(pairsData,totals, window=window)
      ma.indiv <- computeMovingAverage(indivData,totals, window=window)
      d<-addDynamicAssociationToRelations(ma.joint,ma.indiv,measures = valueCol)
      d <- addRelationName(d, staticData, excludeConceptsFromName)
      computeTrend(d,indicator,valueCol)
      surges <- detectSurges(d, globalThreshold=NA)
      surges[,indicator:=indicator,]
      surges[,window:=window,]
      l[[length(l)+1]]<-surges
    }
  }
  r<-rbindlist(l)
  g <- ggplot(r,aes_string('year',valueCol,fill='surge'))+geom_col()+facet_grid(window~indicator)+scale_fill_manual(values = c("#AAAAAA","#FF766D"))+ theme(legend.position="none")
  if (withTitle) {
#    d <- merge(pairsData, staticData,by=c('c1','c2'))
#    d[,n1 := paste0(term.c1,' (',c1,")"),]
#    d[,n2 := paste0(term.c2,' (',c2,")"),]
#    d[,fulldescr := paste0(n1,' - ',n2),]
#    name <- unique(d[,fulldescr])
#    print(name)
    name <- d[1,relation]
    g <- g +ggtitle(name)
  }
  g
}


displayAvgPerfByParam <- function(evalDT,fontsize=14, evalWindow=3) {
  d <- evalDT[is.na(eval.at) & eval.window==evalWindow,]
  # window
  x<-perfByParameter(d,'ma.window')
  x[,ma.window := as.factor(ma.window)]
  g1 <- ggplot(x,aes(ma.window,V1,fill=mode))+geom_col(position='dodge')+theme(text=element_text(size=fontsize),legend.position="none")+ylab('Mean recall')
  # measure
  x<-perfByParameter(d,'measure')
  x$measure[x$measure=='prob.joint'] <- 'prob'
  g2<-ggplot(x,aes(measure,V1,fill=mode))+geom_col(position='dodge')+theme(text=element_text(size=fontsize),legend.position="none")+ylab(NULL)
  # indicator
  x<-perfByParameter(d,'indicator')
  g3<-ggplot(x,aes(indicator,V1,fill=mode))+geom_col(position='dodge')+theme(text=element_text(size=fontsize))+ylab(NULL)
  plot_grid(g1,
            g2,
            g3,
            labels = NULL,
            label_x = 0.2,
            ncol = 3)
  
}



heatMapCommonMatrix <- function(m, varname='overlap',asPercentage=TRUE) {
  d<-reshape2::melt(m,value.name=varname) 
  d$Var1 <- factor(d$Var1, levels = unique(d$Var1[order(d$Var1, d$Var1)]))
  d$Var2 <- factor(d$Var2, levels = unique(d$Var2[order(d$Var2, d$Var2)]))
  ggplot(d, aes(Var1, Var2)) +
    geom_tile(aes_string(fill = varname)) + geom_text(aes(label = roundIfNotZero(overlap,asPercentage))) + scale_fill_gradient(low = "white", high = "red",labels = percent) + xlab(NULL)+ylab(NULL)+theme(axis.text.x=element_text(angle = -90, hjust = 0),text=element_text(size=20),legend.position = c(.15, .75))
}


#plotSurgesAcrossTime <- function(dir='data/input',window=3,measure='scp',indicator='diff',bins=30,fontsize=14) {
#  d<-loadSurgesData(dir, ma_window = window,measure=measure,indicator=indicator)
#  d[,first.surge:=(year==min(year)),by=key(d)]
#  ggplot(d,aes(year,fill=first.surge))+geom_histogram(position='stack',bins=bins)+theme(text=element_text(size=fontsize))
#}


doublePlotAcrossTime <- function(dynamicJointDT, surgesDT,bins=71,fontsize=14,marginAdjustMm=2,withLegend=TRUE) {
  firstcooc <- dynamicJointDT[,.SD[year==min(year),],by=key(dynamicJointDT)]
  plot.cooc <- ggplot(firstcooc,aes(year))+geom_histogram(bins=bins)+theme(text=element_text(size=fontsize),plot.margin = margin(0, marginAdjustMm, 0, 0, "mm"))+xlim(c(1950,2020))+xlab(NULL)+ylab(NULL)+ggtitle('First cooccurrences')
  surgesDT[,first.surge:=(year==min(year)),by=key(surgesDT)]
  plot.surges <- ggplot(surgesDT,aes(year,fill=first.surge))+geom_histogram(position='stack',bins=bins)+xlim(c(1950,2020))+xlab(NULL)+ylab(NULL)+ggtitle('Surges')
  if (withLegend) {
    plot.surges <- plot.surges +theme(text=element_text(size=fontsize),legend.position = c(.88, .75),plot.margin = margin(0, 0, 0, marginAdjustMm, "mm"))
  } else {
    plot.surges <- plot.surges +theme(text=element_text(size=fontsize),legend.position = 'none',plot.margin = margin(0, 0, 0, marginAdjustMm, "mm"))
  }
  plot_grid(plot.cooc,
            plot.surges,
            labels = NULL,
            label_x = 0.2,
            nrow = 2)
}

# calculates also for joint if not null
calculateDiffYears <- function(surgesDT, indivDT,jointDT=NULL) {
  firstocc <- indivDT[,.SD[year==min(year),],by=key(indivDT)]
  firstocc[,freq:=NULL]
  d<-merge(surgesDT,firstocc,by.x='c1',by.y='concept',suffixes=c('','.first.c1'))
  d<-merge(d,firstocc,by.x='c2',by.y='concept',suffixes=c('','.first.c2'))
  d[,year.first.both:=pmax(year.first.c1,year.first.c2)]
  d[,duration:=year-year.first.both]
  if (!is.null(jointDT)) {
    firstjoint <- jointDT[,.SD[year==min(year),],by=key(jointDT)]
    firstjoint[,freq:=NULL]
    d <- merge(d, firstjoint,by=c('c1','c2'),suffixes=c('','.first.joint')) 
    d[,duration.joint:=year-year.first.joint]
  }
  d
}

plotDiffYears <- function(surgesDT, indivDT,bins=71,fontsize=14) {
  surgesDT[,first.surge:=(year==min(year)),by=key(surgesDT)]
  d <- calculateDiffYears(surgesDT, indivDT)
  setnames(d,'duration','years')
  ggplot(d,aes(years,fill=first.surge))+geom_histogram(position='stack',bins=bins)+ggtitle('Time before surge')+theme(text=element_text(size=fontsize),legend.position = 'bottom',legend.direction = "vertical")+ylab(NULL)
}

threePlotsAboutYears <- function(surgesDT,indivDT,jointDT,bins=71,fontsize=14,marginAdjustMm=2) {
  g1 <- doublePlotAcrossTime(jointDT, surgesDT,bins,fontsize,marginAdjustMm,withLegend = FALSE)
  # the -1 is a hack, for some reason the graph looks better like this
  g2 <- plotDiffYears(surgesDT,indivDT, bins-1,fontsize)
  plot_grid(g1,
            g2,
            labels = NULL,
            rel_widths = c(2,1),
            ncol = 2)
}


# surgesDT <- loadSurgesData(...)
surgesAmongRelationsByFreq <- function(surgesDT, static_data, bins=100,freqmax=100000,xLabel='joint frequency', fontsize=11, asProportion=FALSE, logY=TRUE) {
  d<-surgesDT
  d[,first.surge:=(year==min(year)),by=key(d)]
  d<-d[first.surge==TRUE,]
  x<-merge(d, static_data,suffixes=c('.dynamic','.static'),all.y=TRUE)
  x[,surge:=!is.na(year),]
  if (asProportion) {
    mypos ='fill'
    yLabel='proportion'
  } else {
    mypos = 'stack'
    yLabel='count'
  }
  g <- ggplot(x, aes(freq.joint.static,fill=surge))+geom_histogram(position=mypos,bins=bins)+xlim(0,freqmax)+ylab(yLabel)+xlab(xLabel)+theme(text=element_text(size=fontsize),legend.position=c(.8,.8))
  if (logY) {
    g <- g +scale_y_log10(labels = scales::comma_format(accuracy=1,big.mark = ",", decimal.mark = "."))
  }
  g
}


###### EVALUATION

countSurgesByRelation <- function(surgesDT) {
  surges.by.key <- surgesDT[,.(n.surges=length(surge[surge])),by=key(surgesDT)]
  r<-surges.by.key[,.(n=.N,prop=.N/nrow(surges.by.key)),by=n.surges]
  setorder(r,n.surges)
  r
}


statsSurges <- function(total.pairs,total.rel, dir='data/input', ma_windows=c(1,3,5),measures=c('prob.joint','pmi','npmi','mi','nmi','scp'), indicators=c('diff','rate')) {
  l <- list()
  for (w in ma_windows) {
    for (m in measures) {
      for (i in indicators) {
        #        print(paste('w=',w,'m=',m,'i=',i))
        d<-loadSurgesData(dir,w,m,i)
        prop.pairs <- nrow(d) / total.pairs
        prop.rel <- nrow(unique(d,by=key(d)))/total.rel
        #        print(paste('pairs=',nrow(d),'rel=',nrow(unique(d,by=key(d)))))
        l[[length(l)+1]]<-data.table(window=w,measure=m,indicator=i,prop.pairs=prop.pairs,prop.rel=prop.rel)
      }
    }
  }
  r<-rbindlist(l)
  print('stats pairs:')
  print(summary(r$prop.pairs))
  print('stats relations:')
  print(summary(r$prop.rel))
  r
}



#
# matches the gold discoveries with a surges DT and returns a DT which has exactly one row for every
# gold discovery, and the closest year found in the surges (NA if not found at all).
# In order to match using the first year surge by relation, the surges DT should be filtered first.
#
matchSurgesWithGold <- function(surgesDT, goldDT, selectedCols=c('c1','c2','year.gold','year.pred')) {
  merged <- merge(goldDT, surgesDT, by=c('c1','c2'),suffixes=c('.gold','.pred'),all.x=TRUE)
  merged[,diff := abs(year.gold - year.pred)]
  matched <- merged[,.SD[is.na(year.pred) | diff==min(diff)],by=c('c1','c2')]
  matched <- matched[,.SD[1,],by=c('c1','c2')] # this is in case two rows are selected because e.g. both Y-1 and Y+1 have diff=1
  matched[,..selectedCols]
}


collectEvalDataSurgesAgainstGold <- function(goldDT, dir='data/input',suffix='.min100.ND',evalAt=NA,ma_windows=c(1,3,5),measures=c('prob.joint','pmi','npmi','mi','nmi','scp'), indicators=c('diff','rate')) {
  l <- list()
  for (w in ma_windows) {
    for (m in measures) {
      for (i in indicators) {
#        print(paste('w=',w,'m=',m,'i=',i))
        d0<-loadSurgesData(dir,suffix,w,m,i)
        setkey(d0,c1,c2)
        originalSize <- nrow(d0)
        for (maxSize in evalAt) {
          d <- copy(d0)
          # first year version
          fy <- d[,.SD[year==min(year)],by=key(d)]
          if (!is.na(maxSize)) {
            d <- head(d[order(-trend),],maxSize)
            fy <- head(fy[order(-trend),],maxSize)
          }
          # regular version
          e1 <- matchSurgesWithGold(d,goldDT)
          e1[,mode:='any.year',]
          e2 <- matchSurgesWithGold(fy,goldDT)
          e2[,mode:='first.year',]
          this <- rbind(e1,e2)
          this[,ma.window := w]
          this[,measure := m]
          this[,indicator:=i]
          this[,original.size:=originalSize]
          this[,eval.at:=maxSize]
          l[[length(l)+1]]<-this
        }
      }
    }
  }
  rbindlist(l)
}


evalSingleCase <- function(dt, eval_windows,proportion=TRUE) {
  l <- list()
  for (w in eval_windows) {
    n <- nrow(dt[!is.na(year.pred) & abs(year.gold-year.pred)<=w,])
    if (proportion) {
      n <- n / nrow(dt)
    }
    l[[length(l)+1]]<-data.table(perf=n,eval.window=w)
  }
  rbindlist(l)
}

#   d <- collectEvalDataSurgesAgainstGold(..)
evalSurgessAgainstGold <- function(d,eval_windows) {
  setkey(d,ma.window,measure,indicator,mode,eval.at)
  d[,evalSingleCase(.SD, eval_windows),by=key(d)]
}

# perfDT <- evalSurgessAgainstGold(...)
perfByParameter <- function(perfDT, param='ma.window') {
  perfDT[,mean(perf),by=c('mode','eval.at',param)]
}


# multiParamSurgesDT <- readMultipleSurgesFiles(...)
correlationMatrixMultiParam <- function(multiParamSurgesDT, yearWindow=NA) {
  if (length(unique(multiParamSurgesDT$indicator))==1) {
    multiParamSurgesDT[,config:=paste(measure,ma.window,sep='.')]
  } else {
    multiParamSurgesDT[,config:=paste(measure,ma.window,indicator,sep='.')]
  }
  buildCommonMatrix(multiParamSurgesDT,yearWindow)
}



calculateOverlapTwoConfigs <- function(multiSurgesDT, config1, config2,yearWindow=NA) {
  d1 <- multiSurgesDT[config==config1,]
  d2 <- multiSurgesDT[config==config2,]
  merged <- merge(d1,d2,by=c('c1','c2'))
  if (!is.na(yearWindow)) {
    merged <- merged[abs(year.x-year.y)<=yearWindow,]
  }
  nrow(merged)/min(nrow(d1),nrow(d2))
}

# requires 'config' column
# note: the matrix is long to compute, it can be saved with write.table(m,filename)
buildCommonMatrix <-function(d, yearWindow=NA) {
  configs <- sort(unique(d[,config]))
  l <- list()
  for (i in 1:(length(configs)-1)) {
    cat(i)
    cat(': ')
    v <- rep(0,i)
    for (j in (i+1):length(configs)) {
      cat(j)
      cat(' ')
      x <- calculateOverlapTwoConfigs(d, configs[i],configs[j],yearWindow)
      v <- c(v,x)
    }
    cat('\n')
    l[[i]] <- v
  }
  l[[length(configs)]] <- rep(0,length(configs))
  names(l) <- configs[1:(length(configs))]
  x<-as.data.frame(l)
  rownames(x) <- configs
  as.matrix(x)
}


