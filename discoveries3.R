
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
  indiv <- fread(paste(dir, indivFile,sep='/'),col.names=c('concept','freq','term','group'),drop=3)
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
  totalsDT[,ma.total := ma(total, window, padWithNA = TRUE),]
#  mainDT[,c('ma','ma.total') := list(ma(freq, window, padWithNA = TRUE), selectTotalYears(totalsDT,min(year),max(year))),by=key(mainDT)]
  mainDT[,ma := ma(freq, window, padWithNA = TRUE),by=key(mainDT)]
  mainDT[,ma.total := as.double(selectTotalYears(totalsDT,min(year),max(year))),by=key(mainDT)]
}


# returns x * log(y) if y>1, NA otherwise (in order to prevent negative log result)
calculateProdLog <- function(x,y) {
  res <- rep(NA, length(x))
  regularCase <- is.finite(y) & y > 1
  res[regularCase] <- x[regularCase] * log2(y[regularCase])
  res
}


# new version:
# - fixing bug about possible wrong shift previous/next year in previous version
# - input as data.table, modified by reference
# - combines the 2 steps of computing the 'trend' and selecting outliers as in 'computeSurgeYears'
# - possible indicators: 'prob.rate', 'prob.diff, 'prod.log'
# - outliers can be taken either globally (across all the 'trend' values) or locally (among the 'trend' values by key)
#
# d <- computeMovingAverage(...)
#
computeSurgeYears <- function(d, indicator='prob.rate', globalOutliers=FALSE) {
  d[,prob:=ma/ma.total,]
  d[, c('ma.prev', 'prob.prev') := list( c(NA,head(ma,-1)), c(NA,head(prob,-1)) ), by=key(d)]
  if ('prob.rate' == indicator) {
    d[, trend := (prob-prob.prev)/prob.prev,]
  }
  if ('prob.diff' == indicator) {
    d[, trend := prob-prob.prev,]
  }
  if ('prod.log' == indicator) {
    d[, freq.diff := ma-ma.prev,]
    d[, prob.rate := (prob-prob.prev)/prob.prev,]
    d[, trend := calculateProdLog(prob.rate, freq.diff),]
  }
  # always removing any non-finite trend value
  d <- d[is.finite(trend),]
  if (globalOutliers) {
    d[, surge := trend>=calculateThresholdTopOutliers(trend),]
  } else {
    d[, surge := trend>=calculateThresholdTopOutliers(trend), by=key(d)]
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


