
library(ggplot2)
library(reshape2)
library(plyr)

loadTotalFiles <- function(dataPath='data/21-extract-discoveries', by='by-doc') {
  f <- paste(dataPath,'totals',by,'KD',sep='/')
  kd<-read.table(f,sep='\t')
  f <- paste(dataPath,'totals',by,'PTC',sep='/')
  ptc<-read.table(f,sep='\t')
  kd$source <- rep('KD',nrow(kd))
  ptc$source <- rep('PTC',nrow(ptc))
  d<-rbind(kd,ptc)
  colnames(d) <- c('year', 'unique_concepts','nb','concepts_mentions','source')
  d
}

loadDiscoveriesFile <- function(dataPath='data/21-extract-discoveries',f='full.min100/joint/KD.05.stats') {
  f <- paste(dataPath,f,sep='/')
#  print(f)
  d<-read.table(f,sep='\t',quote='',comment.char = '')
  if (ncol(d)==5) {
    colnames(d)<-c('concept','first_year','year','ratio','freq')
  } else {
    if (ncol(d) == 6) {
      colnames(d)<-c('c1','c2','first_year','year','ratio','freq')
    } else {
      stop(paste("Error: unexpected number of columns in ",f,":",ncol(d)))
    }
  }
  d
}

loadAllDiscoveriesFiles <- function(dataPath='data/21-extract-discoveries', dir='full.min100',types=c('indiv','joint'),sources=c('KD','PTC'),windowSizes=c('05','10'),suffix='.stats') {
  ldply(types, function(t) {
    ldply(sources, function(source) {
      ldply(windowSizes, function(windowSize) {
        f <- paste0(source,'.',windowSize,suffix)
        path <- paste(dir,t,f,sep='/')
        print(path)
        data <- loadDiscoveriesFile(dataPath,path)
        data$source <- rep(source,nrow(data))
        data$indivOrJoint <- rep(t,nrow(data))
        data$window_size <- rep(windowSize,nrow(data))
        data$window_size<-as.numeric(data$window_size)
        data$diff <- data$year- data$first_year
        data
      })
    })
  })
}

# d0 <- loadAllDiscoveriesFiles()
# d <- d0[sample(nrow(d0),100000),]
plotDistribRatio <- function(d,maxRatio=NA) {
  if (!is.na(maxRatio)) {
    d <- d[!is.na(d$ratio) & d$ratio<=maxRatio,]
  }
  ggplot(d,aes(ratio))+geom_histogram()+scale_y_log10()+facet_grid(window_size+indivOrJoint~source,scales="free_x")
}

# d <- loadAllDiscoveriesFiles()
statsByGroup <- function(d) {
  d0 <- d[,c('first_year','year','ratio','freq','source','indivOrJoint', 'window_size')]
  ddply(d0,c('source','indivOrJoint','window_size'), function(s) {
    nas = nrow(s[!complete.cases(s),])
    data.frame(rows=nrow(s), nas=nas, propNA=nas/nrow(s))
  })
}

# d0 <- loadAllDiscoveriesFiles()
# d <- d0[sample(nrow(d0),100000),]
plotDistribDiscoveryYears <- function(d, col='year', byGroup=TRUE) {
  if (byGroup) {
    ggplot(d,aes_string(col))+geom_histogram(bins = 34)+facet_grid(source+indivOrJoint~window_size,scales="free_y")
  } else {
    ggplot(d,aes_string(col))+geom_histogram(bins = 34)
  }
}

# 
# d <- loadAllDiscoveriesFiles()
plotDiffFirstYear <- function(d,col='diff',byGroup=TRUE,log=FALSE,colourCol=NULL) {
  if (is.null(colourCol)) {
    g <- ggplot(d,aes_string(col))
  } else {
    g <- ggplot(d,aes_string(col,fill=colourCol))
  }
  if (byGroup) {
    if ('indivOrJoint' %in% colnames(d)) {
      g <- g+geom_histogram()+facet_grid(source+indivOrJoint~window_size,scales="free_y")
    } else {
      g <- g+geom_histogram()+facet_grid(source~window_size,scales="free_y")
      
    }
  } else {
    g <- g+geom_histogram()
  }
  if (log) {
    g +scale_y_log10()
  } else {
    g
  }
}


# d <- loadAllDiscoveriesFiles()
proportionDiffLowerWindowByGroup <- function(d) {
  ddply(d,c('source','indivOrJoint','window_size'),function(s) {
    data.frame(prop=nrow(s[s$diff<s$window_size,])/nrow(s)) 
  })
}

# d <- loadAllDiscoveriesFiles()
addIndivInfoToJoint <- function(d) {
  ddply(d,c('source','window_size'),function(s) {
    print(s[1,c('source','window_size')])
    joint <- s[s$indivOrJoint=='joint',c('c1','c2','first_year','year','ratio','freq','diff')]
    indiv <- s[s$indivOrJoint=='indiv',c('concept','first_year','year','ratio','freq')]
    joint1 <- merge(joint,indiv,by.x='c1',by.y='concept',suffixes=c('','.c1'))
    joint2 <- merge(joint1,indiv,by.x='c2',by.y='concept',suffixes=c('','.c2'))
    # pmax = pairwise max, see https://stackoverflow.com/questions/19994543/how-can-i-take-pairwise-parallel-maximum-between-two-vectors
    joint2$latest_indiv_first_year <- pmax(joint2[,'first_year.c1'],joint2[,'first_year.c2'])
    joint2$diff_both_indiv <- joint2$year - joint2$latest_indiv_first_year
    joint2
  })
}

# d0 <- loadAllDiscoveriesFiles() 
# d <- addIndivInfoToJoint(d0)
proportionDiffLowerWindowByGroupJoint <- function(d,col='diff_both_indiv') {
  ddply(d,c('source','window_size'),function(s) {
    data.frame(prop=nrow(s[s[,col]<s$window_size,])/nrow(s)) 
  })
}


# d0 <- loadAllDiscoveriesFiles() 
# d <- addIndivInfoToJoint(d0)
addWindowLabel <- function(d) {
  d$firstYearsWithinWindow <- rep('no',nrow(d))
  d[is.na(d$diff) | is.na(d$diff_both_indiv),]$firstYearsWithinWindow <- rep('NA',nrow(d[is.na(d$diff) | is.na(d$diff_both_indiv),]))
  d[!is.na(d$diff) & d$diff<d$window_size,]$firstYearsWithinWindow <- rep('joint_first',nrow(d[!is.na(d$diff) &d$diff<d$window_size,]))
  d[!is.na(d$diff_both_indiv) & d$diff_both_indiv<d$window_size,]$firstYearsWithinWindow <- rep('indiv_first',nrow(d[!is.na(d$diff_both_indiv) & d$diff_both_indiv<d$window_size,]))
  d
}

addRatioBins <- function(d,bins=c(0,2,3,5,10,100,99999)) {
  d$ratio_bin <- cut(d$ratio,bins)
  d
}
