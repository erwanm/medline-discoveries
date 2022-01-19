
#
# EM Jan 22 - shiny code for exploring the surges/discoveries extracted from Medline MeSH data
#       
#

library(shiny)
library(shinyWidgets)
source('../discoveries.R')

inputDir <- '../data/input'
staticInputDir <- paste(inputDir,'static',sep='/')
outputDir <- '../data/output'
UMLSSemanticGroupsFile <- '../data/SemGroups.txt'

sourcesList <- c('Full Medline' = '.min100','Medline ND subset' = '.min100.ND')
windowsList <- c(1,3,5)
indicatorsList <- c('diff','rate')
associationMeasureList <- c('proba'='prob.joint', 
                            'PMI'='pmi', 
                            'NPMI'='npmi',
                            'MI'='mi',
                            'NMI'='nmi',
                            'SCP'='scp')
defaultSource <- '.min100.ND'
defaultWindow <- 5
defaultIndicator <- 'diff'
defaultMeasure <- 'scp'

groupsSeparator <- ' '
defaultSemGroups <- c('DISO','CHEM','GENE','ANAT')

defaultCondiThreshold <- 0.7
defaultYearRange <- c(2000,2020)
yearMin<-1950
yearMax<-2020

viewOptionsList <- c('Show first surge only for each relation' = 'firstYear',
                     'Adjust surge year (next non-zero year)' = 'adjustYear',
                     'Sort by year (default: sort by tend)' = 'sortByYear',
                     'Show MeSH descriptors' = 'showDescriptors',
                     'Show conditional probabilities' = 'showConditional',
                     'Show semantic groups' = 'showGroups'
                     )
selected=defaultViewOptions <- c('firstYear','showDescriptors')


ui <- fluidPage(
  
  # Application title
  titlePanel("Medline impactful discoveries"),
  
  sidebarLayout(
    sidebarPanel(
#      helpText(tags$a(href=linkUserGuide,target="_blank", "User Guide")),
    
      # dataset selection
      selectInput("sourceSelectionInput", 
                  label = "Source dataset",
                  choices = sourcesList,
                  selected = defaultSource),
      fluidRow(
        column(7, radioButtons("windowInput", "Moving avg window",
                  windowsList, 
                  selected = defaultWindow, inline=TRUE)),
        column(5, radioButtons("indicatorInput", "Indicator",
                  indicatorsList,
                  selected = defaultIndicator, inline=TRUE))
      ),
        
      radioButtons("measureInput", "Measure",
                   associationMeasureList, 
                   selected = defaultMeasure, inline=TRUE),

      # semantic groups selection
      uiOutput("semanticTypesSelection"),
      # concept selection
      uiOutput("conceptSelection"),
      
      sliderInput("condiThresholdInput", "Conditional prob. maximum", min = 0, max = 1, value = defaultCondiThreshold),
      sliderInput("yearRangeInput", "Years Range", min = yearMin, max=yearMax, value= defaultYearRange),
      
      checkboxGroupInput("viewOptionsInput","View options", 
                         viewOptionsList,
                         selected=defaultViewOptions),
      
      downloadButton("downloadData", "Download")

    ),
    
    mainPanel(
      DT::dataTableOutput("discoveriesTable")
    )
  )
)


server <- function(input, output) {

  getStaticData <- reactive({
    loadStaticData(dir=staticInputDir,suffix=input$sourceSelectionInput)
  })
  
  getSurgesData <- reactive({
    loadSurgesData(dir=outputDir, suffix=input$sourceSelectionInput, ma_window=input$windowInput,measure=input$measureInput, indicator=input$indicatorInput)
  })

  getUMLSSemGroupsFile <- reactive({
    d<-fread(UMLSSemanticGroupsFile,sep='|',header=FALSE,drop = c(3,4))
    colnames(d) <- c('group.short','group.long')
    unique(d)
  })
  
  addNameAndGroup <- reactive({
    d <- getSurgesData()
    static_data <- getStaticData()
    if (!is.null(d)) {
      addRelationName(d,static_data)
    }
  })
  
  availableSemanticGroups <- reactive({
    d <- addNameAndGroup()
    if (!is.null(d)) {
      groups <- unlist(strsplit(c(d[,group.c1],d[,group.c2]),split=groupsSeparator,fixed=TRUE))
      data.table(group=sort(as.character(unique(groups))))
    }
  })
  

  output$semanticTypesSelection <- renderUI({
    g <- availableSemanticGroups()
    umls <- getUMLSSemGroupsFile()
    merged <- merge(g,umls,by.x='group',by.y='group.short')
    choicesAsList <- as.list(merged[,group])
    names(choicesAsList) <- merged[,group.long]
    pickerInput(
      inputId = "semTypesInput",
      label = "Filter by semantic group",
      choices = choicesAsList,
      selected = defaultSemGroups,
      options = list(
        `actions-box` = TRUE,
        size = 10,
        `selected-text-format` = "count > 4"
      ),
      multiple = TRUE
    )
  })

  applyGroupFilter <- reactive({
    d <- addNameAndGroup()
    groups <- input$semTypesInput
    if (!is.null(d)) {
      selectRelationsGroups(d,groups1=groups, groups2=groups)
    }
  })
  
  availableConcepts <- reactive({
    d <- applyGroupFilter()
    if (!is.null(d) & nrow(d)>0) {
      concepts <- unique(data.table(concept=c(d$c1,d$c2),term=c(d$term.c1,d$term.c2)))
      concepts <- concepts[order(term),]
      concepts <- rbind(data.table(concept='NONE', term='(no concept filter)'),concepts)
      concepts
    }
  })
  
  output$conceptSelection <- renderUI({
    concepts <-availableConcepts()
    if(!is.null(concepts)) {
        choices<-as.list(concepts$concept)
        names(choices)<-paste(concepts$concept, concepts$term,sep=' ')
        selectInput("conceptInput", 
                    label = "Filter  concept",
                    choices = choices,
                    selected = choices[[1]])
      }
  })
  
  
  applyConceptFilter <- reactive({
    d <- applyGroupFilter()
    selectedConcept <- input$conceptInput
    if (!is.null(d) & (nrow(d)>0) & !is.null(selectedConcept)) {
      if (selectedConcept == 'NONE') {
        d
      } else {
        filterConcepts(d,selectedConcept)
      }
    }
  })

  
  applyAllFilters <- reactive({
    d <- applyConceptFilter()
    if (!is.null(d) & !is.null(input$condiThresholdInput) & !is.null(input$viewOptionsInput) &!is.null(input$yearRangeInput) ) {
      filterCondProb(d, input$condiThresholdInput, 'firstYear' %in% input$viewOptionsInput, input$yearRangeInput[1], input$yearRangeInput[2])
    }
  })

  
  orderedByAssociation <- reactive({
    orderConceptsByAssociation(freqFiltered(), statsTargetAndView(), input$associationMeasure)
  })
  
  
  
  
  finalData <- reactive({
    df <- copy(applyAllFilters()) # copy because the year column is modified... not sure if needed
    if (!is.null(df)) {
      if ('adjustYear' %in% input$viewOptionsInput) {
        df[, year := next.nonzero.year]
      }
      if ('sortByYear' %in% input$viewOptionsInput) {
        df <- df[order(year),]
      } else {
        df <- df[order(-trend),]
      }
      cols <- 'year'
      if ('showDescriptors' %in% input$viewOptionsInput) {
        cols <- c(cols, c('c1','c2'))
      }
      cols <- c(cols, c('term.c1', 'term.c2'))
      if ('showGroups' %in% input$viewOptionsInput) {
        cols <- c(cols, c('group.c1','group.c2'))
      }
      if ('showConditional' %in% input$viewOptionsInput) {
        cols <- c(cols, c('prob.C1GivenC2','prob.C2GivenC1'))
      }
      cols <- c(cols, 'trend')
      df[,..cols]
    }
    
  })
  

  output$discoveriesTable <- DT::renderDataTable({
    d <- finalData()
    validate(need(d, 'Data not ready yet, please wait...'))
    d
    #    datatable(df) %>% formatPercentage(c("probCuiGivenTarget.high", "probTargetGivenCui.high","probCuiGivenTarget.low", "probTargetGivenCui.low"), 3)
    #    df[,c('CUI','firstTerm','probCuiGivenTarget.high','probTargetGivenCui.high','probCuiGivenTarget.low','probTargetGivenCui.low','MI')]
  },escape=FALSE,rownames= FALSE,options = list(  # escape=FALSE for displaying links as html
    # columnDefs = list(list(className = 'dt-center', targets = 5)),
    pageLength = 15
    #lengthMenu = c(5, 10, 15, 20)
  ))
  
  output$downloadData <- downloadHandler(
    filename = function() {
      paste(input$dataset, ".csv", sep = "")
    },
    content = function(file) {
      write.csv(datasetInput(), file, row.names = FALSE)
    }
  )
  
}

# Run the application 
shinyApp(ui = ui, server = server)

