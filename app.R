library(shiny)
# library(httr)
# library(jsonlite)
# library(xml2)
# library(org.Hs.eg.db)
# library(org.Mm.eg.db)
# library(org.Dm.eg.db)
# library(org.Dr.eg.db)
library(parallel)





# ui <- tagList(
#     tags$head(
#       tags$title('windowTitle'),
#       tags$link(rel="stylesheet", type="text/css",
#                 href="app.css"),
#       tags$h1(a(href="http://www.ohsu.edu/xd/health/services/heart-vascular/"))
#     )
#   )

# collects all of the tab UIs
#shinyUI(
#

ui <- tagList(
  tags$head(
    tags$style(HTML(" .shiny-output-error-validation {color: darkred; } ")),
    tags$style(".mybuttonclass{background-color:#CD0000;} .mybuttonclass{color: #fff;} .mybuttonclass{border-color: #9E0000;}")
  ),
  navbarPage(
    
    theme = "https://bootswatch.com/3/spacelab/bootstrap.min.css",
    inverse = TRUE,
    title = "Merge Counts Files",
    tabPanel("", 
             fluidRow(column(4,wellPanel(
               h4("Upload all files containing counts data"),
               h4("(select multiple .CSV)"),
               fileInput('datafile', '',
                         accept=c('text/csv', 
                                  'text/comma-separated-values,text/plain', 
                                  '.csv'),multiple = TRUE
               ),
               checkboxInput("addOne", "Add +1 to counts", FALSE),
               checkboxInput("addGeneNames", "Retrieve gene names from ensembl ids", FALSE),
               conditionalPanel("input.addGeneNames",
                                radioButtons('refGenome','',
                                             c('Homo_sapiens.GRCh38.81',
                                               'Homo_sapiens.GRCh38.84',
                                               'Mus_musculus.GRCm38.82',
                                               'Danio_rerio.GRCz10.84',
                                               'Drosophila_melanogaster.BDGP6.81'
                                             ),selected = "Homo_sapiens.GRCh38.81"),
                                radioButtons('geneNameColumn','',
                                             c('Add gene.names column after gene ids'="add",
                                               'Replace gene ids column by gene names'="replace"
                                             ),selected = "add")),
               conditionalPanel("output.filesUploaded",
                                actionButton("upload_data","Merge Files",
                                             style="color: #fff; background-color: #CD0000; border-color: #9E0000")),
               conditionalPanel("output.filesMerged",
                                downloadLink('downloadData', 'Download Merged File',class = "btn btn-primary", style="color: #fff; background-color: #9E0000; border-color: #9E0000"))
             )
             ),#column
             column(8,h2("Merged Table"),hr(),
                    dataTableOutput("contents")
             )#column
             )#fluidrow
    ),#tabpanel
    
    
    
    ## ==================================================================================== ##
    ## FOOTER
    ## ==================================================================================== ##              
    footer=p(hr(),p("ShinyApp created by ", strong("{Ayman Yousif}")," of ",align="center",width=4),
             p(("Center for Genomics and Systems Biology, NYU Abu Dhabi"),align="center",width=4),
             p(("Copyright (C) 2017, code licensed under GPLv3"),align="center",width=4)
    )
  ) #end navbarpage
) #end taglist




options(shiny.maxRequestSize = 60*1024^2)
# Define server logic required to draw a histogram ----
server <- function(input, output,session) {
  
  output$contents <- renderDataTable({
    tmp <- analyzeDataReactive()
    if(!is.null(tmp)) tmp$data
  })
  
  output$downloadData <- downloadHandler(
    filename = function() {
      paste('data-', Sys.Date(), '.csv', sep='')
    },
    content = function(con) {
      write.csv(analyzeDataReactive()$data, con,row.names=FALSE)
    }
  )
  
  inputDataReactive <- reactive({
    
    inFile <- input$datafile
    if (is.null(inFile))
      return(NULL)
    
    
    return(inFile)
  })
  
  
  analyzeDataReactive <- 
    eventReactive(input$upload_data,
                  ignoreNULL = FALSE, {
                    
                    inFile <- inputDataReactive()
                    if(is.null(inFile))
                      return(NULL)
                    
                    progress <- Progress$new(session, min=0, max=1)
                    on.exit(progress$close())
                    
                    progress$set(message = 'Merging files ...')
                    
                    files <- list();
                    
                    sep = '\t'
                    if(length(inFile$datapath) > 0 ){
                      testSep = read.csv(inFile$datapath[1], header = FALSE, sep = '\t')
                      if(ncol(testSep) < 2)
                        sep = ','
                    }
                    else
                      return(NULL)
                    
                    #remove zero size files
                    inFile <- inFile[inFile$size != 0,]
                    
                    fileContent = read.csv(inFile$datapath[1], header = FALSE, sep = sep)
                    fileContent = fileContent[!grepl("__", fileContent$V1),] #remove rows containing underscores
                    
                    #Sort by gene_id incase they are not sorted
                    fileContent = fileContent[order(fileContent[,1]),]
                    
                    #Create data frame table for merged counts
                    total = data.frame(matrix(ncol = length(inFile$datapath) + 1, nrow = nrow(fileContent)))
                    
                    #Extract sample name from file name
                    samplename = tools::file_path_sans_ext(inFile$name[1])
                    
                    #Create first column "gene.ids", and second column for first sample Counts
                    total[,1] = fileContent[,1]
                    total[,2] = fileContent[,2]
                    colnames(total)[1] = "gene.ids"
                    colnames(total)[2] = samplename
                    
                    
                    #Merge the remaining files and store in 'total'
                    for(i in 2:length(inFile$datapath))
                    {
                      
                      if(length(inFile$datapath) == 1)
                        break
                      
                      fileContent = read.csv(inFile$datapath[i], header = FALSE, sep = sep)
                      fileContent = fileContent[!grepl("__", fileContent$V1),] #remove rows containing underscores
                      fileContent = fileContent[order(fileContent[,1]),]
                      
                      samplename = tools::file_path_sans_ext(inFile$name[i])
                      
                      #add new column to total
                      total[,i+1] = fileContent[,2]
                      colnames(total)[i+1] = samplename
                      
                      progress$set(value = i/length(inFile$datapath))
                    }
                    
                    
                    
                    if(input$addGeneNames)
                    {
                      geneNames <- getNamesFromEnsembl(total$gene.ids, progress)
                      #browser()
                      if(input$geneNameColumn == "add")
                        total = as.data.frame(append(total, list(gene.names= geneNames), after = 1))
                      else{
                        #total[,1] = list(gene.names= geneNames)
                        
                        total[,1] = make.names(geneNames, unique=TRUE)
                        colnames(total)[1] = "gene.names"
                      }
                        
                    }
                    
                    if(input$addOne)
                      total[,!(names(total) %in% c("gene.ids","gene.names"))] = total[,!(names(total) %in% c("gene.ids","gene.names"))] + 1
                    
                    browser()
                    
                    return(list('data'=total))
                    
                  })
  
  getNamesFromEnsembl <- function(ensNames, progress)
  {
    # <- Progress$new(session, min=0, max=1)
    progress$set(value = 0.3)
    progress$set(message = 'Adding gene names ...')
    
                 
    #load("geneid2name.Rda")
    
    if(input$refGenome == "Homo_sapiens.GRCh38.81")
      load("Homo_sapiens.GRCh38.81.Rda")
    else if(input$refGenome == "Homo_sapiens.GRCh38.84")
      load("Homo_sapiens.GRCh38.84.Rda")
    else if(input$refGenome == "Mus_musculus.GRCm38.82")
      load("Mus_musculus.GRCm38.82.Rda")
    else if(input$refGenome == "Danio_rerio.GRCz10.84")
      load("Danio_rerio.GRCz10.84.Rda")
    else if(input$refGenome == "Drosophila_melanogaster.BDGP6.81")
      load("Drosophila_melanogaster.BDGP6.81.Rda")
    
    # geneStartStr = as.character(ensNames[1])
    # 
    # annoDb <- NULL
    # if(gdata::startsWith(geneStartStr, "ENSDAR",ignore.case=TRUE))
    #   annoDb= org.Dr.eg.db
    # else if(gdata::startsWith(geneStartStr, "ENSMUS",ignore.case=TRUE))
    #   annoDb <- org.Mm.eg.db
    # else if(gdata::startsWith(geneStartStr, "FB",ignore.case=TRUE))
    #   annoDb <- org.Dm.eg.db
    # else
    #   annoDb <- org.Hs.eg.db
    
    
    # Calculate the number of cores
    no_cores <- detectCores() - 1
    
    # Initiate cluster
    cl <- makeCluster(no_cores)
    
    print(paste(format(Sys.time(), "%H:%M:%OS3"),": Started Renaming ",length(ensNames), " genes"))
    #levelsList = character(length(ensNames))
    levelsList = parallel::parLapply(cl,ensNames, function(x){
      return(geneid2name[geneid2name$gene_id == as.character(x),]$gene_name)
    })

    #browser()
    print(paste(format(Sys.time(), "%H:%M:%OS3"),": Finished renaming"))
    stopCluster(cl)
    
    # for (i in 1:length(ensNames))
    # {
    #   #genename = geneid2name[geneid2name$gene_id == as.character(ensNames[i]),]$gene_name
    #   
    #   
    #   tryCatch({
    #     genename = select(annoDb, keys=as.character(ensNames[i]), keytype="ENSEMBL", column="SYMBOL")$SYMBOL
    #   }, error = function(e){
    #     print(as.character(ensNames[i]))
    #     genename = geneid2name[geneid2name$gene_id == as.character(ensNames[i]),]$gene_name
    #   })
    #   
    # 
    #   levelsList[i] = genename
    #   if(i %% 100)
    #     progress$set(value = i/length(ensNames))
    # }
    
    progress$set(value = 0.8)
    
    flatList = unlist(levelsList)
    #browser()
    progress$set(value = 1)
    return(flatList)
    
    #return(annoDb$SYMBOL)
  }
  
  
  output$filesUploaded <- reactive({
    return(!is.null(inputDataReactive()))
  })
  outputOptions(output, 'filesUploaded', suspendWhenHidden=FALSE)
  
  output$filesMerged <- reactive({
    return(!is.null(analyzeDataReactive()))
  })
  outputOptions(output, 'filesMerged', suspendWhenHidden=FALSE)
  
  
  observe({
    # Check if example selected, or if not then ask to upload a file.
    shiny:: validate(
      need((input$data_file_type=="examplecounts")|((!is.null(input$rdatafile))|(!is.null(input$datafile))), 
           message = "Please select a file")
    )
    inFile <- input$datafile
    
  })
  
}

shinyApp(ui = ui, server = server)