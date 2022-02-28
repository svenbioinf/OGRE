#' SHREC SHiny interface for REgion Comparison
#'
#' `SHREC()` is a graphical user interface for OGRE
#' @import shiny 
#' @importFrom DT renderDT datatable JS
#' @importFrom shinyFiles shinyDirButton shinyDirChoose
#' @import shinydashboard
#' @export
SHREC <- function(){
  runApp(shinyApp(
    ui = dashboardPage(
      dashboardHeader(title="OGRE Overlapping Genomic REgions",titleWidth = 500),
      dashboardSidebar(
        sidebarMenu(
          menuItem("Information", tabName = "information", icon = icon("exclamation")),
          menuItem("Preparations", tabName = "preparations", icon = icon("cogs")),
          menuItem("Charts", tabName = "charts", icon = icon("chart-bar")),
          menuItem("Tables", tabName = "tables", icon = icon("table")),
          menuItem("Plots", tabName = "plots", icon = icon("dna")),
          actionButton("runOGRE", "Start analysis",icon("play"),style=
              "color: #fff; background-color: #ff0e00; border-color: #ff0e00"),
          textAreaInput("datasets", "Datasets", rows = 4,value = "none"),
          br(),br(),
          h4("Status:"),
          textOutput("status1"),
          textOutput("status2"),
          HTML('<script type="text/javascript">$(document).ready(function() {
          $("#addHardDrive").click(function() {$("#status1").text("Loading...");
          });});</script>'),
          HTML('<script type="text/javascript">$(document).ready(function() {
          $("#addAnnotationHub").click(function() {$("#status2").text("Loading...");
          });});</script>')
        )
      ),
      dashboardBody(
        tags$head(tags$script(HTML('Shiny.addCustomMessageHandler("jsCode",
                                function(message) {eval(message.value);});'))),
        tabItems(
          tabItem(tabName = "information",
            fluidRow(
            box(width = 3,height = 200,title="Welcome!",column(12,align="center",
             img(src='extFolder/logo.png',align="center",width="249",hight="108"))),
            box(width = 3,height = 200,title =" Before you start...",
                  h4("Please make sure your datasets strand information uses format ('+' '-' '*')
                Not (1,-1,minus,plus), originate from a common genome version 
                and build. HG19 and HG38 coordinates can differ slightly."))),
            fluidRow(
            box(width = 3,height = 200,title="Welcome!"),
            box(width = 3,height = 200,title="Overlap!",column(12,align="center",
              img(src='extFolder/overlap.png',align="center",width="336",hight="140"))))
          ),
          tabItem(tabName = "preparations",
            box(title = "Datasets from hard drive",
            shinyFiles::shinyDirButton('queryFolder', 'Select query folder', 'Please select a query folder', FALSE),
            textInput("queryFolderText",label=NULL,
                      value=file.path(system.file('extdata', package = 'OGRE'),"query")),
            shinyFiles::shinyDirButton("subjFolder","Select subject folder","Please select a subject folder",FALSE),
            textInput("subjFolderText",label=NULL,
                      value=file.path(system.file('extdata', package = 'OGRE'),"subject")),
            actionButton("addHardDrive","Add datasets")
            ),
            box(title="Datasets from AnnotationHub",
            radioButtons("checkboxQuery", "Query",selected = character(0),
                             choices = listPredefinedDataSets()),                
            checkboxGroupInput("checkboxSubjects", label = ("Subjects"), 
                               choices = listPredefinedDataSets()),
            actionButton("addAnnotationHub","Add datasets")
            ),
            box(title = "Manipulate datasets",
                solidHeader=TRUE,splitLayout(
                  textAreaInput("subsetIdentifier", "Subset dataset by ID", rows = 3,value = "ID1\nID2\n..."),
                  textInput("subsetName","Name of dataset to subset:","myDataset"),
                ),
                splitLayout(title="Manipulate datasets",
                            textInput("extendIdentifier","Name of dataset to extend:","myDataset"),
                            textInput("extendUpstream","upstream(bp)",0),
                            textInput("extendDownstream","downstream(bp)",0),
                ),
                splitLayout(
                  textInput("promotersIdentifier","Name of dataset to extract promoters:","myDataset"),
                  textInput("promotersUpstream","upstream(bp)",0),
                  textInput("promotersDownstream","downstream(bp)",0),
                )),
            box(title = "GViz plot settings",
            radioButtons("queriesToPlot", h5("Queries to plot"),
                         choices = list("All queries" = 1, 
                                        "First 10 queries" = 2,
                                        "User defined"=3)
                         ,selected = 3,),
            textAreaInput("queriesToPlotCustom", "Status", rows = 4,
                          value="ENSG00000269011\nENSG00000142168",resize="none")
            )
          ),
          tabItem(tabName = "charts",
                  box(title = "Summary",
                  plotOutput(outputId = "barplot_summary")),
                  box(title = "Histogram",
                      plotOutput(outputId = "histogram"),
                      uiOutput("selectHist")
                  )
          ),
          tabItem(tabName = "tables",
                  box(title="Overlap numbers", DT::DTOutput("quickDT")),
                  tabBox(
                    title = "Overlap statistics",
                    id = "tabset1",
                    tabPanel("Summary",DT::DTOutput("summary")),
                    tabPanel("Summary(only overlapping)",DT::DTOutput("summary2")),
                  ),
                  box(title="Overlap details (wide-format)",width=12,DT::DTOutput("BestHits"),downloadButton("downloadb", "")),
                  hr(),
                  hr(),
                  box("Overlap details (long-format)",width=12,DT::DTOutput("summtD"),  downloadButton("download1", "")),# no label: this button will be hidden
          ),
          tabItem(tabName = "plots",
                  sidebarLayout(
                    sidebarPanel(
                      #DT::DTOutput("quickDT")
                    ),
                    mainPanel(
                      plotOutput("distPlot")
                    )
                  )
          ), 
          tabItem(tabName = "widgets",
                  h2("Widgets tab content")
          )
        )
      )#end dashboardbody
    ),#end dashboardpage
    server =  function(input, output,session) {#--------------------------------
      addResourcePath("extFolder", system.file('extdata', package = 'OGRE'))
      v = reactiveValues(myOGRE=OGREDataSet(),
                         queryFolder=NULL, 
                         subjFolder=NULL,
                         queriesToPlot=NULL,
                         status="Ready")
      
      shinyFiles::shinyDirChoose(input,id='queryFolder',roots = c(root = '/'))
      shinyFiles::shinyDirChoose(input,id='subjFolder',roots = c(root = '/'))
      observeEvent(input$queryFolder,{updateTextInput(session,"queryFolderText",
                    value=parseDirPath(roots = c(root = '/'),input$queryFolder))
        v$queryFolder <- parseDirPath(roots = c(root = '/'),input$queryFolder)
      })
      observeEvent(input$subjFolder,{updateTextInput(session,"subjFolderText",
                    value=parseDirPath(roots = c(root = '/'),input$subjFolder))
        v$subjFolder <- parseDirPath(roots = c(root = '/'),input$subjFolder)
      })
      observeEvent(input$queryFolderText,{
        v$queryFolder <- input$queryFolderText
      })
      observeEvent(input$subjFolderText,{
        v$subjFolder <- input$subjFolderText
      })
      observeEvent(input$queriesToPlot, {
        v$queriesToPlot <- input$queriesToPlot
      })
      observeEvent(input$addHardDrive,{#Add data from hardDrive
        if(v$queryFolder!=""){
          v$myOGRE <- addGRanges(v$myOGRE,readQuery(v$myOGRE)[[1]],"query")}
        if(v$subjectFolder!=""){
          v$myOGRE <- readSubject(v$myOGRE)}
        updateTextAreaInput(session,"datasets",value = paste0(names(v$myOGRE),"\n"))
        output$status1 <- renderText({"Dataset added. Ready"})
      })
      observeEvent(input$addAnnotationHub,{#Add data from AnnotationHub
        if(!is.null(input$checkboxQuery)){
          v$myOGRE <- addDataSetFromHub(v$myOGRE,input$checkboxQuery ,"query")
        }
        if(!is.null(input$checkboxSubjects)){
          for(x in input$checkboxSubjects){
            v$myOGRE <- addDataSetFromHub(v$myOGRE,x ,"subject")
        }}
        updateTextAreaInput(session,"datasets",value = paste0(names(v$myOGRE),"\n"))
        output$status2 <- renderText({"Dataset added. Ready"})
        
      })
      
      

      observeEvent(input$runOGRE, { #start main processing----------------------
        session$sendCustomMessage(type='jsCode', list(value = 'alert("Started calculation");'))
        
        #check if input is as file/annotationHub
        if(!is.null(v$queryFolder)&!is.null(v$queryFolder)){#if dir supplied
          #read folders
          v$myOGRE <- OGREDataSetFromDir(queryFolder = v$queryFolder,
                                         subjectFolder = v$subjFolder)
          v$myOGRE <- loadAnnotations(v$myOGRE)
          
          
        }
        else{
          
        }
        #check if user required subsetting
        if(input$subsetIdentifier!="ID1\nID2\n..."){
          subsetIdentifier <- strsplit(input$subsetIdentifier,split = "\n")[[1]]
          v$myOGRE <- subsetGRanges(v$myOGRE,subsetIdentifier,input$subsetName)
        }
        #check if user requires extending
        if(input$extendIdentifier!="myDataset"){
          v$myOGRE <- extendGRanges(v$myOGRE,input$extendIdentifier,
                                    input$extendUpstream,input$extendDownstream)
        }
        #check if user requires extending
        if(input$promotersIdentifier!="myDataset"){
          v$myOGRE <- extractPromoters(v$myOGRE,input$promotersIdentifier,
                                       input$promotersUpstream,input$promotersDownstream)
        }
        #number of queries to plot after processing
        getQueriesToPlot <- function(OGREDataSet,queriesToPlot){
          if(queriesToPlot==1){
            return(mcols(OGREDataSet[[1]])$ID)
          }else if(queriesToPlot==2){
            return(mcols(OGREDataSet[[1]])$ID[seq(10)])
          }else if(queriesToPlot==3){
            return(strsplit(input$queriesToPlotCustom,split = "\n")[[1]])
          }
        }
        v$myOGRE <- fOverlaps(v$myOGRE)
        v$myOGRE <- sumPlot(v$myOGRE)
        v$myOGRE <- plotHist(v$myOGRE)
        v$myOGRE <- summarizeOverlap(v$myOGRE)
        v$myOGRE <- gvizPlot(v$myOGRE,getQueriesToPlot(v$myOGRE,v$queriesToPlot),
                             showPlot = FALSE,
                             trackRegionLabels = setNames(c("name","name"),c("genes","CGI")))
        addResourcePath("gvizPlotsFolder", metadata(v$myOGRE)$gvizPlotsFolder)
        #link plots to QueryIDs
        metadata(v$myOGRE)$sumDT[,queryID:=paste0("<a target='_blank' href='",
                    "gvizPlotsFolder","/",queryID,".pdf","'>",queryID,"</a>")]
        session$sendCustomMessage(type='jsCode', list(value = 'alert("Finished calculation");'))
        ###Tables
        callback <- JS( #for custom download button
          "var a = document.createElement('a');",
          "$(a).addClass('dt-button');",
          "a.href = document.getElementById('download1').href;",
          "a.download = '';",
          "$(a).attr('target', '_blank');",
          "$(a).text('Download(Full table)');",
          "$('div.dwnld').append(a);",
          "$('#download1').hide();"
        )
        
        output$BestHits =DT::renderDT(server=TRUE,{
          datatable(metadata(v$myOGRE)$sumDT,extensions = 'Buttons',callback = callback2,
                    escape=FALSE,
                    options = list(
                      dom = 'B<"dwnldb">frtip',
                      autoWidth=FALSE,
                      scrollX = TRUE,
                      pageLength=5,
                      buttons = list(list(extend='csv',filename="BestHits"),list(extend='excel',filename="BestHits")))
          )})
        output$downloadb <- downloadHandler(
          filename = function() {
            paste("BestHitsFull.csv")
          },
          content = function(file) {
            fwrite(metadata(v$myOGRE)$sumDT, file)
          })
        callback2 <- JS( #for custom download button
          "var b = document.createElement('a');",
          "$(b).addClass('dt-button');",
          "b.href = document.getElementById('downloadb').href;",
          "b.download = '';",
          "$(b).attr('target', '_blank');",
          "$(b).text('Download(Full table)');",
          "$('div.dwnldb').append(b);",
          "$('#downloadb').hide();"
        )
        output$summtD = DT::renderDT(server=TRUE,{
          datatable(metadata(v$myOGRE)$detailDT,extensions = 'Buttons',selection=list(mode = 'multiple', selected = c(1)),callback = callback,
                    options = list(
                      dom = 'B<"dwnld">frtip',
                      autoWidth=FALSE,
                      scrollX = TRUE,
                      pageLength=5,
                      buttons = list(list(extend='csv',filename="DetailedTable"),list(extend='excel',filename="DetailedTable")))
          )})
        #download button
        output$download1 <- downloadHandler(
          filename = function() {
            paste("DetailedTableFull.csv")
          },
          content = function(file) {
            fwrite(metadata(v$myOGRE)$detailDT, file)
          }
        )
        ###Charts
        output$barplot_summary <- renderPlot({metadata(v$myOGRE)$barplot_summary})
        output$selectHist = renderUI({selectInput("plot", "Choose plot:", 
                    choices=names(metadata(v$myOGRE)$hist))})
        output$histogram <- renderPlot({metadata(v$myOGRE)$hist[[input$plot]]})
        ###Tables
        output$quickDT = DT::renderDT({
          metadata(v$myOGRE)$quickDT[,queryID:=paste0("<a target='_blank' href='",
            "gvizPlotsFolder","/",queryID,".pdf","'>",queryID,"</a>")]
          datatable(metadata(v$myOGRE)$quickDT,extensions = 'Buttons',escape = FALSE,
                    options = list(
                      dom = 'Bfrtip',
                      autoWidth=FALSE,
                      scrollX = TRUE,
                      scrollY="40vh",
                      paging = FALSE,
                      #columnDefs = list(list(width = '50', targets = c(1,2,3))),
                      buttons = c('copy', 'csv', 'excel', 'pdf', 'print')
                    ) #options = list(scrollY="300px",scrollX="300px", pageLength = 100, autoWidth = TRUE))
          )
        },server=FALSE)
        output$summary <- DT::renderDT({datatable(metadata(myOGRE)$summaryDT[["includes0"]])})
        output$summary2 <- DT::renderDT({datatable(metadata(myOGRE)$summaryDT[["excludes0"]])})
        
      })#end run
    }
  )
  )
}#end SHREC
