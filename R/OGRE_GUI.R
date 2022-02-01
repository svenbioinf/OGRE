#' SHREC SHiny interface for REgion Comparison
#'
#' `SHREC()` is a graphical user interface for OGRE
#' @import shiny 
#' @importFrom DT renderDT datatable JS
#' @importFrom shinyFiles shinyDirButton shinyDirChoose
#' @examples
#' SHREC()
#' @export
SHREC <- function(){

runApp(shinyApp(
  ui = shinyUI(
    navbarPage("OGRE",tags$head(tags$script(HTML('Shiny.addCustomMessageHandler("jsCode",function(message) {eval(message.value);});'))),
      tabPanel("Notes",

      ),
      tabPanel("1) Preparations",
               h3("Input from hard drive"),
               hr(),
               shinyFiles::shinyDirButton('queryFolder', 'Select query folder', 'Please select a query folder', FALSE),
               textInput("queryFolderText",label=NULL,
                value=file.path(system.file('extdata', package = 'OGRE'),"query")),
               shinyFiles::shinyDirButton("subjFolder","Select subject folder","Please select a subject folder",FALSE),
               textInput("subjFolderText",label=NULL,
                value=file.path(system.file('extdata', package = 'OGRE'),"subject")),
               h3("GViz plot settings"),
               hr(),
               radioButtons("queriesToPlot", h5("Queries to plot"),
                            choices = list("All queries" = 1, 
                                           "First 10 queries" = 2,
                                           "User defined"=3)
                            ,selected = 3,),
               textAreaInput("queriesToPlotCustom", "Status", rows = 4,
                 value="ENSG00000269011\nENSG00000142168",resize="none")
               
      ),
      tabPanel("3) Run",
               actionButton("runOGRE", "Run OGRE")
      ),
      tabPanel("Graphs",
       plotOutput(outputId = "barplot_summary")
               
      ),
      navbarMenu("Tables",
       tabPanel("Hits(wide-format)",DT::DTOutput("BestHits"),downloadButton("downloadb", "")),
       tabPanel("Hits(long-format)",DT::DTOutput("summtD"),  downloadButton("download1", "")),# no label: this button will be hidden
               
      ),
      tabPanel("Plots",
      )
    )
  ),
  server =  function(input, output,session) {
    v = reactiveValues(myOGRE=OGREDataSet(),
                       queryFolder=NULL, 
                       subjFolder=NULL,
                       queriesToPlot=NULL)
    
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
    observeEvent(input$runOGRE, {
      session$sendCustomMessage(type='jsCode', list(value = 'alert("Started calculation");'))
      getQueriesToPlot <- function(OGREDataSet,queriesToPlot){
        if(queriesToPlot==1){
          return(mcols(OGREDataSet[[1]])$ID)
        }else if(queriesToPlot==2){
          return(mcols(OGREDataSet[[1]])$ID[seq(10)])
        }else if(queriesToPlot==3){
          return(strsplit(input$queriesToPlotCustom,split = "\n")[[1]])
        }
      }
      if(!is.null(v$queryFolder)&!is.null(v$queryFolder)){#if dir supplied
        #read folders
        v$myOGRE <- OGREDataSetFromDir(queryFolder = v$queryFolder,
                                       subjectFolder = v$subjFolder)
        v$myOGRE <- loadAnnotations(v$myOGRE)
        
        v$myOGRE <- fOverlaps(v$myOGRE)
        v$myOGRE <- sumPlot(v$myOGRE)
        v$myOGRE <- gvizPlot(v$myOGRE,getQueriesToPlot(v$myOGRE,v$queriesToPlot),
                     showPlot = FALSE,
                     trackRegionLabels = setNames(c("name","name"),c("genes","CGI")))
        
      }
      else{
        
      }
      #link plots to QueryIDs
      metadata(v$myOGRE)$sumDT[,queryID:=paste0("<a href='","file:///",
       metadata(v$myOGRE)$gvizPlotsFolder,"/",queryID,".pdf","'>",
       queryID,"</a>")]
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
    ###Plots
    output$barplot_summary <- renderPlot({metadata(v$myOGRE)$barplot_summary})
    })#end run
  }
)
)
}#end SHREC
