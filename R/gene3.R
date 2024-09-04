library(GENIE3)
library(shinyjs)

ui.genie3 <- function(id) {
  ge <- NS(id)
  
  fluidPage(
    useShinyjs(),
    titlePanel("Inferring Regulatory Networks by GENIE3"),
    
    sidebarLayout(
      
      sidebarPanel(
        # 左侧参数调整区域
        tags$style(type="text/css",
                   ".shiny-output-error { visibility: hidden; }",
                   ".shiny-output-error:before { visibility: hidden; }",
                   ".shiny-output-error:after { visibility: hidden; }"),
        
        fileInput(ge("Exper"), "Upload exprMatr file(csv):"),
        
        selectInput(ge("Candidat"), "Candidate Regulator",
                    choices = c("All the genes in exprMatr", "Selected genes in list"),
                    selected = "All the genes in exprMatr"),
        fileInput(ge("geneList"), "Upload Candidate Regulator file(one line one gene in csv file):"),
        selectInput(ge("method"), "tree-based method(Random Forests(RF) OR Extra-Trees(ET))",
                    choices = c("RF", "ET"),selected = "RF"),
        sliderInput(ge("Core"), "select Run Core Number.", value =1, min = 1, max = 20,sep=1),
       
        actionButton(ge("run"), " Run GENIE3!!!"),
        downloadButton(ge("download"), "Download Result",disabled = TRUE),
      ),

        mainPanel(
        # 右侧图形展示区域
        wellPanel(
          DTOutput(ge("NOG"))
        )
      )
    )
  )
}
server.genie3 <- function(input, output, session) {
  Network <- reactiveValues(
    data = NA
    
  )
  observe({
    tryCatch({
    if(input$Candidat == "All the genes in exprMatr"){
      shinyjs::disable("geneList")
    }
    else{
      shinyjs::enable("geneList")
    }
    }, error = function(e) {
      # 捕获并显示错误
      shiny::showNotification(paste("Error:", e$message), type = "error")
    })
  })
  
  observeEvent(input$run, {
    tryCatch({
      req(input$Exper)
      set.seed(123)
      shinyjs::disable("run")
      #input$Exper$datapath
      exprMatr <- as.matrix(read.table(input$Exper$datapath, sep = ',', header = TRUE, row.names = 1, quote = "", comment = "", check.names = FALSE))
      
      if(input$Candidat == "All the genes in exprMatr"){#从这里选择是创建全部的基因网络，还是指定候选基因的调控网络
        weightMat <- GENIE3(exprMatr,treeMethod=input$method,nCores=input$Core)
        
      }
      else{
        req(input$geneList)
        regulators <- read.csv(input$geneList$datapath, header = FALSE)[[1]]
        weightMat <- GENIE3(exprMatr,treeMethod=input$method,nCores=input$Core, regulators=regulators)
        
      }
      
      linkList <-   as.data.frame(getLinkList(weightMat))
      
      Network$data <- linkList
      ## after runing finsh, the download button is activate, user can download the result!!
      shinyjs::enable("download")
      shinyjs::enable("run") 
      
      
    }, error = function(e) {
      # 捕获并显示错误
      shinyjs::enable("run") # if ERROR has merge, we can run again!!
      shiny::showNotification(paste("Error:", e$message), type = "error")
    })
  })
  
  
  output$download <- downloadHandler(
    filename = function() {
      paste("Regulatory Networks by GENIE3", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      tryCatch({
        # 写入 CSV 文件
        write.csv(unique(Network$data), file, col.names = NA, row.names = FALSE, quote = FALSE)
      }, error = function(e) {
        # 捕获并显示错误
        shinyjs::enable("run") # 如果发生错误，启用按钮以便可以再次运行
        shiny::showNotification(paste("Error:", e$message), type = "error")
      })
    }
  )
}