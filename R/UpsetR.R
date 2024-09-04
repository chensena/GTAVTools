#UpsetR
library(UpSetR)

ui.upset <- function(id) {
  us <- NS(id)
  
  fluidPage(
    useShinyjs(),
    titlePanel("UpsetR"),
    
    sidebarLayout(
      sidebarPanel(
        style = "overflow-y: scroll; height: 600px;",
        # 左侧参数调整区域
        tags$style(type="text/css",
                   ".shiny-output-error { visibility: hidden; }",
                   ".shiny-output-error:before { visibility: hidden; }",
                   ".shiny-output-error:after { visibility: hidden; }"),
        
        fileInput(us("file1"), "Upload data set withe csv Format:"),
        
        sliderInput(us("Point_size"), "Point Size", value = 4, min = 1, max = 10,sep=1),
        sliderInput(us("Line_size"), "Line Size", value = 1, min = 0.1, max = 5,sep=0.1),
        numericInput(us("Number"),"Show the set Number:",value = 30, min = 1, max = 100),
        sliderInput(us("inter_title"), "intersection size title", value = 1, min = 0.1, max = 5,sep=0.1),
        sliderInput(us("inter_tick"), "intersection size tick labels ", value = 1, min = 0.1, max = 5,sep=0.1),
        sliderInput(us("set_title"), "set size title", value = 1, min = 0.1, max = 5,sep=0.1),
        sliderInput(us("set_tick"), "set size tick labels", value = 1, min = 0.1, max = 5,sep=0.1),
        sliderInput(us("set_name"), "set names", value = 1, min = 0.1, max = 5,sep=0.1),
        sliderInput(us("above_bar"), "numbers above bars", value = 1, min = 0.1, max = 5,sep=0.1),
        
        actionButton(us("Run"), " Run VennDiagram!!!"),
        downloadButton(us("download_plot"), "Download Plot", disabled = TRUE)
      ),
      
      mainPanel(
        # 右侧图形展示区域
        wellPanel(
          plotOutput(us("Upset_plot"), height = "600px")
        )
      )
    )
  )
}



server.upset <- function(input, output, session) {
  
  upset_Result <- reactiveValues(
    plot = NULL
    
  )
  
  #定义一个响应用户选择的逻辑函数
  plot_data <- reactive({
    req(input$file1)
    List <- read.table(input$file1$datapath, header=TRUE, sep=",")
   
  })
  
  observeEvent(input$Run, {
    req(input$file1)
    

    output$Upset_plot <- renderPlot({
      
      upset_Result$plot <- upset(
        fromList(plot_data()),
        nsets=length(plot_data()),
        order.by = "freq",
        nintersects = input$Number,
        point.size = input$Point_size,
        line.size = input$Line_size,
        text.scale = c(input$inter_title,input$inter_tick,input$set_title,input$set_tick,input$set_name,input$above_bar)
        
      )
      upset_Result$plot
      })

    shinyjs::enable("download_plot")
    })
  

  
  output$download_plot <- downloadHandler(
    filename = function() {
      paste("Upset_plot", Sys.time(), ".pdf", sep = "")
    },
    content = function(file) {
      pdf(file,width = 11, height = 8.5)  # A4纸大小 (单位: 英寸)
      upset_Result$plot
      dev.off()  # 关闭PDF设备0
      
    }
  )
    
    
  }
  

  

