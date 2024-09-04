#Venn
library(VennDiagram)
library(venn)
library(shinyjs)

ui.Venn <- function(id) {
  vd <- NS(id)
  
  fluidPage(
    useShinyjs(),
    titlePanel("Venn Diagram"),
    
    sidebarLayout(
      sidebarPanel(
        style = "overflow-y: scroll; height: 600px;",
        # 左侧参数调整区域
        tags$style(type="text/css",
                   ".shiny-output-error { visibility: hidden; }",
                   ".shiny-output-error:before { visibility: hidden; }",
                   ".shiny-output-error:after { visibility: hidden; }"),
        
        fileInput(vd("file1"), "Upload data set withe csv Format:"),
        
        sliderInput(vd("Number"), "Number font size:", min = 1, max = 10, value = 5, step = 0.1),
        sliderInput(vd("Group"), "Group Name font size:", min = 1, max = 5, value = 1, step = 0.1),
        
        actionButton(vd("Run"), " Run VennDiagram!!!"),
        downloadButton(vd("download_result"), "Download result", disabled = TRUE),
        downloadButton(vd("download_plot"), "Download Plot", disabled = TRUE)
      ),
      
      mainPanel(
        # 右侧图形展示区域
        wellPanel(
          plotOutput(vd("Venn_plot"), height = "600px")
        )
      )
    )
  )
}



server.Venn <- function(input, output, session) {
  
  Venn_Result <- reactiveValues(
    plot = NULL,
    result = NULL,
    venn = NULL
    
  )

  #定义一个响应用户选择的逻辑函数
  plot_data <- reactive({
    req(input$file1)
    DEG <- read.delim(input$file1$datapath,sep=",")

    
  })
  
  observeEvent(input$Run, {
    req(input$file1)
    venn_dat <- plot_data()
    venn_list <- lapply(1:ncol(venn_dat), function(i) venn_dat[, i])
    names(venn_list) <- colnames(plot_data())
    
    #print(head(venn_list))
    #View(venn_list)
    
    output$Venn_plot <- renderPlot({
      
      Venn_Result$plot <- venn(venn_list,
                                     zcolor='style', # 调整颜色，style是默认颜色，bw是无颜色，当然也可以自定义颜???
                                     opacity = 0.3,  # 调整颜色透明???
                                     box = F,        # 是否添加边框
                                     ilcs= input$Number,     # 数字大小
                                     sncs = input$Group        # 组名字体大小
      )
      
    })
    
    shinyjs::enable("download_plot")
    shinyjs::enable("download_result")
    
  })
  
  output$download_plot <- downloadHandler(
    filename = function() {
      paste("Venn_plot", Sys.time(), ".pdf", sep = "")
    },
    content = function(file) {
      venn_dat <- plot_data()
      venn_list <- lapply(1:ncol(venn_dat), function(i) venn_dat[, i])
      names(venn_list) <- colnames(plot_data())
      
      pdf(file)  # A4纸大小 (单位: 英寸)
      venn(venn_list,
           zcolor='style', # 调整颜色，style是默认颜色，bw是无颜色，当然也可以自定义颜???
           opacity = 0.3,  # 调整颜色透明???
           box = F,        # 是否添加边框
           ilcs= 0.6,     # 数字大小
           sncs = input$Group        # 组名字体大小
      )
      dev.off()  # 关闭PDF设备0
    
    }
  )
  
  
  output$download_result <- downloadHandler(
    filename = function() {
      paste("Venn_Result", Sys.time(), ".csv", sep = "")
    },
    content = function(file) {
      venn_dat <- plot_data()
      venn_list <- lapply(1:ncol(venn_dat), function(i) venn_dat[, i])
      names(venn_list) <- colnames(plot_data())
      inter <- get.venn.partitions(venn_list)
      for (i in 1:nrow(inter)) inter[i,'values'] <- paste(inter[[i,'..values..']], collapse = '.')
      inter <- subset(inter, select = -..values.. )
      inter <- subset(inter, select = -..set.. )
      
      
      
      write.csv(unique(inter),file, col.names = NA, row.names = FALSE,quote = FALSE)
      
      
    }
  )
  
  
}
