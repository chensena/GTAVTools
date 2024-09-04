library(shiny)
library(ggplot2)

ui.volcano <- function(id) {
  ns <- NS(id)
  
  fluidPage(
    titlePanel("Volcano Plot Generator"),
    
    sidebarLayout(
      sidebarPanel(
        # 左侧参数调整区域
        tags$style(type="text/css",
                   ".shiny-output-error { visibility: hidden; }",
                   ".shiny-output-error:before { visibility: hidden; }",
                   ".shiny-output-error:after { visibility: hidden; }"),
        
        fileInput(ns("file1"), "Upload data file:"),
        
        sliderInput(ns("axis_font_size"), "Axis text font size:", min = 1, max = 50, value = 14, step = 1),
        sliderInput(ns("title_font_size"), "Axis title font size:", min = 1, max = 50, value = 10, step = 1),
        numericInput(ns("axis_line"), 
                     "axis line size:", 
                     value = 1, 
                     min = 0.5, 
                     max = 5, 
                     step = 0.5),
        
        selectInput(ns("legend_position"), "Legend position:",
                    choices = c("none", "left", "right", "top", "bottom", "inset"),
                    selected = "right"),
        
        sliderInput(ns("legend_title_font_size"), "Legend title font size:", min = 1, max = 50, value = 10, step = 1),
        sliderInput(ns("legend_text_font_size"), "Legend text font size:", min = 1, max = 50, value = 10, step = 1),
        actionButton(ns("draw_volcano"), "Draw Volcano"),
        
        downloadButton(ns("download_plot"), "Download Plot")
      ),
      
      mainPanel(
        # 右侧图形展示区域
        wellPanel(
          plotOutput(ns("volcano_plot"), height = "600px")
        )
      )
    )
  )
}



server.volcano <- function(input, output, session) {
  #定义一个响应用户选择的逻辑函数
  plot_data <- reactive({
    req(input$file1)
    cat(input$file1$datapath)
    DEG <- read.csv(input$file1$datapath, header = TRUE)
    
   
  })
  
  observeEvent(input$draw_volcano, {
    output$volcano_plot <- renderPlot({
      DEG <- plot_data()
      ggplot(data = DEG, aes(x = log2FoldChange, y = -log10(padj), color = change)) + 
        geom_point(shape = 16, size = 2) + 
        theme_classic() + 
        xlab("log2 Fold Change") + 
        ylab("-log10 Adjusted p-value") +
        #theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold")) + 
        scale_colour_manual(values = c('green', 'grey', 'red')) +
        theme(
          text = element_text(),
          axis.text = element_text(size = input$axis_font_size),
          axis.title = element_text(size = input$title_font_size),
          legend.position = input$legend_position,
          legend.text = element_text(size = input$legend_text_font_size),
          legend.title = element_text(size = input$legend_title_font_size),
          axis.line = element_line(size=input$axis_line),
          axis.ticks = element_line(size=input$axis_line),
          aspect.ratio = 1  # 设置长宽比为1:1
        )
    })
  })
  
  # 下载绘图按钮的响应函数
  output$download_plot <- downloadHandler(
    filename = function() {
      paste("volcano_plot", Sys.time(), ".pdf", sep = "")
    },
    content = function(file) {
      DEG <- plot_data()
      ggplot(data = DEG, aes(x = log2FoldChange, y = -log10(padj), color = change)) + 
        geom_point(shape = 16, size = 2) + 
        theme_classic() + 
        xlab("log2 Fold Change") + 
        ylab("-log10 Adjusted p-value") +
        #theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold")) + 
        scale_colour_manual(values = c('green', 'grey', 'red')) +
        theme(
          text = element_text(),
          axis.text = element_text(size = input$axis_font_size),
          axis.title = element_text(size = input$title_font_size),
          legend.position = input$legend_position,
          legend.text = element_text(size = input$legend_text_font_size),
          legend.title = element_text(size = input$legend_title_font_size),
          axis.line = element_line(size=input$axis_line),
          axis.ticks = element_line(size=input$axis_line),
          aspect.ratio = 1  # 设置长宽比为1:1
        )
      ggsave(file, width = 8, height = 6, units = "in",device = "pdf")
    }
  )
  
}
