#PCA分析

library(ggbiplot)

ui.pca <- function(id) {
  pca <- NS(id)
  
  fluidPage(
    useShinyjs(),
    titlePanel("Principal Component Analysis (PCA)"),
    
    sidebarLayout(
      sidebarPanel(
        
        useShinyjs(),
        style = "overflow-y: scroll; height: 600px;",
        # 左侧参数调整区域
        tags$style(type="text/css",
                   ".shiny-output-error { visibility: hidden; }",
                   ".shiny-output-error:before { visibility: hidden; }",
                   ".shiny-output-error:after { visibility: hidden; }"),
        
        fileInput(pca("Exper"), "Upload exprMatr file (csv):"),
        
        fileInput(pca("Group"), "Upload Group file (csv)"),
        
        checkboxInput(pca("label"), "Show Labels", value = FALSE),
        sliderInput(pca("label_size"), "Labels Size", value = 12, min = 6, max = 40,sep=1),
        sliderInput(pca("point_size"), "Point Size", value = 1.5, min = 0.5, max = 4,sep=0.1),
        
        
        checkboxInput(pca("ellipse"), "Add Ellipse", value = TRUE),
        checkboxInput(pca("ellipse_fill"), "Ellipse fill", value = TRUE),
        sliderInput(pca("ellipse_linewidth"), " Ellipse Linewidth", value = 1.3, min = 0.5, max = 4,sep=0.1),
        numericInput(pca("ellipse_alpha"), "Ellipse Alpha", value = 0.25, min = 0, max = 1),
        sliderInput(pca("axis_font_size"), "Axis text font size:", min = 1, max = 50, value = 14, step = 1),
        sliderInput(pca("title_font_size"), "Axis title font size:", min = 1, max = 50, value = 10, step = 1),
        numericInput(pca("axis_line"), 
                     "axis line size:", 
                     value = 1, 
                     min = 0.5, 
                     max = 5, 
                     step = 0.5),
        
        selectInput(pca("legend_position"), "Legend position:",
                    choices = c("none", "left", "right", "top", "bottom", "inset"),
                    selected = "right"),
        
        sliderInput(pca("legend_title_font_size"), "Legend title font size:", min = 1, max = 50, value = 10, step = 1),
        sliderInput(pca("legend_text_font_size"), "Legend text font size:", min = 1, max = 50, value = 10, step = 1),
       
        actionButton(pca("run"), " Click to Run PCA!!!"),
        actionButton(pca("draw"), " Click to Draw!!!",disabled = T),
        
        downloadButton(pca("download"), "Download",disabled = TRUE)
        
      ),
      
      mainPanel(
        # 右侧图形展示区域
        wellPanel(
          plotOutput(pca("NOG"), height = "600px")
        )
      )
    )
  )
}
server.pca <- function(input, output, session) {
  #定义一个响应用户选择的逻辑函数
  pca_plot <- reactiveValues(
    plot = NULL,
    pca =NULL,
    label =NULL,
    Group = NULL,
    Exper =NULL,
    Design = NULL
    
  )

  
  observeEvent(input$run, {
    shinyjs::disable("run")
    shinyjs::disable("draw")
    req(input$Exper)
    req(input$Group)
    
    Exper <- read.delim(input$Exper$datapath, row.names= 1,  header=T, sep=",")
    m.mad <- apply(Exper, 1, mad)
    Exper <- Exper[which(m.mad > max(quantile(m.mad, probs = seq(0, 1, 0.25))[2], 0.01)), ]
    Exper <- as.data.frame(Exper)
    
    Group <- read.table(input$Group$datapath, header=T, row.names= 1, sep=",") 
    colnames(Group) <- c("Group")
    
    pca_plot$pca <- prcomp(t(Exper), scale. = TRUE)
    
    pca_plot$Group <- Group$Group
      
    pca_plot$Design <- Group
    shinyjs::enable("run")
    shinyjs::enable("draw")
 
  })
  
  observe({
    if(input$label == T){
      
      pca_plot$label <- row.names(pca_plot$Design)
      
    }else
    {
      pca_plot$label <- NULL
      
    }
    
  })
  
    observeEvent(input$draw, {
      req(pca_plot$pca)
      pca_plot$plot <- ggbiplot(pca_plot$pca, 
                             obs.scale = 1, 
                             var.scale = 1,
                             groups = pca_plot$Group, 
                             ellipse = input$ellipse,
                             ellipse.fill = input$ellipse_fill,
                             ellipse.alpha = input$ellipse_alpha,
                             ellipse.linewidth = input$ellipse_linewidth,
                             var.axes = F,
                             labels =  pca_plot$label,
                             labels.size = input$label_size,
                             point.size = input$point_size
                             )+
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
    
    output$NOG <- renderPlot({
      
      pca_plot$plot
      
    })
  })
}
