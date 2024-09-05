# Correlation heatmap

# UpsetR
library(pheatmap)


ui.corheat <- function(id) {
  ch <- NS(id)

  fluidPage(
    useShinyjs(),
    titlePanel("Correlation Heatmap"),
    sidebarLayout(
      sidebarPanel(
        style = "overflow-y: scroll; height: 600px;",
        # 左侧参数调整区域
        tags$style(
          type = "text/css",
          ".shiny-output-error { visibility: hidden; }",
          ".shiny-output-error:before { visibility: hidden; }",
          ".shiny-output-error:after { visibility: hidden; }"
        ),
        fileInput(ch("file1"), "Upload data set withe csv Format:"),
        sliderInput(ch("fontsize_row"), "Row Font Size", value = 14, min = 6, max = 20, sep = 1),
        sliderInput(ch("fontsize_col"), "Column Font Size", value = 14, min = 6, max = 20, sep = 1),
        selectInput(ch("angle"), "Select angle of the column labels:",
          choices = c(45, 0, 90, 270, 315)
        ),
        actionButton(ch("Run"), " Run Correlation Heatmap!!!"),
        downloadButton(ch("download_plot"), "Download Plot", disabled = TRUE)
      ),
      mainPanel(
        # 右侧图形展示区域
        wellPanel(
          plotOutput(ch("heat_plot"), height = "700px")
        )
      )
    )
  )
}



server.corheat <- function(input, output, session) {
  upset_Result <- reactiveValues(
    plot = NULL,
    M = NULL
  )


  observeEvent(input$Run, {
    req(input$file1)
    dataExpr <- read.table(input$file1$datapath, header = TRUE, row.names = 1, sep = ",")
    m.mad <- apply(dataExpr, 1, mad)
    dataExprVar <- dataExpr[which(m.mad > max(quantile(m.mad, probs = seq(0, 1, 0.25))[2], 0.01)), ]
    dataExpr <- as.data.frame(dataExprVar)
    # View(dataExpr)
    M <- cor(log2(dataExpr + 1)) # 相关性分析
    upset_Result$M <- M

    output$heat_plot <- renderPlot({
      upset_Result$plot <- pheatmap(M,
        fontsize_row = input$fontsize_row,
        fontsize_col = input$fontsize_col,
        angle_col = input$angle
      )
      upset_Result$plot
    })

    shinyjs::enable("download_plot")
  })

  output$download_plot <- downloadHandler(
    filename = function() {
      paste("Correlation_plot", Sys.time(), ".pdf", sep = "")
    },
    content = function(file) {
      pdf(file, width = 11.69, height = 11.69) # A4纸大小 (单位: 英寸)
      pheatmap(upset_Result$M,
        fontsize_row = input$fontsize_row,
        fontsize_col = input$fontsize_col,
        angle_col = input$angle
      )
      dev.off() # 关闭PDF设备0
    }
  )
}
