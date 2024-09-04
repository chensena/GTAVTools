# KEGG Erichment
# Go Erichment

# 读入TERM2NAME
library(DT)
library(clusterProfiler)


ui.KeggErich <- function(id) {
  KE <- NS(id)

  fluidPage(
    titlePanel("kEGG Pathway Erichment Analysis"),
    sidebarLayout(
      sidebarPanel(
        # 左侧参数调整区域
        tags$style(
          type = "text/css",
          ".shiny-output-error { visibility: hidden; }",
          ".shiny-output-error:before { visibility: hidden; }",
          ".shiny-output-error:after { visibility: hidden; }"
        ),
        fileInput(KE("file1"), "Upload Go set file:"),
        fileInput(KE("geneList"), "Upload gene list file:"),
        selectInput(KE("cutoff"), "Pvalue cutoff:",
          choices = c(1, 0.05, 0.01),
          selected = 0.05
        ),
        sliderInput(KE("font_size"), "Label font size:", min = 1, max = 50, value = 14, step = 1),
        actionButton(KE("Run"), "Erich!!"),
        downloadButton(KE("download_table"), "downTable", disable = TRUE),
        downloadButton(KE("download_plot"), "downPlot", disable = TRUE)
      ),
      mainPanel(
        # 右侧图形展示区域
        wellPanel(
          plotOutput(KE("plot"), height = "600px")
        )
      )
    )
  )
}


server.KeggErich <- function(input, output, session) {
  Go_Erich <- reactiveValues(
    TERM2NAME = NULL, ## 储存TERM2NAME
    TERM2GENE = NULL, ## 储存TERM2GENE
    geneList = NULL, ## 储存
    result = NULL, # 存储计算结果的变量
    plot = NULL
  )

  # 定义一个响应用户选择的逻辑函数
  observe({ # 读取文件
    req(input$file1) # 检查两个文件都被读入
    req(input$geneList) # 检查文件是否读入
    Go_Erich$TERM2GENE <- read.table(input$file1$datapath, header = TRUE, sep = ",") # TERM2Gene
    Go_Erich$TERM2NAME <- read.csv("data/PATHWAY2NAME.csv", header = T) # TERM2NAME
    Go_Erich$geneList <- read.table(input$geneList$datapath, header = TRUE, sep = ",")[, 1] # geneList
    print("文件全部读入！！")
  }) # 输入文件必须为csv
  ### 进行富集分析

  ### 成功获取了 TERN2NAME,TERM2GENE,GENELIST


  ## Erinchment,先把结果计算了储存起来，然后后面就可以调用
  observeEvent(input$Run, {
    tryCatch(
      {
        shinyjs::disable("Run")
        Go_Erich$result <- enricher(
          Go_Erich$geneList,
          pvalueCutoff = as.numeric(input$cutoff),
          pAdjustMethod = "BH",
          minGSSize = 10,
          maxGSSize = 500,
          qvalueCutoff = 0.2,
          TERM2GENE = Go_Erich$TERM2GENE,
          TERM2NAME = Go_Erich$TERM2NAME
        )

        # 画图
        Go_Erich$plot <- barplot(Go_Erich$result,
          showCategory = 20, color = "pvalue", label_format = 100,
          font.size = input$font_size
        )
        output$plot <- renderPlot(
          # fuzzy c-means 聚类，需手动定义聚类个数，比方说设置 12 个簇
          barplot(Go_Erich$result,
            showCategory = 20, color = "pvalue", label_format = 100,
            font.size = input$font_size
          )
        )



        ### 更新图表和图片
        ### 更新图表和图片
        shinyjs::enable("Run")
        shinyjs::enable("download_table")
        shinyjs::enable("download_plot")
      },
      error = function(e) {
        # 捕获并显示错误
        shinyjs::enable("Run") # if ERROR has merge, we can run again!!
        shiny::showNotification(paste("Error:", e$message), type = "error")
      }
    )
  })

  output$download_table <- downloadHandler(
    filename = function() {
      paste("KEGG Erichment Result", Sys.time(), ".csv", sep = "")
    },
    content = function(file) {
      table <- as.data.frame(Go_Erich$result)

      write.csv(table, file, row.names = FALSE, quote = FALSE)
    }
  )



  # 下载图片
  output$download_plot <- downloadHandler(
    filename = function() {
      paste("KEGG Erichment plot", Sys.time(), ".pdf", sep = "")
    },
    content = function(file) {
      pdf(file, width = 11, height = 8.5)
      print(Go_Erich$plot)
      dev.off()
    }
  )
}
