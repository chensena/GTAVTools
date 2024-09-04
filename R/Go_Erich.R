# Go Erichment

# 读入TERM2NAME
library(DT)
library(clusterProfiler)


ui.GoErich <- function(id) {
  GE <- NS(id)

  fluidPage(
    titlePanel("Go Erichment Analysis"),
    sidebarLayout(
      sidebarPanel(
        # 左侧参数调整区域
        tags$style(
          type = "text/css",
          ".shiny-output-error { visibility: hidden; }",
          ".shiny-output-error:before { visibility: hidden; }",
          ".shiny-output-error:after { visibility: hidden; }"
        ),
        fileInput(GE("file1"), "Upload Go set file:"),
        fileInput(GE("geneList"), "Upload gene list file:"),
        actionButton(GE("Run"), "Erich!!"),
        downloadButton(GE("download_table"), "downTable"),
        downloadButton(GE("download_plot"), "downPlot")
      ),
      mainPanel(
        # 右侧图形展示区域
        wellPanel(
          plotOutput(GE("plot"), height = "600px"),
          DTOutput(GE("table"), height = "500px")
        )
      )
    )
  )
}


server.GoErich <- function(input, output, session) {
  Go_Erich <- reactiveValues(
    TERM2NAME = NULL, ## 储存TERM2NAME
    TERM2GENE = NULL, ## 储存TERM2GENE
    geneList = NULL, ## 储存
    result = NULL # 存储计算结果的变量
  )

  # 定义一个响应用户选择的逻辑函数
  observe({ # 读取文件
    req(input$file1) # 检查两个文件都被读入
    req(input$geneList) # 检查文件是否读入
    Go_Erich$TERM2GENE <- read.table(input$file1$datapath, header = TRUE, sep = ",") # TERM2Gene
    Go_Erich$TERM2NAME <- read.csv("data/TERM2NAME.csv", header = T) # TERM2NAME
    Go_Erich$geneList <- read.table(input$geneList$datapath, header = TRUE, sep = ",")[, 1] # geneList
    print("文件全部读入！！")
  }) # 输入文件必须为csv
  ### 进行富集分析

  ### 成功获取了 TERN2NAME,TERM2GENE,GENELIST


  ## Erinchment,先把结果计算了储存起来，然后后面就可以调用
  observeEvent(input$Run, {
    Go_Erich$result <- enricher(
      Go_Erich$geneList,
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH",
      minGSSize = 10,
      maxGSSize = 500,
      qvalueCutoff = 0.2,
      TERM2GENE = Go_Erich$TERM2GENE,
      TERM2NAME = Go_Erich$TERM2NAME
    )

    # 画图
    output$plot <- renderPlot(
      # fuzzy c-means 聚类，需手动定义聚类个数，比方说设置 12 个簇
      dotplot(Go_Erich$result,
        showCategory = 20, color = "pvalue", label_format = 100,
        font.size = 10
      )
    )

    output$table <- renderDT(
      {
        table <- as.data.frame(Go_Erich$result)
        table[, 1:5]
      },
      options = list(
        pageLength = 10,
        # autoWidth = TRUE,  # 自适应列宽
        columnDefs = list(list(width = "auto", targets = "_all")), # 自适应列宽
        paging = TRUE # 启用分页
      )
    )




    ### 更新图表和图片
    ### 更新图表和图片
  })
}
