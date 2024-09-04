# Go Erichment

# 读入TERM2NAME
library(DT)
library(clusterProfiler)


ui.gsea <- function(id) {
  GA <- NS(id)

  fluidPage(
    titlePanel("GSEA Erichment Analysis"),
    sidebarLayout(
      sidebarPanel(
        # 左侧参数调整区域
        tags$style(
          type = "text/css",
          ".shiny-output-error { visibility: hidden; }",
          ".shiny-output-error:before { visibility: hidden; }",
          ".shiny-output-error:after { visibility: hidden; }"
        ),
        fileInput(GA("file1"), "Upload Gene set file:"),
        fileInput(GA("geneList"), "Upload gene list file:"),
        selectInput(GA("GeneSetType"), "Select GeneSet Type:",
          choices = c("none", "KeggSet", "GoSet"),
          selected = "none"
        ),
        selectInput(GA("cutoff"), "Pvalue cutoff:",
          choices = c(1, 0.05, 0.01),
          selected = 0.05
        ),
        sliderInput(GA("font_size"), "Label font size:", min = 1, max = 50, value = 14, step = 1),
        sliderInput(GA("strip"), "Strip tiltle font size:", min = 1, max = 50, value = 14, step = 1),
        actionButton(GA("Run"), "Erich!!"),
        downloadButton(GA("download_table"), "downTable", disable = TRUE),
        downloadButton(GA("download_plot"), "downPlot", disable = TRUE)
      ),
      mainPanel(
        # 右侧图形展示区域
        wellPanel(
          plotOutput(GA("plot"), height = "600px"),
          DTOutput(GA("table"), width = "600px")
        )
      )
    )
  )
}


server.gsea <- function(input, output, session) {
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
    req(input$GeneSetType != "none") # 检查是否选择了数据集

    Go_Erich$TERM2GENE <- read.table(input$file1$datapath, header = TRUE, sep = ",") # TERM2Gene
    if (input$GeneSetType == "GoSet") {
      Go_Erich$TERM2NAME <- read.csv("data/TERM2NAME.csv", header = T) # TERM2NAME
    } else if (input$GeneSetType == "KeggSet") {
      Go_Erich$TERM2NAME <- read.csv("data/PATHWAY2NAME.csv", header = T)
    }
    DEGs <- read.table(input$geneList$datapath, header = TRUE, sep = ",")

    new_DEGs <- DEGs[order(-DEGs$log2FoldChange), ] # 按着log2FoldChange的值排序
    gene_list <- new_DEGs$log2FoldChange # 提取logfold2的值
    names(gene_list) <- new_DEGs$GID # 给值命名，构建好genelist
    Go_Erich$geneList <- gene_list # 数据赋值
    print("文件全部读入！！")
  }) # 输入文件必须为csv
  ### 进行富集分析

  ### 成功获取了 TERN2NAME,TERM2GENE,GENELIST


  ## Erinchment,先把结果计算了储存起来，然后后面就可以调用
  observeEvent(input$Run, {
    tryCatch(
      {
        shinyjs::disable("Run")
        isolate({
          Go_Erich$result <- GSEA(
            Go_Erich$geneList,
            exponent = 1,
            minGSSize = 10,
            maxGSSize = 500,
            eps = 1e-10,
            pvalueCutoff = 0.05,
            pAdjustMethod = "BH",
            TERM2GENE = Go_Erich$TERM2GENE,
            TERM2NAME = Go_Erich$TERM2NAME,
            verbose = TRUE,
            seed = FALSE,
            by = "fgsea"
          )
        })
        # 画图
        Go_Erich$plot <- enrichplot::dotplot(Go_Erich$result,
          split = ".sign",
          showCategory = 20,
          color = "pvalue",
          label_format = 100,
          font.size = input$font_size
        ) + facet_wrap(~.sign, scales = "free") +
          theme(strip.text = element_text(size = input$strip))
        output$plot <- renderPlot(
          # 点图可视化
          enrichplot::dotplot(Go_Erich$result,
            split = ".sign",
            showCategory = 20,
            color = "pvalue",
            label_format = 100,
            font.size = input$font_size
          ) + facet_wrap(~.sign, scales = "free") +
            theme(strip.text = element_text(size = input$strip))
        )

        shinyjs::enable("Run")
        shinyjs::enable("download_table")
        shinyjs::enable("download_plot")

        # 下载表格
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
      paste("GSEA Result", Sys.time(), ".csv", sep = "")
    },
    content = function(file) {
      table <- as.data.frame(Go_Erich$result)

      write.csv(table, file, row.names = FALSE, quote = FALSE)
    }
  )



  # 下载图片
  output$download_plot <- downloadHandler(
    filename = function() {
      paste("GSEA dot plot", Sys.time(), ".pdf", sep = "")
    },
    content = function(file) {
      pdf(file, width = 11, height = 8.5)
      print(Go_Erich$plot)
      dev.off()
    }
  )
}
