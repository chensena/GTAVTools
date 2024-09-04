# mfuzz
library(Mfuzz) # biomanger


ui.mfuzz <- function(id) {
  mf <- NS(id)

  fluidPage(
    titlePanel("mfuzz analysis"),
    sidebarLayout(
      sidebarPanel(
        # 左侧参数调整区域
        tags$style(
          type = "text/css",
          ".shiny-output-error { visibility: hidden; }",
          ".shiny-output-error:before { visibility: hidden; }",
          ".shiny-output-error:after { visibility: hidden; }"
        ),
        fileInput(mf("file1"), "Upload data csv file:"),
        sliderInput(mf("Cluster_number"), "Set Cluster Number:", min = 1, max = 10, value = 2, step = 1),
        sliderInput(mf("row_graphic"), "Set graphic row Number:", min = 1, max = 10, value = 2, step = 1),
        sliderInput(mf("col_graphic"), "Set graphic col Number:", min = 1, max = 10, value = 2, step = 1),
        actionButton(mf("mfuzz_analysis"), " Run mfuzz!!!"),
        downloadButton(mf("download_result"), "Download result"),
        downloadButton(mf("download_plot"), "Download Plot")
      ),
      mainPanel(
        # 右侧图形展示区域
        wellPanel(
          plotOutput(mf("mfuzz_plot"), height = "600px")
        )
      )
    )
  )
}



server.mfuzz <- function(input, output, session) {
  # 定义一个响应用户选择的逻辑函数
  plot_data <- reactive({
    req(input$file1)

    dat <- read.csv(input$file1$datapath, header = TRUE, row.names = 1)
    dat <- as.matrix(dat)
    # creat object
    dat <- new("ExpressionSet", exprs = dat)

    # rm na
    dat <- filter.NA(dat, thres = 1)
    dat <- fill.NA(dat, mode = "mean")

    # filter low experssion gene
    dat <- filter.std(dat, min.std = 0.3)

    # Data normalization
    standardise(dat)
  })

  observeEvent(input$mfuzz_analysis, {
    # fuzzy c-means 聚类，需手动定义聚类个数，比方说设置 12 个簇
    dat <- plot_data()
    n <- input$Cluster_number
    print(n)
    # 评估出最佳的 m 值，防止随机数据聚类
    m <- mestimate(dat)
    # 聚类
    set.seed(123)
    cl <- mfuzz(dat, c = n, m = m)


    output$mfuzz_plot <- renderPlot(
      # fuzzy c-means 聚类，需手动定义聚类个数，比方说设置 12 个簇

      mfuzz.plot(dat, cl = cl, mfrow = c(input$row_graphic, input$col_graphic), new.window = F)
    )
  })

  # 下载绘图按钮的响应函数
  output$download_plot <- downloadHandler(
    filename = function() {
      paste("mfuzz_plot", Sys.Date(), ".pdf", sep = "")
    },
    content = function(file) {
      dat <- plot_data()

      n <- input$Cluster_number
      print(n)
      # 评估出最佳的 m 值，防止随机数据聚类
      m <- mestimate(dat)
      # 聚类
      set.seed(123)
      cl <- mfuzz(dat, c = n, m = m)

      # fuzzy c-means 聚类，需手动定义聚类个数，比方说设置 12 个簇


      pdf(file)

      # 生成 mfuzz.plot() 图像
      mfuzz.plot(dat, cl = cl, mfrow = c(input$row_graphic, input$col_graphic), new.window = FALSE)

      # 关闭 PDF 设备
      dev.off()
    }
  )


  ### 下载结果数据


  output$download_result <- downloadHandler(
    filename = function() {
      paste("mfuzz_result", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      dat <- plot_data()

      n <- input$Cluster_number
      print(n)
      # 评估出最佳的 m 值，防止随机数据聚类
      m <- mestimate(dat)
      # 聚类
      set.seed(123)
      cl <- mfuzz(dat, c = n, m = m)
      # 在这里生成要下载的数据
      gene_cluster <- cbind(cl$cluster, cl$membership)
      colnames(gene_cluster)[1] <- "cluster"
      write.csv(gene_cluster, file, col.names = NA, quote = FALSE)
    }
  )
}
