# GSVA分析

# Go Erichment

# 读入TERM2NAME
library(DT)
library(clusterProfiler)
library(GSVA)


ui.gsva <- function(id) {
  GV <- NS(id)

  fluidPage(
    titlePanel("GSVA Erichment Analysis"),
    sidebarLayout(
      sidebarPanel(

        # 左侧参数调整区域
        tags$style(
          type = "text/css",
          ".shiny-output-error { visibility: hidden; }",
          ".shiny-output-error:before { visibility: hidden; }",
          ".shiny-output-error:after { visibility: hidden; }"
        ),
        fileInput(GV("file1"), "Upload Gene set file:"),
        fileInput(GV("geneExp"), "Upload Gene Experssion Matrix file:"),
        selectInput(GV("GeneSetType"), "Select GeneSet Type:",
          choices = c("none", "KeggSet", "GoSet"),
          selected = "none"
        ),
        actionButton(GV("checkFile"), "check file"),
        actionButton(GV("Run"), "Erich!!"),
        downloadButton(GV("download_table"), "downTable", disable = TRUE) # 设计一个下载GSVA表达矩阵的按钮
      ),
      mainPanel(
        # 右侧图形展示区域
        wellPanel(
          DTOutput(GV("InputTable"), height = "500px"),
          DTOutput(GV("OutputTable"), height = "500px")
        )
      )
    )
  )
}


server.gsva <- function(input, output, session) {
  Go_Erich <- reactiveValues(
    TERM2NAME = NULL, ## 储存TERM2NAME
    TERM2GENE = NULL, ## 储存TERM2GENE
    GeneExp = NULL, ## 储存
    result = NULL # 存储计算结果的变量
  )

  # 读取文件并处理数据的 reactive 表达式
  observeEvent(input$checkFile, {
    req(input$file1)
    req(input$geneExp)
    req(input$GeneSetType != "none")

    # 读取并处理 TERM2GENE 文件
    Go_Erich$TERM2GENE <- read.csv(input$file1$datapath, header = TRUE)
    print(head(Go_Erich$TERM2GENE))


    # 读取 TERM2NAME 或 PATHWAY2NAME 文件
    if (input$GeneSetType == "GoSet") {
      Go_Erich$TERM2NAME <- read.csv("data/TERM2NAME.csv", header = TRUE)
    } else if (input$GeneSetType == "KeggSet") {
      Go_Erich$TERM2NAME <- read.csv("data/PATHWAY2NAME.csv", header = TRUE)
    }

    # 读取和处理表达矩阵文件
    datExpr <- read.table(input$geneExp$datapath, header = TRUE, row.names = 1, sep = ",")
    datExpr <- log2(datExpr + 1)
    Go_Erich$GeneExp <- as.matrix(datExpr)

    print("文件全部读入！！")
  })

  # 进行富集分析的 observeEvent
  observeEvent(input$Run, {
    tryCatch(
      {
        req(Go_Erich$TERM2GENE)
        req(Go_Erich$TERM2NAME)
        req(input$GeneSetType != "none")

        if (input$GeneSetType == "GoSet") {
          print(head(Go_Erich$TERM2GENE))
          getset <- split(Go_Erich$TERM2GENE$GID, Go_Erich$TERM2GENE$GO)

          gsva.es <- gsva(Go_Erich$GeneExp, getset,
            method = "gsva", kcdf = "Gaussian", verbose = TRUE,
            parallel.sz = parallel::detectCores()
          )
          gsva.df <- as.data.frame(gsva.es)
          print(dim(gsva.df))
          gsva2 <- cbind(rownames(gsva.df), gsva.df)
          colnames(gsva2) <- c("GO", colnames(gsva.df))
          print(head(gsva.df)[, 1:3])
          print(head(Go_Erich$TERM2NAME))
          colnames(Go_Erich$TERM2NAME) <- c("GO", "NAME")
          Go_Erich$result <- merge(Go_Erich$TERM2NAME, gsva2, by = "GO")
        } else if (input$GeneSetType == "KeggSet") {
          getset <- split(Go_Erich$TERM2GENE$GID, Go_Erich$TERM2GENE$Pathway)

          print(head(getset))

          gsva.es <- gsva(Go_Erich$GeneExp, getset,
            method = "gsva", kcdf = "Gaussian", verbose = TRUE,
            parallel.sz = parallel::detectCores()
          )
          gsva.df <- as.data.frame(gsva.es)
          print(dim(gsva.df))
          gsva2 <- cbind(rownames(gsva.df), gsva.df)
          colnames(gsva2) <- c("Pathway", colnames(gsva.df))
          Go_Erich$result <- merge(Go_Erich$TERM2NAME, gsva2, by = "Pathway")
        }
      },
      error = function(e) {
        # 捕获并显示错误
        shiny::showNotification(paste("Error:", e$message), type = "error")
      }
    )
  })

  # 渲染结果表格
  output$OutputTable <- renderDT(
    {
      req(Go_Erich$result)
      Go_Erich$result[1:20, 1:3, drop = FALSE]
    },
    options = list(
      pageLength = 10,
      scrollX = TRUE,
      columnDefs = list(list(width = "auto", targets = "_all")),
      paging = TRUE
    )
  )

  # 渲染输入表格
  output$InputTable <- renderDT(
    {
      req(Go_Erich$GeneExp)
      Go_Erich$GeneExp[1:20, 1:3, drop = FALSE]
    },
    options = list(
      pageLength = 10,
      columnDefs = list(list(width = "auto", targets = "_all")),
      paging = TRUE
    )
  )

  output$download_table <- downloadHandler(
    filename = function() {
      paste("GSVA  Result", Sys.time(), ".csv", sep = "")
    },
    content = function(file) {
      table <- as.data.frame(Go_Erich$result)

      write.csv(table, file, row.names = FALSE, quote = FALSE)
    }
  )
}
