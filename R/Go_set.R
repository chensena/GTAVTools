# Go set create

library(tidyr)

library(dplyr)
library(AnnotationDbi) # biomanger
library(DT)

ui.Goset <- function(id) {
  GS <- NS(id)

  fluidPage(
    titlePanel("Generate Gene Set"),
    sidebarLayout(
      sidebarPanel(
        # 左侧参数调整区域
        tags$style(
          type = "text/css",
          ".shiny-output-error { visibility: hidden; }",
          ".shiny-output-error:before { visibility: hidden; }",
          ".shiny-output-error:after { visibility: hidden; }"
        ),
        fileInput(GS("file1"), "Upload data file:"),

        # actionButton(GS("Create_Go_Set"), "Generate"),

        downloadButton(GS("download_TERM2GENE"), "TERM2GENE"),
        downloadButton(GS("download_pathway2GENE"), "pathway2GENE")
      ),
      mainPanel(
        # 右侧图形展示区域
        wellPanel(
          DTOutput(GS("NOG"), height = "800px")
        )
      )
    )
  )
}



server.Goset <- function(input, output, session) {
  # 定义一个响应用户选择的逻辑函数
  plot_data <- reactive({ # 读取文件
    req(input$file1)

    emapper <- read.delim(input$file1$datapath, comment.char = "#")
  })

  output$NOG <- renderDT(
    {
      data <- plot_data() # 返回要显示的数据框
      # data <- subset(data,select=c(query,KEGG_ko,KEGG_rclass,KEGG_TC))[1:20,]
      data
    },
    options = list(
      pageLength = 5,
      # autoWidth = TRUE,  # 自适应列宽
      scrollX = TRUE,
      scrolly = TRUE,
      columnDefs = list(list(width = "auto", targets = "_all")), # 自适应列宽
      paging = TRUE # 启用分页
    )
  )

  ### 生成并下载TERM2GENE
  output$download_TERM2GENE <- downloadHandler(
    filename = function() {
      paste("TERM2GENE", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      emapper <- plot_data()
      emapper[emapper == "-"] <- NA
      emapper <- dplyr::select(emapper, GID = query, GO = GOs)
      TERM2GENE <- dplyr::select(emapper, GO, GID) %>%
        separate_rows(GO, sep = ",", convert = F) %>%
        dplyr::filter(GO != "-", !is.na(GO))

      write.csv(unique(TERM2GENE), file, col.names = NA, row.names = FALSE, quote = FALSE)
    }
  )

  ### 生成并下载pathway2GENE
  output$download_pathway2GENE <- downloadHandler(
    filename = function() {
      paste("pathway2GENE", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      emapper <- plot_data()
      emapper[emapper == "-"] <- NA
      emapper <- dplyr::select(emapper, GID = query, Pathway = KEGG_Pathway)
      pathway2GENE <- dplyr::select(emapper, Pathway, GID) %>%
        separate_rows(Pathway, sep = ",", convert = F) %>%
        dplyr::filter(!is.na(Pathway)) %>%
        dplyr::filter(!grepl("map", Pathway, ignore.case = TRUE))

      write.csv(pathway2GENE, file, row.names = FALSE, quote = FALSE)
    }
  )
}
