library(GENIE3)
library(shinyjs)
library(minet)
library(dplyr)
library(DT)




ui.aracne <- function(id) {
  ar <- NS(id)

  fluidPage(
    useShinyjs(),
    titlePanel("Inferring Regulatory Networks by ARACNe Method"),
    sidebarLayout(
      sidebarPanel(
        useShinyjs(),
        # 左侧参数调整区域
        tags$style(
          type = "text/css",
          ".shiny-output-error { visibility: hidden; }",
          ".shiny-output-error:before { visibility: hidden; }",
          ".shiny-output-error:after { visibility: hidden; }"
        ),
        fileInput(ar("Exper"), "Upload exprMatr file(csv):"),
        fileInput(ar("geneList"), "Extract the Regulatory by geneList(one line one gene in csv file):(Option)"),
        actionButton(ar("run"), " Click to Run!!!"),
        downloadButton(ar("download"), "Download Full Regulatory Networks", disabled = TRUE),
        downloadButton(ar("download_list"), "Download the geneList Regulatory Networks ", disabled = TRUE)
      ),
      mainPanel(
        # 右侧图形展示区域
        wellPanel(
          DTOutput(ar("NOG"))
        )
      )
    )
  )
}
server.aracne <- function(input, output, session) {
  Network <- reactiveValues(
    data = NA,
    filter = NA
  )

  observeEvent(input$run, {
    tryCatch(
      {
        req(input$Exper)
        shinyjs::disable("run")
        data <- read.csv(input$Exper$datapath, header = T, sep = ",", row.names = 1)
        correlation_matrix <- cor(t(data), method = "pearson")
        network <- aracne(correlation_matrix)
        network_list <- as.data.frame(as.table(network))
        network_list <- na.omit(network_list[network_list$Freq != 0, ])
        colnames(network_list) <- c("node1", "node2", "Freq")

        print(head(network_list))

        if (!is.null(req(input$geneList))) {
          regulators <- read.csv(input$geneList$datapath, header = FALSE)[[1]]
          filter_Network <- subset(network_list, node1 %in% regulators | node2 %in% regulators)
          Network$filter <- filter_Network
        }


        shinyjs::enable("download")
        shinyjs::enable("run")
        # download_list
        if (!is.null(req(input$geneList))) {
          shinyjs::enable("download_list")
        }
      },
      error = function(e) {
        # 捕获并显示错误
        shinyjs::enable("run") # if ERROR has merge, we can run again!!
        shiny::showNotification(paste("Error:", e$message), type = "error")
      }
    )
  })


  output$download <- downloadHandler(
    filename = function() {
      paste("Full Regulatory Networks Construct by ARACNe", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      tryCatch(
        {
          # 写入 CSV 文件
          write.csv(unique(Network$data), file, col.names = NA, row.names = FALSE, quote = FALSE)
        },
        error = function(e) {
          # 捕获并显示错误
          shinyjs::enable("run") # 如果发生错误，启用按钮以便可以再次运行
          shiny::showNotification(paste("Error:", e$message), type = "error")
        }
      )
    }
  )

  output$download_list <- downloadHandler(
    filename = function() {
      paste("Select gene's  Regulatory Networks Construct by ARACNe", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      tryCatch(
        {
          # 写入 CSV 文件
          write.csv(unique(Network$filter), file, col.names = NA, row.names = FALSE, quote = FALSE)
        },
        error = function(e) {
          # 捕获并显示错误
          shinyjs::enable("run") # 如果发生错误，启用按钮以便可以再次运行
          shiny::showNotification(paste("Error:", e$message), type = "error")
        }
      )
    }
  )
}
