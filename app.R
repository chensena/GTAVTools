library(shiny)
library(shinydashboard)
library(shinythemes)
library(bslib)
options(shiny.maxRequestSize = 10 * 1024^3)

# 手动加载 R/ 目录下的所有 R 文件
sourceDir <- function(directory) {
  files <- list.files(directory, pattern = "\\.R$", full.names = TRUE)
  for (file in files) {
    source(file)
  }
}

# 调用 sourceDir() 加载 R/ 目录下的所有文件
sourceDir("R/")  # 假设 R/ 目录与 app.R 在同一目录下

# Define UI for application that draws a histogram
ui <- dashboardPage(
  skin = "green",
  dashboardHeader(title = "GTAVTools"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Home", tabName = "home", icon = icon("home")),
      
      menuItem("Data Analysis", icon = icon("table"),
               menuSubItem("Multi-group DEG Analysis (Pipeline)", tabName = "DEG"),
               menuSubItem("O2PLS Analysis (Pipeline)", tabName = "O2PLS"),
               menuSubItem("mfuzz", tabName = "mfuzz"),
               menuSubItem("Differential Analysis for GSVA", tabName = "DAGSVA")
               
               ),
      
      menuItem("Enrichment Analysis", icon = icon("chart-line"),#这些全都要注释文件，得开始做。
               menuSubItem("Go/KEGG Set Create", tabName = "Erich_Set"),
               menuSubItem("GO Erichment", tabName = "GO"),
               menuSubItem("KEGG Erichment", tabName = "KEGG"),
               menuSubItem("GSVA", tabName = "GSVA"),
               menuSubItem("GSEA", tabName = "GSEA"),
               menuSubItem("Promotor Element Statistics", tabName = "promotor")
      ),
      
      menuItem("Plot", icon = icon("chart-pie"),
               menuSubItem("Volcano Plot", tabName = "Volcano"),
               menuSubItem("Upset Plot", tabName = "Upset"),
               menuSubItem("Heatmap", tabName = "Heatmap"),
               menuSubItem("Venn Diagram", tabName = "Venn"),
               menuSubItem("Sample Correlation Heatmap", tabName = "Correlation"),
               menuSubItem("PCA", tabName = "PCA")
      ),
      
      menuItem("Machine Learning", icon = icon("cog"),
               
               menuSubItem("Lasso (Pipeline)", tabName = "Lasso"),
               
               menuSubItem("RandomForest (Pipeline)",tabName="RandomForest")
      ),
      
      menuItem("Regulatory Networks", icon = icon("cog"),
               menuSubItem("One Step WGCNA (Pipeline)", tabName = "WGCNA"),
               menuSubItem("GENIE3", tabName = "GENIE3"),
               menuSubItem("ARACNe", tabName = "ARACNe")
               
      ),
      menuItem("Little Tools", tabName = "tools", icon = icon("sign-out-alt"))
    )
  ),
  
  
  dashboardBody(
    tabItems(
      tabItem(tabName = "home", "Home tab content"),
      tabItem(tabName = "Volcano", ui.volcano("volcano")),  # 调用火山图模块的 UI 函数
      tabItem(tabName = "mfuzz", ui.mfuzz("mfuzz")),
      tabItem(tabName = "Erich_Set", ui.Goset("Erich_Set")),
      tabItem(tabName = "GO", ui.GoErich("Erich_Go")),
      tabItem(tabName = "KEGG", ui.KeggErich("Erich_KEGG")),
      tabItem(tabName = "GSEA", ui.gsea("GSEA")),
      tabItem(tabName = "GSVA", ui.gsva("gsva")),
      tabItem(tabName = "Heatmap", ui.heatmap("heatmap")),
      tabItem(tabName = "promotor", ui.promotor("promotor")),
      tabItem(tabName = "WGCNA", ui.wgcna("wgcna")),
      tabItem(tabName = "GENIE3", ui.genie3("genie3")),
      tabItem(tabName = "ARACNe", ui.aracne("aracne")),
      tabItem(tabName = "Lasso", ui.lasso("Lasso")),
      tabItem(tabName = "RandomForest", ui.randomforest("RandomForest")),
      tabItem(tabName = "DEG", ui.deg("DEG")),
      tabItem(tabName = "O2PLS", ui.o2pls("O2PLS")),
      tabItem(tabName = "PCA", ui.pca("PCA")),
      tabItem(tabName = "Venn", ui.Venn("Venn")),
      tabItem(tabName = "Upset", ui.upset("upset")),
      tabItem(tabName = "Correlation", ui.corheat("Correlation"))
    )
  )
  
)

# Define server logic
server <- function(input, output, session) {
  callModule(server.volcano, "volcano")  # 调用火山图模块的服务器逻辑
  callModule(server.mfuzz, "mfuzz")  #调用mfuzz分析逻辑
  callModule(server.Goset, "Erich_Set")  #调用Go set 分析逻辑
  callModule(server.GoErich, "Erich_Go")  #Go Erichment
  callModule(server.KeggErich, "Erich_KEGG")
  callModule(server.gsea, "GSEA")
  callModule(server.gsva, "gsva")
  callModule(server.heatmap, "heatmap")
  callModule(server.promotor, "promotor")
  callModule(server.wgcna, "wgcna")
  callModule(server.genie3, "genie3")
  callModule(server.aracne, "aracne")
  callModule(server.lasso, "Lasso")
  callModule(server.randomforest, "RandomForest")
  callModule(server.deg, "DEG")
  callModule(server.o2pls, "O2PLS")
  callModule(server.pca, "PCA")
  callModule(server.Venn, "Venn")
  callModule(server.upset, "upset")
  callModule(server.corheat, "Correlation")
}


# Run the application
shinyApp(ui = ui, server = server)
