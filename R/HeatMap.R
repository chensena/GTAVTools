#Heatmap
library(pheatmap)


ui.heatmap <- function(id) {
  hp <- NS(id)
  
  fluidPage(
    titlePanel("Heatmap"),
    
    sidebarLayout(
      sidebarPanel(
        style = "overflow-y: scroll; height: 600px;",
        # 左侧参数调整区域
        tags$style(type="text/css",
                   ".shiny-output-error { visibility: hidden; }",
                   ".shiny-output-error:before { visibility: hidden; }",
                   ".shiny-output-error:after { visibility: hidden; }"),
        fileInput(hp("expr_file"), "Choose Expression Data File (CSV)", accept = ".csv"),
        # fileInput(hp("gene_file"), "Choose Gene Annotation File (CSV)(Option)", accept = ".csv"),##col_Annotation
        #fileInput(hp("sample_file"), "Choose Sample Annotation File (CSV)(Option)", accept = ".csv"),##row_Annotation
        checkboxInput(hp("log2_transform"), "Log Transform Data", value = FALSE),
        #checkboxInput(hp("auto_gaps_row"), "Auto Generate Row Gaps", value = FALSE),
        #checkboxInput(hp("auto_gaps_col"), "Auto Generate Column Gaps", value = FALSE),
        sliderInput(hp("fontsize_row"), "Row Font Size", value = 12, min = 6, max = 20,sep=1),
        sliderInput(hp("fontsize_col"), "Column Font Size", value = 14, min = 6, max = 20,sep=1),
        selectInput(hp("clustering_method"), "Clustering Method", 
                    choices = c("ward.D","ward.D2","single","complete","average","mcquitty","median","centroid"),
                    selected = "complete"
        ),
        
        
        
        fluidRow(
          column(3, # 每列占据的宽度
                 checkboxInput(hp("cluster_rows"), "Cluster Rows", value = FALSE)
          ),
          column(3,
                 checkboxInput(hp("cluster_cols"), "Cluster Cols", value = FALSE)
          )
        ),
        fluidRow(
          column(3,
                 checkboxInput(hp("show_rownames"), "Show Row Names", value = FALSE)
          ),
          column(3,
                 checkboxInput(hp("show_colnames"), "Show Col Names", value = FALSE)
          )
        ),
        fluidRow(
          column(3,
                 checkboxInput(hp("SetCellheight"), "Set Cellheight", value = FALSE)
          ),
          column(3,
                 checkboxInput(hp("SetCellwidth"), "Set Cellwidth", value = FALSE)
          )
        ),
        radioButtons(hp("Scale"), 
                     label = "Scale method",
                     choices = c("none", "row", "column"),
                     selected = "none"),
        numericInput(hp("cellwidth"), "cellwidth", value = 14, min = 1, max = 100),
        numericInput(hp("cellheight"), "cellheight", value = 14, min = 1, max = 100), 
        #checkboxInput(hp("cluster_rows"), "Cluster Rows", value = FALSE),##row cluater
        #checkboxInput(hp("cluster_cols"), "Cluster Cols", value = FALSE),#col cluster 
        #checkboxInput(hp("show_rownames"), "Show Row Names", value = FALSE),#
        #checkboxInput(hp("show_colnames"), "Show Col Names", value = FALSE),#
        #checkboxInput(hp("show_legend"), "Annotation Legend", value = FALSE),
        actionButton(hp("update"), "Draw Heatmap"),
        downloadButton(hp("download_plot"), "downPlot"),
        
        
        numericInput(hp("width"),"OutPut Plot Width",value = 6, min = 1, max = 10),
        numericInput(hp("height"),"OutPut Plot Height",value = 8, min = 1, max = 10),
        
        
      ),
      
      mainPanel(
        # 右侧图形展示区域
        wellPanel(
          
          plotOutput(hp("plot"), height = "800px",width="800px")
          
        )
      )
    )
  )
}


server.heatmap <- function(input, output, session) {
  
  
  PlotData <- reactiveValues(
    expr_data = NULL,
    col_annotation = NULL,
    row_annotation = NULL,
    gap_col = NULL,
    gap_row = NULL,
    cellwidth = NA,
    cellheight = NA
  )
  
  
  
  observe({#Gene Group Annotation
    req(input$gene_file)
    gene_data <- read.csv(input$gene_file$datapath, header = TRUE, stringsAsFactors = TRUE)
    
    # 确保文件中包含正确的列
    if (!all(c("Gene", "Group") %in% names(gene_data))) {
      stop("CSV Must Contain Cols 'Gene' and 'Group' ")
    }
    
    # 按照样本的顺序提取分组信息
    group_info <- gene_data$Group
    #print(group_info)
    # 计算每个分组的样本数
    gene_counts <- table(group_info)
    
    # 生成 annotation_col 数据框，保持原始分组顺序
    annotation_row <- data.frame(
      Group = factor(group_info, levels = names(gene_counts))
    )
    print(annotation_row)
    if (input$auto_gaps_row == TRUE){PlotData$gap_row <- as.vector(gene_counts)}else{PlotData$gap_row <- NULL}
    
    PlotData$row_annotation <-  annotation_row
    
    
  })
  
  observe({#Sample Group Annotation
    req(input$sample_file)
    sample_data <- read.csv(input$sample_file$datapath, header = TRUE, stringsAsFactors = TRUE)
    
    # 确保文件中包含正确的列
    if (!all(c("Sample", "Group") %in% names(sample_data))) {
      stop("CSV Must Contain Cols 'Sample' and 'Group' ")
    }
    
    # 按照样本的顺序提取分组信息
    group_info <- sample_data$Group
    #print(group_info)
    # 计算每个分组的样本数
    sample_counts <- table(group_info)
    print("sample_count")
    print(sample_counts)
    
    # 生成 annotation_col 数据框，保持原始分组顺序
    annotation_col <- data.frame(
      Group = factor(group_info, levels = names(sample_counts))
    )
    print(annotation_col)
    if (input$auto_gaps_col == TRUE){PlotData$gap_col <- as.vector(sample_counts)}else{PlotData$gap_col<- NULL}
    PlotData$col_annotation <-  annotation_col
  })
  
  observe({#experssion data
    if(input$SetCellwidth == TRUE){
      
      PlotData$cellwidth <- input$cellwidth
      
    }else{
      PlotData$cellwidth <- NA
      
    }
    
    if(input$SetCellheight == TRUE){
      
      PlotData$cellheight <- input$cellheight
      
    }else{
      PlotData$cellheight <- NA
      
    }    
    
    
    
  })
  
  observeEvent(input$update, {
    req(input$expr_file)
    if (input$log2_transform == FALSE){
      
      PlotData$expr_data <- read.csv(input$expr_file$datapath, header = TRUE, row.names = 1)
      
    }else{
      
      PlotData$expr_data <-  log2(read.csv(input$expr_file$datapath, header = TRUE, row.names = 1) + 1)
      
    }
    output$plot <- renderPlot({
      
      pheatmap(
        PlotData$expr_data,
        #annotation_row = PlotData$row_annotation,
        #annotation_col = PlotData$col_annotation,
        fontsize_row = input$fontsize_row,
        fontsize_col = input$fontsize_col,
        cluster_rows = input$cluster_rows,
        cluster_cols = input$cluster_cols,
        scale = input$Scale,
        #gaps_row = PlotData$gaps_row,
        #gaps_col = PlotData$gaps_col,
        clustering_method = input$clustering_method,
        angle_col = 45,
        cellwidth = PlotData$cellwidth,
        cellheight = PlotData$cellheight,
        show_rownames = input$show_rownames,
        show_colnames = input$show_colnames
        #color = colorRampPalette(colors = c("blue", "white", "red"))(50)
      )
      
    })
  })
  
  
}