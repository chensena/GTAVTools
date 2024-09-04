#启动子分析


ui.promotor <- function(id) {
  pc<- NS(id)
  
  fluidPage(
    titlePanel("Promotor Element Statistics"),
    
    sidebarLayout(
      sidebarPanel(
        # 左侧参数调整区域
        tags$style(type="text/css",
        ".shiny-output-error { visibility: hidden; }",
        ".shiny-output-error:before { visibility: hidden; }",
        ".shiny-output-error:after { visibility: hidden; }"),
        fileInput(pc("file1"), "Upload Experssion file:"),
        actionButton(pc("Check"), "Read The File"),
        actionButton(pc("Run"), "Element Statistics"),
        downloadButton(pc("download"), "downLoad Result")
      ),
      
      mainPanel(
        # 右侧图形展示区域
        
        wellPanel(
          DTOutput(pc("NOG"))
        )
        
      )
    )
  )
}


server.promotor <- function(input, output, session) {
  
  data <- reactiveValues(
    input = NULL,
    Result = NULL
  )

  observeEvent(input$Check,{
    req(input$file1)
    data$input <-  read.table(input$file1$datapath,na.strings = " ",fill=TRUE,header=F,sep=",")
    
    output$NOG <- renderDT({
      data$input
      
    },
    options = list(pageLength = 10,
                   #autoWidth = TRUE,  # 自适应列宽
                   columnDefs = list(list(width = 'auto', targets = "_all")),  # 自适应列宽
                   paging = TRUE  # 启用分页
                   
    )
    )
  })
  
  observeEvent(input$Run, {
    if (!is.null(data$input)){
    promotor_list <- read.table(input$file1$datapath,na.strings = " ",fill=TRUE,header=T,sep="\t")
    #promotor_list[is.na(promotor_list)] <- 0
    colnames(promotor_list) <- c("ID","motif","motif_seq","Location","Length","Strand","Species","function")
    
    ID_table <- as.data.frame(table(promotor_list$ID))  #获得所有的ID
    motif_table <- as.data.frame(table(promotor_list$motif_seq))  #获得所有的motifseq
    names(ID_table) <- c("ID","num")
    names(motif_table) <-c("motif_seq","total_num")
    ID <- ID_table$ID
    motif_function <-  unique(dplyr::select(promotor_list,c("motif_seq","motif","function")))
    
    
    
    
    for (i in 1:length(ID_table$ID)){
      #print(promotor_list$ID[ID[1]])#获得一个geneID
      tmp <- (subset(promotor_list, ID == ID_table$ID[i])$motif_seq) #获得每个基因预测到的MOtif
      #print(tmp)#得到这个序列的所有motif
      tmp_ID_motif_table <- as.data.frame(table(tmp))
      names(tmp_ID_motif_table) <- c("motif_seq",as.character(ID_table$ID[i][[1]]))
      
      #print(tmp_ID_motif_table)
      #合并
      motif_table <- merge(motif_table,tmp_ID_motif_table,by ="motif_seq",all=TRUE)
      #names(motif_table[i+2]) <- ID_table$ID[i]
      #print(dim(motif_table))
      #x <- ID_table$ID[i]
    }
    
    #此时的motif_tabla就是一个motif_统计表
    #print(4)
    motif_table[is.na(motif_table)] <- 0
    #print(1)
    motif_table<-merge(motif_table,motif_function,by="motif_seq",all=TRUE)
    #print(2)
    data$Result <-  motif_table
    #print(3)
    output$NOG <- renderDT({
    data$Result[,1:5]
      
    },
    options = list(pageLength = 10,
                   #autoWidth = TRUE,  # 自适应列宽
                   columnDefs = list(list(width = 'auto', targets = "_all")),  # 自适应列宽
                   paging = TRUE  # 启用分页
                   
    )
    )
    }
  })
  
  output$download <- downloadHandler(
    filename = function() {
      paste("Promotor Element Statistics", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      
      write.csv(unique(data$Result),file, col.names = NA, row.names = FALSE,quote = FALSE)
      
    }
  )
 
}