#DEG Analysis

library(DESeq2)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(apeglm)
#RandomForest 
ui.deg <- function(id) {
  deg <- NS(id)
  
  fluidPage(
    useShinyjs(),
    titlePanel("Multi-group Differential Gene Expression Analysis By DESeq2"),
    
    sidebarLayout(
      sidebarPanel(
        style = "overflow-y: scroll; height: 600px;",
        useShinyjs(),
        # 左侧参数调整区域
        tags$style(type="text/css",
                   ".shiny-output-error { visibility: hidden; }",
                   ".shiny-output-error:before { visibility: hidden; }",
                   ".shiny-output-error:after { visibility: hidden; }"),
        
        
        fileInput(deg("Exper"), "Upload exprMatr file(csv):"),
        fileInput(deg("Group"), "Upload Group file(CSV):"),
        selectInput(deg("padjust"), "Adjusted p-value cutoff:",
                    choices = c(0.05, 0.01)),
        selectInput(deg("LFC"), "Log2 Fold Change (log2FC) Values Cutoff:",
                    choices = c(0.5, 1,1.5)),
        actionButton(deg("run"), " Click to Run!!!"),
        downloadButton(deg("download"), "Download All Result",disabled = TRUE)
      ),
      mainPanel(
        # 右侧图形展示区域
        wellPanel(
          DTOutput(deg("NOG"))
        )
      )
    )
  )
} 

server.deg <- function(input, output, session) {
  
  create_dir_if_not_exists <- function(dir_path) {
    if (!dir.exists(dir_path)) {
      dir.create(dir_path, showWarnings = FALSE, recursive = TRUE)
      cat("Created directory at:", dir_path, "\n")
    } else {
      cat("Directory already exists at:", dir_path, "\n")
    }
  }
  
  # Initialize paths
  
  temp_folder <- gsub("\\\\", "/", tempdir())
  result_folder <- file.path(temp_folder, "DESeq2")
  
  observeEvent(input$run, {
    
    tryCatch({
      req(input$Exper)
      req(input$Group)
      shinyjs::disable("run")
      shinyjs::disable("download")
      
      
      
      # Create directories if they don't exist
      create_dir_if_not_exists(result_folder)
      files <- list.files(path = result_folder, full.names = TRUE)
      
      # 删除所有文件
      file.remove(files)
      
      # 读取数据
      expr_file <- input$Exper$datapath
      group_file <- input$Group$datapath
      
      count <- read.table(expr_file,header = T,row.names = 1,sep=",")
      data <- read.table(group_file,header = T,row.names = 1,sep=",")
      
      colnames(data) <- c("Group")
      
      group <- factor(data$Group, levels = unique(data$Group))#分组信息
      group_summary <- table(group) 
      # 如果需要，获取每组的样本数目
      num <- as.numeric(group_summary)#每个组包含的样本数量
      
      group <-  levels(group)
    
      
      print(num)
      Batch_Deseq_differnece<-function(exprSet,group,num){
        try({
          ##create a folder 
          save_dir<-file.path(result_folder,"DEG")
          dir.create(save_dir)
          ## creat a group
          print(group)
          print(num)
          replicated_groups <- unlist(mapply(rep, group, times=num))
          group_list= factor(replicated_groups,levels=group)
          group_list
          print(group_list)
          colData=data.frame(row.names = colnames(exprSet),
                             group=group_list)
          
          print(colData)
          
          print(group_list)
          loopcount = 1
          #dat<-data.frame()
          ## use the Deseq2 to have Diffence analyse
          for (i in 1:length(group)){
            name=unique(group)[i]
            print(name)
            colData$group<-relevel(colData$group,ref=name)
            dds=DESeq2::DESeqDataSetFromMatrix(countData = round(exprSet),#整数化#如果是
                                               colData = colData,
                                               design = ~group) 
            dds <- dds[ rowSums(DESeq2::counts(dds)) > 10, ]
            dds <- DESeq2::DESeq(dds)
            print(dds)
            
            for (j in 2:length(DESeq2::resultsNames(dds))){
              
              try({
                resname=DESeq2::resultsNames(dds)[j]
                print(resname)
                
                res=DESeq2::results(dds, name=resname)
                
                res_lfc <- lfcShrink(dds, coef=j, res=res, type="apeglm")
                res_lfc
                #res=res_lfc
                
                summary(res_lfc)
                summary(res)
                
                dir.create(save_dir)
                
                #write.csv(res,paste0(save_dir,resname,".csv"))
                
                save_dir2=file.path(result_folder,"MA")
                dir.create(save_dir2)
                
                
                
                save_dir_MA=paste0(save_dir2,"/",resname)
                dir.create(save_dir_MA)
                
                write.csv(res,paste0(save_dir_MA,"/",resname,"_res.csv"))
                write.csv(res_lfc,paste0(save_dir_MA,"/",resname,"_reslfc.csv"))
                
                
                DEG <- res
   
                DEG <- na.omit(DEG)
                
                # 确定差异表达倍数,abs表示绝对值
                logFC_cutoff <- with(DEG,mean(abs(log2FoldChange)) + 2*sd(abs(log2FoldChange)))
                # 取前两位小数
                
                logFC_cutoff <- round(logFC_cutoff, 2)
                DEG$change = as.factor(ifelse(DEG$padj < input$padjust & abs(DEG$log2FoldChange) > input$LFC,
                                              ifelse(DEG$log2FoldChange >  1,'UP','DOWN'),'STABLE'))
                
                DEG <- data.frame(DEG)
                
                write.csv(DEG,paste0(save_dir,"/",resname,"_DEG.csv"))
                
              })
            }
            
          }
          
          
        })
      }
      
      Batch_Deseq_differnece(count,group=group,num=num)
      
      
      
      shinyjs::enable("run")
      shinyjs::enable("download")
    }, error = function(e) {
      #捕获并显示错误
      shinyjs::enable("run") # if ERROR has merge, we can run again!!
      shiny::showNotification(paste("Error:", e$message), type = "error")
    })
  })
  
  output$download <- downloadHandler(
    filename = function() {
      paste("DESeq2_Result_", Sys.Date(), ".zip", sep = "")
    },
    content = function(file) {
      temp_folder <- gsub("\\\\", "/", tempdir())
      result_folder <- file.path(temp_folder, "DESeq2")
      
      zip_file_path <- file.path(temp_folder, "DESeq2_Result.zip")
      
      if (dir.exists(result_folder)) {
        # Zip the result folder contents
        files_to_zip <- list.files(result_folder, full.names = TRUE)
        zip(zipfile = zip_file_path, files = files_to_zip)
        file.copy(zip_file_path, file)
      } else {
        showNotification("No files available for download.", type = "error")
      }
    }
  )
  
}