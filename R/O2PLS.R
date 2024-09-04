#O2PLS
library(OmicsPLS)
library(magrittr)
library(ggplot2)
library(ggrepel)


ui.o2pls <- function(id) {
  o2pls <- NS(id)
  
  fluidPage(
    useShinyjs(),
    titlePanel("Orthogonal Partial Least Squares (O2PLS)"),
    
    sidebarLayout(
      sidebarPanel(
        useShinyjs(),
        # 左侧参数调整区域
        tags$style(type="text/css",
                   ".shiny-output-error { visibility: hidden; }",
                   ".shiny-output-error:before { visibility: hidden; }",
                   ".shiny-output-error:after { visibility: hidden; }"),
        
        
        fileInput(o2pls("Exper1"), "Upload  Matrix file 1 (csv):"),
        fileInput(o2pls("Exper2"), "Upload  Matrix file 2 (CSV):"),
        numericInput(o2pls("Load"), "Set Output Loading Number:", value = 20, min = 1, max = 100),
        actionButton(o2pls("run"), " Click to Run!!!"),
        downloadButton(o2pls("download"), "Download All Result",disabled = TRUE)
      ),
      mainPanel(
        # 右侧图形展示区域
        wellPanel(
          DTOutput(o2pls("NOG"))
        )
      )
    )
  )
}
server.o2pls <- function(input, output, session) {
  
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
  result_folder <- file.path(temp_folder, "O2PLS")
  
  observeEvent(input$run, {
    
    tryCatch({
    req(input$Exper1)
    req(input$Exper2)
    shinyjs::disable("run")
    shinyjs::disable("download")
    
    
    
    # Create directories if they don't exist
    create_dir_if_not_exists(result_folder)
    files <- list.files(path = result_folder, full.names = TRUE)
    
    # 删除所有文件
    file.remove(files)
    
    print(1)
    # 读取数据
    expr_file1 <- input$Exper1$datapath
    expr_file2 <- input$Exper2$datapath
    
    gene <- read.table(expr_file1,sep=",",row.names = 1,header = T)
    met <- read.table(expr_file2,sep=",",row.names=1,header=T)
    print(2)
    gene <- na.omit(t(gene))
    
    met <- na.omit(t(met))
    print(3)
    gene <- scale(gene,scale=F)
    met <- scale(met,scale=F)
    
    print(4)
    set.seed(123)
    c <- crossval_o2m(gene,met,2:5,1:3,1:3,nr_folds = 10) #10折交叉验证
    output <- capture.output(c)
    
    
    ax <- str_extract(output, "(?<=ax=)\\d+")
    ax <- na.omit(ax)
    ax <- if(length(ax) > 0) ax[1] else NA
    ax <- as.numeric(ax)
    
    
    ay <- str_extract(output, "(?<=ay=)\\d+")
    ay <- na.omit(ay)
    ay <- if(length(ay) > 0) ay[1] else NA
    ay <- as.numeric(ay)
    
    a <- str_extract(output, "(?<=a=)\\d+")
    a <- na.omit(a)
    a <- if(length(a) > 0) a[1] else NA
    a <- as.numeric(a)
    
    print(5)
    # 使用最佳参数构建最终模型
    modelfit <- o2m(gene, met, a, ax, ay)
    
    print(6)
    xj<- loadings(modelfit, "Xjoint", 1:2) %>% abs %>% rowSums
    xj[-(order(xj,decreasing=T)[1:input$Load])] = 0
    xj <- sign(xj)
    
    print(7)
    write.csv(as.data.frame(xj),file.path(result_folder,"Matrix1_loading.csv"))
    gene1 <- plot(modelfit, loading_name="Xj", i=1, j=2, label = "c", use_ggplot2 = TRUE,
                  alpha = xj,
                  aes(label = stringr::str_sub(colnames(gene), start = 1)),size=3.5,col='red')+
      theme_bw() +
      coord_fixed(1, c(-1,1),c(-1,1)) +
      geom_point(alpha = 0.5+0.5*xj, col = 'blue',size=1.5) +
      labs(title = "Matrix1 joint loadings",
           x = "First Joint Loadings", y = "Second Joint Loadings") +
      theme(plot.title = element_text(face='bold'))
    print(8)
    #data1 <- data.frame(
    #  First_Loading = modelfit$loadings[,1],  # 第一主成分的加载
    #  Second_Loading = modelfit$loadings[,2], # 第二主成分的加载
    #  Label = stringr::str_sub(colnames(gene), start = 1),  # 标签
    #  Alpha = xj  # 透明度
    #)
    
    # 创建 ggplot2 图
    #gene1 <- ggplot(data1, aes(x = First_Loading, y = Second_Loading)) +
    #  geom_point(aes(alpha = 0.5 + 0.5 * Alpha), color = 'blue', size = 1.5) +
    #  geom_text_repel(aes(label = Label), size = 3.5, color = 'red', vjust = 1.5, hjust = 1.5) +
    #  theme_bw() +
    #  coord_fixed(ratio = 1) +
    #  labs(title = "Matrix1 joint loadings",
    #       x = "First Joint Loadings", 
    #       y = "Second Joint Loadings") +
    #  theme(plot.title = element_text(face = 'bold'))
    
    
    Matrix1_path <- file.path(result_folder, "Matrix1 Loading filter.pdf")
    pdf(file = Matrix1_path, width = 11, height = 8.5)
    print(gene1)
    dev.off()
    
    
    print(9)
    #因变量筛选
    yj<- loadings(modelfit, "Yjoint", 1:2) %>% abs %>% rowSums
    yj[-(order(yj,decreasing=T)[1:input$Load])] = 0
    yj <- sign(yj)
    print (yj)
    write.csv(as.data.frame(yj),file.path(result_folder,"Matrix2_loading.csv"))
    
    
    met1 <- plot(modelfit, loading_name="Yj", i=1, j=2, label = "c", use_ggplot2 = TRUE,
                 alpha = yj,
                 aes(label = stringr::str_sub(colnames(met), start = 1)),size=3.5,col='red')+
      theme_bw() +
      coord_fixed(1, c(-1,1),c(-1,1)) +
      geom_point(alpha = 0.5+0.5*yj, col = 'blue',size=1.5) +
      labs(title ="Matrix2 joint loadings",
           x = "First Joint Loadings", y = "Second Joint Loadings") +
      theme(plot.title = element_text(face='bold'))
    print(10)
    #data2 <- data.frame(
    #  First_Loading = modelfit$loadings[,1],  # 第一主成分的加载
    #  Second_Loading = modelfit$loadings[,2], # 第二主成分的加载
    #  Label = stringr::str_sub(colnames(met), start = 1),  # 标签
    #  Alpha = xj  # 透明度
    #)
    
    # 创建 ggplot2 图
    #met1 <- ggplot(data2, aes(x = First_Loading, y = Second_Loading)) +
    #  geom_point(aes(alpha = 0.5 + 0.5 * Alpha), color = 'blue', size = 1.5) +
    #  geom_text_repel(aes(label = Label), size = 3.5, color = 'red', vjust = 1.5, hjust = 1.5) +
    #  theme_bw() +
    #  coord_fixed(ratio = 1) +
    #  labs(title = "Matrix2 joint loadings",
    #       x = "First Joint Loadings", 
    #       y = "Second Joint Loadings") +
    #  theme(plot.title = element_text(face = 'bold'))
    
    
    Matrix2_path <- file.path(result_folder, "Matrix2 Loading filter.pdf")
    pdf(file = Matrix2_path, width = 11, height = 8.5)
    print(met1)
    dev.off()
    
    print(11)
    #前20载荷输出
    gene_looding <- as.data.frame(modelfit$W.)
    meta_looding <- as.data.frame(modelfit$C.)
    print(12)
    abs <- NA
    gene_looding <- cbind(gene_looding,abs)
    meta_looding <- cbind(meta_looding,abs)
    gene_looding[,3] <- abs(gene_looding[,1])
    meta_looding[,3] <- abs(meta_looding[,1])
    #按第3列abs降序#
    gene_looding <- gene_looding[order(gene_looding$abs, decreasing= T), ]
    meta_looding <- meta_looding[order(meta_looding$abs, decreasing= T), ]
    colnames(gene_looding) <- c("pq1","pq2","abs")
    colnames(meta_looding) <- c("pq1","pq2","abs")
    
    print(13)
   #write.table(gene_looding,"Matrix1_looding.list",sep="\t")
    #write.table(meta_looding,"Matrix2_looding.list",sep="\t")  
    write.csv(gene_looding, file.path(result_folder,"Matrix1_looding.csv"))
    write.csv(meta_looding, file.path(result_folder,"Matrix2_looding.csv"))
    
    
    xj<- loadings(modelfit, "Xjoint", 1:2) %>% abs %>% rowSums
    
    yj<- loadings(modelfit, "Yjoint", 1:2) %>% abs %>% rowSums
    
    print(14)
    #write.table(xj,"Matrix1_score.list",sep="\t")
    #write.table(yj,"Matrix2_score.list",sep="\t")
    
    write.csv(xj, file.path(result_folder,"Matrix1_score.csv"))
    write.csv(yj, file.path(result_folder,"Matrix2_score.csv"))
    
  print(15)
    
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
      paste("O2PLS_", Sys.Date(), ".zip", sep = "")
    },
    content = function(file) {
      temp_folder <- gsub("\\\\", "/", tempdir())
      result_folder <- file.path(temp_folder, "O2PLS")
      
      zip_file_path <- file.path(temp_folder, "O2PLS_Result.zip")
      
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