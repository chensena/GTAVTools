library(randomForest)
library(caret)
library(Boruta)
library(dplyr)
library(reshape2)
library(pROC)
library(ggplot2)
library(ggrepel)
library(MLmetrics)
library(ImageGP)
#RandomForest 
ui.randomforest <- function(id) {
  rf <- NS(id)
  
  fluidPage(
    useShinyjs(),
    titlePanel("RandomForest"),
    
    sidebarLayout(
      sidebarPanel(
        useShinyjs(),
        # 左侧参数调整区域
        tags$style(type="text/css",
                   ".shiny-output-error { visibility: hidden; }",
                   ".shiny-output-error:before { visibility: hidden; }",
                   ".shiny-output-error:after { visibility: hidden; }"),
        
        fileInput(rf("Exper"), "Upload exprMatr file(csv):"),
        fileInput(rf("Trait"), "Upload Trait file(CSV):"),
        selectInput(rf("ModelType"), "Select Model Type:",
                    choices = c("multinomial", "binomial")),
        
        actionButton(rf("run"), " Click to Run!!!"),
        downloadButton(rf("download"), "Download All Result",disabled = TRUE)
      ),
      mainPanel(
        # 右侧图形展示区域
        wellPanel(
          DTOutput(rf("NOG"))
        )
      )
    )
  )
}
server.randomforest <- function(input, output, session) {

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
  result_folder <- file.path(temp_folder, "RandomForest_Result")
    
  observeEvent(input$run, {
    tryCatch({
      req(input$Exper)
      req(input$Trait)
      shinyjs::disable("run")
      shinyjs::disable("download")
      
      
      
      # Create directories if they don't exist
      create_dir_if_not_exists(result_folder)
      files <- list.files(path = result_folder, full.names = TRUE)
      
      # 删除所有文件
      file.remove(files)
      
      
      # 读取数据
      expr_file <- input$Exper$datapath
      metadata_file <- input$Trait$datapath
      
      if(input$ModelType == "binomial")
      {
        
        # 读取数据
        expr <- read.table(expr_file, row.names = 1, header = TRUE, sep = ",", stringsAsFactors = TRUE, check.names = FALSE)
        metadata <- read.table(metadata_file, row.names = 1, header = TRUE, sep = ",", stringsAsFactors = TRUE)
        
        # 处理表达数据
        m.mad <- apply(expr, 1, mad)
        expr_mat <- expr[which(m.mad > max(quantile(m.mad, probs = seq(0, 1, 0.25))[2], 0.01)), ]
        rownames(expr_mat) <- make.names(rownames(expr_mat))
        
        # 转置表达数据
        expr_mat <- t(expr_mat)
        
        # 确保表达数据和metadata样品一致
        expr_mat_sampleL <- rownames(expr_mat)
        metadata_sampleL <- rownames(metadata)
        common_sampleL <- intersect(expr_mat_sampleL, metadata_sampleL)
        expr_mat <- expr_mat[common_sampleL, , drop = FALSE]
        metadata <- metadata[common_sampleL, , drop = FALSE]
        colnames(metadata) <- c("class")
        # 判断回归还是分类
        group <- "class"
        metadata[[group]] <- as.factor(metadata[[group]])
        
        # 根据group的类型进行转换
        if (is.numeric(metadata[[group]])) {
          metadata[[group]] <- as.numeric(metadata[[group]])
        } else {
          metadata[[group]] <- as.factor(metadata[[group]])
        }
        
        # 设置随机数种子
        set.seed(304)
        
        # 随机森林模型训练
        rf <- randomForest(expr_mat, metadata[[group]], importance = TRUE)
        
        # 拆分训练集和测试集
        set.seed(1)
        train_index <- createDataPartition(metadata[[group]], p = 0.75, list = FALSE)
        train_data <- expr_mat[train_index, ]
        train_data_group <- metadata[[group]][train_index]
        test_data <- expr_mat[-train_index, ]
        test_data_group <- metadata[[group]][-train_index]
        
        # 特征选择
        set.seed(1)
        boruta <- Boruta(x = train_data, y = train_data_group, pValue = 0.01, mcAdj = TRUE, maxRuns = 300)
        
        # 处理Boruta结果
        boruta_imp <- function(x) {
          imp <- melt(x$ImpHistory, na.rm = TRUE)[, -1]
          colnames(imp) <- c("Variable", "Importance")
          imp <- imp[is.finite(imp$Importance), ]
          
          variableGrp <- data.frame(Variable = names(x$finalDecision), finalDecision = x$finalDecision)
          showGrp <- data.frame(Variable = c("shadowMax", "shadowMean", "shadowMin"),
                                finalDecision = c("shadowMax", "shadowMean", "shadowMin"))
          
          variableGrp <- rbind(variableGrp, showGrp)
          boruta_variable_imp <- merge(imp, variableGrp, all.x = TRUE)
          
          sortedVariable <- boruta_variable_imp %>% group_by(Variable) %>% 
            summarise(median = median(Importance)) %>% arrange(median)
          sortedVariable <- as.vector(sortedVariable$Variable)
          
          boruta_variable_imp$Variable <- factor(boruta_variable_imp$Variable, levels = sortedVariable)
          
          invisible(boruta_variable_imp)
        }
        
        boruta_variable_imp <- boruta_imp(boruta)
        
        # 绘制特征重要性图plot output Feature Importance Distribution Boxplot

        p1 <- sp_boxplot(boruta_variable_imp, melted = TRUE, xvariable = "Variable", yvariable = "Importance",
                   legend_variable = "finalDecision", legend_variable_order = c("shadowMax", "shadowMean", "shadowMin", "Confirmed"),
                   xtics_angle = 90)
        Boxplot_path <- file.path(result_folder, "Feature Importance Distribution Boxplot.pdf")
        pdf(file = Boxplot_path, width = 11, height = 8.5)
        print(p1)
        dev.off()
        
        
        # 特征选择的最终结果
        boruta_finalVarsWithTentative <- data.frame(Item = getSelectedAttributes(boruta, withTentative = TRUE), Type = "Boruta_with_tentative")
        
        # 创建模型并训练
        
        generateTestVariableSet <- function(num_toal_variable){
          max_power <- ceiling(log10(num_toal_variable))
          tmp_subset <- c(unlist(sapply(1:max_power, function(x) (1:10)^x, simplify = F)), ceiling(max_power/3))
          #return(tmp_subset)
          base::unique(sort(tmp_subset[tmp_subset<num_toal_variable]))
        }
        boruta_train_data <- train_data[, boruta_finalVarsWithTentative$Item]
        boruta_mtry <- generateTestVariableSet(ncol(boruta_train_data))
        
        # 模型训练
        trControl <- trainControl(method = "repeatedcv", number = 10, repeats = 5)
        tuneGrid <- expand.grid(mtry = boruta_mtry)
        
        set.seed(1)
        borutaConfirmed_rf_default <- train(x = boruta_train_data, y = train_data_group, method = "rf", 
                                            tuneGrid = tuneGrid, metric = "Accuracy", trControl = trControl)
        #plot output 模型性能图 Model Performance Plot

        p2 <- plot(borutaConfirmed_rf_default)
        Performance_path <- file.path(result_folder, "Model Performance Plot.pdf")
        pdf(file = Performance_path, width = 11, height = 8.5)
        print(p2)
        dev.off()
        
        
        #plot output 特征主要性点图 Feature Importance Dot Plot

        p3 <- dotPlot(varImp(borutaConfirmed_rf_default))
        Importance_path <- file.path(result_folder, "Feature Importance Dot Plot.pdf")
        pdf(file = Importance_path, width = 11, height = 8.5)
        print(p3)
        dev.off()
        
        
        importance <- varImp(borutaConfirmed_rf_default)
        imp_results <- as.data.frame(importance$importance)
        
        #import gene output
        write.csv(imp_results ,file.path(result_folder,"IMPORTANCE.csv"))
        
        
        # 模型评估
        predictions_train <- predict(borutaConfirmed_rf_default$finalModel, newdata = train_data)
        confusionMatrix(predictions_train, train_data_group)
        
        prediction_prob <- predict(borutaConfirmed_rf_default$finalModel, newdata = test_data, type = "prob")
        roc_curve <- roc(test_data_group, prediction_prob[, 1])
        
        # 结果评估y
        predictions <- predict(borutaConfirmed_rf_default$finalModel, newdata = test_data)
        confusionMatrix(predictions, test_data_group)
        
        best_thresh <- data.frame(coords(roc = roc_curve, x = "best", input = "threshold", 
                                         transpose = FALSE, best.method = "youden"))
        best_thresh$best <- apply(best_thresh, 1, function(x) 
          paste0('threshold: ', x[1], ' (', round(1 - x[2], 3), ", ", round(x[3], 3), ")"))
        
        # ROC曲线绘制
        ROC_data <- data.frame(FPR = 1 - roc_curve$specificities, TPR = roc_curve$sensitivities)
        ROC_data <- ROC_data[with(ROC_data, order(FPR, TPR)), ]
        
        p <- ggplot(data = ROC_data, mapping = aes(x = FPR, y = TPR)) +
          geom_step(color = "red", size = 1, direction = "vh") +
          geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1)) +
          theme_classic() +
          xlab("False positive rate") + 
          ylab("True positive rate") +
          coord_fixed(1) +
          
          xlim(0, 1) +
          
          ylim(0, 1) +
          
          annotate('text', x = 0.5, y = 0.25, label = paste('AUC =', round(roc_curve$auc, 2))) +
          geom_point(data = best_thresh, mapping = aes(x = 1 - specificity, y = sensitivity), color = 'blue', size = 2) +
          geom_text_repel(data = best_thresh, mapping = aes(x = 1.05 - specificity, y = sensitivity, label = best))
        
        
        ROC_path <- file.path(result_folder, "ROC Plot.pdf")
        pdf(file = ROC_path, width = 11, height = 8.5)  
        print(p)
        dev.off()
        
      }
      else{

        #dev.off()              # 读取数据
        expr <- read.table(expr_file, row.names = 1, header = TRUE, sep = "\t", stringsAsFactors = TRUE, check.names = FALSE)
        metadata <- read.table(metadata_file, row.names = 1, header = TRUE, sep = "\t", stringsAsFactors = TRUE)
        
        # 处理表达数据
        m.mad <- apply(expr, 1, mad)
        expr_mat <- expr[which(m.mad > max(quantile(m.mad, probs = seq(0, 1, 0.25))[2], 0.01)), ]
        rownames(expr_mat) <- make.names(rownames(expr_mat))
        
        # 转置表达数据
        expr_mat <- t(expr_mat)
        
        # 确保表达数据和metadata样品一致
        expr_mat_sampleL <- rownames(expr_mat)
        metadata_sampleL <- rownames(metadata)
        common_sampleL <- intersect(expr_mat_sampleL, metadata_sampleL)
        expr_mat <- expr_mat[common_sampleL, , drop = FALSE]
        metadata <- metadata[common_sampleL, , drop = FALSE]
        colnames(metadata) <- c("class")
        # 判断回归还是分类
        group <- "class"
        
        # 根据group的类型进行转换
        metadata[[group]] <- as.factor(metadata[[group]])
        
        # 设置随机数种子
        set.seed(304)
        
        # 随机森林模型训练（多分类）
        rf <- randomForest(expr_mat, metadata[[group]], importance = TRUE, ntree = 500)
        
        # 拆分训练集和测试集
        set.seed(1)
        train_index <- createDataPartition(metadata[[group]], p = 0.75, list = FALSE)
        train_data <- expr_mat[train_index, ]
        train_data_group <- metadata[[group]][train_index]
        test_data <- expr_mat[-train_index, ]
        test_data_group <- metadata[[group]][-train_index]
        
        # 特征选择
        set.seed(1)
        boruta <- Boruta(x = train_data, y = train_data_group, pValue = 0.01, mcAdj = TRUE, maxRuns = 300)
        
        # 处理Boruta结果
        boruta_imp <- function(x) {
          imp <- melt(x$ImpHistory, na.rm = TRUE)[, -1]
          colnames(imp) <- c("Variable", "Importance")
          imp <- imp[is.finite(imp$Importance), ]
          
          variableGrp <- data.frame(Variable = names(x$finalDecision), finalDecision = x$finalDecision)
          
          showGrp <- data.frame(Variable = c("shadowMax", "shadowMean", "shadowMin"),
                                finalDecision = c("shadowMax", "shadowMean", "shadowMin"))
          
          variableGrp <- rbind(variableGrp, showGrp)
          boruta_variable_imp <- merge(imp, variableGrp, all.x = TRUE)
          
          sortedVariable <- boruta_variable_imp %>% group_by(Variable) %>% 
            summarise(median = median(Importance)) %>% arrange(median)
          sortedVariable <- as.vector(sortedVariable$Variable)
          
          boruta_variable_imp$Variable <- factor(boruta_variable_imp$Variable, levels = sortedVariable)
          
          invisible(boruta_variable_imp)
        }
        
        
        boruta_variable_imp <- boruta_imp(boruta)
        
        # 绘制特征重要性图

        p1 <- sp_boxplot(boruta_variable_imp, melted = TRUE, xvariable = "Variable", yvariable = "Importance",
                   legend_variable = "finalDecision", legend_variable_order = c("shadowMax", "shadowMean", "shadowMin", "Confirmed"),
                   xtics_angle = 90)
        Boxplot_path <- file.path(result_folder, "Feature Importance Distribution Boxplot.pdf")
        pdf(file = Boxplot_path, width = 11, height = 8.5)
        print(p1)
        dev.off()
        
        # 特征选择的最终结果
        boruta_finalVarsWithTentative <- data.frame(Item = getSelectedAttributes(boruta, withTentative = TRUE), Type = "Boruta_with_tentative")
        
        # 创建模型并训练
        boruta_train_data <- train_data[, boruta_finalVarsWithTentative$Item]
        generateTestVariableSet <- function(num_toal_variable){
          max_power <- ceiling(log10(num_toal_variable))
          tmp_subset <- c(unlist(sapply(1:max_power, function(x) (1:10)^x, simplify = F)), ceiling(max_power/3))
          #return(tmp_subset)
          base::unique(sort(tmp_subset[tmp_subset<num_toal_variable]))
        }
        boruta_mtry <- generateTestVariableSet(ncol(boruta_train_data))
        
        # 模型训练（多分类）
        trControl <- trainControl(method = "repeatedcv", number = 10, repeats = 5, 
                                  classProbs = TRUE, summaryFunction = multiClassSummary)
        tuneGrid <- expand.grid(mtry = boruta_mtry)
        
        set.seed(1)
        borutaConfirmed_rf_default <- train(x = boruta_train_data, y = train_data_group, method = "rf", 
                                            tuneGrid = tuneGrid, metric = "Accuracy", trControl = trControl)
        
        # 绘制模型性能图

        p2 <- plot(borutaConfirmed_rf_default)
        Performance_path <- file.path(result_folder, "Model Performance Plot.pdf")
        pdf(file = Performance_path, width = 11, height = 8.5)
        print(p2)
        dev.off()
        
        # 绘制特征重要性点图

        p3 <- dotPlot(varImp(borutaConfirmed_rf_default))
        Importance_path <- file.path(result_folder, "Feature Importance Dot Plot.pdf")
        pdf(file = Importance_path, width = 11, height = 8.5)
        print(p3)
        dev.off()
        
        importance <- varImp(borutaConfirmed_rf_default)
        imp_results <- as.data.frame(importance$importance)
        
        # 导出重要性结果
        write.csv(imp_results ,file.path(result_folder,"IMPORTANCE.csv"))
        
        #print("importance!!")
        # 模型评估
        predictions_train <- predict(borutaConfirmed_rf_default$finalModel, newdata = train_data)
        confusionMatrix(predictions_train, train_data_group)
        print(1)
        prediction_prob <- predict(borutaConfirmed_rf_default$finalModel, newdata = test_data, type = "prob")
        # 创建一个函数来计算 ROC 数据
        calculate_roc_data <- function(true_labels, prob_matrix) {
          roc_data_list <- lapply(levels(true_labels), function(cls) {
            roc_curve <- roc(true_labels == cls, prob_matrix[, cls])
            data.frame(
              FPR = 1 - roc_curve$specificities,
              TPR = roc_curve$sensitivities,
              Class = cls
            )
          })
          do.call(rbind, roc_data_list)
        }
        
        # 计算 ROC 数据
        roc_data <- calculate_roc_data(test_data_group, prediction_prob)
        
        # 计算每个类别的 AUC 值
        roc_list <- lapply(levels(test_data_group), function(cls) {
          roc(test_data_group == cls, prediction_prob[, cls])
        })
        auc_values <- sapply(roc_list, function(r) auc(r))
        mean_auc <- mean(auc_values)  # 计算所有类别的平均 AUC
        
        # 绘制 ROC 曲线
        p <- ggplot(data = roc_data, aes(x = FPR, y = TPR, color = Class)) +
          geom_step(size = 1) +
          geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
          theme_classic() +
          xlab("False Positive Rate") +
          ylab("True Positive Rate") +
          ggtitle(paste("Multiclass ROC Curve (Mean AUC:", round(mean_auc, 2), ")"))

        ROC_path <- file.path(result_folder, "ROC Plot.pdf")
        pdf(file = ROC_path, width = 11, height = 8.5)  
        print(p)
        dev.off()
      }
      
      
      
      
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
      paste("RandomForest_Result_", Sys.Date(), ".zip", sep = "")
    },
    content = function(file) {
      temp_folder <- gsub("\\\\", "/", tempdir())
      result_folder <- file.path(temp_folder, "RandomForest_Result")
      
      zip_file_path <- file.path(temp_folder, "RandomForest_Result.zip")
      
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