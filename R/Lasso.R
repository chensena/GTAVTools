library(glmnet)
# devtools::install_github("Tong-Chen/ImageGP")
library(caret)
library(dplyr)
library(purrr)
library(ggplot2)
library(patchwork)
library(pROC) # 处理ROC曲线
library(DT)


# Lasso 回归分析
ui.lasso <- function(id) {
  la <- NS(id)

  fluidPage(
    useShinyjs(),
    titlePanel("Least absolute shrinkage and selection operator (Lasso)"),
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
        fileInput(la("Exper"), "Upload exprMatr file(csv):"),
        fileInput(la("Trait"), "Upload Trait file(CSV):"),
        selectInput(la("ModelType"), "Select Model Type:",
          choices = c("multinomial", "binomial")
        ),
        actionButton(la("run"), " Click to Run!!!"),
        downloadButton(la("download"), "Download All Result", disabled = TRUE)
      ),
      mainPanel(
        # 右侧图形展示区域
        wellPanel(
          DTOutput(la("NOG"))
        )
      )
    )
  )
}
server.lasso <- function(input, output, session) {
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
  result_folder <- file.path(temp_folder, "Lasso_Result")

  observeEvent(input$run, {
    tryCatch(
      {
        req(input$Exper)
        req(input$Trait)
        shinyjs::disable("run")



        # Create directories if they don't exist
        create_dir_if_not_exists(result_folder)
        files <- list.files(path = result_folder, full.names = TRUE)

        # 删除所有文件
        file.remove(files)


        # 读取数据
        expr_file <- input$Exper$datapath
        metadata_file <- input$Trait$datapath
        expr_mat <- read.table(expr_file, row.names = 1, header = TRUE, sep = ",", stringsAsFactors = TRUE)
        metadata <- read.table(metadata_file, row.names = 1, header = TRUE, sep = ",")
        colnames(metadata) <- c("class")

        # 处理基因名字
        rownames(expr_mat) <- make.names(rownames(expr_mat))

        # 转置数据矩阵
        expr_mat <- t(expr_mat)

        if (input$ModelType == "binomial") { # 二分类数据
          metadata$class <- as.factor(metadata$class)

          set.seed(1)
          train_index <- createDataPartition(metadata$class, p = 0.75, list = FALSE)
          train_data <- expr_mat[train_index, ]
          train_data_group <- metadata$class[train_index]
          test_data <- expr_mat[-train_index, ]
          test_data_group <- metadata$class[-train_index]

          lambda <- 0.01

          la.eq <- glmnet(as.matrix(train_data), as.matrix(train_data_group),
            lambda = lambda, family = "binomial",
            type.measure = "class", intercept = FALSE, alpha = 1
          )

          mod_cv <- cv.glmnet(as.matrix(train_data), as.matrix(train_data_group),
            family = "binomial", type.measure = "class",
            intercept = FALSE, alpha = 1
          )
          cv_model_path <- file.path(result_folder, "Cross-Validation Error for Lasso Model.pdf")
          pdf(file = cv_model_path, width = 11, height = 8.5)
          plot(mod_cv)
          dev.off()

          best_lambda <- mod_cv$lambda.min

          alphalist <- seq(0, 1, by = 0.1)
          elasticnet <- lapply(alphalist, function(a) {
            cv.glmnet(as.matrix(train_data), as.matrix(train_data_group),
              alpha = a, family = "binomial", lambda.min.ratio = 0.001
            )
          })

          cvm_values <- sapply(elasticnet, function(model) min(model$cvm))
          best_alpha_index <- which.min(cvm_values)
          best_alpha <- alphalist[best_alpha_index]

          cv_model_ls %>%
            map2(names(.), function(model, alpha) {
              lambda_min <- model$lambda[model$index_min] %>% log()
              lambda_1se <- model$lambda[model$index_1se] %>% log()

              data.frame(
                x = log(model$lambda),
                y = model$cvm,
                std = model$cvsd
              ) %>%
                ggplot(aes(x = x, y = y)) +
                geom_point(color = "Red") +
                geom_errorbar(aes(ymin = y - std, ymax = y + std), color = "gray") +
                geom_vline(xintercept = c(lambda_min, lambda_1se), linetype = 3) +
                geom_text(aes(x = lambda_min + 3, y = 0.5, label = paste("lambda_Min=", round(lambda_min, 2)))) +
                labs(x = expression(paste("Log(", lambda, ")")), y = "Binomial Deviance") +
                theme_bw() +
                ggtitle(paste0("alpha: ", alpha))
            }) %>%
            wrap_plots(ncol = 4)

          # Cross-Validation Error vs. Log(lambda) for Different Alpha Values
          Alpha_path <- file.path(result_folder, "Cross-Validation Error vs. Log(lambda) for Different Alpha Values.pdf")
          pdf(file = Alpha_path, width = 11, height = 8.5)
          cv_model_ls
          dev.off()

          p <- cv_model_ls %>%
            map(function(x) list(cvm = x$cvm, lambda = x$lambda)) %>%
            map(as_tibble) %>%
            enframe(name = "alpha") %>%
            unnest(value) %>%
            mutate(lambda = log(lambda)) %>%
            ggplot(aes(x = lambda, y = cvm)) +
            geom_line(aes(color = alpha), linewidth = 1) +
            labs(x = expression(paste("Log(", lambda, ")")), y = "Binomial Deviance") +
            ggthemes::theme_base()

          Alpha_path2 <- file.path(result_folder, "Cross-Validation Error vs. Log(lambda) for Different Alpha Values(Total).pdf")
          pdf(file = Alpha_path2, width = 11, height = 8.5)
          p + (p + xlim(c(-5.5, -2)) + ylim(c(0, 0.3))) +
            guide_area() +
            plot_layout(guides = "collect", ncol = 3, widths = c(5, 5, 1))
          dev.off()

          best_model <- glmnet(as.matrix(train_data), as.matrix(train_data_group),
            alpha = best_alpha, lambda = best_lambda,
            family = "binomial", type.measure = "class",
            intercept = FALSE
          )

          coef.min <- coef(best_model)
          active.min <- which(coef.min@i != 0)
          lasso_geneids <- coef.min@Dimnames[[1]][coef.min@i + 1]

          write(lasso_geneids, file.path(result_folder, "lasso_genes.csv"))


          pre_train <- predict(best_model, s = best_lambda, newx = as.matrix(train_data), type = "response")
          act_train <- ifelse(train_data_group == "NOT", 1, 0)

          pre_test <- predict(best_model, s = best_lambda, newx = as.matrix(test_data), type = "response")
          act_test <- ifelse(test_data_group == "NOT", 1, 0)

          # 开始绘图！
          ROC_plot_path <- file.path(result_folder, "ROC_Plot.pdf")
          pdf(file = AROC_plot_path, width = 11, height = 8.5)

          train <- plot.roc(act_train, pre_train, smooth = FALSE, lwd = 2, ylim = c(0, 1), xlim = c(0.95, 0), legacy.axes = TRUE, main = "", col = "red")
          test <- plot.roc(act_test, pre_test, smooth = FALSE, add = TRUE, lwd = 2, ylim = c(0, 1), xlim = c(0.95, 0), legacy.axes = TRUE, main = "", col = "green")

          train[["auc"]]
          test[["auc"]]

          legend.name <- c(
            paste("Train:AUC", sprintf("%.4f", train[["auc"]])),
            paste("Test:AUC", sprintf("%.4f", test[["auc"]]))
          )

          legend("bottomright",
            legend = legend.name,
            lwd = 2,
            col = c("red", "green"),
            bty = "n"
          )

          dev.off()

          # 绘图结束！
        } else if (input$ModelType == "multinomial") {
          metadata$class <- as.factor(metadata$class) # 转换类别列为因子类型

          set.seed(1)
          train_index <- createDataPartition(metadata$class, p = 0.75, list = FALSE)
          train_data <- expr_mat[train_index, ]
          train_data_group <- metadata$class[train_index]
          test_data <- expr_mat[-train_index, ]
          test_data_group <- metadata$class[-train_index]

          lambda <- 0.01

          la.eq <- glmnet(as.matrix(train_data), as.matrix(train_data_group),
            lambda = lambda, family = "multinomial",
            type.measure = "class", intercept = FALSE, alpha = 1
          )

          mod_cv <- cv.glmnet(as.matrix(train_data), as.matrix(train_data_group),
            family = "multinomial", type.measure = "class",
            intercept = FALSE, alpha = 1
          )

          cv_model_path <- file.path(result_folder, "Cross-Validation Error for Lasso Model.pdf")
          pdf(file = cv_model_path, width = 11, height = 8.5)
          plot(mod_cv)
          dev.off()


          best_lambda <- mod_cv$lambda.min

          alphalist <- seq(0, 1, by = 0.1)
          elasticnet <- lapply(alphalist, function(a) {
            cv.glmnet(as.matrix(train_data), as.matrix(train_data_group),
              alpha = a, family = "multinomial", lambda.min.ratio = 0.001
            )
          })

          cvm_values <- sapply(elasticnet, function(model) min(model$cvm))
          best_alpha_index <- which.min(cvm_values)
          best_alpha <- alphalist[best_alpha_index]


          cv_model_ls <- set_names(alphalist, alphalist) %>%
            map(function(alpha) {
              cv.model <- cv.glmnet(as.matrix(train_data), as.matrix(train_data_group),
                family = "multinomial", alpha = alpha,
                lambda.min.ratio = 0.001
              )
              list(
                cvm = cv.model$cvm,
                cvsd = cv.model$cvsd,
                lambda = cv.model$lambda,
                index_min = cv.model$index["min", 1],
                index_1se = cv.model$index["1se", 1]
              )
            })


          CV_plot <- cv_model_ls %>%
            map2(names(.), function(model, alpha) {
              lambda_min <- model$lambda[model$index_min] %>% log()
              lambda_1se <- model$lambda[model$index_1se] %>% log()

              data.frame(
                x = log(model$lambda),
                y = model$cvm,
                std = model$cvsd
              ) %>%
                ggplot(aes(x = x, y = y)) +
                geom_point(color = "Red") +
                geom_errorbar(aes(ymin = y - std, ymax = y + std), color = "gray") +
                geom_vline(xintercept = c(lambda_min, lambda_1se), linetype = 3) +
                geom_text(aes(x = lambda_min + 3, y = 0.5, label = paste("lambda_Min=", round(lambda_min, 2)))) +
                labs(x = expression(paste("Log(", lambda, ")")), y = "Multinomial Deviance") +
                theme_bw() +
                ggtitle(paste0("alpha: ", alpha))
            }) %>%
            wrap_plots(ncol = 4)
          Alpha_path <- file.path(result_folder, "Cross-Validation Error vs. Log(lambda) for Different Alpha Values.pdf")
          pdf(file = Alpha_path, width = 11, height = 8.5)
          print(CV_plot)
          dev.off()



          best_model <- glmnet(as.matrix(train_data), as.matrix(train_data_group),
            alpha = best_alpha, lambda = best_lambda,
            family = "multinomial", type.measure = "class",
            intercept = FALSE
          )

          coef_min <- coef(best_model)

          # 提取系数的维度名称
          coef_matrix <- coef_min[[1]]
          nonzero_indices <- which(coef_matrix != 0, arr.ind = TRUE)

          # 提取非零系数对应的特征名称
          feature_names <- rownames(coef_matrix)
          nonzero_features <- feature_names[nonzero_indices[, 1]]
          nonzero_features_per_lambda <- lapply(coef_min, function(coef_matrix) {
            nonzero_indices <- which(coef_matrix != 0, arr.ind = TRUE)
            feature_names <- rownames(coef_matrix)
            feature_names[nonzero_indices[, 1]]
          })
          print(6)
          my_list <- nonzero_features_per_lambda
          max_length <- max(sapply(my_list, length))
          pad_vector <- function(vec, length) {
            length_out <- length - length(vec)
            if (length_out > 0) {
              vec <- c(vec, rep(NA, length_out))
            }
            return(vec)
          }
          print(7)
          # 应用填充函数
          padded_list <- lapply(my_list, pad_vector, length = max_length)
          df <- as.data.frame(padded_list)

          # 保存为 CSV 文件
          write.csv(df, file.path(result_folder, "lasso_genes.csv"))




          # 预测
          pre_train <- predict(best_model, s = best_lambda, newx = as.matrix(train_data), type = "class")
          pre_test <- predict(best_model, s = best_lambda, newx = as.matrix(test_data), type = "class")

          # 计算准确率
          # pre_train <- pre_train[, levels(train_data_group), drop = FALSE]

          # 将预测标签从矩阵转换为向量
          predicted_labels_train <- as.vector(pre_train)
          predicted_labels_test <- as.vector(pre_test)

          # 确保实际标签是因子
          train_data_group <- factor(train_data_group)
          test_data_group <- factor(test_data_group)

          # 计算训练集和测试集的准确度
          accuracy_train <- sum(predicted_labels_train == train_data_group) / length(train_data_group)
          accuracy_test <- sum(predicted_labels_test == test_data_group) / length(test_data_group)

          # 输出准确度

          output_text <- paste(
            "Training Accuracy: ", accuracy_train, "\n",
            "Testing Accuracy: ", accuracy_test, "\n"
          )

          # 将结果写入文本文件
          writeLines(output_text, file.path(result_folder, "accuracy_results.txt"))
        }







        shinyjs::enable("run")
        shinyjs::enable("download")
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
      paste("Lasso_Result_", Sys.Date(), ".zip", sep = "")
    },
    content = function(file) {
      temp_folder <- gsub("\\\\", "/", tempdir())
      result_folder <- file.path(temp_folder, "Lasso_Result")

      zip_file_path <- file.path(temp_folder, "Lasso_Result.zip")

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
