library(WGCNA) # biomanger
library(shinyjs)
library(reshape2)
library(stringr)
library(zip)
# one step WGCNA
# 输入表型文件
# 输入表达矩阵
# 模式选择 sign unsign
# run
# 得到最终结果，ZIP压缩包
# pipline
ui.wgcna <- function(id) {
  wg <- NS(id)

  fluidPage(
    useShinyjs(),
    titlePanel("One Step WGCNA"),
    sidebarLayout(
      sidebarPanel(
        # 左侧参数调整区域
        tags$style(
          type = "text/css",
          ".shiny-output-error { visibility: hidden; }",
          ".shiny-output-error:before { visibility: hidden; }",
          ".shiny-output-error:after { visibility: hidden; }"
        ),
        fileInput(wg("Exper"), "Upload Experssion Data:"),
        fileInput(wg("Trait"), "Upload Trait Data:(Option)"),
        selectInput(wg("NetworkType"), "Select NetWork Type:",
          choices = c("signed", "unsigned"),
          selected = "signed"
        ),
        selectInput(wg("corType"), "Selecting Correlation Methods:",
          choices = c("pearson", "sperman"),
          selected = "pearson"
        ),
        actionButton(wg("run"), " Run WGCNA!!!"),
        downloadButton(wg("download"), "Download Result", disabled = TRUE)
      ),
      mainPanel(
        # 右侧图形展示区域,这里可以展示运行进度
        wellPanel(
          plotOutput(wg("WGCNA_plot"), height = "600px")
        )
      )
    )
  )
}


library(shiny)
library(shinyjs)
library(zip)

server.wgcna <- function(input, output, session) {
  # 使用 reactiveVal 存储临时文件夹路径
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
  result_folder <- file.path(temp_folder, "WGCNA_Result")

  # Create directories if they don't exist


  # 确保结果文件夹存在



  print(temp_folder)
  observeEvent(input$run, {
    tryCatch(
      {
        req(input$Exper) # 确保输入文件存在
        print("runing!!!")
        shinyjs::disable("run")
        shinyjs::disable("run")
        # 设置选项和初始化
        create_dir_if_not_exists(result_folder)
        options(stringsAsFactors = FALSE)
        enableWGCNAThreads(nThreads = 10)

        exprMat <- input$Exper$datapath

        # 选择网络类型
        type <- if (input$NetworkType == "signed") "signed" else "unsigned"
        # 读取数据
        dataExpr <- read.table(exprMat, sep = ",", header = TRUE, row.names = 1, quote = "", comment = "", check.names = FALSE)
        dataExpr <- log(dataExpr + 1)
        # 过滤数据
        m.mad <- apply(dataExpr, 1, mad)
        dataExprVar <- dataExpr[which(m.mad > max(quantile(m.mad, probs = seq(0, 1, 0.25))[2], 0.01)), ]
        dataExpr <- as.data.frame(t(dataExprVar))
        # 检查样本和基因
        gsg <- goodSamplesGenes(dataExpr, verbose = 3)
        if (!gsg$allOK) {
          if (sum(!gsg$goodGenes) > 0) printFlush(paste("Removing genes:", paste(names(dataExpr)[!gsg$goodGenes], collapse = ",")))
          if (sum(!gsg$goodSamples) > 0) printFlush(paste("Removing samples:", paste(rownames(dataExpr)[!gsg$goodSamples], collapse = ",")))
          dataExpr <- dataExpr[gsg$goodSamples, gsg$goodGenes]
        }
        # 计算样本树
        nGenes <- ncol(dataExpr)
        nSamples <- nrow(dataExpr)
        sampleTree <- hclust(dist(dataExpr), method = "average")

        # 保存图形
        tree_file_path <- file.path(result_folder, "sample_tree_plot.pdf")
        pdf(file = tree_file_path, width = 11, height = 8.5)
        plot(sampleTree, main = "Sample clustering to detect outliers", sub = "", xlab = "")
        dev.off()
        print("样本聚类结束！！")
        print(tree_file_path)
        # 选择软阈值
        powers <- c(c(1:10), seq(from = 12, to = 30, by = 2))
        sft <- pickSoftThreshold(dataExpr, powerVector = powers, networkType = type, verbose = 5)
        #
        # # 绘制图形
        power_file_path <- file.path(result_folder, "Soft_Threshold_plot.pdf")
        pdf(file = power_file_path, width = 11, height = 8.5)
        par(mfrow = c(1, 2))
        cex1 <- 0.9
        plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2], xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit, signed R^2", type = "n", main = "Scale independence")
        text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2], labels = powers, cex = cex1, col = "red")
        abline(h = 0.85, col = "red")
        plot(sft$fitIndices[, 1], sft$fitIndices[, 5], xlab = "Soft Threshold (power)", ylab = "Mean Connectivity", type = "n", main = "Mean connectivity")
        text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, cex = cex1, col = "red")
        dev.off()
        print("软阈值筛选结束")
        # 更新下载按钮状态
        print(power_file_path)
        power <- sft$powerEstimate


        if (is.na(power)) {
          power <- ifelse(nSamples < 20, ifelse(type == "unsigned", 9, 18),
            ifelse(nSamples < 30, ifelse(type == "unsigned", 8, 16),
              ifelse(nSamples < 40, ifelse(type == "unsigned", 7, 14),
                ifelse(type == "unsigned", 6, 12)
              )
            )
          )
        }
        ### step 1 end

        ### step 2 begin
        net <- blockwiseModules(dataExpr,
          power = power, maxBlockSize = nGenes,
          TOMType = type, minModuleSize = 30,
          reassignThreshold = 0, mergeCutHeight = 0.25,
          numericLabels = TRUE, pamRespectsDendro = FALSE,
          saveTOMs = F,
          loadTOMs = F,
          # saveTOMFileBase = paste0(exprMat, ".tom"),
          verbose = 3
        )
        table(net$colors)
        ## 绘画结果展示### open a graphics window
        # sizeGrWindow(12, 9)
        # Convert labels to colors for plotting
        mergedColors <- labels2colors(net$colors)
        # Plot the dendrogram and the module colors underneath
        Dendro_file_path <- file.path(result_folder, "DendroAndColor_plot.pdf")
        pdf(file = tree_file_path, width = 11, height = 8.5)
        plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]], "Module colors",
          dendroLabels = FALSE, hang = 0.03,
          addGuide = TRUE, guideHang = 0.05
        )
        dev.off()
        moduleLabels <- net$colors
        moduleColors <- labels2colors(net$colors)
        table(moduleColors)
        x <- table(moduleColors)
        write.csv(x, file.path(result_folder, "Module_Gene_Number.csv"), row.names = FALSE)

        MEs <- net$MEs
        geneTree <- net$dendrograms[[1]]
        MEs <- net$MEs
        MEs_col <- MEs
        colnames(MEs_col) <- paste0("ME", labels2colors(
          as.numeric(str_replace_all(colnames(MEs), "ME", ""))
        ))
        MEs_col <- orderMEs(MEs_col)

        pdf(file = file.path(result_folder, "Eigengene adjacency heatmap.pdf"), width = 11, height = 8.5)
        plotEigengeneNetworks(MEs_col, "Eigengene adjacency heatmap",
          marDendro = c(3, 3, 2, 4),
          marHeatmap = c(3, 4, 2, 2), plotDendrograms = T,
          xLabelsAngle = 90
        )
        dev.off()

        TOM <- TOMsimilarityFromExpr(dataExpr, power = power)
        ## 结果输出

        y <- as.data.frame(table(moduleColors)) # table转化数据框
        names(y) <- c("color", "num") # 数据框重新命名
        color <- as.character(y$color) # 获得分组

        for (i in color) {
          modules <- i
          probes <- names(dataExpr)
          inModule <- is.finite(match(moduleColors, modules))
          modProbes <- probes[inModule]
          # modGenes = annot$gene_symbol[match(modProbes, annot$substanceBXH)];
          # Select the corresponding Topological Overlap
          modTOM <- TOM[inModule, inModule]
          dimnames(modTOM) <- list(modProbes, modProbes)
          # Export the network into edge and node list files Cytoscape can read
          cyt <- exportNetworkToCytoscape(modTOM,
            edgeFile = file.path(result_folder, paste("edges-", paste(modules, collapse = "-"), ".txt", sep = "")),
            nodeFile = file.path(result_folder, paste("nodes-", paste(modules, collapse = "-"), ".txt", sep = "")),
            weighted = TRUE,
            threshold = 0.02,
            nodeNames = modProbes,
            # altNodeNames = modGenes,
            nodeAttr = moduleColors[inModule]
          )
        }


        #######################  模块可视化热图
        person <- cor(dataExpr, use = "p")
        corr <- TOM
        Colors <- mergedColors
        colnames(corr) <- colnames(dataExpr)
        rownames(corr) <- colnames(dataExpr)
        names(Colors) <- colnames(dataExpr)
        colnames(person) <- colnames(dataExpr)
        rownames(person) <- colnames(dataExpr)
        umc <- unique(mergedColors)
        lumc <- length(umc)


        for (i in c(1:lumc))
        {
          if (umc[i] == "grey") {
            next
          }
          ME <- MEs_col[, paste("ME", umc[i], sep = "")]
          pdf(file.path(result_folder, paste0(umc[i], "-expression.pdf")))
          par(mfrow = c(2, 1), mar = c(0.3, 5.5, 3, 2))
          plotMat(t(scale(dataExpr[, Colors == umc[i]])),
            nrgcols = 30,
            rlabels = F, rcols = umc[i], main = umc[i], cex.main = 2
          )
          par(mar = c(5, 4.2, 0, 0.7))

          barplot(ME,
            col = umc[i], main = "", cex.main = 2, ylab = "eigengene expression",
            xlab = "array sample"
          )
          dev.off()


          png(file.path(result_folder, paste0(umc[i], "-expression.png")))
          par(mfrow = c(2, 1), mar = c(0.3, 5.5, 3, 2))
          plotMat(t(scale(dataExpr[, Colors == umc[i]])),
            nrgcols = 30,
            rlabels = F, rcols = umc[i], main = umc[i], cex.main = 2
          )
          par(mar = c(5, 4.2, 0, 0.7))

          barplot(ME,
            col = umc[i], main = "", cex.main = 2, ylab = "eigengene expression",
            xlab = "array sample"
          )
          dev.off()
        }
        ## 性状表型关联

        if (!is.null(data$Trait)) {
          trait <- data$Trait$datapath
          traitData <- read.table(trait,
            sep = ",", header = T, row.names = 1,
            check.names = FALSE, comment = "", quote = ""
          )

          nGenes <- ncol(dataExpr)
          nSamples <- nrow(dataExpr)
          # 重新计算带有颜色标签的模块
          MEs0 <- moduleEigengenes(dataExpr, moduleColors)$eigengenes
          MEs <- orderMEs(MEs0)
          moduleTraitCor <- cor(MEs_col, traitData, use = "p", method = input$corType)
          moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)
          # 通过相关值对每个关联进行颜色编码

          # 展示模块与表型数据的相关系数和 P值
          textMatrix <- paste(signif(moduleTraitCor, 2), "(",
            signif(moduleTraitPvalue, 1), ")",
            sep = ""
          )

          dim(textMatrix) <- dim(moduleTraitCor)
          Module_Trait_path <- file.path(result_folder, "Module Trait relationships.pdf")
          pdf(file = tree_file_path, width = 11, height = 8.5)
          sizeGrWindow(10, 10)
          par(mar = c(3, 10, 3, 3)) ## 下 左 上 右
          labeledHeatmap(
            Matrix = moduleTraitCor,
            xLabels = names(traitData),
            xLabelsAngle = 35,
            xLabelsAdj = 1,
            yLabels = names(MEs_col),
            ySymbols = names(MEs_col),
            colorLabels = FALSE,
            colors = blueWhiteRed(50),
            textMatrix = textMatrix,
            setStdMargins = FALSE,
            cex.text = 1,
            zlim = c(-1, 1),
            main = paste("Module-Trait relationships")
          )
          dev.off()
        }



        ####
        # shinyjs::toggleState("download")
        shinyjs::enable("download")
        shinyjs::enable("run")
      },
      error = function(e) {
        # 捕获并显示错误
        shinyjs::enable("run")
        shiny::showNotification(paste("Error:", e$message), type = "error")
      }
    )
  })

  output$download <- downloadHandler(
    filename = function() {
      paste("WGCNA_Result_", Sys.Date(), ".zip", sep = "")
    },
    content = function(file) {
      zip_file_path <- file.path(temp_folder, "WGCNA_Result.zip")

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
