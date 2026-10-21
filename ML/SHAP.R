rm(list=ls())
setwd("\\SHAP")

library(glmnet)
library(survival)
library(randomForestSRC)
library(DALEX)
library(shapper)
library(ggplot2)
library(dplyr)
library(shapviz)

candidate_genes <- read.csv("MLgene.csv", stringsAsFactors = FALSE)
candidate_genes <- unique(trimws(candidate_genes$Gene))
candidate_genes <- candidate_genes[
  !is.na(candidate_genes) & nzchar(candidate_genes)
]

common_feature <- c("ID", "OS.time", "OS", candidate_genes)
pre_var <- candidate_genes

library(data.table)
train_data<-as.data.frame(load("list_train_vali_Data.RData"))
est_dd <- as.data.frame(list_train_vali_Data$TCGA)[, common_feature[-1]] 
#####LassoRSF#########
message( paste0('--- 10.Lasso + ','RSF ---'))

x1 <- as.matrix(est_dd[, pre_var])
x2 <- as.matrix(Surv(est_dd$OS.time, est_dd$OS))
set.seed(seed)
fit = cv.glmnet(x1, x2,
                nfold = 10, #10-fold cross-validation
                family = "cox", alpha = 1,
                type.measure = "class")
fit$lambda.min

myCoefs <- coef(fit, s = "lambda.min");
rid <- myCoefs@Dimnames[[1]][Matrix::which(myCoefs != 0 )]
rf_nodesize <- 5

est_dd2 <- as.data.frame(list_train_vali_Data$TCGA)[,c('OS.time', 'OS', rid)]

fit <- rfsrc(Surv(OS.time,OS)~., data = est_dd2,
             ntree = 1000, nodesize = rf_nodesize, 
             splitrule = 'logrank',
             importance = T,
             proximity = T,
             forest = T,
             seed = seed)

pred_wrapper <- function(model, newdata) {
  surv_mat <- predict(model, newdata = newdata)$survival
  t0_idx <- which.min(abs(model$time.interest - 1825))
  surv_mat[, t0_idx]
}

X_train <- est_dd2[, - (1:2)]

#background
set.seed(seed)
bg <- X_train[sample(nrow(X_train), 100), ]
#####shap#########
shap_values <- kernelshap(
  object   = fit,
  X        = X_train,
  pred_fun = pred_wrapper,
  bg_X     = bg
)

####fig###
head(shap_values$S)
shap_vis <- shapviz(shap_values, X_train)
feature_importance <- colMeans(abs(shap_values$S))
sorted_features <- names(sort(feature_importance, decreasing = TRUE))
sv_importance(shap_vis, kind="bar", show_numbers=TRUE)
visualization_theme <-theme_minimal(base_size = 26) +
  theme(plot.title =element_text(face ="bold", colour = "black", size =18),
        axis.title =element_text(size =20, colour = "black"))
pdf("SHAP_Feature_Importance_Barplot.pdf", width = 8.2, height = 6)
?sv_importance
sv_importance(shap_vis, kind = "bar", show_numbers = TRUE,number_size = 5.8) +
  visualization_theme +
  labs(
    #title    = "Feature Importance (Mean Absolute SHAP)",
    #subtitle = "Average impact magnitude of each feature on model predictions",
    x        = "Mean |SHAP value|",
    y        = "Feature",
    #caption  = "Bar height indicates feature importance, with value showing mean absolute SHAP"
  )
dev.off()
pdf("SHAP_BeeSwarm_Plot.pdf", width = 8, height = 6)
sv_importance(shap_vis, kind = "bee", show_numbers = FALSE) +
  visualization_theme +
  labs(
    #title    = "SHAP Value Distribution (Bee Swarm)",
    #subtitle = "Each point represents one sample, color indicates feature value",
    x        = "SHAP value (impact on model output)",
    y        = "Feature",
    #caption  = "Red: high feature values | Blue: low feature values | Horizontal spread: direction of effect"
  )
dev.off()
visualization_theme <-theme(plot.title =element_text(face ="plain", size =14),
                            axis.title =element_text(size =12))
pdf("SHAP_Feature_Dependence.pdf", width = 12, height = 12)
sv_dependence(shap_vis, sorted_features[1:11])

dev.off()
?sv_dependence

visualization_theme=theme_minimal(base_size = 18) +theme(plot.title =element_text(face ="bold", colour = "black", size =18),
                                                         axis.title =element_text(size =18, colour = "black"))
train_data<-X_train
?sv_waterfall
dir.create('Sample_Waterfall', showWarnings = FALSE)
for (i in 1:min(10, nrow(train_data))) {
  pdf(paste0("5Sample_Waterfall/sample_", i, "_waterfall.pdf"), width = 9, height = 5)
  p <- sv_waterfall(shap_vis, row_id = i,  show_annotation = TRUE,
                    annotation_size = 5.2,max_display=Inf) +
    visualization_theme
  print(p)
  dev.off()
}
dir.create('Force_Plot', showWarnings = FALSE)
visualization_theme=theme_minimal(base_size = 18) +theme(plot.title =element_text(face ="bold", colour = "black", size =18),
                                                         axis.title =element_text(size =18, colour = "black"))
for (i in 1:min(10, nrow(train_data))) {
  pdf(paste0("Force_Plot/sample_", i, "_force.pdf"), width = 9, height =4)
  p <- sv_force(shap_vis, row_id = i,  bar_label_size = 5.2,
                
                show_annotation = TRUE,
                annotation_size = 5.2,) +visualization_theme
  
  print(p)
  dev.off()
}

