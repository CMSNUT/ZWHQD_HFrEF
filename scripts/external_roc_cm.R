plot_external_roc <- function(
    model, 
    external_expr, 
    external_group, 
    dataset_name = "External", 
    plot_p = TRUE,
    save_path = NULL) {
  # prob
  pred_prob <- predict(model, newdata = external_expr, type = "response")
  
  # ROC & AUC
  roc_obj <- roc(external_group, pred_prob, direction = "auto", quiet = TRUE)
  auc_val <- round(auc(roc_obj), 3)
  
  roc_obj_smooth <- smooth(roc_obj, method = "binormal")
  
  roc_df <- data.frame(
    FPR = 1 - roc_obj_smooth$specificities,
    TPR = roc_obj_smooth$sensitivities
  )
  
  # plot
  p <- ggplot(roc_df, aes(x = FPR, y = TPR)) +
    geom_line(color = "red", size = 1.2) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
    labs(x = "1 - Specificity", y = "Sensitivity",
         title = paste0(dataset_name, " (AUC = ", auc_val, ")")) +
    theme_bw() +
    theme(legend.position = "none")
  
  if (plot_p) print(p)
  
  if (!is.null(save_path)) {
    ggsave(save_path, p, width = 5, height = 5)
    cat("ROC plot saved to", save_path, "\n")
  }
  
  invisible(list(roc = roc_obj, auc = auc_val, plot = p))
}


plot_external_cm <- function(
    model, 
    external_expr, 
    external_group,
    dataset_name = "External", 
    threshold = NULL, 
    save_path = NULL
  ) {

  pred_prob <- predict(model, newdata = external_expr, type = "response")
  
  if (is.null(threshold)) {
    roc_obj <- roc(external_group, pred_prob, direction = "auto", quiet = TRUE)
    threshold <- coords(roc_obj, "best", ret = "threshold", best.method = "youden")$threshold
    cat("Optimal threshold (Youden index):", round(threshold, 3), "\n")
  }

  pred_class <- ifelse(pred_prob >= threshold, 1, 0)

  cm <- confusionMatrix(as.factor(pred_class), as.factor(external_group), positive = "1")
  
  cm_table <- as.data.frame(cm$table)
  colnames(cm_table) <- c("Predicted", "Actual", "Count")
  cm_table <- cm_table %>%
    group_by(Actual) %>%
    mutate(Percent = Count / sum(Count) * 100)
  
  # plot
  p <- ggplot(cm_table, aes(x = Actual, y = Predicted, fill = Count)) +
    geom_tile(color = "white") +
    geom_text(aes(label = paste0(Count, "\n(", round(Percent, 1), "%)")), size = 4) +
    scale_fill_gradient(low = "white", high = "steelblue", name = "Count") +
    labs(x = "Actual Class", y = "Predicted Class",
         title = paste0(dataset_name, " Confusion Matrix (threshold = ", round(threshold, 3), ")")) +
    theme_minimal() +
    theme(legend.position = "none")
  
  print(p)
  
  if (!is.null(save_path)) {
    ggsave(save_path, p, width = 5, height = 4)
    cat("Confusion matrix saved to", save_path, "\n")
  }
  
  cat("\n=== Performance on", dataset_name, "===\n")
  metrics <- cm$byClass[c("Sensitivity", "Specificity", "Precision", "F1", "Balanced Accuracy")]
  print(round(metrics, 3))
  
  invisible(list(confusion_matrix = cm, threshold = threshold, plot = p))
}