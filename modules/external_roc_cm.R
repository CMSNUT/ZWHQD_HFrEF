# Z-score
normalize_genes <- function(expr_mat) {
  scaled_mat <- t(scale(t(expr_mat)))
  scaled_mat[is.na(scaled_mat)] <- 0
  return(scaled_mat)
}

# Batch Calculation of ROC and Prediction Probability
calc_roc <- function(model, data, dataset_name){
  pred <- predict(model, data, type = "response")
  roc_obj <- roc(data$group, pred)
  auc_val <- round(auc(roc_obj),3)
  roc_smooth <- smooth(roc_obj, method="binormal")
  
  df <- data.frame(
    Specificity = 1 - roc_smooth$specificities,
    Sensitivity = roc_smooth$sensitivities,
    Group = paste0(dataset_name, " AUC=", auc_val),
    Dataset = dataset_name
  )
  return(list(df=df, auc=auc_val, pred=pred, roc=roc_obj))
}


# Correction curve data calculation
calc_calibration <- function(model, data, dataset_name){
  pred <- predict(model, data, type="response")
  obs <- data$group
  cal <- ResourceSelection::hoslem.test(obs, pred, g=10)
  data.frame(
    pred_prob = pred,
    observed = obs,
    Dataset = dataset_name
  )
}

# Confusion matrix plotting function
plot_cm <- function(model, data, title){
  pred <- ifelse(predict(model, data, type="response")>0.5, 1, 0)
  cm <- table(Actual = data$group, Predicted = pred)
  cm_df <- as.data.frame(cm)
  
  ggplot(cm_df, aes(x=Actual, y=Predicted, fill=Freq)) +
    geom_tile(color="black") +
    geom_text(aes(label=Freq), size=3) +
    scale_fill_gradient(low="white", high="#2E86AB") +
    labs(title=title) +
    theme(legend.position = "none")
}



