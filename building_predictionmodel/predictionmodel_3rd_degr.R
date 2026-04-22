# Rscript to test various log regression models 
# Call with Rscript /hpc/hers_en/pbeele/simPed/building_predictionmodel/training_famsize.R > predictionmodel.txt 2>&1

install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

packages <- c("dplyr", "glmnet", "DescTools", "ROCR", "rcompanion", "car", "pROC", "boot", "caret",
              "ggplot2", "reshape2", "gridExtra", "grid", "cowplot", "PRROC", "patchwork", "wesanderson", "purrr")
options(repos = c(CRAN = "https://cloud.r-project.org/"))  # Set CRAN mirror

invisible(lapply(packages, install_if_missing))
options(crayon.enabled = FALSE)
options(pillar.max_footer_lines = 999)

## Helper function to easily get all required model outputs:
evaluate_glm_model <- function(model, data, model_name = "GLM Model") {
  # --- 1. Model Fit Statistics ---
  nk_out <- rcompanion::nagelkerke(model)
  nagelkerke_r2 <- nk_out$Pseudo.R.squared.for.model.vs.null["Nagelkerke (Cragg and Uhler)", "Pseudo.R.squared"]

  null_deviance <- model$null.deviance
  residual_deviance <- model$deviance
  aic <- AIC(model)

  # --- 2. Predictions and Classification Metrics ---
  predicted_probabilities <- predict(model, newdata = data, type = "response")
  predicted_classes <- ifelse(predicted_probabilities > 0.5, "mendel", "polygen")
  actual_classes <- ifelse(data$mendel_ALS_Y == 1, "mendel", "polygen")
  conf_matrix <- table(Predicted = predicted_classes, Actual = actual_classes)

  # Robust extraction of confusion matrix components
  TP <- ifelse("mendel" %in% rownames(conf_matrix) && "mendel" %in% colnames(conf_matrix), conf_matrix["mendel", "mendel"], 0)
  TN <- ifelse("polygen" %in% rownames(conf_matrix) && "polygen" %in% colnames(conf_matrix), conf_matrix["polygen", "polygen"], 0)
  FP <- ifelse("mendel" %in% rownames(conf_matrix) && "polygen" %in% colnames(conf_matrix), conf_matrix["mendel", "polygen"], 0)
  FN <- ifelse("polygen" %in% rownames(conf_matrix) && "mendel" %in% colnames(conf_matrix), conf_matrix["polygen", "mendel"], 0)

  sensitivity <- ifelse((TP + FN) > 0, TP / (TP + FN), NA)
  specificity <- ifelse((TN + FP) > 0, TN / (TN + FP), NA)
  accuracy <- sum(predicted_classes == actual_classes) / length(predicted_classes)
  PPV <- ifelse((TP + FP) > 0, TP / (TP + FP), NA)
  NPV <- ifelse((TN + FN) > 0, TN / (TN + FN), NA)

  pred <- ROCR::prediction(predicted_probabilities, data$mendel_ALS_Y)
  perf <- ROCR::performance(pred, measure = "auc")
  auc <- perf@y.values[[1]]

  # --- 3. Coefficients, Odds Ratios, CIs, and P-values ---
  summary_model <- summary(model)
  coefs <- coef(model)
  odds_ratios <- exp(coefs)
  se <- coef(summary_model)[, "Std. Error"]
  lower_ci <- exp(coefs - 1.96 * se)
  upper_ci <- exp(coefs + 1.96 * se)
  p_values <- coef(summary_model)[, "Pr(>|z|)"]

  coef_table <- data.frame(
    Coefficient = coefs,
    Odds_Ratio = odds_ratios,
    OR_2.5_CI = lower_ci,
    OR_97.5_CI = upper_ci,
    P_value = p_values
  )

  # --- 4. Multicollinearity Diagnostics (VIF) ---
  # predictor_names <- names(coef(model))[-1]
  # if (exists("vif", mode = "function") && length(predictor_names) > 1) {
  #   vif_values <- vif(model)
  #   vif_table <- as.data.frame(vif_values)
  #   colnames(vif_table) <- "VIF"
  #   vif_output <- vif_table
  #   vif_message <- NULL
  # } else if (length(predictor_names) > 1) {
  #   vif_output <- NULL
  #   vif_message <- "VIF not calculated: 'vif' function not found. Load the 'car' package if you want VIFs."
  # } else {
  #   vif_output <- NULL
  #   vif_message <- "VIF not applicable: only one predictor in the model."
  # }

  # # --- 5. Correlation Matrix for Predictors ---
  # if (length(predictor_names) > 1) {
  #   predictors_df <- model.matrix(model)[, -1, drop = FALSE] # remove intercept
  #   cor_matrix <- cor(predictors_df)
  # } else {
  #   cor_matrix <- NULL
  # }

  # --- 6. Print Results ---
  cat(paste0("\n--- ", model_name, " ---\n"))
  cat("\n--- Model Fit ---\n")
  cat(paste("Null deviance:", round(null_deviance, 3), "\n"))
  cat(paste("Residual deviance:", round(residual_deviance, 3), "\n"))
  cat(paste("Nagelkerke R-squared:", round(nagelkerke_r2, 3), "\n"))
  cat(paste("AIC:", round(aic, 3), "\n"))

  cat("\n--- Classification Metrics ---\n")
  cat(paste("Sensitivity:", round(sensitivity, 3), "\n"))
  cat(paste("Specificity:", round(specificity, 3), "\n"))
  cat(paste("PPV:", round(PPV, 3), "\n"))
  cat(paste("NPV:", round(NPV, 3), "\n"))
  cat(paste("Accuracy:", round(accuracy, 3), "\n"))
  cat(paste("AUC:", round(auc, 3), "\n"))

  cat("\n--- Confusion Matrix (Predicted as rows, Actual as columns) ---\n")
  # Custom print for annotated confusion matrix
  #           Actual
  # Predicted   mendel   polygen
  # mendel      TP       FP     | PPV
  # polygen     FN       TN     | NPV
  #             |Sens    |Spec

  # Extract values for display (ensure all are present)
  TP_disp <- TP
  FP_disp <- FP
  TN_disp <- TN
  FN_disp <- FN

  # Print header
  cat(sprintf("%-12s %-10s %-10s\n", "", "mendel", "polygen"))
  cat(sprintf("%-12s %-10s %-10s %s\n", 
              "mendel", 
              format(TP_disp, justify="right"), 
              format(FP_disp, justify="right"), 
              paste0("| PPV: ", ifelse(is.na(PPV), "NA", round(PPV, 2)))))
  cat(sprintf("%-12s %-10s %-10s %s\n", 
              "polygen", 
              format(FN_disp, justify="right"), 
              format(TN_disp, justify="right"), 
              paste0("| NPV: ", ifelse(is.na(NPV), "NA", round(NPV, 2)))))
  cat(sprintf("%-12s %-10s %-10s\n", 
              "", 
              paste0("Sens: ", ifelse(is.na(sensitivity), "NA", round(sensitivity, 2))), 
              paste0("Spec: ", ifelse(is.na(specificity), "NA", round(specificity, 2)))))
  cat("\n")

  cat("\n--- Coefficients, Odds Ratios, 95% CI, and P-values ---\n")
  print(coef_table)

  # cat("\n--- Multicollinearity Diagnostics (VIF) ---\n")
  # if (!is.null(vif_output)) {
  #   print(vif_output)
  # } else {
  #   cat(vif_message, "\n")
  # }

#  cat("\n--- Correlation Matrix for Predictors ---\n")
#  if (!is.null(cor_matrix)) {
#    print(round(cor_matrix, 3))
#  } else {
#    cat("Correlation matrix not applicable: only one predictor in the model.\n")
#  }
#
  cat("\n-----------------------------------------\n")

  # --- 7. Return Results (Optional) ---
  results <- list(
    Null_Deviance = null_deviance,
    Residual_Deviance = residual_deviance,
    Nagelkerke_R2 = nagelkerke_r2,
    AIC = aic,
    Sensitivity = sensitivity,
    Specificity = specificity,
    PPV = PPV,
    NPV = NPV,
    Accuracy = accuracy,
    AUC = auc,
    TP = TP,
    FP = FP,
    TN = TN,
    FN = FN,
    Confusion_Matrix = conf_matrix,
    Coefficient_Table = coef_table
    # VIF_Table = vif_output,
    # VIF_Message = vif_message,
    # Correlation_Matrix = cor_matrix
  )
  return(invisible(results))
}

#load file from last full simulation 
results_df <- read.csv("../results/trainingset/combined_simulations.csv")

## Step 0: preparation of data, so they fit nicely in the model
## Calculated unaffected individuals per degree, per affected disease
results_df$relatives_1st_unaffected <- results_df$relatives_1st - results_df$relatives_1st_als
results_df$relatives_2nd_unaffected <- results_df$relatives_2nd - results_df$relatives_2nd_als
results_df$relatives_3rd_unaffected <- results_df$relatives_3rd - results_df$relatives_3rd_als

results_df$relatives_1st_ALSFTD_unaffected <- results_df$relatives_1st - results_df$relatives_1st_als - results_df$relatives_1st_ftd_unique
results_df$relatives_2nd_ALSFTD_unaffected <- results_df$relatives_2nd - results_df$relatives_2nd_als - results_df$relatives_2nd_ftd_unique
results_df$relatives_3rd_ALSFTD_unaffected <- results_df$relatives_3rd - results_df$relatives_3rd_als - results_df$relatives_3rd_ftd_unique

results_df$relatives_1st_ALSdem_unaffected <- results_df$relatives_1st - results_df$relatives_1st_als - results_df$relatives_1st_dementia_unique
results_df$relatives_2nd_ALSdem_unaffected <- results_df$relatives_2nd - results_df$relatives_2nd_als - results_df$relatives_2nd_dementia_unique
results_df$relatives_3rd_ALSdem_unaffected <- results_df$relatives_3rd - results_df$relatives_3rd_als - results_df$relatives_3rd_dementia_unique

# Step 1: demonstrate effect of looking beyond 1st degree, to 2nd and 3rd degree

model_1st_degree <- glm(mendel_ALS_Y ~ relatives_1st_als,
              data = results_df,
              family = binomial(link = "logit"))
evaluate_glm_model(model_1st_degree, results_df, model_name = "Model 1st degree relatives with ALS")

model_2nd_degree <- glm(mendel_ALS_Y ~  relatives_1st_als + relatives_2nd_als,
                data = results_df,
                family = binomial(link = "logit"))
evaluate_glm_model(model_2nd_degree, results_df, model_name = "Model 1st or 2nd degree relatives with ALS")

model_3rd_degree <- glm(mendel_ALS_Y ~  relatives_1st_als + relatives_2nd_als + relatives_3rd_als,
                data = results_df,
                family = binomial(link = "logit"))
evaluate_glm_model(model_3rd_degree, results_df, model_name = "Model 1st, 2nd, or 3rd degree relatives with ALS")

# Step 2: demonstrate effect of looking at FTD or other dementias, beyond only ALS affected relatives

model_ALS <- glm(mendel_ALS_Y ~  relatives_1st_als + relatives_2nd_als + relatives_3rd_als,
                data = results_df,
                family = binomial(link = "logit"))
evaluate_glm_model(model_ALS, results_df, model_name = "Model 1st, 2nd, or 3rd degree relatives with ALS")

model_FTD <- glm(mendel_ALS_Y ~  relatives_1st_als + relatives_2nd_als + relatives_3rd_als +
                relatives_1st_ftd_unique + relatives_2nd_ftd_unique + relatives_3rd_ftd_unique,
                data = results_df,
                family = binomial(link = "logit"))
evaluate_glm_model(model_FTD, results_df, model_name = "Model 1st, 2nd, or 3rd degree relatives with ALS or FTD")

model_dementia <- glm(mendel_ALS_Y ~  relatives_1st_als + relatives_2nd_als + relatives_3rd_als +
                relatives_1st_dementia_unique + relatives_2nd_dementia_unique + relatives_3rd_dementia_unique,
                data = results_df,
                family = binomial(link = "logit"))
evaluate_glm_model(model_dementia, results_df, model_name = "Model 1st, 2nd, or 3rd degree relatives with ALS or dementia")

# Calculate ROC data for detailed analysis
# Create ROC objects for your models
roc_ftd <- roc(results_df$mendel_ALS_Y, predict(model_FTD, type="response"))
roc_dementia <- roc(results_df$mendel_ALS_Y, predict(model_dementia, type="response"))

# Verify AUC using different approaches
auc(roc_ftd)
auc(roc_dementia)

# Step 3: compare effect of just looking at affecteds to looking at affecteds + unaffecteds

model_ALS_unaffected <- glm(
  mendel_ALS_Y ~ 
    relatives_1st_als + relatives_2nd_als + relatives_3rd_als + 
    relatives_1st_unaffected + relatives_2nd_unaffected + relatives_3rd_unaffected,
  data = results_df,
  family = binomial(link = "logit")  
)
evaluate_glm_model(model_ALS_unaffected, results_df, model_name = "Model looking at 1st, 2nd, or 3rd degree ALS-affected and unaffected relatives")

model_FTD_unaffected <- glm(
  mendel_ALS_Y ~ 
    relatives_1st_als + relatives_2nd_als + relatives_3rd_als + 
    relatives_1st_ftd_unique + relatives_2nd_ftd_unique + relatives_3rd_ftd_unique +
    relatives_1st_ALSFTD_unaffected + relatives_2nd_ALSFTD_unaffected + relatives_3rd_ALSFTD_unaffected,
  data = results_df,
  family = binomial(link = "logit")  
)
evaluate_glm_model(model_FTD_unaffected, results_df, model_name = "Model looking at 1st, 2nd, or 3rd degree ALS or FTD-affected and unaffected relatives")

model_dementia_unaffected <- glm(
  mendel_ALS_Y ~ 
    relatives_1st_als + relatives_2nd_als + relatives_3rd_als + 
    relatives_1st_dementia_unique + relatives_2nd_dementia_unique + relatives_3rd_dementia_unique +
    relatives_1st_ALSdem_unaffected + relatives_2nd_ALSdem_unaffected + relatives_3rd_ALSdem_unaffected,
  data = results_df,
  family = binomial(link = "logit")  
)
evaluate_glm_model(model_dementia_unaffected, results_df, model_name = "Model looking at 1st, 2nd, or 3rd degree ALS or dementia-affected and unaffected relatives")


# interim evaluation through plotted roc curves
# Step 1: ALS, 1st, 1st+2nd, 1st+2nd+3rd degree
roc_list1 <- list(
  "ALS 1st deg" = roc(results_df$mendel_ALS_Y, predict(model_1st_degree, type="response")),
  "ALS 1st+2nd deg" = roc(results_df$mendel_ALS_Y, predict(model_2nd_degree, type="response")),
  "ALS 1st+2nd+3rd deg" = roc(results_df$mendel_ALS_Y, predict(model_3rd_degree, type="response"))
)

# Step 2: Adding FTD and dementia
roc_list2 <- list(
  "ALS only" = roc(results_df$mendel_ALS_Y, predict(model_ALS, type="response")),
  "ALS+FTD" = roc(results_df$mendel_ALS_Y, predict(model_FTD, type="response")),
  "ALS+dementia" = roc(results_df$mendel_ALS_Y, predict(model_dementia, type="response"))
)

# Step 3: Adding unaffected individuals
roc_list3 <- list(
  "ALS only" = roc(results_df$mendel_ALS_Y, predict(model_ALS, type="response")),
  "ALS+unaffected" = roc(results_df$mendel_ALS_Y, predict(model_ALS_unaffected, type="response")),
  "ALS+FTD+unaffected" = roc(results_df$mendel_ALS_Y, predict(model_FTD_unaffected, type="response")),
  "ALS+dementia+unaffected" = roc(results_df$mendel_ALS_Y, predict(model_dementia_unaffected, type="response"))
)

# select colors
darjeeling_cols <- wesanderson::wes_palette("Darjeeling1", 4, type = "continuous")

plot_roc_curves <- function(roc_list, title = "ROC Curves", subtitle = NULL, model_labels = NULL) {
  plot_df <- do.call(rbind, lapply(names(roc_list), function(n) {
    data.frame(
      model = n,
      specificity = rev(roc_list[[n]]$specificities),
      sensitivity = rev(roc_list[[n]]$sensitivities),
      auc = round(auc(roc_list[[n]]), 3)
    )
  }))
  plot_df$model <- factor(plot_df$model, levels = names(roc_list))
  
  # Create AUC labels
  auc_df <- plot_df %>%
    group_by(model) %>%
    summarise(auc = unique(auc), .groups = "drop")
  
  if (!is.null(model_labels)) {
    auc_df$display_name <- ifelse(names(model_labels) %in% auc_df$model,
                                 model_labels[match(auc_df$model, names(model_labels))],
                                 as.character(auc_df$model))
  } else {
    auc_df$display_name <- as.character(auc_df$model)
  }
  auc_df$label <- paste0(auc_df$display_name, "\n(AUC=", auc_df$auc, ")")
  
  p <- ggplot(plot_df, aes(x = specificity, y = sensitivity, color = model)) +
    geom_line(size = 1.1) +
    geom_abline(slope = 1, intercept = 1, linetype = "dashed", color = "black", linewidth = 0.8) +
    scale_x_reverse(limits = c(1, 0), labels = scales::percent_format(accuracy = 1)) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
    labs(title = title, subtitle = subtitle, 
         x = "Specificity", 
         y = "Sensitivity") +
    scale_color_manual(
      values = darjeeling_cols[1:length(unique(plot_df$model))],
      labels = auc_df$label,
      guide = guide_legend(nrow = length(unique(plot_df$model)), byrow = TRUE)
    ) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      plot.title = element_text(size = 10, hjust = 0.5, margin = margin(b = 2), face = "bold"),
      plot.subtitle = element_text(size = 8, hjust = 0.5, margin = margin(b = 6)),
      axis.title = element_text(size = 8),
      legend.title = element_blank(),
      legend.text = element_text(size = 7, margin = margin(t = 1, b = 1), lineheight = 1),
      legend.key.height = unit(0.35, "cm"),
      legend.spacing.y = unit(0.15, "cm"),
      axis.text.x = element_text(color = "black", size = 8),
      axis.text.y = element_text(color = "black", size = 8),
      axis.ticks = element_line(color = "black"),
      axis.title.x = element_text(size = 8),
      axis.title.y = element_text(size = 8),
      plot.margin = unit(c(5, 5, 5, 5), "pt")
    ) +
    coord_fixed(ratio = 1)
  
  return(p)
}


# Create output directories
if (!dir.exists("roccurves_trainingset")) dir.create("roccurves_trainingset")

# Step 1: Effect of more distant relatives
roc_plot1 <- plot_roc_curves(
  roc_list1, 
  title = "(A) Monogenic ALS prediction model",
  subtitle = "Accuracy including more distant relatives",
  model_labels = c(
    "ALS 1st deg" = "ALS within 1st° relatives", 
    "ALS 1st+2nd deg" = "ALS within 2nd° relatives", 
    "ALS 1st+2nd+3rd deg" = "ALS within 3rd° relatives"
  )
) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE))

ggsave("roccurves_trainingset/roc_step1_ALS.pdf", roc_plot1, width = 4.5, height = 4, units = "in", dpi = 300)

# Step 2: Effect of FTD and Dementia
roc_plot2 <- plot_roc_curves(
  roc_list2, 
  title = "(B) Monogenic ALS prediction model among 3rd° relatives",
  subtitle = "Accuracy also including FTD and Dementia",
  model_labels = c(
    "ALS only" = "ALS-affected relatives",
    "ALS+FTD" = "ALS/FTD-affected relatives",
    "ALS+dementia" = "ALS/Dementia-affected relatives"
  )
) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE))

ggsave("roccurves_trainingset/roc_step2_FTD_dementia.pdf", roc_plot2, width = 4.5, height = 4, units = "in", dpi = 300)


# Step 3: Effect of unaffected individuals
roc_plot3 <- plot_roc_curves(
  roc_list3, 
  title = "(C) Monogenic ALS prediction model among 3rd° relatives",
  subtitle = "Accuracy also including unaffected individuals",
  model_labels = c(
    "ALS only" = "ALS-affected relatives",
    "ALS+unaffected" = "ALS-affected + unaffected relatives",
    "ALS+FTD+unaffected" = "ALS/FTD-affected + unaffected relatives",
    "ALS+dementia+unaffected" = "ALS/Dementia-affected + unaffected relatives"
  )
) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE))

ggsave("roccurves_trainingset/roc_step3_unaffected.pdf", roc_plot3, width = 4.5, height = 4, units = "in", dpi = 300)

cat("ROC plots saved as 4x4 PDFs with Darjeeling1 Wes Anderson palette\n")

# Combine into one figure
combined_roc_plot <- wrap_plots(roc_plot1, roc_plot2, roc_plot3, ncol = 3) +
  plot_layout(widths = c(0.8, 0.8, 1)) +
  plot_annotation(
    title = "Supplementary Figure 6: ROC curves for monogenic ALS prediction",
    subtitle = "Stepwise effects of including more distant relatives, FTD/dementia, and unaffected individuals",
    theme = theme(
      plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5)
    )
  )

ggsave(
  "roccurves_trainingset/roc_combined.pdf",
  combined_roc_plot,
  width = 14,
  height = 5,
  units = "in",
  dpi = 300
)

cat("ROC plots saved as one combined PDF with Darjeeling1 Wes Anderson palette\n")

## Part 2: crossvalidation of the trainingset

## 5) Final set-up for the ALS, ALS+FTD and ALS+dementia models, including crossvalidation 
# Create helper function
evaluate_cv_glm_model <- function(caret_model, positive_class = "Mendelian", threshold = 0.5, model_name = NULL) {
  # Print model title if provided
  if (!is.null(model_name)) {
    cat("\n==============================\n")
    cat("Model:", model_name, "\n")
    cat("==============================\n\n")
  }  # Extract final glm model
  final_glm <- caret_model$finalModel
  
  # Calculate Nagelkerke R-squared
  n <- length(final_glm$y)
  LL_null <- -final_glm$null.deviance / 2
  LL_full <- -final_glm$deviance / 2
  nagelkerke_r2 <- (1 - exp((2 / n) * (LL_null - LL_full))) / (1 - exp((2 / n) * LL_null))  

  # Extract cross-validated predictions (ensure savePredictions = "final" in trainControl)
  preds <- caret_model$pred
  
  # Filter predictions to best tuning parameter if needed
  if (ncol(caret_model$bestTune) > 0) {
    for (param in names(caret_model$bestTune)) {
      preds <- preds[preds[[param]] == caret_model$bestTune[[param]], ]
    }
  }
  
  preds <- na.omit(preds)
  
  # Define consistent factor levels order: negative first, positive second
  levels_order <- c("Mendelian", "Polygenic")

  # True labels and predicted probabilities for positive class, with explicit levels
  true_labels <- factor(preds$obs, levels = levels_order)
  pred_probs <- preds[[positive_class]]

  # Predicted classes using threshold, with the same levels order
  pred_classes <- factor(ifelse(pred_probs >= threshold, positive_class,
                             setdiff(levels_order, positive_class)),
                       levels = levels_order)

  # Create confusion matrix with consistent factor levels
  cm <- table(Predicted = pred_classes, Actual = true_labels)
  
  TP <- cm[positive_class, positive_class]
  FP <- sum(cm[positive_class, ]) - TP
  TN <- cm[setdiff(levels(true_labels), positive_class), setdiff(levels(true_labels), positive_class)]
  FN <- sum(cm[, positive_class]) - TP
  
  # Calculate metrics
  PPV <- TP / (TP + FP)
  NPV <- TN / (TN + FN)
  Sens <- TP / (TP + FN)
  Spec <- TN / (TN + FP)
  Acc <- sum(diag(cm)) / sum(cm)
  
  # Ensure positive class is the second level
  true_labels <- factor(preds$obs, levels = c(setdiff(levels(preds$obs), positive_class), positive_class))
  pred_probs <- preds[[positive_class]]

  pred_obj <- prediction(pred_probs, true_labels)
  auc <- performance(pred_obj, "auc")@y.values[[1]]
  # Correct inverted AUC if necessary
  if (auc < 0.5) {
    auc <- 1 - auc
  }
  
  # Print deviance and AIC from final glm
  cat("\nModel Fit Statistics:\n")
  cat(sprintf("Null Deviance: %.3f\n", final_glm$null.deviance))
  cat(sprintf("Residual Deviance: %.3f\n", final_glm$deviance))
  cat(sprintf("Nagelkerke R-squared: %.3f\n\n", nagelkerke_r2))
  cat(sprintf("AIC: %.3f\n\n", final_glm$aic))
  
  # Print metrics
  cat("Performance Metrics (threshold =", threshold, "):\n")
  cat(sprintf("Sensitivity: %.3f\n", Sens))
  cat(sprintf("Specificity: %.3f\n", Spec))
  cat(sprintf("PPV: %.3f\n", PPV))
  cat(sprintf("NPV: %.3f\n", NPV))
  cat(sprintf("Accuracy: %.3f\n", Acc))
  cat(sprintf("AUC: %.3f\n\n", auc))
  
  # --- Print simple confusion matrix ---
  cat("--- Confusion Matrix (Predicted as rows, Actual as columns) ---\n")
  cat(sprintf("%-12s", ""))
  for (col_label in levels_order) {
    cat(sprintf("%-10s", col_label))
  }
  cat("\n")

  for (row_label in levels_order) {
    cat(sprintf("%-12s", row_label))
    for (col_label in levels_order) {
      cat(sprintf("%-10d", cm[row_label, col_label]))
    }
    if (row_label == positive_class) {
      cat(sprintf(" | PPV: %.3f", PPV))
    } else {
      cat(sprintf(" | NPV: %.3f", NPV))
    }
    cat("\n")
  }

  cat(sprintf("%-12s %-10s %-10s\n", "", paste0("Sens: ", round(Sens, 3)), paste0("Spec: ", round(Spec, 3))))
  cat("\n")
  
  # Extract coefficients, OR, 95% CI, p-values
  coefs <- summary(final_glm)$coefficients
  OR <- exp(coefs[, "Estimate"])
  CI_lower <- exp(coefs[, "Estimate"] - 1.96 * coefs[, "Std. Error"])
  CI_upper <- exp(coefs[, "Estimate"] + 1.96 * coefs[, "Std. Error"])
  p_vals <- coefs[, "Pr(>|z|)"]
  
  coef_df <- data.frame(
    Term = rownames(coefs),
    Estimate = coefs[, "Estimate"],
    OR = OR,
    CI_lower = CI_lower,
    CI_upper = CI_upper,
    p_value = p_vals,
    stringsAsFactors = FALSE
  )
  
  cat("Model Coefficients:\n")
  print(coef_df, row.names = FALSE)
  
  invisible(list(
    confusion_matrix = cm,
    performance = list(PPV = PPV, NPV = NPV, Sensitivity = Sens, Specificity = Spec, Accuracy = Acc, AUC = auc),
    coefficients = coef_df
  ))
}

# Set seed for reproducibility
set.seed(123)

# Define cross-validation method: 10-fold CV
train_control <- trainControl(method = "cv", number = 10, classProbs = TRUE, summaryFunction =twoClassSummary, savePredictions = "final")
# Convert outcome to factor with meaningful levels
results_df$mendel_ALS_Y <- factor(results_df$mendel_ALS_Y, levels = c(0, 1), labels = c("Polygenic", "Mendelian"))

# Fit logistic regression model with cross-validation
model_cv_ALS <- train(
  mendel_ALS_Y ~ 
    relatives_1st_als + relatives_2nd_als + relatives_3rd_als + 
    relatives_1st_unaffected + relatives_2nd_unaffected + relatives_3rd_unaffected,
  data = results_df,
  method = "glm",
  family = binomial(link = "logit"),
  trControl = train_control,
  metric = "ROC"
)

model_cv_FTD <- train(
  mendel_ALS_Y ~ 
    relatives_1st_als + relatives_2nd_als + relatives_3rd_als + 
    relatives_1st_ftd_unique + relatives_2nd_ftd_unique + relatives_3rd_ftd_unique +
    relatives_1st_ALSFTD_unaffected + relatives_2nd_ALSFTD_unaffected + relatives_3rd_ALSFTD_unaffected,
  data = results_df,
  method = "glm",
  family = binomial(link = "logit"),
  trControl = train_control,
  metric = "ROC"
)

model_cv_dementia <- train(
  mendel_ALS_Y ~ 
    relatives_1st_als + relatives_2nd_als + relatives_3rd_als + 
    relatives_1st_dementia_unique + relatives_2nd_dementia_unique + relatives_3rd_dementia_unique +
    relatives_1st_ALSdem_unaffected + relatives_2nd_ALSdem_unaffected + relatives_3rd_ALSdem_unaffected,
  data = results_df,
  method = "glm",
  family = binomial(link = "logit"),
  trControl = train_control,
  metric = "ROC"
)

evaluate_cv_glm_model(model_cv_ALS, model_name = "ALS model including unaffected individuals")
evaluate_cv_glm_model(model_cv_FTD, model_name = "ALS + FTD model including unaffected individuals")
evaluate_cv_glm_model(model_cv_dementia, model_name = "ALS + Dementia model including unaffected individuals")

# ## PART 3; testing the cross-validated trained models on newly simulated testing data

# # Load new data
new_data <- read.csv("../results/testset/combined_simulations.csv")

new_data$relatives_1st_unaffected <- new_data$relatives_1st - new_data$relatives_1st_als
new_data$relatives_2nd_unaffected <- new_data$relatives_2nd - new_data$relatives_2nd_als
new_data$relatives_3rd_unaffected <- new_data$relatives_3rd - new_data$relatives_3rd_als

new_data$relatives_1st_ALSFTD_unaffected <- new_data$relatives_1st - new_data$relatives_1st_als - new_data$relatives_1st_ftd_unique
new_data$relatives_2nd_ALSFTD_unaffected <- new_data$relatives_2nd - new_data$relatives_2nd_als - new_data$relatives_2nd_ftd_unique
new_data$relatives_3rd_ALSFTD_unaffected <- new_data$relatives_3rd - new_data$relatives_3rd_als - new_data$relatives_3rd_ftd_unique

new_data$relatives_1st_ALSdem_unaffected <- new_data$relatives_1st - new_data$relatives_1st_als - new_data$relatives_1st_dementia_unique
new_data$relatives_2nd_ALSdem_unaffected <- new_data$relatives_2nd - new_data$relatives_2nd_als - new_data$relatives_2nd_dementia_unique
new_data$relatives_3rd_ALSdem_unaffected <- new_data$relatives_3rd - new_data$relatives_3rd_als - new_data$relatives_3rd_dementia_unique

# Ensure outcome column is a factor with correct levels (matching training)
new_data$mendel_ALS_Y <- factor(new_data$mendel_ALS_Y, levels = c(0, 1), labels = c("Polygenic", "Mendelian"))

# List of cross-validated caret models
model_list <- list(
  ALS = model_cv_ALS,
  FTD = model_cv_FTD,
  Dementia = model_cv_dementia
)

model_labels <- c(
  ALS = "ALS-affected + unaffected relatives",
  FTD = "ALS/FTD-affected + unaffected relatives",
  Dementia = "ALS/Dementia-affected + unaffected relatives"
)

# Function to get metrics from caret model on new data
get_caret_model_metrics <- function(caret_model, new_data, outcome_col = "mendel_ALS_Y", 
                                    positive_class = "Mendelian", threshold = 0.5, plot_roc = FALSE) {
  # Ensure outcome is factor with correct levels: negative first, positive second
  new_data[[outcome_col]] <- factor(new_data[[outcome_col]], 
                                   levels = c(setdiff(levels(new_data[[outcome_col]]), positive_class), positive_class))
  
  # Predicted probabilities for positive class
  probs <- predict(caret_model, newdata = new_data, type = "prob")[, positive_class]
  
  # Predicted classes using threshold
  preds <- ifelse(probs > threshold, positive_class, setdiff(levels(new_data[[outcome_col]]), positive_class))
  preds <- factor(preds, levels = levels(new_data[[outcome_col]]))
  
  truth <- new_data[[outcome_col]]
  
  # Confusion matrix and metrics
  tab <- table(Predicted = preds, Actual = truth)
  
  TP <- ifelse(positive_class %in% rownames(tab) && positive_class %in% colnames(tab), tab[positive_class, positive_class], 0)
  TN <- ifelse(setdiff(levels(truth), positive_class) %in% rownames(tab) && setdiff(levels(truth), positive_class) %in% colnames(tab), 
               tab[setdiff(levels(truth), positive_class), setdiff(levels(truth), positive_class)], 0)
  FP <- ifelse(positive_class %in% rownames(tab) && setdiff(levels(truth), positive_class) %in% colnames(tab), 
               tab[positive_class, setdiff(levels(truth), positive_class)], 0)
  FN <- ifelse(setdiff(levels(truth), positive_class) %in% rownames(tab) && positive_class %in% colnames(tab), 
               tab[setdiff(levels(truth), positive_class), positive_class], 0)
  
  sens <- ifelse((TP + FN) > 0, TP / (TP + FN), NA)
  spec <- ifelse((TN + FP) > 0, TN / (TN + FP), NA)
  ppv <- ifelse((TP + FP) > 0, TP / (TP + FP), NA)
  npv <- ifelse((TN + FN) > 0, TN / (TN + FN), NA)
  acc <- mean(preds == truth)
  
  # Calculate AUC and plot ROC if requested
  auc <- NA
  if (requireNamespace("pROC", quietly = TRUE)) {
    # Explicitly specify levels to ensure correct positive class
    roc_obj <- pROC::roc(response = truth, predictor = probs, levels = levels(truth), direction = "<")
    auc <- as.numeric(pROC::auc(roc_obj))
    if (plot_roc) {
      plot(roc_obj, main = "ROC Curve", col = "#2C3E50", lwd = 2)
      abline(a = 0, b = 1, lty = 2, col = "gray")
      legend("bottomright", legend = sprintf("AUC = %.2f", auc), bty = "n")
    }
  }
  
  # Calculate AUPRC (PR-AUC)
  auprc <- NA
  if (requireNamespace("PRROC", quietly = TRUE)) {
    truth_binary <- as.numeric(truth == positive_class)
    pr <- PRROC::pr.curve(scores.class0 = probs[truth_binary == 1],
                          scores.class1 = probs[truth_binary == 0],
                          curve = FALSE)
    auprc <- as.numeric(pr$auc.integral)
  }
  
  list(
    accuracy = acc,
    sensitivity = sens,
    specificity = spec,
    ppv = ppv,
    npv = npv,
    auc = auc,
    auprc = auprc
  )
}

# Evaluate all models on new data (assuming model_list and new_data exist)
thresholds <- seq(0.10, 0.85, by = 0.025)

metrics_list <- lapply(model_list, function(mod) {
  # For each model, evaluate at all thresholds
  do.call(rbind, lapply(thresholds, function(thresh) {
    metrics <- get_caret_model_metrics(mod, new_data, positive_class = "Mendelian", threshold = thresh)
    metrics$threshold <- thresh
    return(as.data.frame(metrics)[1, ])  # Ensure data.frame row format
  }))
})

metrics_df <- do.call(rbind, metrics_list)

metrics_df$model <- rep(model_labels[names(model_list)], each = length(thresholds))
rownames(metrics_df) <- NULL
metrics_df <- metrics_df[, c("model", "threshold", setdiff(names(metrics_df), c("model", "threshold")))]

metrics_df <- metrics_df %>%
  mutate(F1 = 2 * (ppv * sensitivity) / (ppv + sensitivity))
metrics_df <- metrics_df %>%
  mutate(
    sensitivity = round(sensitivity, 3),
    ppv         = round(ppv, 3),
    auc         = round(auc, 3),
    auprc       = round(auprc, 3),
    F1          = round(2 * (ppv * sensitivity) / (ppv + sensitivity), 3),
    youdens_j = round(sensitivity + specificity - 1, 3)
  )

print(metrics_df)

best_f1 <- metrics_df %>%
  group_by(model) %>%
  filter(F1 == max(F1, na.rm = TRUE))
print("Best f1 statistic values:")
print(best_f1)

# Best AUPRC per model
best_auprc <- metrics_df %>%
  group_by(model) %>%
  filter(auprc == max(auprc, na.rm = TRUE))
print("Best AUPRC values:")
print(best_auprc)

# Best Youden's J per model
best_youdens_j <- metrics_df %>%
  group_by(model) %>%
  filter(youdens_j == max(youdens_j, na.rm = TRUE))
print("Best youdens values:")
print(best_youdens_j)



library(ggplot2)
f1_plot <- ggplot(metrics_df, aes(x = threshold, y = F1, color = model)) +
  geom_line(size = 1.1) +
  scale_color_manual(values = darjeeling_cols[1:length(unique(metrics_df$model))]) +
  labs(title = "Optimal F1 Score Threshold", x = "Threshold", y = "F1 Score") +
  theme_bw(base_size = 10) +
  theme(
    plot.title = element_text(size = 10, hjust = 0.5),
    axis.title = element_text(size = 8),
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = 8),
    axis.text = element_text(color = "black", size = 8)
  )
ggsave("optimal_F1_model.pdf", f1_plot, width = 6, height = 4, units = "in", dpi = 300)


# Compute ROC objects for each model on the test set
roc_list <- lapply(model_list, function(mod) {
  probs <- predict(mod, newdata = new_data, type = "prob")[, "Mendelian"]
  truth <- new_data$mendel_ALS_Y
  pROC::roc(response = truth, predictor = probs, levels = c("Polygenic", "Mendelian"), direction = "<", quiet = TRUE)
})
names(roc_list) <- names(model_list)

# Updated ROC plotting function (matching trainingset style)
plot_roc_curves <- function(roc_list, title = "ROC Curves", subtitle = NULL, model_labels = NULL) {
  plot_df <- do.call(rbind, lapply(names(roc_list), function(n) {
    data.frame(
      model = n,
      specificity = rev(roc_list[[n]]$specificities),
      sensitivity = rev(roc_list[[n]]$sensitivities),
      auc = round(auc(roc_list[[n]]), 3)
    )
  }))
  plot_df$model <- factor(plot_df$model, levels = names(roc_list))
  
  # Create AUC labels
  auc_df <- plot_df %>%
    group_by(model) %>%
    summarise(auc = unique(auc), .groups = "drop")
  
  if (!is.null(model_labels)) {
    auc_df$display_name <- ifelse(names(model_labels) %in% auc_df$model,
                                 model_labels[match(auc_df$model, names(model_labels))],
                                 as.character(auc_df$model))
  } else {
    auc_df$display_name <- as.character(auc_df$model)
  }
  auc_df$label <- paste0(auc_df$display_name, "\n(AUC=", auc_df$auc, ")")
  
  p <- ggplot(plot_df, aes(x = specificity, y = sensitivity, color = model)) +
    geom_line(size = 1.1) +
    geom_abline(slope = 1, intercept = 1, linetype = "dashed", color = "black", linewidth = 0.8) +
    scale_x_reverse(limits = c(1, 0), labels = scales::percent_format(accuracy = 1)) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
    labs(title = title, subtitle = subtitle, 
         x = "Specificity", 
         y = "Sensitivity") +
    scale_color_manual(
      values = darjeeling_cols[1:length(unique(plot_df$model))],
      labels = auc_df$label,
      guide = guide_legend(nrow = length(unique(plot_df$model)), byrow = TRUE)
    ) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      plot.title = element_text(size = 10, hjust = 0.5, margin = margin(b = 2), face = "bold"),
      plot.subtitle = element_text(size = 8, hjust = 0.5, margin = margin(b = 6)),
      axis.title = element_text(size = 8),
      legend.title = element_blank(),
      legend.text = element_text(size = 7, margin = margin(t = 1, b = 1), lineheight = 1),
      legend.key.height = unit(0.35, "cm"),
      legend.spacing.y = unit(0.15, "cm"),
      axis.text.x = element_text(color = "black", size = 8),
      axis.text.y = element_text(color = "black", size = 8),
      axis.ticks = element_line(color = "black"),
      axis.title.x = element_text(size = 8),
      axis.title.y = element_text(size = 8),
      plot.margin = unit(c(5, 5, 5, 5), "pt")
    ) +
    coord_fixed(ratio = 1)
  
  return(p)
}

if (!dir.exists("roccurves_testset")) dir.create("roccurves_testset")

roc_plot_testset <- plot_roc_curves(
  roc_list, 
  title = "Supplementary Figure 7: Monogenic ALS prediction model",
  subtitle = "Accuracy performance on main simulation",
  model_labels = model_labels) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE))

ggsave("roccurves_testset/roc_trained_models_on_testset.pdf", 
       roc_plot_testset, width = 4.8, height = 5, units = "in", dpi = 300)


## Save crossvalidated models for use in other scripts
if (!dir.exists("models")) dir.create("models")
saveRDS(model_cv_ALS,      file = "models/model_cv_ALS.rds")
saveRDS(model_cv_FTD,      file = "models/model_cv_FTD.rds")
saveRDS(model_cv_dementia, file = "models/model_cv_dementia.rds")