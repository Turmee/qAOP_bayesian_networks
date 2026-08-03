# Reusable k-fold and bootstrap validation for Gaussian Bayesian networks.
# This file defines functions only; it does not load data or run an analysis.

safe_r2 <- function(observed, predicted) {
  keep <- complete.cases(observed, predicted)
  if (sum(keep) < 2 || var(observed[keep]) == 0 || var(predicted[keep]) == 0) return(NA_real_)
  cor(observed[keep], predicted[keep])^2
}

safe_nmse <- function(observed, predicted) {
  keep <- complete.cases(observed, predicted)
  if (sum(keep) < 2 || var(observed[keep]) == 0) return(NA_real_)
  mean((observed[keep] - predicted[keep])^2) / var(observed[keep])
}

resolve_thresholds <- function(target_nodes, thresholds) {
  if (length(thresholds) == 1 && is.null(names(thresholds))) {
    return(setNames(rep(as.numeric(thresholds), length(target_nodes)), target_nodes))
  }
  if (is.null(names(thresholds)) || !all(target_nodes %in% names(thresholds))) {
    stop("Provide one threshold or a named threshold for every target node.")
  }
  setNames(as.numeric(thresholds[target_nodes]), target_nodes)
}

calculate_bn_metrics <- function(results, target_nodes) {
  values <- list()
  for (node in target_nodes) {
    observed <- results[[paste0("observed_", node)]]
    predicted <- results[[paste0("predicted_", node)]]
    truth <- results[[paste0(node, "_high")]]
    probability <- results[[paste0("probability_", node, "_high")]]
    values[[paste0("R2_", node)]] <- safe_r2(observed, predicted)
    values[[paste0("NMSE_", node)]] <- safe_nmse(observed, predicted)
    keep <- complete.cases(truth, probability)
    values[[paste0("AUC_", node)]] <- if (sum(keep) > 1 && length(unique(truth[keep])) == 2) {
      as.numeric(pROC::auc(pROC::roc(truth[keep], probability[keep], quiet = TRUE)))
    } else NA_real_
  }
  data.frame(metric = names(values), value = unlist(values), row.names = NULL)
}

predict_bn_rows <- function(fit, test, evidence_nodes, target_nodes, thresholds,
                            method = "lw", n_samples = 50000) {
  output <- vector("list", nrow(test))
  for (i in seq_len(nrow(test))) {
    evidence <- as.list(test[i, evidence_nodes, drop = FALSE])
    posterior <- bnlearn::cpdist(fit, nodes = target_nodes, evidence = evidence,
                                 method = method, n = n_samples)
    row <- list()
    for (node in target_nodes) {
      row[[paste0("predicted_", node)]] <- mean(posterior[[node]])
      row[[paste0("prediction_sd_", node)]] <- sd(posterior[[node]])
      row[[paste0("probability_", node, "_high")]] <-
        mean(posterior[[node]] >= thresholds[[node]])
    }
    output[[i]] <- as.data.frame(row, check.names = FALSE)
  }
  do.call(rbind, output)
}

bn_kfold_cv <- function(data, dag, target_nodes, thresholds, k_folds = 5,
                        method = "lw", n_samples = 50000, metadata = NULL,
                        seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  dag_nodes <- bnlearn::nodes(dag)
  missing <- setdiff(c(dag_nodes, metadata), names(data))
  if (length(missing)) stop("Missing columns: ", paste(missing, collapse = ", "))
  thresholds <- resolve_thresholds(target_nodes, thresholds)
  evidence_nodes <- setdiff(dag_nodes, target_nodes)
  folds <- sample(rep(seq_len(k_folds), length.out = nrow(data)))
  output <- vector("list", k_folds)
  for (fold in seq_len(k_folds)) {
    train <- data[folds != fold, , drop = FALSE]
    test <- data[folds == fold, , drop = FALSE]
    fit <- bnlearn::bn.fit(dag, train[, dag_nodes, drop = FALSE], method = "mle-g")
    observed <- data.frame(fold = rep(fold, nrow(test)))
    if (length(metadata)) observed <- cbind(observed, test[, metadata, drop = FALSE])
    for (node in target_nodes) {
      observed[[paste0("observed_", node)]] <- test[[node]]
      observed[[paste0(node, "_high")]] <- as.integer(test[[node]] >= thresholds[[node]])
    }
    output[[fold]] <- cbind(observed, predict_bn_rows(
      fit, test, evidence_nodes, target_nodes, thresholds, method, n_samples
    ))
  }
  do.call(rbind, output)
}

bootstrap_bn_validation <- function(data, dag, target_nodes, thresholds,
                                    repetitions = 200, k_folds = 5,
                                    n_samples = 50000, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  runs <- replicate(repetitions, {
    sampled <- data[sample(seq_len(nrow(data)), replace = TRUE), , drop = FALSE]
    metrics <- calculate_bn_metrics(bn_kfold_cv(
      sampled, dag, target_nodes, thresholds, k_folds, n_samples = n_samples
    ), target_nodes)
    setNames(metrics$value, metrics$metric)
  })
  data.frame(metric = rownames(runs), mean = apply(runs, 1, mean, na.rm = TRUE),
             sd = apply(runs, 1, sd, na.rm = TRUE),
             lower_95 = apply(runs, 1, quantile, 0.025, na.rm = TRUE),
             upper_95 = apply(runs, 1, quantile, 0.975, na.rm = TRUE),
             row.names = NULL)
}

# Example (not executed):
if (FALSE) {
  dag <- bnlearn::model2network("[Exposure][KeyEvent1|Exposure][KeyEvent2|KeyEvent1]")
  targets <- c("KeyEvent1", "KeyEvent2")
  thresholds <- c(KeyEvent1 = 1.5, KeyEvent2 = 1.5)
  cv <- bn_kfold_cv(my_data, dag, targets, thresholds, metadata = c("Study", "Chemical"), seed = 123)
  calculate_bn_metrics(cv, targets)
  bootstrap_bn_validation(my_data, dag, targets, thresholds, seed = 123)
}
