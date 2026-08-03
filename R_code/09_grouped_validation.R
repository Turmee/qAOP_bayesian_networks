# Leave-one-group-out validation and permutation testing.
# Source 08_resampling_validation.R before using these functions.

leave_one_group_out_cv <- function(data, dag, target_nodes, group_column,
                                   thresholds, method = "lw", n_samples = 50000,
                                   metadata = NULL) {
  dag_nodes <- bnlearn::nodes(dag)
  missing <- setdiff(c(dag_nodes, group_column, metadata), names(data))
  if (length(missing)) stop("Missing columns: ", paste(missing, collapse = ", "))
  thresholds <- resolve_thresholds(target_nodes, thresholds)
  evidence_nodes <- setdiff(dag_nodes, target_nodes)
  groups <- unique(na.omit(data[[group_column]]))
  output <- vector("list", length(groups))
  for (index in seq_along(groups)) {
    held_out <- groups[[index]]
    test_rows <- !is.na(data[[group_column]]) & data[[group_column]] == held_out
    train <- data[!test_rows, , drop = FALSE]
    test <- data[test_rows, , drop = FALSE]
    fit <- bnlearn::bn.fit(dag, train[, dag_nodes, drop = FALSE], method = "mle-g")
    observed <- data.frame(held_out_group = rep(as.character(held_out), nrow(test)))
    if (length(metadata)) observed <- cbind(observed, test[, metadata, drop = FALSE])
    for (node in target_nodes) {
      observed[[paste0("observed_", node)]] <- test[[node]]
      observed[[paste0(node, "_high")]] <- as.integer(test[[node]] >= thresholds[[node]])
    }
    output[[index]] <- cbind(observed, predict_bn_rows(
      fit, test, evidence_nodes, target_nodes, thresholds, method, n_samples
    ))
  }
  do.call(rbind, output)
}

group_permutation_test <- function(data, dag, target_nodes, group_column,
                                   thresholds, repetitions = 200,
                                   n_samples = 50000, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  observed_cv <- leave_one_group_out_cv(
    data, dag, target_nodes, group_column, thresholds, n_samples = n_samples
  )
  observed <- calculate_bn_metrics(observed_cv, target_nodes)
  observed_values <- setNames(observed$value, observed$metric)
  null <- replicate(repetitions, {
    permuted <- data
    permuted[[group_column]] <- sample(data[[group_column]], replace = FALSE)
    metrics <- calculate_bn_metrics(leave_one_group_out_cv(
      permuted, dag, target_nodes, group_column, thresholds, n_samples = n_samples
    ), target_nodes)
    setNames(metrics$value, metrics$metric)
  })
  p_values <- vapply(seq_along(observed_values), function(i) {
    values <- null[i, is.finite(null[i, ])]
    if (!length(values) || !is.finite(observed_values[i])) return(NA_real_)
    higher_is_better <- grepl("^(R2|AUC)_", names(observed_values)[i])
    extreme <- if (higher_is_better) values >= observed_values[i] else values <= observed_values[i]
    (1 + sum(extreme)) / (1 + length(values))
  }, numeric(1))
  list(
    observed = data.frame(metric = names(observed_values),
                          value = unname(observed_values),
                          empirical_p = p_values, row.names = NULL),
    null_distribution = as.data.frame(t(null))
  )
}

# Example (not executed):
if (FALSE) {
  source("08_resampling_validation.R")
  grouped_cv <- leave_one_group_out_cv(
    my_data, dag, targets, "Study", thresholds, metadata = "Chemical"
  )
  calculate_bn_metrics(grouped_cv, targets)
  group_permutation_test(my_data, dag, targets, "Study", thresholds,
                         repetitions = 200, seed = 123)
}
