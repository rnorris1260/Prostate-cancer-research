compute_cn_deletion_threshold <- function(
    df,
    cn_col       = "CN",
    expr_col     = "expr",       # ideally log2(TPM + 1)
    purity_col   = NULL,         # optional: column name with tumour purity
    log_expr     = TRUE,         # set FALSE if expr is raw TPM/CPM
    pct_drop     = 0.5,          # X% drop (0.5 = 50% drop)
    diploid_range = c(1.8, 2.2), # CN range considered "diploid"
    min_n        = 20,           # minimum samples needed
    plot         = FALSE,        # set TRUE to see a diagnostic plot
    gene_name    = NULL          # optional, for plot title
) {
  if (!requireNamespace("mgcv", quietly = TRUE)) {
    stop("Please install 'mgcv' for smoothing: install.packages('mgcv')")
  }
  if (plot && !requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Please install 'ggplot2' for plotting: install.packages('ggplot2')")
  }
  
  # subset & rename
  cols <- c(cn_col, expr_col, purity_col)
  cols <- cols[!is.na(cols)]
  d <- df[, cols, drop = FALSE]
  names(d)[names(d) == cn_col]   <- "CN"
  names(d)[names(d) == expr_col] <- "expr"
  if (!is.null(purity_col)) {
    names(d)[names(d) == purity_col] <- "purity"
  }
  
  # drop NAs
  d <- d[complete.cases(d[, c("CN", "expr", if (!is.null(purity_col)) "purity")]), ]
  d$CN <- as.numeric(d$CN)
  if (!is.numeric(d$expr)) d$expr <- as.numeric(d$expr)
  
  if (nrow(d) < min_n) {
    stop("Not enough complete samples (n = ", nrow(d), "). Need at least ", min_n, ".")
  }
  
  # ensure log2 scale
  if (!log_expr) {
    d$expr <- log2(d$expr + 1)
  }
  
  # purity adjustment (optional)
  if (!is.null(purity_col)) {
    fit_pur <- lm(expr ~ purity, data = d)
    # residuals + mean to keep scale interpretable
    d$expr_adj <- resid(fit_pur) + mean(d$expr, na.rm = TRUE)
  } else {
    d$expr_adj <- d$expr
  }
  
  # fit smooth curve expr_adj ~ s(CN)
  library(mgcv)
  gam_fit <- mgcv::gam(expr_adj ~ s(CN, k = 5), data = d, method = "REML")
  
  # prediction grid across observed CN range
  cn_grid <- seq(min(d$CN, na.rm = TRUE),
                 max(d$CN, na.rm = TRUE),
                 length.out = 200)
  pred_log <- predict(gam_fit, newdata = data.frame(CN = cn_grid))
  
  # 1) diploid reference expression (in log scale)
  dip_idx <- which(cn_grid >= diploid_range[1] & cn_grid <= diploid_range[2])
  if (length(dip_idx) < 5) {
    warning("Few points in diploid_range; threshold may be unstable.")
  }
  expr_dip_log <- mean(pred_log[dip_idx], na.rm = TRUE)
  
  # convert to linear TPM-like scale, apply % drop, convert back to log
  expr_dip_lin <- 2^expr_dip_log - 1
  expr_target_lin <- expr_dip_lin * (1 - pct_drop)      # e.g. 50% of diploid
  expr_target_log <- log2(expr_target_lin + 1)
  
  # 2) find CN* where predicted expression falls below that target
  # for deletions we care about CN <= diploid upper bound
  loss_region_idx <- which(cn_grid <= diploid_range[2])
  idx_below <- loss_region_idx[which(pred_log[loss_region_idx] <= expr_target_log)]
  
  if (length(idx_below) == 0) {
    cn_star <- NA_real_
    warning("No CN where predicted expression is at least ",
            round(pct_drop * 100), "% below diploid; CN* set to NA.")
  } else {
    # we want the *highest* CN in the loss region where it's below target
    cn_star <- max(cn_grid[idx_below])
  }
  
  # optional plot
  p <- NULL
  if (plot) {
    library(ggplot2)
    plot_df <- data.frame(
      CN          = cn_grid,
      expr_pred   = pred_log
    )
    p <- ggplot(d, aes(x = CN, y = expr_adj)) +
      geom_point(alpha = 0.5) +
      geom_line(data = plot_df, aes(y = expr_pred), linewidth = 1) +
      geom_vline(xintercept = cn_star, linetype = "dashed") +
      geom_hline(yintercept = expr_target_log, linetype = "dotted") +
      annotate("text",
               x = cn_star,
               y = expr_target_log,
               label = paste0("CN* ≈ ", round(cn_star, 2)),
               vjust = -1, hjust = 0.5, size = 3) +
      labs(
        title = paste0("Deletion threshold CN* (", round(pct_drop * 100), "% drop)",
                       if (!is.null(gene_name)) paste0(" for ", gene_name) else ""),
        x = "Copy number",
        y = "Purity-adjusted expression (log2 scale)"
      ) +
      theme_bw()
  }
  
  return(list(
    CN_star              = cn_star,
    diploid_log_expr     = expr_dip_log,
    target_log_expr      = expr_target_log,
    pct_drop             = pct_drop,
    diploid_range        = diploid_range,
    gam_fit              = gam_fit,
    plot                 = p
  ))
}


res_RB1 <- compute_cn_deletion_threshold(
  RB1_cnv_rna_purity,
  cn_col        = "CN",
  expr_col      = "expr",   # log2(TPM+1)
  purity_col    = "tumorPurity",
  log_expr      = TRUE,
  pct_drop      = 0.5,          # 50% drop threshold
  diploid_range = c(1.8, 2.2),
  plot          = TRUE,
  gene_name     = "RB1"
)

res_RB1$CN_star      # this is your biologically interpretable CN*
res_RB1$plot         # see the diagnostic plot
