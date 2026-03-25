library(dplyr)

otus <- read.csv("MAG_covem_oil.csv", row.names = 1)

num_samples <- nrow(otus)
cat(paste("Total number of samples:", num_samples, "\n"))

calculate_typeI_ple <- function(sample_data) {
  m_values <- numeric(nrow(sample_data))
  V_values <- numeric(nrow(sample_data))
  
  for (i in 1:nrow(sample_data)) {
    row_data <- sample_data[i, ]
    non_zero_values <- row_data[row_data > 0]
    
    if (length(non_zero_values) >= 2) {
      m_values[i] <- mean(non_zero_values)
      V_values[i] <- var(non_zero_values)
    } else {
      m_values[i] <- NA
      V_values[i] <- NA
    }
  }
  
  valid_mask <- !is.na(m_values) & !is.na(V_values) & m_values > 0
  m_valid <- m_values[valid_mask]
  V_valid <- V_values[valid_mask]
  
  ln_m <- log(m_valid)
  ln_V <- log(V_valid)
  
  if (length(ln_m) < 2) {
    return(list(
      b = NA,
      ln_a = NA,
      R = NA,
      R_squared = NA,
      p_value = NA,
      lm_fit = NA,
      m_values = m_values,
      V_values = V_values
    ))
  }
  
  lm_fit <- lm(ln_V ~ ln_m)
  b <- coef(lm_fit)[2]
  ln_a <- coef(lm_fit)[1]
  r_squared <- summary(lm_fit)$r.squared
  r <- sqrt(r_squared)
  p_value <- summary(lm_fit)$coefficients[2, 4]
  
  return(list(
    b = b,
    ln_a = ln_a,
    R = r,
    R_squared = r_squared,
    p_value = p_value,
    lm_fit = lm_fit,
    m_values = m_values,
    V_values = V_values
  ))
}

overall_results <- calculate_typeI_ple(otus)
overall_results$N <- num_samples

cat("PLE results for oil samples:\n")
cat(paste("  Number of samples =", overall_results$N, "\n"))
cat(paste("  b =", round(overall_results$b, 4), "\n"))
cat(paste("  ln(a) =", round(overall_results$ln_a, 4), "\n"))
cat(paste("  R =", round(overall_results$R, 4), "\n"))
cat(paste("  R-squared =", round(overall_results$R_squared, 4), "\n"))
cat(paste("  p-value =", format(overall_results$p_value, scientific = TRUE, digits = 5), "\n\n"))

sample_mv <- data.frame(
  SampleID = rownames(otus),
  m = overall_results$m_values,
  V = overall_results$V_values
)
write.csv(sample_mv, "mv_values_oil.csv", row.names = FALSE)

final_results <- data.frame(
  Group = "Oil",
  N = overall_results$N,
  b = overall_results$b,
  ln_a = overall_results$ln_a,
  R = overall_results$R,
  R_squared = overall_results$R_squared,
  p_value = overall_results$p_value
)
write.csv(final_results, "results_oil.csv", row.names = FALSE)

out_file <- "PLE_oil.pdf"

if (!is.na(overall_results$b) && !is.null(overall_results$lm_fit)) {
  tryCatch({
    pdf(out_file, width = 6, height = 5)
    
    valid_data <- sample_mv[complete.cases(sample_mv), ]
    
    plot(
      log(valid_data$m),
      log(valid_data$V),
      main = "Type-I Taylor's Power Law in oil samples",
      xlab = "ln(m)",
      ylab = "ln(V)",
      pch = 21,
      col = "#e67e22",
      bg = adjustcolor("#e67e22", alpha.f = 0.35),
      cex = 1,
      lwd = 0.7
    )
    
    abline(overall_results$lm_fit, col = "#e67e22", lwd = 2)
    
    legend(
      "topleft",
      legend = c(
        paste("R² =", round(overall_results$R_squared, 2)),
        paste("b =", round(overall_results$b, 3)),
        paste("p =", format(overall_results$p_value, scientific = TRUE, digits = 2))
      ),
      bty = "n",
      cex = 0.9
    )
    
    dev.off()
    cat("Plot saved to:", out_file, "\n")
    
  }, error = function(e) {
    cat("Plotting error:", e$message, "\n")
    if (dev.cur() != 1) dev.off()
  })
} else {
  cat("Insufficient data for plotting\n")
}
