# Helper Functions - Rui Tang
library(tidyverse)
library(data.table)
library(boot)
library(ggplot2)
library(diptest)
library(mixtools)
library(philentropy)

find_peak_mode_original <- function(test_Cell_Num) {
  # Create the density plot
  p <- ggplot() +
    geom_density(aes(x = test_Cell_Num)) +
    scale_x_continuous(trans = "log2")
  # Find the peak point of the density plot
  mode_value <- 2^(ggplot_build(p)$data[[1]]$x[which.max(ggplot_build(p)$data[[1]]$y)])
  return(mode_value)
}

find_peak_mode <- function(cell_num, min_threshold = 1000) {
  # Filter data by threshold only if there are more than 20 data points above the threshold
  above_threshold <- cell_num[cell_num > min_threshold]

  # Return NA if insufficient data points
  if (length(above_threshold) < 20) {
    min_threshold <- 0
  }
  cell_num <- cell_num[cell_num > min_threshold]

  # Calculate density directly
  d <- density(log2(cell_num))

  # Find peak point
  mode_value <- 2^(d$x[which.max(d$y)])

  return(mode_value)
}

find_1st_mode <- function(df) {
  df %>%
    group_by(gene, sgID) %>%
    summarise(
      # Calculate density for each group and find primary peak mode
      mode_x = if (n() > 1) 2^(density(log2(Cell_Num), n = 512)$x[which.max(density(log2(Cell_Num), n = 512)$y)]) else NA,
      mode_y = if (n() > 1) max(density(log2(Cell_Num), n = 512)$y) else NA,
      .groups = "drop"
    )
}

find_2nd_mode <- function(df) {
  df %>%
    filter(Cell_Num < 2000) %>%
    group_by(gene, sgID) %>%
    summarise(
      # Calculate density for each group and find secondary peak mode
      mode_x = if (n() > 1) 2^(density(log2(Cell_Num), n = 512)$x[which.max(density(log2(Cell_Num), n = 512)$y)]) else NA,
      mode_y = if (n() > 1) max(density(log2(Cell_Num), n = 512)$y) else NA,
      .groups = "drop"
    )
}

# GeneList == target creates logical vector (TRUE/FALSE)
# Example: c("gene1", "gene2", "gene1") == "gene1"
#          -> c(TRUE, FALSE, TRUE)
# sum() converts TRUE to 1, FALSE to 0 and adds them
calculate_MetSeed <- function(GeneList, target) {
  MetSeed <- sum(GeneList == target)
  return(MetSeed)
}

calculate_dormancy <- function(cell_num_vector, cutoff) {
  below_cutoff <- sum(cell_num_vector < cutoff, na.rm = TRUE)
  total_count <- sum(!is.na(cell_num_vector))
  percentage_below_cutoff <- (below_cutoff / total_count) * 100
  return(percentage_below_cutoff)
}

calculate_supermet <- function(ctc_vector) {
  yes_count <- sum(ctc_vector == "Yes", na.rm = TRUE)
  total_count <- sum(!is.na(ctc_vector))
  percentage_yes <- (yes_count / total_count) * 100
  return(percentage_yes)
}

# functions related to statistic test and p-value
getLower <- function(Data){
  return (quantile(Data, 0.025, na.rm = T))
}

getUpper <- function(Data){
  return (quantile(Data, 0.975, na.rm = T))
}

getEmpiricalP <- function(Data, Baseline=1, NoAdjust = FALSE, Side=2){
  # Side can be adjusted to return one sided or two sided p-values
  # NoAdjust can either be false or the value you want to put in
  # which corresponds to the observed value of that metric before bootstrapping
  Data <- Data[is.finite(Data)]
  Counts = (sum(Data > Baseline) + sum(Data == Baseline)/2)/length(Data)
  
  if (Side == 2){ # Two sided test
    return(min(Counts,1-Counts)*2)
  }else{
    
    if (!NoAdjust){ # One sided test
      return(min(Counts,1-Counts))
    }else{
      if (NoAdjust == Baseline){
        return (0.5)
      }else if(NoAdjust > Baseline){
        return (1-Counts)
      }else{
        return (Counts)
      }
    }
  }
}

# Function to calculate MetBurden given the gene_df, safe_df, target gene, preTran_cell_num and avg_MetBurden_safe
calc_MetBurden <- function(gene_df, safe_df, target, preTran_cell_num) {
  filtered_gene_df <- gene_df %>% filter(sgID == target)
  PreTran_Cell_Num <- preTran_cell_num$cell_num[preTran_cell_num$sgID == target]

  if (nrow(filtered_gene_df) == 0) {
    return(data.frame(target = target, MetBurden = NA, CI_lower = NA, CI_upper = NA, pval = NA, FDR_pval = NA))
  }

  preTran_cutoff <- 5 * quantile(preTran_cell_num$cell_num, 0.01, na.rm = TRUE)
  sgAll <- unique(preTran_cell_num$sgID[preTran_cell_num$cell_num > preTran_cutoff])

  if (!(target %in% sgAll)) {
    return(data.frame(target = target, MetBurden = NA, CI_lower = NA, CI_upper = NA, pval = NA, FDR_pval = NA))
  }

  # Determining the average MetBurden for all of the sgSafe
  avg_MetBurden_safe <- safe_df %>%
    filter(sgID %in% sgAll) %>%
    group_by(sgID) %>%
    summarize(Total_Cell_Num = sum(cell_num)) %>%
    left_join(preTran_cell_num %>% select(sgID, gene, cell_num) %>% dplyr::rename(PreTran_Cell_Num = cell_num), by = "sgID") %>%
    mutate(MetBurden = Total_Cell_Num / PreTran_Cell_Num) %>%
    summarize(avg_MetBurden_safe = mean(MetBurden, na.rm = TRUE))

  # Bootstrap test
  boot_result <- boot(filtered_gene_df, function(d, i) {
    sum(d$cell_num[i]) / PreTran_Cell_Num / avg_MetBurden_safe$avg_MetBurden_safe
  }, R = 100)

  # Get p-value
  pval <- getEmpiricalP(boot_result$t, Baseline = 1, Side = 2)
  FDR_pval <- p.adjust(pval, method = "BH")

  return(data.frame(target = target, MetBurden = boot_result$t0, CI_lower = getLower(boot_result$t), CI_upper = getUpper(boot_result$t), pval = pval, FDR_pval = FDR_pval, row.names = NULL))
}

# SizePercentile Analysis
calc_SizePercentile <- function(gene_df, safe_df, target, percentile = 0.9) {
  gene_filtered_df <- gene_df %>% filter(sgID == target) #%>% filter(n() >= 20)

  if (nrow(gene_filtered_df) == 0) {
    return(data.frame(target = target, SizePercentile = NA, CI_lower = NA, CI_upper = NA, pval = NA, FDR_pval = NA))
  }

  safe_percentile_size <- safe_df %>%
    summarize(cell_num = quantile(cell_num, percentile, na.rm = TRUE))

  # Bootstrap test to identify 90th percentile size
  boot_result <- boot(gene_filtered_df, function(d, i) {
    quantile(d$cell_num[i], percentile) / safe_percentile_size$cell_num
  }, R = 1000)

  # Get p-value
  pval <- getEmpiricalP(boot_result$t, Baseline = 1, NoAdjust = FALSE, Side = 2)
  FDR_pval <- p.adjust(pval, method = "BH")

  return(data.frame(target = target, SizePercentile = boot_result$t0, CI_lower = getLower(boot_result$t), CI_upper = getUpper(boot_result$t), pval = pval, FDR_pval = FDR_pval, row.names = NULL))
}

calc_MetSeeding <- function(gene_df, safe_df, target, preTran_cell_num) {
  gene_filtered_df <- gene_df %>% filter(sgID == target)
  PreTran_Cell_Num <- preTran_cell_num$cell_num[preTran_cell_num$sgID == target]

  if (nrow(gene_filtered_df) == 0) {
    return(data.frame(target = target, MetSeeding = NA, CI_lower = NA, CI_upper = NA, pval = NA, FDR_pval = NA))
  }

  preTran_cutoff <- 5 * quantile(preTran_cell_num$cell_num, 0.01, na.rm = TRUE)
  sgAll <- unique(preTran_cell_num$sgID[preTran_cell_num$cell_num > preTran_cutoff])

  if (!(target %in% sgAll)) {
    return(data.frame(target = target, MetSeeding = NA, CI_lower = NA, CI_upper = NA, pval = NA, FDR_pval = NA))
  }

  avg_MetSeed_safe <- safe_df %>%
    filter(sgID %in% sgAll) %>%
    group_by(sgID) %>%
    summarize(MetSeed = n()) %>%
    left_join(preTran_cell_num %>% select(sgID, gene, cell_num) %>% dplyr::rename(PreTran_Cell_Num = cell_num), by = "sgID") %>%
    mutate(Ratio = MetSeed / PreTran_Cell_Num) %>%
    summarize(avg_MetSeed_safe = mean(Ratio, na.rm = TRUE))

  # Bootstrap test
  boot_result <- boot(gene_filtered_df, function(d, i) {
    nrow(d[i, ]) / PreTran_Cell_Num / avg_MetSeed_safe$avg_MetSeed_safe
  }, R = 100)

  # Get p-value
  pval <- getEmpiricalP(boot_result$t, Baseline = 1, Side = 2)
  FDR_pval <- p.adjust(pval, method = "BH")

  return(data.frame(target = target, MetSeeding = boot_result$t0, CI_lower = getLower(boot_result$t), CI_upper = getUpper(boot_result$t), pval = pval, FDR_pval = FDR_pval, row.names = NULL))
}

calc_PeakMode <- function(gene_df, safe_df, target) {
  gene_filtered_df <- gene_df %>% filter(sgID == target) %>% filter(n() >= 20)

  if (nrow(gene_filtered_df) == 0) {
    return(data.frame(target = target, PeakMode = NA, CI_lower = NA, CI_upper = NA, pval = NA, FDR_pval = NA))
  }

  safe_peak_mode <- safe_df %>%
    summarize(pm_safe = find_peak_mode(cell_num))

  # Bootstrap test to identify peak mode
  boot_result <- boot(gene_filtered_df, function(d, i) {
    find_peak_mode(d$cell_num[i]) / safe_peak_mode$pm_safe
  }, R = 1000)

  # Get p-value
  pval <- getEmpiricalP(boot_result$t, Baseline = 1, NoAdjust = FALSE, Side = 2)
  FDR_pval <- p.adjust(pval, method = "BH")

  return(data.frame(target = target, PeakMode = boot_result$t0, CI_lower = getLower(boot_result$t), CI_upper = getUpper(boot_result$t), pval = pval, FDR_pval = FDR_pval, row.names = NULL))
}

calc_MetDormancy <- function(gene_df, safe_df, target) {
  gene_filtered_df <- gene_df %>% filter(sgID == target)

  if (nrow(gene_filtered_df) == 0) {
    return(data.frame(target = target, MetDormancy = NA, CI_lower = NA, CI_upper = NA, pval = NA, FDR_pval = NA))
  }

  safe_dormancy <- safe_df %>%
    group_by(sgID) %>%
    summarize(dormancy = calculate_dormancy(cell_num, 2000)) %>%
    ungroup() %>%
    summarize(avg_md_safe = mean(dormancy, na.rm = TRUE))

  # Bootstrap test
  boot_result <- boot(gene_filtered_df, function(d, i) {
    calculate_dormancy(d$cell_num[i], 2000) / safe_dormancy$avg_md_safe
  }, R = 1000)

  # Get p-value
  pval <- getEmpiricalP(boot_result$t, Baseline = 1, NoAdjust = FALSE, Side = 2)
  FDR_pval <- p.adjust(pval, method = "BH")

  return(data.frame(target = target, MetDormancy = boot_result$t0, CI_lower = getLower(boot_result$t), CI_upper = getUpper(boot_result$t), pval = pval, FDR_pval = FDR_pval, row.names = NULL))
}

plot_dormancy <- function(gene_df, target) {
  # Skip sgID if it has less than 100 rows
  if (nrow(gene_df) < 20) {
    return(NULL)
  }

  # Create log2 transformed Cell_Num column
  gene_df$Cell_Num_log2 <- log2(gene_df$Cell_Num)

  # Between Peak1 and Peak2, find the lowest peak


  # Select the data vector with log2 Cell_Num
  data_vector <- df_gene$Cell_Num_log2

  # Set the fixed means for Components 1 and 2
  Mode1 <- find_1st_mode(df_gene)
  Mode2 <- find_2nd_mode(df_gene)

  fixed_mean1 <- log2(Mode1$mode_x[Mode1$sgID == target])
  fixed_mean2 <- log2(512) # log2(Mode2$mode_x[Mode2$sgID == sg])

  # Fit a mixture model with constraints (fix means for both components)
  custom_mix_model <- normalmixEM(data_vector, k = 2, mean.constr = c(fixed_mean1, fixed_mean2))

  # Get the posterior probabilities for each component
  posterior_probs <- custom_mix_model$posterior

  # Calculate the percentage of data in each component
  component1_percentage <- sum(posterior_probs[, 1]) / length(data_vector) * 100
  component2_percentage <- sum(posterior_probs[, 2]) / length(data_vector) * 100


  # set dormancy cutoff
  DorCutOff <- 1000
  df_dorm <- df_gene %>% filter(Cell_Num < DorCutOff)

  # Ensure single values for gene, Tissue, and Mouse_Genotype; set NA if none found
  gene_value <- ifelse(length(unique(df_gene$gene)) > 0, unique(df_gene$gene)[1], NA)
  tissue_value <- ifelse(length(unique(df_gene$Tissue)) > 0, unique(df_gene$Tissue)[1], NA)
  genotype_value <- ifelse(length(unique(df_gene$Mouse_Genotype)) > 0, unique(df_gene$Mouse_Genotype)[1], NA)
  timepoint_value <- ifelse(length(unique(df_gene$Time_Point)) > 0, unique(df_gene$Time_Point)[1], NA)

  # Add a new row to the summary dataframe
  summary_df <- rbind(summary_df, data.frame(
    sgID = sg,
    gene = gene_value,
    Tissue = tissue_value,
    Mouse_Genotype = genotype_value,
    Time_Point = timepoint_value,
    Dormancy_Percent = 100 * nrow(df_dorm) / nrow(df_gene),
    FirstPeak_Mode = 2^fixed_mean1, # Convert back from log2 scale
    FirstPeak_Percent = component1_percentage,
    SecondPeak_Mode = 2^fixed_mean2, # Convert back from log2 scale
    SecondPeak_Percent = component2_percentage
  ))
}

check_bimodal <- function(data) {
  data <- log2(data)

  # Hartigan's dip test
  require(diptest)
  dip_result <- dip.test(data)

  # Kernel density estimation
  d <- density(data)

  # Find peaks (local maxima)
  peaks <- which(diff(sign(diff(d$y))) == -2) + 1
  peak_values <- d$x[peaks]
  peak_densities <- d$y[peaks]

  # Filter where peaks are less than log2(50)
  #peaks <- peaks[peak_values > log2(50)]
  #peak_densities <- peak_densities[peak_values > log2(50)]
  #peak_values <- peak_values[peak_values > log2(50)]

  # Find valley (local minimum between peaks)
  if (length(peaks) >= 2) {
    valley <- which(diff(sign(diff(d$y))) == 2) + 1
    valley_values <- d$x[valley]
    valley_densities <- d$y[valley]

    # Filter where valley is between the two peaks
    #valley <- valley[valley_values > peak_values[1] & valley_values < peak_values[2]]
    #valley_densities <- valley_densities[valley_values > peak_values[1] & valley_values < peak_values[2]]
    #valley_values <- valley_values[valley_values > peak_values[1] & valley_values < peak_values[2]]
  } else {
    valley_values <- 0
    valley_densities <- 0
  }

  return(list(
    is_bimodal = dip_result$p.value < 0.05,
    peaks = peak_values,
    peak_heights = peak_densities,
    valley = valley_values,
    valley_heights = valley_densities,
    density = d
  ))
}

fit_normal_sd <- function(data, peak, peak_height, valley) {
  # Filter data by valley
  if (peak < valley) {
    filtered_data <- data[log2(data) < valley] # First peak
    filtered_data <- filtered_data[log2(filtered_data) >= log2(50)]
  } else {
    filtered_data <- data[log2(data) > valley] # Second peak
  }

  # Get density of filtered data
  d <- density(log2(filtered_data))

  # Function to minimize
  fit_func <- function(sd) {
    norm_density <- dnorm(d$x, mean = peak, sd = sd) *
      (peak_height / dnorm(peak, mean = peak, sd = sd))
    sum((d$y - norm_density)^2)
  }

  opt <- optimize(fit_func, interval = c(0.1, 5))
  return(opt$minimum)
}

plot_bimodal <- function(data, peaks_data, show_fit = TRUE, show_peaks = TRUE, show_valleys = TRUE, show_percentage = TRUE, return_data = FALSE) {
  peaks <- peaks_data$peaks
  peak_heights <- peaks_data$peak_heights
  #peak_heights <- peak_heights[peaks > log2(50)]
  #peaks <- peaks[peaks > log2(50)]
  valleys <- peaks_data$valley
  valley_heights <- peaks_data$valley_heights

  # If there are more than 2 peaks, we only consider the highest 2
  if (length(peaks) > 2){
    # Get indices of top 2 peaks by height
    top_indices <- order(peak_heights, decreasing = TRUE)[1:2]
    peaks_subset <- peaks[top_indices]
    heights_subset <- peak_heights[top_indices]

    # Sort by position (left to right)
    position_order <- order(peaks_subset)
    peaks <- peaks_subset[position_order]
    peak_heights <- heights_subset[position_order]
  }

  # Handle single peak case
  if (length(peaks) < 2) {
    # Duplicate single peak for second component
    peaks[2] <- 0
    peak_heights[2] <- 0
  }

  # Calculate fit parameters if requested
  if (show_fit) {
    # Always fit second peak
    sd2 <- fit_normal_sd(data, peaks[2], peak_heights[2], valleys[1])

    # Fit first peak only if two peaks exist
    if (length(peaks) >= 2) {
      sd1 <- fit_normal_sd(data, peaks[1], peak_heights[1], valleys[1])
    }
  }
  custom_mix_model <- normalmixEM(log2(data), k = 2, mean.constr = c(peaks[1], peaks[2]))
  # Get the posterior probabilities for each component
  posterior_probs <- custom_mix_model$posterior
  # Calculate the percentage of data in each component
  component1_percentage <- sum(posterior_probs[, 1]) / length(data) * 100
  component2_percentage <- sum(posterior_probs[, 2]) / length(data) * 100

  plot_data <- data.frame(x = log2(data))
  p <- ggplot() +
    geom_density(data = plot_data, aes(x = x), alpha = 0.3) +
     theme_minimal() +
       scale_x_continuous(
         breaks = seq(0, 30, 1),
         limits = c(0, 30)
       ) +
       labs(
         x = "log2(cell_num)",
         y = "Density"
       )

  fit_line1 <- NA
  fit_line2 <- NA
  if (show_fit) {
    fit_line1 <- data.frame(
      x = seq(peaks[1] - 4 * sd1, peaks[1] + 4 * sd1, length.out = length(data))
    ) %>% mutate(
      y = dnorm(x, mean = peaks[1], sd = sd1) *
        (peak_heights[1] / dnorm(peaks[1], mean = peaks[1], sd = sd1))
    )
    fit_line2 <- data.frame(
      x = seq(peaks[2] - 4 * sd2, peaks[2] + 4 * sd2, length.out = length(data))
    ) %>% mutate(
      y = dnorm(x, mean = peaks[2], sd = sd2) *
        (peak_heights[2] / dnorm(peaks[2], mean = peaks[2], sd = sd2))
    )
    p <- p +
      geom_line(
        data = fit_line1,
        aes(x = x, y = y),
        color = "red"
      ) +
      geom_line(
        data = fit_line2,
        aes(x = x, y = y),
        color = "blue"
      )
  }

  if (show_peaks) {
    p <- p +
      geom_vline(xintercept = peaks, color = "red", linetype = "dashed") +
      geom_text(aes(
        x = peaks[1], 
        y = peak_heights[1] + 0.005, 
        label = format(round(2^peaks[1], 2), big.mark = ",")
      )) +
      geom_text(aes(
        x = peaks[2],
        y = peak_heights[2] + 0.005,
        label = format(round(2^peaks[2], 2), big.mark = ",")
      ))
  }

  if (show_valleys) {
    p <- p +
      geom_vline(xintercept = valleys, color = "blue", linetype = "dashed") +
      geom_text(aes(
        x = valleys[1],
        y = valley_heights[1] + 0.005,
        label = format(round(2^valleys[1], 2), big.mark = ",")
      ))
  }

  if (show_percentage) {
    p <- p +
      geom_text(aes(x = peaks[1], y = peak_heights[1] + 0.01, label = paste0(round(component1_percentage, 2), "%"))) +
      geom_text(aes(x = peaks[2], y = peak_heights[2] + 0.01, label = paste0(round(component2_percentage, 2), "%")))
  }
  res <- data.frame(
    plot_data = plot_data,
    peak1 = peaks[1],
    peak2 = peaks[2],
    valley = valleys[1],
    peak1_height = peak_heights[1],
    peak2_height = peak_heights[2],
    valley_height = valley_heights[1],
    fit_line1 = fit_line1,
    fit_line2 = fit_line2,
    comp1_perc = component1_percentage,
    comp2_perc = component2_percentage
  )
  if (return_data) {
    return(res)
  } else {
    p
  }
}

compare_distributions <- function(dist1, dist2) {
  # KL divergence
  require(philentropy)

  # Create sequence of points
  x <- seq(min(dist1, dist2), max(dist1, dist2), length.out = 1000)

  # Get densities
  d1 <- density(dist1)
  d2 <- density(dist2)

  # Statistical tests
  t_test <- t.test(dist1, dist2)
  var_test <- var.test(dist1, dist2)

  # Calculate overlap coefficient
  overlap <- sum(pmin(d1$y, d2$y)) / sum(d1$y)

  return(list(
    t_test = t_test$p.value,
    var_test = var_test$p.value,
    overlap = overlap
  ))
}
