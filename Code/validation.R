# This is the code used to generate the BL6 vs NSG vs Rag2 KO

# Author: Irenaeus Chan
# Date: 03/20/2025

# Load the necessary libraries
library(data.table)
library(tidyverse)
library(readxl)
library(ggplot2)

# Compare MobaV_FullSet AndyXu vs MobaSeq Pipeline
andy_res <- fread("/Volumes/ruit/Active/AX_Working_Folder/2024_MobaV_Fullset/2024_MobaV_Fullset_FilteredSampleInfo.csv")
moba_res <- fread("/Volumes/ruit/Active/IrenaeusChan/MobaV_FullSet/outputs_og/FilteredSampleInfo.csv")
moba_new_regression <- fread("/Volumes/ruit/Active/IrenaeusChan/MobaV_FullSet/testing/combined_results/Mobaseq_FilteredSampleInfo.csv")

andy_df <- andy_res %>%
    filter(dist == 0) %>%
    filter(sgID != "sgDummy") %>%
    filter(Sample_ID != "preTran") %>%
    filter(Cell_Num >= 1) %>%
    mutate(
        unique_id = paste(Sample_ID, sgID, BC, sep = "_")
    ) %>%
    dplyr::rename(
        AX_cell_num = Cell_Num
    ) %>%
    select(unique_id, Sample_ID, Tissue, AX_cell_num)

moba_new_df <- moba_new_regression %>%
    filter(distance == 0) %>%
    filter(sgID != "sgDummy") %>%
    filter(Sample_ID != "preTran") %>%
    mutate(
        unique_id = paste(Sample_ID, sgID, barcode, sep = "_")
    ) %>%
    dplyr::rename(
        IC_cell_num = cell_num
    ) %>%
    select(unique_id, Sample_ID, Tissue, IC_cell_num)

moba_df <- moba_res %>%
    filter(distance == 0) %>%
    filter(sgID != "sgDummy") %>%
    filter(Sample_ID != "preTran") %>%
    mutate(
        unique_id = paste(Sample_ID, sgID, barcode, sep = "_")
    ) %>%
    dplyr::rename(
        IC_cell_num = cell_num
    ) %>%
    select(unique_id, Sample_ID, Tissue, IC_cell_num)

# Full join to keep all rows
data_to_plot <- full_join(andy_df, moba_df, by = c("unique_id", "Sample_ID", "Tissue"))
data_to_plot <- data_to_plot %>%
    mutate(across(where(is.numeric), ~ replace_na(., 0))) %>%
    mutate(
        Missing = ifelse(AX_cell_num == 0 | IC_cell_num == 0, "Missing", "Present"),
        Missing = as.factor(Missing)
    )
# Calculate regression on filtered data
reg_data <- data_to_plot %>%
    filter(AX_cell_num > 0, IC_cell_num > 0)

lm_fit <- lm(log10(IC_cell_num) ~ log10(AX_cell_num), data = reg_data)
r2 <- summary(lm_fit)$r.squared

# Plot
data_to_plot %>%
    mutate(
        #sample_label = ifelse(abs(AX_cell_num - IC_cell_num) > 100000, Sample_ID, "")
        # sample_label = case_when(
        #     Sample_ID == "73210-Liver" ~ Sample_ID,
        #     Sample_ID == "73213-Liver" ~ Sample_ID,
        #     Sample_ID == "73400-Liver" ~ Sample_ID,
        #     Sample_ID == "73207-Liver" ~ Sample_ID,
        #     TRUE ~ "Other"
        # )
    ) %>%
    ggplot(
        aes(
            x = AX_cell_num,
            y = IC_cell_num,
            #color = sample_label
        )
    ) +
    # geom_point(aes(color = Missing)) +
    geom_point() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    geom_smooth(data = reg_data, method = "lm", formula = y ~ x, color = "blue", se = FALSE) +
    geom_text(
        data = data.frame(
            x = 100000000,
            y = max(reg_data$IC_cell_num),
            text_label = sprintf("R² = %.3f", r2)
        ),
        aes(x = x, y = y, label = text_label),
        hjust = 0, vjust = 1
    ) +
    geom_text(data = data_to_plot %>% filter(AX_cell_num == 0) %>% summarize(total = n()), aes(x = 0.001, y = 10, label = total), hjust = 1, vjust = 1) +
    geom_text(data = data_to_plot %>% filter(IC_cell_num == 0) %>% summarize(total = n()), aes(x = 10, y = 0.001, label = total), hjust = 1, vjust = 1) +
    scale_x_log10(
        limits = c(0.001, NA),
        breaks = 10^seq(-1, 10, 1), # Creates 0.1, 1, 10, 100, 1000
        labels = scales::label_comma()
    ) +
    scale_y_log10(
        limits = c(0.001, NA),
        breaks = 10^seq(-1, 10, 1),
        labels = scales::label_comma()
    ) +
    theme_bw() +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        # legend.position = "none"
    ) +
    labs(
        x = "AndyXu Cell Number",
        y = "MobaSeq Cell Number"
    )

# There was an issue with the Spike-In Regression Algorithm
# In some cases, the regression line was not calculated correctly because there were some Spike-Ins that had too low counts
# This caused the regression line to be skewed and the Z-score was not a good reprensentation of the data
spike_in_reads <- fread("/Volumes/ruit/Active/IrenaeusChan/MobaV_FullSet/outputs/MobaV_Full_SpikeInInfo.csv")
spike_in_reads <- fread("/Volumes/ruit/Active/IrenaeusChan/MobaV_FullSet/testing/combined_results/MobaSeq_SpikeInInfo.csv")
spike_in_info <- fread("/Volumes/ruit/Active/IrenaeusChan/MobaV_FullSet/input_files/B16keyAX01.csv")

# This sample only have 7 spike-ins, we know we have 8
spike_in_reads %>% filter(Sample_ID == "73326-WB")
sum(spike_in_reads %>% filter(Sample_ID == "73326-WB") %>% pull(count))
# Round to the nearest 100000
round(502230, -6)
round(502230)

spike_in_perc_count <- spike_in_reads %>%
    group_by(Sample_ID) %>%
    mutate(
        perc_count = count / sum(count)
    )

spike_in_quant <- spike_in_perc_count %>%
    ungroup() %>%
    group_by(name) %>%
    summarise(
        quant_5 = quantile(perc_count, 0.025),
        quant_95 = quantile(perc_count, 0.975)
    )

pdf("SpikeIn_DensityPlot_Test.pdf")
for (spikein in spike_in_info$name) {
    print(spikein)
    print(spike_in_perc_count %>%
        filter(name == spikein) %>%
        ggplot(aes(x = perc_count)) +
        geom_density() +
        geom_vline(data = spike_in_info %>% filter(name == spikein), aes(xintercept = expected), linetype = "dashed", color = "blue") +
        geom_text(data = spike_in_info %>% filter(name == spikein), aes(x = expected, y = 0.8, label = sprintf("Expected (%.3f)", expected)), color = "blue", angle = 90, vjust = -0.5) +
        geom_vline(data = spike_in_quant %>% filter(name == spikein), aes(xintercept = quant_5), linetype = "dashed", color = "red") +
        geom_text(data = spike_in_quant %>% filter(name == spikein), aes(x = quant_5, y = 0.8, label = sprintf("2.5%% (%.3f)", quant_5)), color = "red", angle = 90, vjust = -0.5) +
        geom_vline(data = spike_in_quant %>% filter(name == spikein), aes(xintercept = quant_95), linetype = "dashed", color = "red") +
        geom_text(data = spike_in_quant %>% filter(name == spikein), aes(x = quant_95, y = 0.8, label = sprintf("97.5%% (%.3f)", quant_95)), color = "red", angle = 90, vjust = -0.5) +
        theme_bw() +
        labs(
            title = paste0("Density Plot of ", spikein, " Spike-In"),
            x = "Spike-In Relative Read Counts",
            y = "Density"
        ))
}
dev.off()

# After changing the Spike-In Regression Line
data_to_plot <- full_join(andy_df, moba_new_df, by = c("unique_id", "Sample_ID", "Tissue"))
data_to_plot <- data_to_plot %>%
    mutate(across(where(is.numeric), ~ replace_na(., 0))) %>%
    mutate(
        Missing = ifelse(AX_cell_num == 0 | IC_cell_num == 0, "Missing", "Present"),
        Missing = as.factor(Missing)
    )
# Calculate regression on filtered data
reg_data <- data_to_plot %>%
    filter(AX_cell_num > 0, IC_cell_num > 0)

lm_fit <- lm(log10(IC_cell_num) ~ log10(AX_cell_num), data = reg_data)
r2 <- summary(lm_fit)$r.squared

# Plot
data_to_plot %>%
    mutate(
        # I'm highlighting these ones because they were outliers during the calculation of the regression (see above)
        sample_label = case_when(
            Sample_ID == "73210-Liver" ~ Sample_ID,
            Sample_ID == "73213-Liver" ~ Sample_ID,
            Sample_ID == "73400-Liver" ~ Sample_ID,
            Sample_ID == "73207-Liver" ~ Sample_ID,
            TRUE ~ "Other"
        )
    ) %>%
    ggplot(
        aes(
            x = AX_cell_num,
            y = IC_cell_num,
            color = sample_label
        )
    ) +
    geom_point() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    geom_smooth(data = reg_data, method = "lm", formula = y ~ x, color = "blue", se = FALSE) +
    scale_x_log10(
        limits = c(0.001, NA),
        breaks = 10^seq(-1, 10, 1), # Creates 0.1, 1, 10, 100, 1000
        labels = scales::label_comma()
    ) +
    scale_y_log10(
        limits = c(0.001, NA),
        breaks = 10^seq(-1, 10, 1),
        labels = scales::label_comma()
    ) +
    theme_bw() +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        # legend.position = "none"
    ) +
    labs(
        x = "AndyXu Cell Number",
        y = "MobaSeq Cell Number"
    )

# Calculate spike-in ratio as in either:
#   1) spike-in reads/total reads
#   2) spike-in cell number/ total gDNA (concentration * volumn from Isabella’s gDNA info

total_reads <- fread("/Volumes/ruit/Active/IrenaeusChan/MobaV_FullSet/outputs/total_reads_per_sgid.csv")
spike_in_reads <- fread("/Volumes/ruit/Active/IrenaeusChan/MobaV_FullSet/outputs/MobaV_Full_SpikeInInfo.csv")
total_dna <- fread("/Volumes/ruit/Active/IrenaeusChan/MobaV_FullSet/outputs/total_dna_per_mouse.csv")

total_reads_per_sample <- total_reads %>%
    select(-sgID) %>%
    summarise(across(everything(), sum, na.rm = TRUE)) %>%
    pivot_longer(everything(), names_to = "Sample_ID", values_to = "total_reads")

data_df <- spike_in_reads %>%
    select(Sample_ID, name, count, expected_cellnum) %>%
    left_join(
        total_reads_per_sample,
        by = "Sample_ID"
    ) %>%
    left_join(
        rsq_df,
        by = "Sample_ID"
    ) %>%
    left_join(
        counting_table,
        by = c("Sample_ID" = "Sample")
    ) %>%
    filter(Sample_ID != "preTran") %>%
    left_join(
        total_dna,
        by = c("Mouse_ID" = "Mouse_ID")
    ) %>%
    mutate(
        spike_in_ratio = count / total_reads,
        spike_in_ratio_cellnum = expected_cellnum / Total_DNA
    )


select(Sample_ID, name, count, expected_cellnum, total_reads, Spike_Rsq, spike_in_ratio, spike_in_ratio_cellnum) %>%
    # filter(Sample_ID == "73320-WB")
    filter(Sample_ID %in% c("73210-Liver", "73400-Liver", "73213-Liver", "73207-Liver"))

library(ggExtra)

quant_df <- data_df %>%
    group_by(Tissue) %>%
    summarise(
        quant_rsq = quantile(Spike_Rsq, 0.05),
        quant_spike = quantile(spike_in_ratio, 0.05)
    )

# Create main plot
p <- data_df %>%
    filter(!(Sample_ID %in% c("73210-Liver", "73400-Liver", "73213-Liver", "73207-Liver"))) %>%
    ggplot(
        aes(
            x = Spike_Rsq,
            y = spike_in_ratio,
            label = Sample_ID,
            color = Tissue
        )
    ) +
    geom_point() +
    geom_text_repel(data = quant_df, aes(x = 0.85, y = quant_spike, label = sprintf("%.5f", quant_spike)), hjust = 0, vjust = -1) +
    # geom_text_repel(size = 3, max.overlaps = 20) +
    # geom_vline(data = quant_df, aes(xintercept = quant_rsq, color = Tissue), linetype = "dashed") +
    geom_hline(data = quant_df, aes(yintercept = quant_spike, color = Tissue), linetype = "dashed") +
    scale_y_log10(labels = scales::comma) +
    theme_bw() +
    theme(legend.position = "bottom") +
    labs(
        x = "Spike-In R²",
        y = "Spike-In Ratio"
    )

# Add marginal plots with tissue-specific percentiles
ggMarginal(p,
    type = "density",
    groupColour = TRUE,
    groupFill = TRUE,
    yparams = list(vline = quant_df$quant_spike)
)


ggplot(
    aes(
        x = Spike_Rsq,
        y = spike_in_ratio_cellnum,
        label = Sample_ID
    )
) +
    geom_point() +
    geom_text_repel() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    theme_bw()