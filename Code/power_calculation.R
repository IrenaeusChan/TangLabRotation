# This is the code used to attempt to compute necessary power for the Met Burden analysis

# Author: Irenaeus Chan
# Date: 03/20/2025

# Load the necessary libraries
library(data.table)
library(tidyverse)
library(readxl)
library(boot)
library(ggrepel)

source("TangLabRotation/Code/help_functions.R")

# Read in the data
df <- fread("/Volumes/bolton/Active/Users/IrenaeusChan/Moba500/out_dir/FilteredSampleInfo.csv")
df <- df[nchar(df$sgID) <= 10, ]
df <- df %>% filter(!is.na(distance))
df <- df %>% filter(distance == 0)
df <- df %>% filter(sgID != "sgDummy")

# define sgSafe and other sgIDs
sgSafe <- c(
    "sgNeo1", "sgNeo2", "sgNeo3", "sgNT1", "sgNT2", "sgNT3",
    "sgSafe1", "sgNT4", "sgSafe10", "sgSafe11", "sgSafe13",
    "sgSafe15", "sgSafe16", "sgSafe17", "sgSafe18", "sgSafe19",
    "sgSafe27", "sgSafe28", "sgSafe32", "sgSafe33", "sgSafe34",
    "sgSafe36", "sgSafe37", "sgSafe38", "sgSafe39", "sgSafe40",
    "sgSafe41", "sgSafe42", "sgSafe2", "sgSafe23", "sgSafe24",
    "sgSafe4", "sgSafe5", "sgSafe6", "sgSafe7", "sgSafe8",
    "sgSafe9", "sgSafe25", "sgSafe26"
)

sgAll <- c(unique(df$sgID))
sg_notSafe <- sgAll[!(sgAll %in% sgSafe)]
sgNfib <- c("sgNf1b-2", "sgNf1b-4", "sgNf1b-5")

# add gene name, safe and nt as control
df$gene <- ifelse(df$sgID %in% sgSafe, "Ctrl", df$sgID)
df$gene <- ifelse(df$gene %in% sgNfib, "Nfib", df$gene)

# 1. Get the data for Crebbp in Liver
filtered_df <- df %>% filter(Tissue == "Liver", Mouse_Genotype == "C57")
safe_df <- df %>% filter(gene == "Ctrl", Tissue == "Liver", Mouse_Genotype == "C57")
message("Total colonies for Crebbp: ", nrow(filtered_df))
preTran_df <- df %>% filter(Tissue == "preTran")

# 2. Identify which sgIDs were in the Pre-Transplantation Sample, if it wasn't in Pre-Transplantation, you can't do analysis
preTran_cell_num <- aggregate(cell_num ~ sgID + Sample_ID + gene, data = preTran_df, FUN = sum)
# Calculate the 1st and 99th percentiles
ci_l <- quantile(preTran_cell_num$cell_num, 0.01, na.rm = TRUE)
ci_r <- quantile(preTran_cell_num$cell_num, 0.99, na.rm = TRUE)
preTran_cutoff <- 5 * ci_l
sgAll <- unique(preTran_cell_num$sgID[preTran_cell_num$cell_num > preTran_cutoff])
message("There are ", length(sgAll), " sgRNAs that are present in the Pre-Transplantation Sample above the 1st percentile")

# 4. Calculate the average Met Burden for the safe sgRNAs from the Met Burden data frame
avg_MetBurden_safe <- safe_df %>%
    filter(sgID %in% sgAll) %>%
    filter(gene == "Ctrl") %>%
    group_by(sgID) %>%
    summarize(Total_Cell_Num = sum(cell_num)) %>%
    left_join(
        preTran_cell_num %>% select(sgID, gene, cell_num) %>% dplyr::rename(PreTran_Cell_Num = cell_num),
        by = "sgID"
    ) %>%
    mutate(MetBurden = Total_Cell_Num / PreTran_Cell_Num) %>%
    summarize(avg_MetBurden_safe = mean(MetBurden, na.rm = TRUE))

crebbp_df <- filtered_df %>% filter(sgID == "Crebbp") # Using this as a testing dataframe
calc_MetBurden(crebbp_df, safe_df, "Crebbp", preTran_cell_num)

# Function to test single sample size
test_sample_size <- function(n_colonies, gene_df, safe_df, target, preTran_cell_num, avg_MetBurden_safe) {
    sampled_df <- gene_df %>% slice_sample(n = n_colonies)
    PreTran_Cell_Num <- preTran_cell_num$cell_num[preTran_cell_num$sgID == target]

    # Bootstrap test
    boot_result <- boot(sampled_df, function(d, i) {
        sum(d$cell_num[i]) / PreTran_Cell_Num / avg_MetBurden_safe$avg_MetBurden_safe
    }, R = 10)

    # Get p-value
    pval <- getEmpiricalP(boot_result$t, Baseline = 1, Side = 2)
    FDR_pval <- p.adjust(pval, method = "BH")

    cat(sprintf("n_colonies: %d, bootstrap t0: %.3f, and FDR p-value: %.3f\n", n_colonies, boot_result$t0, FDR_pval))

    return(log2(boot_result$t0) > 0 & FDR_pval < 0.05)
}
# Example to how to use code:
test_sample_size(323, crebbp_df, safe_df, "Crebbp", preTran_cell_num, avg_MetBurden_safe)

# Function to test single sample mixed with safe sgRNAs
test_sample_size_mixed <- function(n_colonies, gene_df, safe_df, target, preTran_cell_num, avg_MetBurden_safe) {
    # Sample colonies
    #sampled_df <- gene_df %>%
    #    slice_sample(n = n_colonies) %>%
    #    group_by(sgID) %>%
    #    summarize(Total_Cell_Num = sum(cell_num))

    sampled_df <- rbind(
        safe_df %>%
            slice_sample(n = n_colonies),
        gene_df
    )

    PreTran_Cell_Num <- preTran_cell_num %>%
        filter(sgID %in% sampled_df$sgID, sgID != target) %>%
        group_by(sgID) %>%
        summarize(pre_tran_cell_num = cell_num) %>%
        left_join(
            safe_df %>%
                filter(sgID %in% sampled_df$sgID) %>%
                group_by(sgID) %>%
                summarize(n_colonies = n())
        ) %>%
        mutate(pre_tran_cell_num = pre_tran_cell_num / n_colonies) %>%
        summarise(cell_num = sum(pre_tran_cell_num, na.rm = TRUE))
    
    PreTran_Cell_Num = sum(PreTran_Cell_Num$cell_num, preTran_cell_num$cell_num[preTran_cell_num$sgID == target])

    # Bootstrap test
    boot_result <- boot(sampled_df, function(d, i) {
        sum(d$cell_num[i]) / PreTran_Cell_Num / avg_MetBurden_safe$avg_MetBurden_safe
    }, R = 10)

    # Get p-value
    pval <- getEmpiricalP(boot_result$t, Baseline = 1, Side = 2)
    FDR_pval <- p.adjust(pval, method = "BH")
    
    cat(sprintf("%s -- n_colonies: %d, bootstrap t0: %.3f, and FDR p-value: %.3f\n", target, n_colonies, boot_result$t0, FDR_pval))

    res <- c(sum(boot_result$data$cell_num), PreTran_Cell_Num, log2(boot_result$t0))
    return(res)
    #return(log2(boot_result$t0))
    if (FDR_pval < 0.05) {
        #res <- c(met_burden$Total_Cell_Num, PreTran_Cell_Num, log2(boot_result$t0))
        return(log2(boot_result$t0))
        #return(res)
    } else {
        return(NA)
    }
    # return(log2(boot_result$t0) < 0 & p.adjust(pval, method = "BH") < 0.05)
}
#results <- test_sample_size_mixed(299, filtered_df %>% filter(sgID == "Crebbp"), safe_df, "Crebbp", preTran_cell_num, avg_MetBurden_safe)

# Function to run power analysis for one gene
run_power_analysis <- function(target_gene, df, safe_df, preTran_cell_num, avg_MetBurden_safe) {
    print(target_gene)
    # Get gene specific data
    gene_df <- df %>%
        filter(sgID == target_gene)

    # Skip if no data
    if (nrow(gene_df) == 0) {
        return(NULL)
    }

    # Calculate power across sample sizes
    results <- data.frame(
        gene = target_gene,
        n_colonies = seq(1, nrow(gene_df), by = 1),
        power = NA
    )

    for (i in 1:nrow(results)) {
        n <- results$n_colonies[i]
        success <- replicate(100, test_sample_size(n, gene_df, safe_df, target_gene, preTran_cell_num, avg_MetBurden_safe))
        results$power[i] <- mean(success)
    }

    return(results)
}

run_power_analysis_met_burden <- function(target_gene, df, safe_df, preTran_cell_num, avg_MetBurden_safe) {
    print(target_gene)
    gene_df <- df %>%
        filter(sgID == target_gene)

    if (nrow(gene_df) == 0) {
        return(NULL)
    }

    # n_colonies <- seq(1, nrow(gene_df), by = 5)
    n_colonies <- seq(1, nrow(gene_df), by = 5)
    all_results <- map_dfr(n_colonies, function(n) {
        success <- replicate(100, test_sample_size(n, gene_df, safe_df, target_gene, preTran_cell_num, avg_MetBurden_safe))
        tibble(
            gene = target_gene,
            n_colonies = n,
            iteration = 1:100,
            # met_burden = list(success) # Nest met_burden as list column
            total_cell_num = list(success[1]),
            pre_tran_cell_num = list(success[2]),
            met_burden = list(success[3])
        )
    })

    return(all_results)
}

run_power_analysis_met_burden_mixed <- function(target_gene, df, safe_df, preTran_cell_num, avg_MetBurden_safe) {
    print(target_gene)
    gene_df <- df %>%
        filter(sgID == target_gene)

    if (nrow(gene_df) == 0) {
        return(NULL)
    }

    #n_colonies <- seq(1, nrow(gene_df), by = 5)
    n_colonies <- seq(1, nrow(gene_df), by = 5)
    all_results <- map_dfr(n_colonies, function(n) {
        success <- replicate(100, test_sample_size_mixed(n, gene_df, safe_df, target_gene, preTran_cell_num, avg_MetBurden_safe))
        tibble(
            gene = target_gene,
            n_colonies = n,
            iteration = 1:100,
            #met_burden = list(success) # Nest met_burden as list column
            total_cell_num = list(success[1]),
            pre_tran_cell_num = list(success[2]),
            met_burden = list(success[3])
        )
    })

    return(all_results)
}
# Example to run MetBurden analysis for mixing safe sgRNAs
test_results <- run_power_analysis_met_burden_mixed("Crebbp", filtered_df, safe_df, preTran_cell_num, avg_MetBurden_safe)
test_results %>%
    unnest(cols = c(total_cell_num, pre_tran_cell_num, met_burden)) %>%
    as.data.frame() %>%
    select(met_burden) %>%
    max()
test_results %>%
    unnest(cols = c(total_cell_num, pre_tran_cell_num, met_burden)) %>%
    #unnest(cols = c(met_burden)) %>%
    as.data.frame() %>%
    mutate(
        total_cell_num = total_cell_num/100000000,
        pre_tran_cell_num = pre_tran_cell_num/100000,
    ) %>%
    ggplot(
        aes(
            x = n_colonies, group = n_colonies
        )
    ) +
    geom_boxplot(aes(y = met_burden)) +
    geom_boxplot(aes(y = total_cell_num), color = "red") +
    geom_boxplot(aes(y = pre_tran_cell_num), color = "blue") +
    theme_bw() +
    labs(
        x = "Number of Colonies",
        y = "Met Burden",
        title = "Met Burden Distribution by Colony Sampling"
    ) +
    #scale_y_continuous(breaks = seq(-30, 10, by = 1)) +
    scale_x_continuous(breaks = seq(0, 5000, by = 100)) +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1)
    )

test_results %>%
    unnest(cols = c(total_cell_num, pre_tran_cell_num, met_burden)) %>%
    as.data.frame() %>%
    group_by(n_colonies) %>%
    summarise(
        mean_met_burden = mean(met_burden),
        mean_total = mean(total_cell_num / 10000000),
        mean_pretran = mean(pre_tran_cell_num / 10000)
    ) %>%
    ggplot(aes(x = n_colonies)) +
    geom_line(aes(y = mean_met_burden), color = "black") +
    geom_line(aes(y = mean_total), color = "red") +
    geom_line(aes(y = mean_pretran), color = "blue") +
    theme_bw() +
    labs(
        x = "Number of Colonies",
        y = "Met Burden",
        title = "Met Burden Distribution by Colony Sampling"
    ) +
    scale_x_continuous(breaks = seq(0, 5000, by = 100)) +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1)
    )

# Run analysis for all genes
all_results <- map_df(
    #sgAll[!(sgAll %in% sgSafe)],
    c("Crebbp", "Pten", "Tsc2", "Ewsr1", "Apc", "Gata4", "Nae1", "Rac1"),
    ~ run_power_analysis(., filtered_df, safe_df, preTran_cell_num, avg_MetBurden_safe)
)

all_results %>%
    group_by(gene) %>%
    summarise(
        max_power = max(power),
        required_colonies = if (max_power >= 0.8) {
            n_colonies[which.min(power >= 0.8)]
        } else {
            n_colonies[which.max(n_colonies)]
        },
        achieved_power = if (max_power >= 0.8) {
            0.8
        } else {
            max_power
        }
    ) %>%
    arrange(desc(required_colonies))

all_results %>%
    group_by(gene) %>%
    summarise(
        max_power = max(power),
        required_colonies = if (max_power >= 0.8) {
            first_pass <- which(power >= 0.8)[1] # Get first position where power >= 0.8
            if (!is.na(first_pass)) n_colonies[first_pass] else max(n_colonies)
        } else {
            max(n_colonies)
        },
        achieved_power = power[which(n_colonies == required_colonies)]
    ) %>%
    mutate(achieved_threshold = achieved_power >= 0.8) %>%
        arrange(
            desc(achieved_threshold),
            if_else(achieved_threshold, desc(required_colonies), -Inf),
            desc(achieved_power)
        ) %>%
        fwrite("power_analysis_results.csv")


# Crebbp, Tsc2, Ewsr1, Apc, Gata4, Nae1, Rac1
all_results %>%
    mutate(
        label = case_when(
            gene == "Crebbp" ~ "Crebbp",
            gene == "Pten" ~ "Pten",
            gene == "Tsc2" ~ "Tsc2",
            gene == "Ewsr1" ~ "Ewsr1",
            gene == "Apc" ~ "Apc",
            gene == "Gata4" ~ "Gata4",
            gene == "Nae1" ~ "Nae1",
            gene == "Rac1" ~ "Rac1",
            TRUE ~ NA
        )
    ) %>%
filter(gene %in% c("Pten", "Crebbp", "Apc", "Tsc2")) %>%
ggplot(
    aes(
        x = n_colonies, y = power, color = gene, #label = label
    )
) +
geom_line() +
geom_hline(yintercept = 0.8, linetype = "dashed") +
xlim(10, 300) +
#scale_x_log10() +
theme_bw() +
labs(
    x = "Number of Colonies",
    y = "Power",
    title = "Power Analysis by Gene"
) +
theme(
    legend.position = "right"
)