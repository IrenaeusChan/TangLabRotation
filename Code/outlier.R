# This is the code used to attempt to identify outliers within the data using PCA and Read Counts

# Author: Irenaeus Chan
# Date: 03/20/2025

# Load the necessary libraries
library(data.table)
library(tidyverse)
library(DESeq2)
library(readxl)
library(ggrepel)

metadata <- read_excel("/Volumes/ruit/Active/IrenaeusChan/MobaV_FullSet/input_files/CountingTable.xlsx")
df <- fread("/Volumes/ruit/Active/IrenaeusChan/MobaV_FullSet/output_new_outlier/MobaSeq_FilteredSampleInfo.csv")
df <- df %>% filter(!is.na(distance))
df <- df %>% filter(distance == 0)
df <- df %>% filter(sgID != "sgDummy")

# There were 2 Cell Lines used in the experiment
df <- df %>%
    filter(
        BC_end == "GGAA" | BC_end == "CTGA"
    ) %>%
    mutate(
        Cell_Line = ifelse(BC_end == "GGAA", "rp48", "rp116")
    )

# define sgSafe and other sgIDs
sgSafe <- c(
    "sgNT2-MMU", "sgNT3-MMU",
    "sgSafe19-MMU", "sgSafe29-MMU", "sgSafe36-MMU", "sgSafe4-MMU", "sgSafe28-MMU", "sgSafe30-MMU", "sgSafe33-MMU", "sgSafe5-MMU"
)
sgAll <- c(unique(df$sgID))
sg_notSafe <- sgAll[!(sgAll %in% sgSafe)]
sgNfib <- c("sgNf1b-2", "sgNf1b-4", "sgNf1b-5")

sgAll <- c(unique(df$sgID))
sg_notSafe <- sgAll[!(sgAll %in% sgSafe)]
sgNfib <- c("Nfibmmu-sg1", "Nfibmmu-sg2", "Nfibmmu-sg3")

# add gene name, safe and nt as control
df$gene <- ifelse(df$sgID %in% sgSafe, "Ctrl", df$sgID)
df$gene <- ifelse(df$gene %in% sgNfib, "Nfib", df$gene)

# In RNASeq, raw_counts is defined by how many reads are mapped to a gene
all_samples_df <- df %>%
    filter(Sample_ID != "preTran") %>%
    filter(Time_Point == "3weeks") %>%
    filter(Tissue == "Liver" | Tissue == "Lung") %>%
    group_by(Mouse_ID, sgID) %>%
    summarize(
        total_counts = sum(count)
    ) %>%
    ungroup() %>%
    arrange(Mouse_ID) %>%
    spread(Mouse_ID, total_counts, fill = 0) %>%
    as.data.frame() %>% # Convert from tibble
    column_to_rownames("sgID")

all_samples_sg_df <- df %>%
    filter(Sample_ID != "preTran") %>%
    filter(Time_Point == "3weeks") %>%
    filter(Tissue == "Liver" | Tissue == "Lung") %>%
    filter(gene == "Ctrl") %>%
    group_by(Mouse_ID, sgID) %>%
    summarize(
        total_counts = sum(count)
    ) %>%
    ungroup() %>%
    arrange(Mouse_ID) %>%
    spread(Mouse_ID, total_counts, fill = 0) %>%
    as.data.frame() %>% # Convert from tibble
    column_to_rownames("sgID")

coldata_all <- metadata %>%
    filter(Sample != "preTran") %>%
    filter(Time_Point == "3weeks") %>%
    filter(Tissue == "Liver" | Tissue == "Lung") %>%
    arrange(Mouse_ID) %>%
    select(Mouse_ID, Mouse_Genotype, Sex) %>%
    distinct() %>%
    as.data.frame() %>%
    column_to_rownames("Mouse_ID") # Use appropriate column name

# 1. Create DESeq2 objects
dds_all <- DESeqDataSetFromMatrix(
    countData = all_samples_df,
    colData = coldata_all,
    design = ~ 1
)
dds_all <- estimateSizeFactors(dds_all)
dds_all <- estimateDispersions(dds_all)

# Filter out rows with no counts
#all_samples_sg_df_filtered <- all_samples_sg_df[rowSums(all_samples_sg_df) > 0, ] + 1

dds_safe <- DESeqDataSetFromMatrix(
    countData = all_samples_sg_df,
    colData = coldata_all,
    design = ~1
)
dds_safe <- estimateSizeFactors(dds_safe)
dds_safe <- estimateDispersions(dds_safe)

# 2. Transform data (VST or rlog)
rlog_all <- rlog(dds_all)
rlog_safe <- rlog(dds_safe)

# 3. Perform PCA
pca_all <- prcomp(t(assay(rlog_all)), scale = TRUE)
pca_safe <- prcomp(t(assay(rlog_safe)), scale = TRUE)

# 4. Calculate percent variance
percentVar_all <- round(100 * pca_all$sdev^2 / sum(pca_all$sdev^2), 1)
percentVar_safe <- round(100 * pca_safe$sdev^2 / sum(pca_safe$sdev^2), 1)

# 5. Create plotting data
pca_data_all <- data.frame(
    PC1 = pca_all$x[, 1],
    PC2 = pca_all$x[, 2],
    mouse_id = colnames(all_samples_df)
) %>%
left_join(
    metadata %>% select(Sample, Mouse_ID, Tissue, Mouse_Genotype) %>% mutate(Mouse_ID = as.character(Mouse_ID)),
    by = c("mouse_id" = "Mouse_ID")
)

pca_data_safe <- data.frame(
    PC1 = pca_safe$x[, 1],
    PC2 = pca_safe$x[, 2],
    mouse_id = colnames(all_samples_sg_df)
) %>%
    left_join(
        metadata %>% select(Sample, Mouse_ID, Tissue, Mouse_Genotype) %>% mutate(Mouse_ID = as.character(Mouse_ID)),
        by = c("mouse_id" = "Mouse_ID")
    )

# 6. Plot the PCA for all sgIDs
# pdf("PCA_sgAll_Genotype.pdf")
print(ggplot(
    pca_data_all %>% 
        distinct(PC1, PC2, Mouse_Genotype, mouse_id),
    aes(
        x = PC1, 
        y = PC2, 
        color = Mouse_Genotype, 
        label = mouse_id
    )) +
    geom_text_repel(hjust = 0, vjust = -1, max.overlap = 100) +
    geom_point(size = 3) +
    xlab(paste0("PC1: ", percentVar_all[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar_all[2], "% variance")) +
    labs(
        title = "PCA of All sgIDs",
        color = "Genotype"
    ) +
    theme_bw())

#pdf("PCA_sgSafes_Genotype.pdf")
print(ggplot(pca_data_safe %>% distinct(PC1, PC2, Mouse_Genotype, mouse_id), aes(x = PC1, y = PC2, color = Mouse_Genotype, label = mouse_id)) +
    geom_text_repel(hjust = 0, vjust = -1) +
    geom_point(size = 3) +
    xlab(paste0("PC1: ", percentVar_safe[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar_safe[2], "% variance")) +
    labs(
        title = "PCA of sgSafes",
        color = "Genotype"
    ) +
    theme_bw())

# 7. For each PCA Component, calculate the z-scores
pc_zscores_all <- apply(pca_all$x, 2, scale)
colnames(pc_zscores_all) <- colnames(pca_all$x)
rownames(pc_zscores_all) <- colnames(all_samples_df)

pc_zscores_safe <- apply(pca_safe$x, 2, scale)
colnames(pc_zscores_safe) <- colnames(pca_safe$x)
rownames(pc_zscores_safe) <- colnames(all_samples_sg_df)

n_pcs_all <- which(cumsum(percentVar_all) >= 70)[1]
n_pcs_safe <- which(cumsum(percentVar_safe) >= 70)[1]

# 8. Calculate the Euclidean distance for each row
calculate_stable_mahalanobis <- function(x, n_pcs) {
    # Filter PCs
    x_filtered <- x[, 1:n_pcs]
    
    # Calculate mean and covariance
    center <- colMeans(x_filtered)
    cov_mat <- cov(x_filtered) + diag(1e-6, n_pcs)
    
    # Return mahalanobis distances
    sqrt(mahalanobis(x_filtered, center, cov_mat))
}

test_idea_1 <- pc_zscores_all %>%
    as.data.frame() %>%
    rownames_to_column("Mouse_ID") %>%
    mutate(
        euc_dist = sqrt(rowSums(pca_all$x^2)),
        euc_dist_scaled = scale(euc_dist),
        mahal_dist = calculate_stable_mahalanobis(pca_all$x, n_pcs_all),
        mahal_dist_scaled = scale(mahal_dist),
        euclidean_zscore = sqrt(rowSums(pc_zscores_all^2)),
        rms_zscore = sqrt(rowMeans(pc_zscores_all^2)),
        weighted_zscore = sqrt(rowSums((pc_zscores_all^2) * (percentVar_all / 100)))
    ) %>%
    left_join(
        metadata %>% select(Sample, Mouse_ID, Tissue, Mouse_Genotype, Sex) %>% mutate(Mouse_ID = as.character(Mouse_ID)),
        by = "Mouse_ID"
    )
test_idea_2 <- pc_zscores_safe %>%
    as.data.frame() %>%
    rownames_to_column("Mouse_ID") %>%
    mutate(
        euc_dist = sqrt(rowSums(pca_safe$x^2)),
        euc_dist_scaled = scale(euc_dist),
        mahal_dist = calculate_stable_mahalanobis(pca_safe$x, n_pcs_safe),
        mahal_dist_scaled = scale(mahal_dist),
        euclidean_zscore = sqrt(rowSums(pc_zscores_safe^2)),
        rms_zscore = sqrt(rowMeans(pc_zscores_safe^2)),
        weighted_zscore = sqrt(rowSums((pc_zscores_safe^2) * (percentVar_safe / 100))),
    ) %>%
    left_join(
        metadata %>% select(Sample, Mouse_ID, Tissue, Mouse_Genotype, Sex) %>% mutate(Mouse_ID = as.character(Mouse_ID)),
        by = "Mouse_ID"
    )

pca_with_centroids_all <- pca_data_all %>%
    group_by(Mouse_Genotype) %>%
    mutate(
        centroid_x = mean(PC1),
        centroid_y = mean(PC2),
        dist_to_centroid = sqrt((PC1 - centroid_x)^2 + (PC2 - centroid_y)^2),
        dist_to_centroid_scaled = scale(dist_to_centroid)
    )
pca_with_centroids_safe <- pca_data_safe %>%
    group_by(Mouse_Genotype) %>%
    mutate(
        centroid_x = mean(PC1),
        centroid_y = mean(PC2),
        dist_to_centroid = sqrt((PC1 - centroid_x)^2 + (PC2 - centroid_y)^2),
        dist_to_centroid_scaled = scale(dist_to_centroid)
    )

data.frame(
    x = pca_with_centroids_all$dist_to_centroid,
    y = pca_with_centroids_safe$dist_to_centroid,
    label = pca_with_centroids_all$mouse_id
) %>% unique

# We tried a bunch of ideas, in the end we did not decide on a specific one
# However, exploring the euclidean distance between the two PCA components seems to have some potential

data.frame(
    x = test_idea_1$euc_dist_scaled,
    y = test_idea_2$euc_dist_scaled,
    label = test_idea_1$Mouse_ID
) %>%
    ggplot(
        aes(
            x = x,
            y = y,
            label = label,
            #color = genotype
        )
    ) +
    geom_abline(intercept = 0, slope = 1, linetype = "solid") +
    geom_abline(intercept = 1, slope = 1, linetype = "dashed", alpha = 0.5) +
    geom_abline(intercept = -1, slope = 1, linetype = "dashed", alpha = 0.5) +
    geom_smooth(
        method = "lm",
        se = TRUE,         # show confidence interval
        color = "blue",
        linetype = "solid",
        size = 1
    ) +
    geom_point() +
    geom_text(vjust = -1) +
    labs(
        x = "All Genes",
        y = "Safe Genes"
    ) +
    theme_bw()