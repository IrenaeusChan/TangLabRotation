library(tidyverse)
library(data.table)
library(mixtools)

source("RT_Code/help_functions.R")

#metadata <- read_excel("/Volumes/ruit/Active/IrenaeusChan/MobaV_FullSet/input_files/CountingTable.xlsx")
df <- fread("/Users/irenaeuschan/Documents/Irenaeus/TEST_MOBASEQ/2024_Moba1stDraft_CellCount.csv")

df <- df %>%
    rename(
        distance = dist,
        cell_num = Cell_Num,
        barcode = BC,
        count = Count,
        reading_depth = Reading_Depth
    )

# df <- df[nchar(df$sgID) <= 10, ]
df <- df %>% filter(!is.na(distance))
df <- df %>% filter(distance == 0)
df <- df %>% filter(sgID != "sgDummy")
df <- df %>% filter(cell_num > 100)

sgAll <- c(unique(df$sgID))
# define sgSafe and other sgIDs
sgSafe <- c(sgAll[grepl("sgSafe", sgAll)])
sgNfib <- c(sgAll[grepl("sgNf1b", sgAll)])

sg_notSafe <- sgAll[!(sgAll %in% sgSafe)]


# add gene name, safe and nt as control
df$gene <- ifelse(df$sgID %in% sgSafe, "Ctrl", df$sgID)
df$gene <- ifelse(df$gene %in% sgNfib, "Nfib", df$gene)

df <- df %>%
    mutate(
        gene = ifelse(grepl("sgSafe", sgID), "Ctrl", sgID),
        gene = ifelse(grepl("sgNf1b", sgID), "Nfib", gene)
    )

sgSafe_df_2days <- df %>% filter(gene == "Ctrl", Time_Point == "2days", Tissue == "Liver")
sgSafe_df_1week <- df %>% filter(gene == "Ctrl", Time_Point == "1week", Tissue == "Liver")
sgSafe_df_2week <- df %>% filter(gene == "Ctrl", Time_Point == "2weeks", Tissue == "Liver")
sgSafe_df_3week <- df %>% filter(gene == "Ctrl", Time_Point == "3weeks", Tissue == "Liver")

res_2days <- check_bimodal(sgSafe_df_2days$cell_num)
res_1week <- check_bimodal(sgSafe_df_1week$cell_num)
res_2week <- check_bimodal(sgSafe_df_2week$cell_num)
res_3week <- check_bimodal(sgSafe_df_3week$cell_num)

plot_bimodal(sgSafe_df_2days$cell_num, res_2days)
plot_bimodal(sgSafe_df_1week$cell_num, res_1week)
plot_bimodal(sgSafe_df_2week$cell_num, res_2week)
plot_bimodal(sgSafe_df_3week$cell_num, res_3week)

# One big plot
sgSafe_df <- df %>% filter(gene == "Ctrl", Sample_ID != "preTran1", Tissue == "Liver")

peaks <- c(res_2days$peaks[1:2], res_1week$peaks[1:2], res_2week$peaks[1:2], res_3week$peaks[1:2])
peak_heights <- c(res_2days$peak_heights[1:2], res_1week$peak_heights[1:2], res_2week$peak_heights[1:2], res_3week$peak_heights[1:2])
valleys <- c(res_2days$valley[1], res_1week$valley[1], res_2week$valley[1], res_3week$valley[1])
valley_heights <- c(res_2days$valley_heights[1], res_1week$valley_heights[1], res_2week$valley_heights[1], res_3week$valley_heights[1])


sgSafe23_df %>%
    mutate(
        Log2FC = log2(cell_num),
        Time_Point = factor(Time_Point, levels = c("2days", "1week", "2weeks", "3weeks"))
    ) %>%
    ggplot() +
    geom_density(
        aes(
            x = Log2FC, 
            fill = Time_Point
        ), 
        alpha = 0.3
    ) +
    # geom_line(
    #     data = data.frame(
    #         x = seq(peaks[1] - 4 * sd1, peaks[1] + 4 * sd1, length.out = 1000)
    #     ) %>% mutate(
    #         y = dnorm(x, mean = peaks[1], sd = sd1) *
    #             (peak_heights[1] / dnorm(peaks[1], mean = peaks[1], sd = sd1))
    #     ),
    #     aes(x = x, y = y),
    #     color = "red"
    # ) +
    # geom_line(
    #     data = data.frame(
    #         x = seq(peaks[2] - 4 * sd2, peaks[2] + 4 * sd2, length.out = 1000)
    #     ) %>% mutate(
    #         y = dnorm(x, mean = peaks[2], sd = sd2) *
    #             (peak_heights[2] / dnorm(peaks[2], mean = peaks[2], sd = sd2))
    #     ),
    #     aes(x = x, y = y),
    #     color = "blue"
    # ) +
    geom_vline(xintercept = peaks, color = "red", linetype = "dashed") +
    geom_vline(xintercept = valleys, color = "blue", linetype = "dashed") +
    theme_minimal() +
    scale_x_continuous(
        breaks = seq(0, 30, 1),
        limits = c(0, 30)
    ) +
    labs(
        x = "log2(cell_num)",
        y = "Density"
    )




safe_df <- df %>% filter(gene == "Ctrl")
preTran_df <- df %>% filter(Sample_ID == "preTran1")







gene_list <- c("Crebbp", "Pten", "Csnk1a1", "Gpatch8", "Zeb1", "Tsc1", "Tsc2", "Gata6", "Mga", "Zmiz1", "Grhpr", "Hnf4a", "Eif5b", "Soga1", "Nup98", "Ctrl")
for (g in gene_list) {
    # Filter the data for NSG genotype and the specific gene
    df_gene <- df %>%
        filter(Mouse_Genotype == "C57") %>%
        filter(Time_Point == "3weeks") %>%
        filter(cell_num > 50) %>%
        filter(Tissue == "Liver") %>%
        filter(gene == g)


    # Create log2 transformed cell_num column - for fitting the mixture model
    df_gene$cell_num_log2 <- log2(df_gene$cell_num)

    # Select the data vector with log2 cell_num
    data_vector <- df_gene$cell_num_log2

    # Set the fixed means for Components 1 and 2
    Mode1 <- find_1st_mode(df_gene)
    Mode2 <- find_2nd_mode(df_gene)
    # mode_x is the x value (cell_num) of the mode (peak) in the density plot
    # mode_y is the y value (total colonies) of the mode (peak) in the density plot
    print(Mode2)

    # Convert the mode values to log2 scale
    fixed_mean1 <- log2(Mode1$mode_x[Mode1$gene == g])
    fixed_mean2 <- log2(Mode2$mode_x[Mode2$gene == g])

    # Fit a mixture model with constraints (fix means for both components)
    custom_mix_model <- normalmixEM(data_vector, k = 2, mean.constr = c(fixed_mean1, fixed_mean2))
    # Basically, we are fitting a Gaussian mixture model with 2 components and fixing the means of the components to the mode values (calculated above)

    # Get the posterior probabilities for each component
    posterior_probs <- custom_mix_model$posterior

    # Calculate the percentage of data in each component
    component1_percentage <- sum(posterior_probs[, 1]) / length(data_vector) * 100
    component2_percentage <- sum(posterior_probs[, 2]) / length(data_vector) * 100


    # Create a histogram with more bins and a dynamic title
    hist(data_vector,
        breaks = 30, prob = TRUE,
        main = paste("Fitted Gaussian Mixture Model for", g),
        xlab = "log2 cell_num", col = "grey", border = "black", alpha = 0.3
    )

    # Generate sequence of x-values for plotting the density curves
    x_vals <- seq(min(data_vector), max(data_vector), length.out = 1000)

    # Calculate the density for each Gaussian component
    # dnorm creates a normalized density curve for a given mean and standard deviation
    comp1 <- custom_mix_model$lambda[1] * dnorm(x_vals, mean = fixed_mean1, sd = custom_mix_model$sigma[1]) # Fixed mean for Component 1
    comp2 <- custom_mix_model$lambda[2] * dnorm(x_vals, mean = fixed_mean2, sd = custom_mix_model$sigma[2]) # Fixed mean for Component 2

    # Calculate the total density of the mixture model
    total_density <- comp1 + comp2

    # Plot the individual Gaussian components
    lines(x_vals, comp1, col = "darkblue", lwd = 2, lty = 2) # First component
    lines(x_vals, comp2, col = "red", lwd = 2, lty = 2) # Second component

    # Plot the total mixture density
    lines(x_vals, total_density, col = "purple", lwd = 3, lty = 1) # Total mixture

    # Add vertical lines for the mean values
    abline(v = fixed_mean1, col = "darkblue", lwd = 2, lty = 2) # Mean of the first component (fixed)
    abline(v = fixed_mean2, col = "red", lwd = 2, lty = 2) # Mean of the second component (fixed)

    # Add a legend at the top of the plot
    legend("topright",
        legend = c(
            paste("C1 (", round(fixed_mean1, 2), ")", sep = ""),
            paste("C2 (", round(fixed_mean2, 2), ")", sep = "")
        ),
        col = c("darkblue", "red"),
        lwd = 2,
        lty = c(2, 2),
        inset = c(0.02, 0.02)
    ) # Adjust 'inset' to move the legend within the plot area

    # Add percentages to the plot
    text(x = fixed_mean1, y = max(comp1) * 0.9, labels = paste("C1: ", round(component1_percentage, 2), "%", sep = ""), col = "darkblue")
    text(x = fixed_mean2, y = max(comp2) * 0.9, labels = paste("C2: ", round(component2_percentage, 2), "%", sep = ""), col = "red")

    # Output the adjusted percentages
    cat("Adjusted Percentage of data in Component 1:", round(component1_percentage, 2), "%\n")
    cat("Adjusted Percentage of data in Component 2:", round(component2_percentage, 2), "%\n")
}

# Usage:
nsg_safe_1week_df <- nsg_liver_1week_df %>% filter(gene == "Ctrl")
nsg_safe_3week_df <- nsg_liver_3week_df %>% filter(gene == "Ctrl")
nsg_nf1b_1week_df <- nsg_liver_1week_df %>% filter(grepl("Nfib", gene))
nsg_nf1b_3week_df <- nsg_liver_3week_df %>% filter(grepl("Nfib", gene))

res_1 <- check_bimodal(nsg_nf1b_1week_df$cell_num)
res_3 <- check_bimodal(nsg_nf1b_3week_df$cell_num)

plot_bimodal(nsg_nf1b_1week_df$cell_num, res_1)
plot_bimodal(nsg_nf1b_3week_df$cell_num, res_3)

# Fit a mixture model with constraints (fix means for both components)
week1_mix_model <- normalmixEM(nsg_liver_1week_df$log2_cell_num, k = 2, mean.constr = c(res_1$peaks[1], res_1$peaks[2]))
week3_mix_model <- normalmixEM(nsg_liver_3week_df$log2_cell_num, k = 2, mean.constr = c(res_3$peaks[1], res_3$peaks[2]))

# Get the posterior probabilities for each component
posterior_probs <- week1_mix_model$posterior
cell_num_data <- week1_mix_model$x

# Create dataframe with original data and posterior probabilities
classified_data <- data.frame(
    cell_num = cell_num_data,
    prob_peak1 = posterior_probs[, 1],
    prob_peak2 = posterior_probs[, 2]
)

# Split data by dominant component
week1_peak1_data <- classified_data %>%
    filter(prob_peak1 > prob_peak2)
week1_peak2_data <- classified_data %>%
    filter(prob_peak2 >= prob_peak1)


week1_peak1_data %>%
    ggplot(aes(x = cell_num)) +
    geom_density() +
    geom_vline(xintercept = res_1$peaks[1], color = "red", linetype = "dashed") +
    geom_vline(xintercept = res_1$valley[1], color = "blue", linetype = "dashed") +
    theme_minimal() +
    labs(
        x = "log2(cell_num)",
        y = "Density",
    )
week1_peak2_data %>%
    ggplot(aes(x = cell_num)) +
    geom_density() +
    geom_vline(xintercept = res_1$peaks[1], color = "red", linetype = "dashed") +
    geom_vline(xintercept = res_1$valley[1], color = "blue", linetype = "dashed") +
    theme_minimal() +
    labs(
        x = "log2(cell_num)",
        y = "Density",
    )


week3_peak1_data %>%
    ggplot(aes(x = cell_num)) +
    geom_density() +
    geom_vline(xintercept = res_3$peaks[1], color = "red", linetype = "dashed") +
    geom_vline(xintercept = res_3$valley[1], color = "blue", linetype = "dashed") +
    theme_minimal() +
    labs(
        x = "log2(cell_num)",
        y = "Density",
    )
week3_peak2_data %>% ggplot(aes(x = cell_num)) +
    geom_density()


# Calculate the percentage of data in each component
component1_percentage <- sum(posterior_probs[, 1]) / length(nsg_liver_1week_df$cell_num) * 100
component2_percentage <- sum(posterior_probs[, 2]) / length(nsg_liver_1week_df$cell_num) * 100

# Compare distributions
compare_distributions(week1_peak1_data$cell_num, week3_peak1_data$cell_num)
print(comparison)

rbind(
    peak2_data,
    peak1_data
) %>%
    ggplot(
        aes(
            x = log2(cell_num),
            fill = factor(Time_Point)
        )
    ) +
    geom_density(alpha = 0.3) +
    theme_minimal()
