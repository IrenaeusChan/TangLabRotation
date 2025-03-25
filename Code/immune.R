# This is the code used to generate the BL6 vs NSG vs Rag2 KO

# Author: Irenaeus Chan
# Date: 03/20/2025

# Load the necessary libraries
library(data.table)
library(tidyverse)
library(readxl)
library(ggplot2)

source("TangLabRotation/Code/help_functions.R")

# Load the Data
df <- fread("/Volumes/ruit/Active/AX_Working_Folder/2025_MobaV_Immune/outputs/FilteredSampleInfo.csv")
metadata <- read_excel("/Volumes/ruit/Active/AX_Working_Folder/2025_MobaV_Immune/CountingTable.xlsx")

# Set Gene
df <- df %>% filter(!is.na(distance))
df <- df %>% filter(distance == 0)
df <- df %>% filter(sgID != "sgDummy")

# Define sgSafe and other sgIDs
sgSafe <- c(
    "sgNT2-MMU", "sgNT3-MMU",
    "sgSafe19-MMU", "sgSafe29-MMU", "sgSafe36-MMU", "sgSafe4-MMU", "sgSafe28-MMU", "sgSafe30-MMU", "sgSafe33-MMU", "sgSafe5-MMU"
)

sgTSG <- c(
    "Chd7mmu-sg1", "Chd7mmu-sg2", "Chd7mmu-sg3",
    "Crebbpmmu-sg1", "Crebbpmmu-sg2", "Crebbpmmu-sg3",
    "Inpp5emmu-sg1", "Inpp5emmu-sg2", "Inpp5emmu-sg3",
    "Maddmmu-sg1", "Maddmmu-sg2", "Maddmmu-sg3",
    "Ptenmmu-sg1", "Ptenmmu-sg2", "Ptenmmu-sg3",
    "Setd2mmu-sg1", "Setd2mmu-sg2", "Setd2mmu-sg3",
    "Trip12mmu-sg1", "Trip12mmu-sg2", "Trip12mmu-sg3",
    "Tsc1mmu-sg1", "Tsc1mmu-sg2", "Tsc1mmu-sg3",
    "Tsc2mmu-sg1", "Tsc2mmu-sg2", "Tsc2mmu-sg3"
)

sgOnco <- c(
    "Actr3mmu-sg1", "Actr3mmu-sg2", "Actr3mmu-sg3",
    "Arpc4mmu-sg1", "Arpc4mmu-sg2", "Arpc4mmu-sg3",
    "Brk1mmu-sg1", "Brk1mmu-sg2", "Brk1mmu-sg3",
    "Cdh6mmu-sg1", "Cdh6mmu-sg2", "Cdh6mmu-sg3",
    "Gabrr1mmu-sg1", "Gabrr1mmu-sg2", "Gabrr1mmu-sg3",
    "Kdm5ammu-sg1", "Kdm5ammu-sg2", "Kdm5ammu-sg3",
    "KMT5Bmmu-sg1", "KMT5Bmmu-sg2", "KMT5Bmmu-sg3",
    "Lamc1mmu-sg1", "Lamc1mmu-sg2", "Lamc1mmu-sg3",
    "Rac1mmu-sg1", "Rac1mmu-sg2", "Rac1mmu-sg3",
    "Zeb1mmu-sg1", "Zeb1mmu-sg2", "Zeb1mmu-sg3"
)

sgBoth <- c(
    "Bag1mmu-sg1", "Bag1mmu-sg2", "Bag1mmu-sg3",
    "Ctnnd1mmu-sg1", "Ctnnd1mmu-sg2", "Ctnnd1mmu-sg3",
    "Ewsr1mmu-sg1", "Ewsr1mmu-sg2", "Ewsr1mmu-sg3",
    "Gata6mmu-sg1", "Gata6mmu-sg2", "Gata6mmu-sg3",
    "Runx1t1mmu-sg1", "Runx1t1mmu-sg2", "Runx1t1mmu-sg3"
)

sgOther <- c(
    "Gpatch8mmu-sg1", "Gpatch8mmu-sg2", "Gpatch8mmu-sg3",
    "Kctd10mmu-sg1", "Kctd10mmu-sg2", "Kctd10mmu-sg3",
    "Mgammu-sg1", "Mgammu-sg2", "Mgammu-sg3",
    "Zc3h7ammu-sg1", "Zc3h7ammu-sg2", "Zc3h7ammu-sg3",
    "Zfp574mmu-sg1", "Zfp574mmu-sg2", "Zfp574mmu-sg3"
)

sgAll <- c(unique(df$sgID))
sg_notSafe <- sgAll[!(sgAll %in% sgSafe)]
sgNfib <- c("Nfibmmu-sg1", "Nfibmmu-sg2", "Nfibmmu-sg3")

# add gene name, safe and nt as control
df$gene <- ifelse(df$sgID %in% sgSafe, "Ctrl", df$sgID)
df$gene <- ifelse(df$gene %in% sgNfib, "Nfib", df$gene)
df$gene <- ifelse(df$gene %in% sgOnco, "Oncogene", df$gene)
df$gene <- ifelse(df$gene %in% sgTSG, "TSG", df$gene)
df$gene <- ifelse(df$gene %in% sgBoth, "Both", df$gene)
df$gene <- ifelse(df$gene %in% sgOther, "Other", df$gene)

# Creating the constant Dataframes
dfpreTran <- df %>% filter(Sample_ID == "preTran")
preTran_cell_num <- aggregate(cell_num ~ sgID + Sample_ID + gene, data = dfpreTran, FUN = sum)
safe_df <- df %>% filter(gene == "Ctrl")

# Run the calc_MetBurden and calc_MetSeeding Function
all_results <- list()
# This is to generate per Mouse per Tissue MetBurden and MetSeeding
for (mouse in unique(na.omit(df$Mouse_ID))) {
    for (tissue in c("Liver", "Lung", "Brain")) {
        print(paste(mouse, tissue))
        # Filter data
        df_tissue_mouse <- df %>% filter(Tissue == tissue, Mouse_ID == mouse) %>% filter(cell_num > 50)
        safe_df_tissue_mouse <- safe_df %>% filter(Tissue == tissue, Mouse_ID == mouse) %>% filter(cell_num > 50)

        # Calculate statistics
        stats_df <- map_dfr(unique(df_tissue_mouse$sgID), function(sg) {
            print(sg)
            mb <- calc_MetBurden(df_tissue_mouse, safe_df_tissue_mouse, sg, preTran_cell_num)
            ms <- calc_MetSeeding(df_tissue_mouse, safe_df_tissue_mouse, sg, preTran_cell_num)

            data.frame(
                Mouse_ID = mouse,
                Tissue = tissue,
                sgID = sg,
                MetBurden = mb$MetBurden,
                MetSeeding = ms$MetSeeding
            )
        })
        all_results[[length(all_results) + 1]] <- stats_df
    }
}
bind_rows(all_results) %>% write_csv("MobaV_Immune_MetBurden_MetSeeding_MouseTissue.csv")
mobav_immune <- bind_rows(all_results)

#mobav_immune <- fread("TangLabRotation/ProcessedData/MobaV_Immune_MetBurden_MetSeeding_MouseTissue.csv")

mobav_immune <- mobav_immune %>%
    mutate(
        Log2_MetBurden = log2(MetBurden),
        Log2_MetSeeding = log2(MetSeeding)
    )

mobav_immune <- mobav_immune %>% left_join(metadata)

mouse_colors_Vimmune <- c(
    # Red Hues
    "73280" = "#FF4D4D",
    "73281" = "#FF1A1A",
    "73282" = "#E60000",
    "73283" = "#B30000",
    "73286" = "#800000",

    # Blue Hues
    "73390" = "#66B2FF",
    "73391" = "#3399FF",
    "73392" = "#0080FF",
    "73393" = "#0066CC",
    "73394" = "#004C99",

    # Purple Hues
    "73800" = "#F3E6FF",
    "73801" = "#CC99FF",
    "73802" = "#9933FF",
    "73803" = "#6600CC",
    "73804" = "#330066"
)

# Get the slope for each Mouse_ID linear model
summary_stats <- mobav_immune %>%
    filter(!is.na(MetBurden) & !is.na(MetSeeding)) %>%
    group_by(Mouse_ID, Tissue) %>%
    summarise(
        slope = coef(lm(Log2_MetBurden ~ Log2_MetSeeding - 1))[1],
        sd = sd(residuals(lm(Log2_MetBurden ~ Log2_MetSeeding - 1))),
        r.squared = summary(lm(Log2_MetBurden ~ Log2_MetSeeding - 1))$r.squared
    )

text_to_plot <- summary_stats %>%
    left_join(metadata, by = c("Mouse_ID", "Tissue")) %>%
    group_by(Mouse_Genotype, Tissue) %>%
    summarise(
        Mean = mean(slope),
        SD = sd(slope),
        Coeff_Variance = (SD / Mean) * 100,
        x = min(mobav_immune$Log2_MetSeeding), # position for label
        y = max(mobav_immune$Log2_MetBurden) # position for label
    )
p <- mobav_immune %>%
    mutate(
        Mouse_ID = as.factor(Mouse_ID)
    ) %>%
    #filter(Mouse_Genotype != "NSG") %>%
    ggplot(
        aes(
            y = Log2_MetBurden,
            x = Log2_MetSeeding,
            color = Mouse_ID
        )
    ) +
    geom_point(alpha = 0.2) +
    geom_smooth(
        method = "lm",
        formula = y ~ x -1,
        se = TRUE
    ) +
    geom_text(
        data = text_to_plot,
        aes(
            x = -6, y = 12,
            label = sprintf(
                "Mean Slope: %.2f\nSD: %.2f\nCV: %.1f%%",
                Mean, SD, Coeff_Variance
            )
        ),
        inherit.aes = FALSE,
        hjust = 0, vjust = 1
    ) +
    theme_bw() +
    labs(
        y = "Log2FC (Met Burden)",
        x = "Log2FC (Met Seeding)"
    ) +
    facet_wrap(~Mouse_Genotype + Tissue) +
    scale_color_manual(values = mouse_colors_Vimmune)

pdf("MobaV_Immune_MouseGenotype_MetBurden_vs_MetSeeding_byTissue.pdf", width = 10, height = 10)
print(p)
dev.off()

# This is to generate per Mouse MetBurden and MetSeeding
all_results <- list()

for (mouse in unique(na.omit(df$Mouse_ID))) {
    print(paste(mouse))
    # Filter data
    df_mouse <- df %>% filter(Tissue != "Blood", Mouse_ID == mouse) %>% filter(cell_num > 50)
    safe_df_mouse <- safe_df %>% filter(Tissue != "Blood", Mouse_ID == mouse) %>% filter(cell_num > 50)

    # Calculate statistics
    stats_df <- map_dfr(unique(df_tissue_mouse$sgID), function(sg) {
        print(sg)
        mb <- calc_MetBurden(df_mouse, safe_df_mouse, sg, preTran_cell_num)
        ms <- calc_MetSeeding(df_mouse, safe_df_mouse, sg, preTran_cell_num)

        data.frame(
            Mouse_ID = mouse,
            sgID = sg,
            MetBurden = mb$MetBurden,
            MetSeeding = ms$MetSeeding
        )
    })
    all_results[[length(all_results) + 1]] <- stats_df
}
bind_rows(all_results) %>% write_csv("MobaV_Immune_MetBurden_MetSeeding.csv")
mobav_immune <- bind_rows(all_results)

# mobav_immune <- fread("TangLabRotation/ProcessedData/MobaV_Immune_MetBurden_MetSeeding.csv")

mobav_immune <- mobav_immune %>%
    mutate(
        Log2_MetBurden = log2(MetBurden),
        Log2_MetSeeding = log2(MetSeeding)
    )

mobav_immune <- mobav_immune %>% left_join(metadata)

# Get the slope for each Mouse_ID linear model
summary_stats <- mobav_immune %>%
    filter(!is.na(MetBurden) & !is.na(MetSeeding)) %>%
    group_by(Mouse_ID) %>%
    summarise(
        slope = coef(lm(Log2_MetBurden ~ Log2_MetSeeding - 1))[1],
        sd = sd(residuals(lm(Log2_MetBurden ~ Log2_MetSeeding - 1))),
        r.squared = summary(lm(Log2_MetBurden ~ Log2_MetSeeding - 1))$r.squared
    )

text_to_plot <- summary_stats %>%
    left_join(metadata, by = c("Mouse_ID")) %>%
    group_by(Mouse_Genotype) %>%
    summarise(
        Mean = mean(slope),
        SD = sd(slope),
        Coeff_Variance = (SD / Mean) * 100,
        x = min(mobav_immune$Log2_MetSeeding), # position for label
        y = max(mobav_immune$Log2_MetBurden) # position for label
    )

p <- mobav_immune %>%
    mutate(
        Mouse_ID = as.factor(Mouse_ID)
    ) %>%
    # filter(Mouse_Genotype != "NSG") %>%
    ggplot(
        aes(
            y = Log2_MetBurden,
            x = Log2_MetSeeding,
            color = Mouse_ID
        )
    ) +
    geom_point(alpha = 0.2) +
    geom_smooth(
        method = "lm",
        formula = y ~ x - 1,
        se = TRUE
    ) +
    geom_text(
        data = text_to_plot,
        aes(
            x = -4, y = 8,
            label = sprintf(
                "Mean Slope: %.2f\nSD: %.2f\nCV: %.1f%%",
                Mean, SD, Coeff_Variance
            )
        ),
        inherit.aes = FALSE,
        hjust = 0, vjust = 1
    ) +
    theme_bw() +
    labs(
        y = "Log2FC (Met Burden)",
        x = "Log2FC (Met Seeding)"
    ) +
    facet_wrap(~ Mouse_Genotype) +
    scale_color_manual(values = mouse_colors_Vimmune)

pdf("MobaV_Immune_MouseGenotype_MetBurden_vs_MetSeeding.pdf", width = 6, height = 10)
print(p)
dev.off()
