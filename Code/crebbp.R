# This is the code used to generate the Crebbp specific analyses

# Author: Irenaeus Chan
# Date: 03/20/2025

# Load the necessary libraries
library(data.table)
library(tidyverse)
library(readxl)
library(pheatmap)
library(ggrepel)
library(patchwork)
library(ggpubr)
library(sjPlot)

source("TangLabRotation/Code/help_functions.R")

metadata <- read_excel("/Volumes/ruit/Active/IrenaeusChan/MobaV_FullSet/input_files/CountingTable.xlsx")
df <- fread("/Volumes/ruit/Active/IrenaeusChan/MobaV_FullSet/output_new_outlier/MobaSeq_FilteredSampleInfo.csv")
df <- df %>% filter(!is.na(distance))
df <- df %>% filter(distance == 0)
df <- df %>% filter(sgID != "sgDummy")

# There were 2 Cell Lines used in the experiment
df <- df %>% filter(
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
safe_df <- df %>% filter(gene == "Ctrl")
preTran_df <- df %>% filter(Sample_ID == "preTran")
preTran_cell_num <- aggregate(cell_num ~ sgID + Sample_ID + gene + Cell_Line, data = preTran_df, FUN = sum)


# These are all Mouse specific. So let's calculate PER Mouse, all of the below statistics. Then we can compare if any of them are significantly different
all_results <- list()
# This will take a while to run...
for (mouse in unique(na.omit(df$Mouse_ID))) {
    for (tissue in c("Liver", "Lung", "Brain", "BoneMarrow", "WholeBlood", "Cells")) {
        for (cell_line in c("rp48", "rp116")) {
            print(paste(mouse, tissue, cell_line))
            # Filter data
            df_tissue_mouse <- df %>% filter(Tissue == tissue, Mouse_ID == mouse, Cell_Line == cell_line)
            safe_df_tissue_mouse <- safe_df %>% filter(Tissue == tissue, Mouse_ID == mouse, Cell_Line == cell_line)
            preTran_cell_num_mouse <- preTran_cell_num %>% filter(Sample_ID == "preTran", Cell_Line == cell_line)

            # Calculate statistics
            stats_df <- map_dfr(unique(df_tissue_mouse$sgID), function(sg) {
                print(sg)
                mb <- calc_MetBurden(df_tissue_mouse, safe_df_tissue_mouse, sg, preTran_cell_num_mouse)
                ms <- calc_MetSeeding(df_tissue_mouse, safe_df_tissue_mouse, sg, preTran_cell_num_mouse)
                sp_90 <- calc_SizePercentile(df_tissue_mouse, safe_df_tissue_mouse, sg, percentile = 0.9)
                sp_50 <- calc_SizePercentile(df_tissue_mouse, safe_df_tissue_mouse, sg, percentile = 0.5)
                sp_20 <- calc_SizePercentile(df_tissue_mouse, safe_df_tissue_mouse, sg, percentile = 0.2)
                pm <- calc_PeakMode(df_tissue_mouse, safe_df_tissue_mouse, sg)
                md <- calc_MetDormancy(df_tissue_mouse, safe_df_tissue_mouse, sg)

                data.frame(
                    Mouse_ID = mouse,
                    Tissue = tissue,
                    sgID = sg,
                    Cell_Line = cell_line,
                    MetBurden = mb$MetBurden,
                    MetSeeding = ms$MetSeeding,
                    SizePercentile = sp_90$SizePercentile,
                    SizePercentile_50 = sp_50$SizePercentile,
                    SizePercentile_20 = sp_20$SizePercentile,
                    PeakMode = pm$PeakMode,
                    MetDormancy = md$MetDormancy
                )
            })

            all_results[[length(all_results) + 1]] <- stats_df
        }
        all_results[[length(all_results) + 1]] <- stats_df
    }
}
#fwrite(rbindlist(all_results), "all_results_output.csv")
final_df <- fread("TangLabRotation/ProcessedData/all_results_output.csv")
final_df <- final_df %>%
    left_join(
        metadata %>% select(Mouse_ID, Sex, Mouse_Genotype, Tissue, Time_Point),
        by = c("Mouse_ID" = "Mouse_ID", "Tissue" = "Tissue")
    )
final_df <- final_df %>% distinct()
final_df <- final_df %>%
    mutate(
        gene = case_when(
            sgID %in% sgSafe ~ "Ctrl",
            sgID %in% sgNfib ~ "Nfib",
            sgID %in% sgOnco ~ "Oncogene",
            sgID %in% sgTSG ~ "TSG",
            sgID %in% sgBoth ~ "Both",
            sgID %in% sgOther ~ "Other",
            TRUE ~ "Other"
        )
    )
final_df <- final_df %>%
    mutate(
        Log2FC_MetBurden = log2(MetBurden),
        Log2FC_MetSeeding = log2(MetSeeding),
        Log2FC_SizePercentile = log2(SizePercentile),
        Log2FC_SizePercentile_50 = log2(SizePercentile_50),
        Log2FC_SizePercentile_20 = log2(SizePercentile_20),
        Log2FC_PeakMode = log2(PeakMode),
        Log2FC_MetDormancy = log2(MetDormancy)
    )
final_df <- final_df %>%
    mutate(
        target = gsub("mmu-sg[0-9]", "", sgID),
        target = ifelse(grepl("sgSafe", sgID) | grepl("sgNT", sgID), "sgSafe", target)
    )

pdf("Heatmap_byTissue_byCellLine.pdf")
for (statistic in c("Log2FC_MetBurden", "Log2FC_MetSeeding", "Log2FC_SizePercentile")) {
    for (tissue in unique(final_df$Tissue)) {
        for (cellline in c("rp48", "rp116")) {
            print(paste(statistic, tissue, cellline))
            remove_sg <- final_df %>%
                filter(Time_Point == "3weeks", Tissue == tissue, Cell_Line == cellline) %>%
                select(Mouse_ID, target, !!sym(statistic)) %>%
                group_by(Mouse_ID, target) %>%
                summarise(n = sum(!!sym(statistic))) %>%
                filter(is.na(n)) %>%
                ungroup() %>%
                select(target) %>%
                unique()
            heatmap_data <- final_df %>%
                filter(Time_Point == "3weeks" & Tissue == tissue & Cell_Line == cellline) %>%
                filter(!(target %in% remove_sg$target)) %>%
                select(Mouse_ID, target, !!sym(statistic)) %>%
                mutate(
                    across(
                        all_of(statistic),
                        ~if_else(
                            is.na(.x) | is.infinite(.x),
                            0,
                            as.numeric(.x)
                        )
                    )
                ) %>%
                group_by(target, Mouse_ID) %>%
                summarise(
                    average = mean(!!sym(statistic)),
                    .groups = "drop"
                ) %>%
                pivot_wider(
                    names_from = Mouse_ID,
                    values_from = average,
                    values_fill = 0
                ) %>%
                column_to_rownames("target")
            # Remove Gene if all values are 0
            heatmap_data <- heatmap_data %>% filter(rowSums(heatmap_data) != 0)
            if (nrow(heatmap_data) < 2) {
                next
            }
            # Add small noise to break ties
            heatmap_data <- heatmap_data + matrix(
                rnorm(nrow(heatmap_data) * ncol(heatmap_data), mean = 0, sd = 1e-10),
                nrow = nrow(heatmap_data)
            )
            annotations <- final_df %>%
                filter(Time_Point == "3weeks", Tissue == tissue, Cell_Line == cellline) %>%
                select(Mouse_ID, Sex, Mouse_Genotype) %>%
                distinct() %>%
                column_to_rownames("Mouse_ID")
            annotations_row <- final_df %>%
                filter(Time_Point == "3weeks", Tissue == tissue, Cell_Line == cellline) %>%
                select(target, gene) %>%
                distinct() %>%
                column_to_rownames("target")
            # Plot heatmap with additional checks
            print(pheatmap(
                mat = as.matrix(heatmap_data),
                annotation_col = annotations,
                annotation_row = annotations_row,
                scale = "row",
                cluster_rows = TRUE,
                cluster_cols = TRUE,
                show_rownames = TRUE,
                show_colnames = TRUE,
                main = paste0(statistic, " by Mouse and sgID in ", tissue, " - Cell Line: ", cellline),
                color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
                na_col = "grey", # Handle any remaining NAs,
                #fontsize_row = 6
            ))
        }
    }
}
dev.off()

pdf("Heatmap_byCellLine.pdf", width = 14, height = 10)
for (statistic in c("Log2FC_MetBurden", "Log2FC_MetSeeding", "Log2FC_SizePercentile")) {
    for (cellline in c("rp48", "rp116")) {
        print(paste(statistic, cellline))
        remove_sg <- final_df %>%
            filter(Time_Point == "3weeks", Cell_Line == cellline, (Tissue == "Liver" | Tissue == "Lung")) %>%
            mutate(Mouse_ID_Tissue = paste(Mouse_ID, Tissue)) %>%
            select(Mouse_ID_Tissue, target, !!sym(statistic)) %>%
            group_by(Mouse_ID_Tissue, target) %>%
            summarise(n = sum(!!sym(statistic))) %>%
            filter(is.na(n)) %>%
            ungroup() %>%
            select(target) %>%
            unique()
        heatmap_data <- final_df %>%
            filter(Time_Point == "3weeks" & Cell_Line == cellline, (Tissue == "Liver" | Tissue == "Lung")) %>%
            #filter(!(target %in% remove_sg$target)) %>%
            mutate(Mouse_ID_Tissue = paste(Mouse_ID, Tissue)) %>%
            select(Mouse_ID_Tissue, target, !!sym(statistic)) %>%
            mutate(
                across(
                    all_of(statistic),
                    ~ if_else(
                        is.na(.x) | is.infinite(.x),
                        0,
                        as.numeric(.x)
                    )
                )
            ) %>%
            group_by(target, Mouse_ID_Tissue) %>%
            summarise(
                average = mean(!!sym(statistic)),
                .groups = "drop"
            ) %>%
            pivot_wider(
                names_from = Mouse_ID_Tissue,
                values_from = average,
                values_fill = 0
            ) %>%
            column_to_rownames("target")
        # Remove Gene if all values are 0
        heatmap_data <- heatmap_data %>% filter(rowSums(heatmap_data) != 0)
        if (nrow(heatmap_data) < 2) {
            next
        }
        # Add small noise to break ties
        heatmap_data <- heatmap_data + matrix(
            rnorm(nrow(heatmap_data) * ncol(heatmap_data), mean = 0, sd = 1e-10),
            nrow = nrow(heatmap_data)
        )
        annotations <- final_df %>%
            filter(Time_Point == "3weeks", Cell_Line == cellline) %>%
            mutate(Mouse_ID_Tissue = paste(Mouse_ID, Tissue)) %>%
            select(Mouse_ID_Tissue, Sex, Mouse_Genotype) %>%
            distinct() %>%
            column_to_rownames("Mouse_ID_Tissue")
        annotations_row <- final_df %>%
            filter(Time_Point == "3weeks", Cell_Line == cellline) %>%
            select(target, gene) %>%
            distinct() %>%
            column_to_rownames("target")
        # Plot heatmap with additional checks
        print(pheatmap(
            mat = as.matrix(heatmap_data),
            annotation_col = annotations,
            annotation_row = annotations_row,
            scale = "row",
            cluster_rows = TRUE,
            cluster_cols = TRUE,
            show_rownames = TRUE,
            show_colnames = TRUE,
            main = paste0(statistic, " by Mouse and sgID - Cell Line: ", cellline),
            color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
            na_col = "grey", # Handle any remaining NAs,
            # fontsize_row = 6
        ))
    }
}
dev.off()

# Crebbp Specific Analysis - No Mouse Specific
all_results <- list()
for (tissue in c("Liver", "Lung", "Brain", "BoneMarrow", "WholeBlood")) {
    for (cell_line in c("rp48", "rp116")) {
        for (genotype in c("NSG", "BL6")) {
            for (sex in c("F", "M")) {
                print(paste(tissue, cell_line, genotype, sex))
                if (tissue == "Liver") {
                    for (week in c("1week", "3weeks")) {
                        print(paste(tissue, cell_line, genotype, sex, week))
                        df_to_boot <- df %>% filter(Tissue == tissue, Cell_Line == cell_line, Mouse_Genotype == genotype, Sex == sex, Time_Point == week)
                        safe_df_to_boot <- safe_df %>% filter(Tissue == tissue, Cell_Line == cell_line, Mouse_Genotype == genotype, Sex == sex, Time_Point == week)
                        preTran_cell_num_tissue <- preTran_cell_num %>% filter(Sample_ID == "preTran", Cell_Line == cell_line)

                        stats_df <- map_dfr(unique(df_to_boot$sgID), function(sg) {
                            print(sg)
                            print("MetBurden")
                            mb <- calc_MetBurden(df_to_boot, safe_df_to_boot, sg, preTran_cell_num_tissue)
                            print("MetSeeding")
                            ms <- calc_MetSeeding(df_to_boot, safe_df_to_boot, sg, preTran_cell_num_tissue)
                            print("SizePercentile")
                            sp <- calc_SizePercentile(df_to_boot, safe_df_to_boot, sg)
                            print("PeakMode")
                            pm <- calc_PeakMode(df_to_boot, safe_df_to_boot, sg)
                            print("MetDormancy")
                            md <- calc_MetDormancy(df_to_boot, safe_df_to_boot, sg)
                            print("Done")

                            data.frame(
                                Tissue = tissue,
                                sgID = sg,
                                Cell_Line = cell_line,
                                Mouse_Genotype = genotype,
                                Sex = sex,
                                Time_Point = week,
                                MetBurden = mb$MetBurden,
                                MetSeeding = ms$MetSeeding,
                                SizePercentile = sp$SizePercentile,
                                PeakMode = pm$PeakMode,
                                MetDormancy = md$MetDormancy
                            )
                        })
                        all_results[[length(all_results) + 1]] <- stats_df
                    }
                } else {
                    df_to_boot <- df %>% filter(Tissue == tissue, Cell_Line == cell_line, Mouse_Genotype == genotype, Sex == sex, Time_Point == week)
                    safe_df_to_boot <- safe_df %>% filter(Tissue == tissue, Cell_Line == cell_line, Mouse_Genotype == genotype, Sex == sex, Time_Point == week)
                    preTran_cell_num_tissue <- preTran_cell_num %>% filter(Sample_ID == "preTran", Cell_Line == cell_line)

                    stats_df <- map_dfr(unique(df_to_boot$sgID), function(sg) {
                        print(sg)
                        print("MetBurden")
                        mb <- calc_MetBurden(df_to_boot, safe_df_to_boot, sg, preTran_cell_num_tissue)
                        print("MetSeeding")
                        ms <- calc_MetSeeding(df_to_boot, safe_df_to_boot, sg, preTran_cell_num_tissue)
                        print("SizePercentile")
                        sp <- calc_SizePercentile(df_to_boot, safe_df_to_boot, sg)
                        print("PeakMode")
                        pm <- calc_PeakMode(df_to_boot, safe_df_to_boot, sg)
                        print("MetDormancy")
                        md <- calc_MetDormancy(df_to_boot, safe_df_to_boot, sg)
                        print("Done")

                        data.frame(
                            Tissue = tissue,
                            sgID = sg,
                            Cell_Line = cell_line,
                            Mouse_Genotype = genotype,
                            Sex = sex,
                            Time_Point = week,
                            MetBurden = mb$MetBurden,
                            MetSeeding = ms$MetSeeding,
                            SizePercentile = sp$SizePercentile,
                            PeakMode = pm$PeakMode,
                            MetDormancy = md$MetDormancy
                        )
                    })
                    all_results[[length(all_results) + 1]] <- stats_df
                }
            }
            all_results[[length(all_results) + 1]] <- stats_df
        }
        all_results[[length(all_results) + 1]] <- stats_df
    }
    all_results[[length(all_results) + 1]] <- stats_df
}
#fwrite(rbindlist(all_results), "all_results_grouped_output.csv")
grouped_df <- fread("all_results_grouped_output.csv")
grouped_df <- grouped_df %>%
    mutate(
        gene = case_when(
            sgID %in% sgSafe ~ "Ctrl",
            sgID %in% sgNfib ~ "Nfib",
            sgID %in% sgOnco ~ "Oncogene",
            sgID %in% sgTSG ~ "TSG",
            sgID %in% sgBoth ~ "Both",
            sgID %in% sgOther ~ "Other",
            TRUE ~ "Other"
        )
    )
grouped_df <- grouped_df %>%
    mutate(
        Log2FC_MetBurden = log2(MetBurden),
        Log2FC_MetSeeding = log2(MetSeeding),
        Log2FC_SizePercentile = log2(SizePercentile),
        Log2FC_PeakMode = log2(PeakMode),
        Log2FC_MetDormancy = log2(MetDormancy)
    )
grouped_df <- grouped_df %>%
    mutate(
        target = gsub("mmu-sg[0-9]", "", sgID),
        target = ifelse(grepl("sgSafe", sgID) | grepl("sgNT", sgID), "sgSafe", target)
    )

crebbp_grouped_safe_df <- rbind(
        grouped_df %>% filter(target == "Crebbp"),
        grouped_df %>% filter(target == "sgSafe")
    ) %>%
    mutate(
        Mouse_Info = paste(Mouse_Genotype, Sex)
    )
crebbp_safe_df <- rbind(
    final_df %>% filter(target == "Crebbp"),
    final_df %>% filter(target == "sgSafe")
) %>%
    mutate(
        Mouse_Info = paste(Mouse_Genotype, Sex)
    )

# Calculate the statistics using t.test
stat_data_rp116 <- crebbp_safe_df %>%
    filter(Tissue == "Liver", Cell_Line == "rp116") %>%
    group_by(Mouse_Info, Time_Point) %>%
    summarize(
        p_value = t.test(
            Log2FC_MetSeeding[target == "Crebbp"],
            Log2FC_MetSeeding[target == "sgSafe"]
        )$p.value,
        sig = case_when(
            p_value < 0.001 ~ "***",
            p_value < 0.01 ~ "**",
            p_value < 0.05 ~ "*",
            TRUE ~ "ns"
        )
    )
stat_data_rp48 <- crebbp_safe_df %>%
    filter(Tissue == "Liver", Cell_Line == "rp48") %>%
    group_by(Mouse_Info, Time_Point) %>%
    summarize(
        p_value = t.test(
            Log2FC_MetSeeding[target == "Crebbp"],
            Log2FC_MetSeeding[target == "sgSafe"]
        )$p.value,
        sig = case_when(
            p_value < 0.001 ~ "***",
            p_value < 0.01 ~ "**",
            p_value < 0.05 ~ "*",
            TRUE ~ "ns"
        )
    )

# Create function for common plot elements
add_plot_elements <- function(p) {
    p +
    geom_boxplot(
        aes(
            fill = group
        ),
        outlier.shape = NA, alpha = 0.8
    ) +
        geom_point(
            aes(
                fill = group,
            ),
            position = position_dodge(width = 0.8),
            alpha = 0.3, color = "black"
        ) +
        geom_hline(yintercept = 0, linetype = "dashed") +
        facet_wrap(~Cell_Line, ncol = 1) +
        labs(title = "", x = "", y = "MetSeeding (Log2 Fold Change)") +
        scale_fill_manual(
            name = "Time Point",
            values = c(
                "BL6 F.Crebbp.1week" = "#66B2FF",
                "BL6 F.Crebbp.3weeks" = "#0066CC",
                "BL6 M.Crebbp.1week" = "#66B2FF",
                "BL6 M.Crebbp.3weeks" = "#0066CC",
                "NSG F.Crebbp.1week" = "#FF4D4D",
                "NSG F.Crebbp.3weeks" = "#E60000",
                "NSG M.Crebbp.1week" = "#FF4D4D",
                "NSG M.Crebbp.3weeks" = "#E60000",
                "BL6 F.sgSafe.1week" = "#CCCCCC",
                "BL6 F.sgSafe.3weeks" = "#666666",
                "BL6 M.sgSafe.1week" = "#CCCCCC",
                "BL6 M.sgSafe.3weeks" = "#666666",
                "NSG F.sgSafe.1week" = "#CCCCCC",
                "NSG F.sgSafe.3weeks" = "#666666",
                "NSG M.sgSafe.1week" = "#CCCCCC",
                "NSG M.sgSafe.3weeks" = "#666666"
            ),
            labels = c(
                "BL6 F.Crebbp.1week" = "1 Week - Crebbp",
                "BL6 F.Crebbp.3weeks" = "3 Weeks - Crebbp",
                "BL6 M.Crebbp.1week" = "1 Week - Crebbp",
                "BL6 M.Crebbp.3weeks" = "3 Weeks - Crebbp",
                "NSG F.Crebbp.1week" = "1 Week - Crebbp",
                "NSG F.Crebbp.3weeks" = "3 Weeks - Crebbp",
                "NSG M.Crebbp.1week" = "1 Week - Crebbp",
                "NSG M.Crebbp.3weeks" = "3 Weeks - Crebbp",
                "BL6 F.sgSafe.1week" = "1 Week - sgSafe",
                "BL6 F.sgSafe.3weeks" = "3 Weeks - sgSafe",
                "BL6 M.sgSafe.1week" = "1 Week - sgSafe",
                "BL6 M.sgSafe.3weeks" = "3 Weeks - sgSafe",
                "NSG F.sgSafe.1week" = "1 Week - sgSafe",
                "NSG F.sgSafe.3weeks" = "3 Weeks - sgSafe",
                "NSG M.sgSafe.1week" = "1 Week - sgSafe",
                "NSG M.sgSafe.3weeks" = "3 Weeks - sgSafe"
            )
        ) +
        theme_minimal() +
        theme(
            legend.position = "none",
            legend.box = "vertical"
        )
}

p1 <- crebbp_safe_df %>%
    filter(Tissue == "Liver", Cell_Line == "rp116") %>%
    mutate(
        group = interaction(Mouse_Info, target, Time_Point),
        sg_shape = ifelse(target == "sgSafe", target, sgID)
    ) %>%
    ggplot(
        aes(
            x = Mouse_Info,
            y = Log2FC_MetSeeding
        )
    ) +
    geom_signif(y_position = 1.9, xmin = 0.72, xmax = 0.9, annotation = paste0(stat_data_rp116$sig[1]), tip_length = 0.03, size = 0.3) +
    geom_signif(y_position = 1.9, xmin = 1.1, xmax = 1.28, annotation = paste0(stat_data_rp116$sig[2]), tip_length = 0.03, size = 0.3) +
    geom_signif(y_position = 1.9, xmin = 1.72, xmax = 1.9, annotation = paste0(stat_data_rp116$sig[3]), tip_length = 0.03, size = 0.3) +
    geom_signif(y_position = 1.9, xmin = 2.1, xmax = 2.28, annotation = paste0(stat_data_rp116$sig[4]), tip_length = 0.03, size = 0.3) +
    geom_signif(y_position = 1.8, xmin = 2.72, xmax = 2.9, annotation = paste0(stat_data_rp116$sig[5]), tip_length = 0.03, size = 0.3) +
    geom_signif(y_position = 1.8, xmin = 3.1, xmax = 3.28, annotation = paste0(stat_data_rp116$sig[6]), tip_length = 0.03, size = 0.3) +
    geom_signif(y_position = 2.4, xmin = 3.72, xmax = 3.9, annotation = paste0(stat_data_rp116$sig[7]), tip_length = 0.03, size = 0.3) +
    geom_signif(y_position = 2.4, xmin = 4.1, xmax = 4.28, annotation = paste0(stat_data_rp116$sig[8]), tip_length = 0.03, size = 0.3) +
    scale_y_continuous(trans = scales::pseudo_log_trans(base = 2), labels = scales::label_number(), breaks = c(-10, -5, -2, 0, 2, 5, 10), limits = c(NA, 6))
p2 <- crebbp_safe_df %>%
    filter(Tissue == "Liver", Cell_Line == "rp48") %>%
    mutate(
        group = interaction(Mouse_Info, target, Time_Point),
        sg_shape = ifelse(target == "sgSafe", target, sgID)
    ) %>%
    ggplot(
        aes(
            x = Mouse_Info,
            y = Log2FC_MetSeeding
        )
    ) +
    geom_signif(y_position = 1.2, xmin = 0.72, xmax = 0.9, annotation = paste0(stat_data_rp48$sig[1]), tip_length = 0.03, size = 0.3) +
    geom_signif(y_position = 1.2, xmin = 1.1, xmax = 1.28, annotation = paste0(stat_data_rp48$sig[2]), tip_length = 0.03, size = 0.3) +
    geom_signif(y_position = 1.1, xmin = 1.72, xmax = 1.9, annotation = paste0(stat_data_rp48$sig[3]), tip_length = 0.03, size = 0.3) +
    geom_signif(y_position = 1.1, xmin = 2.1, xmax = 2.28, annotation = paste0(stat_data_rp48$sig[4]), tip_length = 0.03, size = 0.3) +
    geom_signif(y_position = 1.1, xmin = 2.72, xmax = 2.9, annotation = paste0(stat_data_rp48$sig[5]), tip_length = 0.03, size = 0.3) +
    geom_signif(y_position = 1.1, xmin = 3.1, xmax = 3.28, annotation = paste0(stat_data_rp48$sig[6]), tip_length = 0.03, size = 0.3) +
    geom_signif(y_position = 1.2, xmin = 3.72, xmax = 3.9, annotation = paste0(stat_data_rp48$sig[7]), tip_length = 0.03, size = 0.3) +
    geom_signif(y_position = 1.2, xmin = 4.1, xmax = 4.28, annotation = paste0(stat_data_rp48$sig[8]), tip_length = 0.03, size = 0.3) +
    scale_y_continuous(trans = scales::pseudo_log_trans(base = 2), labels = scales::label_number(), breaks = c(-10, -5, -2, 0, 2, 5, 10), limits = c(NA, 2.4))
p1 <- add_plot_elements(p1)
p2 <- add_plot_elements(p2)
    
pdf("Crebbp_vs_sgSafe_MetSeeding_Liver_rp116vsrp48_1weekVS3week.pdf", width = 8, height = 12)
print(p1 / p2)
dev.off()

# Add Target to the Dataframe
df <- df %>%
    mutate(
        target = gsub("mmu-sg[0-9]", "", sgID),
        target = ifelse(grepl("sgSafe", sgID) | grepl("sgNT", sgID), "sgSafe", target)
    )

# Look at Dormancy
df_1week_rp48_bl6 <- df %>% filter(target == "Crebbp", Tissue == "Liver", Cell_Line == "rp48", Time_Point == "1week", Mouse_Genotype == "BL6")
df_1week_rp48_nsg <- df %>% filter(target == "Crebbp", Tissue == "Liver", Cell_Line == "rp48", Time_Point == "1week", Mouse_Genotype == "NSG")
df_1week_rp116_bl6 <- df %>% filter(target == "Crebbp", Tissue == "Liver", Cell_Line == "rp116", Time_Point == "1week", Mouse_Genotype == "BL6")
df_1week_rp116_nsg <- df %>% filter(target == "Crebbp", Tissue == "Liver", Cell_Line == "rp116", Time_Point == "1week", Mouse_Genotype == "NSG")
df_3weeks_rp48_bl6 <- df %>% filter(target == "Crebbp", Tissue == "Liver", Cell_Line == "rp48", Time_Point == "3weeks", Mouse_Genotype == "BL6")
df_3weeks_rp48_nsg <- df %>% filter(target == "Crebbp", Tissue == "Liver", Cell_Line == "rp48", Time_Point == "3weeks", Mouse_Genotype == "NSG")
df_3weeks_rp116_bl6 <- df %>% filter(target == "Crebbp", Tissue == "Liver", Cell_Line == "rp116", Time_Point == "3weeks", Mouse_Genotype == "BL6")
df_3weeks_rp116_nsg <- df %>% filter(target == "Crebbp", Tissue == "Liver", Cell_Line == "rp116", Time_Point == "3weeks", Mouse_Genotype == "NSG")

# Check for bimodality
res_1week_rp48_bl6 <- check_bimodal(df_1week_rp48_bl6$cell_num)
res_1week_rp48_nsg <- check_bimodal(df_1week_rp48_nsg$cell_num)
res_1week_rp116_bl6 <- check_bimodal(df_1week_rp116_bl6$cell_num)
res_1week_rp116_nsg <- check_bimodal(df_1week_rp116_nsg$cell_num)
res_3weeks_rp48_bl6 <- check_bimodal(df_3weeks_rp48_bl6$cell_num)
res_3weeks_rp48_nsg <- check_bimodal(df_3weeks_rp48_nsg$cell_num)
res_3weeks_rp116_bl6 <- check_bimodal(df_3weeks_rp116_bl6$cell_num)
res_3weeks_rp116_nsg <- check_bimodal(df_3weeks_rp116_nsg$cell_num)

# Make the plots
den_p1 <- plot_bimodal(df_3weeks_rp116_bl6$cell_num, res_3weeks_rp116_bl6, show_fit = FALSE, show_valleys = FALSE, return_data = TRUE)
den_p2 <- plot_bimodal(df_3weeks_rp48_bl6$cell_num, res_3weeks_rp48_bl6, show_fit = FALSE, show_valleys = FALSE, return_data = TRUE)
den_p3 <- plot_bimodal(df_3weeks_rp116_nsg$cell_num, res_3weeks_rp116_nsg, show_fit = FALSE, show_valleys = FALSE, return_data = TRUE)
den_p4 <- plot_bimodal(df_3weeks_rp48_nsg$cell_num, res_3weeks_rp48_nsg, show_fit = FALSE, show_valleys = FALSE, return_data = TRUE)

to_plot <- rbind(
    den_p1 %>% mutate(Cell_Line = "rp116", Mouse_Genotype = "BL6"),
    den_p2 %>% mutate(Cell_Line = "rp48", Mouse_Genotype = "BL6"),
    den_p3 %>% mutate(Cell_Line = "rp116", Mouse_Genotype = "NSG"),
    den_p4 %>% mutate(Cell_Line = "rp48", Mouse_Genotype = "NSG")
)

to_plot %>%
    mutate(
        group = paste(Mouse_Genotype, Cell_Line)
    ) %>%
    ggplot(
        aes(
            x = x,
            fill = group
        )
    ) +
    geom_density(alpha = 0.3) +
    geom_vline(xintercept = res_3weeks_rp116_bl6$peaks[1], linetype = "dashed") +
    geom_vline(xintercept = res_3weeks_rp116_bl6$peaks[2], linetype = "dashed") +
    geom_vline(xintercept = res_3weeks_rp48_bl6$peaks[1], linetype = "dashed") +
    geom_vline(xintercept = res_3weeks_rp48_bl6$peaks[2], linetype = "dashed") +
    geom_vline(xintercept = res_3weeks_rp116_nsg$peaks[1], linetype = "dashed") +
    geom_vline(xintercept = res_3weeks_rp116_nsg$peaks[2], linetype = "dashed") +
    geom_vline(xintercept = res_3weeks_rp48_nsg$peaks[1], linetype = "dashed") +
    geom_vline(xintercept = res_3weeks_rp48_nsg$peaks[2], linetype = "dashed") +
    # BL6 RP116
    geom_text(x = res_3weeks_rp116_bl6$peaks[1], y = res_3weeks_rp116_bl6$peak_heights[1], label = paste0(round(den_p1$comp1_perc[1], 2), "%"), vjust = -1) +
    geom_text(x = res_3weeks_rp116_bl6$peaks[2], y = res_3weeks_rp116_bl6$peak_heights[2], label = paste0(round(den_p1$comp2_perc[1], 2), "%"), vjust = -1) +
    # BL6 RP48
    geom_text(x = res_3weeks_rp48_bl6$peaks[1], y = res_3weeks_rp48_bl6$peak_heights[1], label = paste0(round(den_p2$comp1_perc[1], 2), "%"), vjust = -1) +
    geom_text(x = res_3weeks_rp48_bl6$peaks[2], y = res_3weeks_rp48_bl6$peak_heights[2], label = paste0(round(den_p2$comp2_perc[1], 2), "%"), vjust = -1) +
    # NSG RP116
    geom_text(x = res_3weeks_rp116_nsg$peaks[1], y = res_3weeks_rp116_nsg$peak_heights[1], label = paste0(round(den_p3$comp1_perc[1], 2), "%"), vjust = -1) +
    geom_text(x = res_3weeks_rp116_nsg$peaks[2], y = res_3weeks_rp116_nsg$peak_heights[2], label = paste0(round(den_p3$comp2_perc[1], 2), "%"), vjust = -1) +
    # NSG RP48
    geom_text(x = res_3weeks_rp48_nsg$peaks[1], y = res_3weeks_rp48_nsg$peak_heights[1], label = paste0(round(den_p4$comp1_perc[1], 2), "%"), vjust = -1) +
    geom_text(x = res_3weeks_rp48_nsg$peaks[2], y = res_3weeks_rp48_nsg$peak_heights[2], label = paste0(round(den_p4$comp2_perc[1], 2), "%"), vjust = -1) +
    labs(
        title = "Density Plot of Cell Number for Crebbp in Liver",
        x = "Log2 (Cell Number)",
        y = "Density"
    ) + 
    theme_minimal()

# Consider Dormancy with Sex
# Male
df_1week_rp48_bl6_m <- df %>% filter(target == "Crebbp", Tissue == "Liver", Cell_Line == "rp48", Time_Point == "1week", Mouse_Genotype == "BL6", Sex == "M")
df_1week_rp48_nsg_m <- df %>% filter(target == "Crebbp", Tissue == "Liver", Cell_Line == "rp48", Time_Point == "1week", Mouse_Genotype == "NSG", Sex == "M")
df_1week_rp116_bl6_m <- df %>% filter(target == "Crebbp", Tissue == "Liver", Cell_Line == "rp116", Time_Point == "1week", Mouse_Genotype == "BL6", Sex == "M")
df_1week_rp116_nsg_m <- df %>% filter(target == "Crebbp", Tissue == "Liver", Cell_Line == "rp116", Time_Point == "1week", Mouse_Genotype == "NSG", Sex == "M")
df_3weeks_rp48_bl6_m <- df %>% filter(target == "Crebbp", Tissue == "Liver", Cell_Line == "rp48", Time_Point == "3weeks", Mouse_Genotype == "BL6", Sex == "M")
df_3weeks_rp48_nsg_m <- df %>% filter(target == "Crebbp", Tissue == "Liver", Cell_Line == "rp48", Time_Point == "3weeks", Mouse_Genotype == "NSG", Sex == "M")
df_3weeks_rp116_bl6_m <- df %>% filter(target == "Crebbp", Tissue == "Liver", Cell_Line == "rp116", Time_Point == "3weeks", Mouse_Genotype == "BL6", Sex == "M")
df_3weeks_rp116_nsg_m <- df %>% filter(target == "Crebbp", Tissue == "Liver", Cell_Line == "rp116", Time_Point == "3weeks", Mouse_Genotype == "NSG", Sex == "M")

# Female
df_1week_rp48_bl6_f <- df %>% filter(target == "Crebbp", Tissue == "Liver", Cell_Line == "rp48", Time_Point == "1week", Mouse_Genotype == "BL6", Sex == "F")
df_1week_rp48_nsg_f <- df %>% filter(target == "Crebbp", Tissue == "Liver", Cell_Line == "rp48", Time_Point == "1week", Mouse_Genotype == "NSG", Sex == "F")
df_1week_rp116_bl6_f <- df %>% filter(target == "Crebbp", Tissue == "Liver", Cell_Line == "rp116", Time_Point == "1week", Mouse_Genotype == "BL6", Sex == "F")
df_1week_rp116_nsg_f <- df %>% filter(target == "Crebbp", Tissue == "Liver", Cell_Line == "rp116", Time_Point == "1week", Mouse_Genotype == "NSG", Sex == "F")
df_3weeks_rp48_bl6_f <- df %>% filter(target == "Crebbp", Tissue == "Liver", Cell_Line == "rp48", Time_Point == "3weeks", Mouse_Genotype == "BL6", Sex == "F")
df_3weeks_rp48_nsg_f <- df %>% filter(target == "Crebbp", Tissue == "Liver", Cell_Line == "rp48", Time_Point == "3weeks", Mouse_Genotype == "NSG", Sex == "F")
df_3weeks_rp116_bl6_f <- df %>% filter(target == "Crebbp", Tissue == "Liver", Cell_Line == "rp116", Time_Point == "3weeks", Mouse_Genotype == "BL6", Sex == "F")
df_3weeks_rp116_nsg_f <- df %>% filter(target == "Crebbp", Tissue == "Liver", Cell_Line == "rp116", Time_Point == "3weeks", Mouse_Genotype == "NSG", Sex == "F")

# Check for bimodality
# Male
res_1week_rp48_bl6_m <- check_bimodal(df_1week_rp48_bl6_m$cell_num)
res_1week_rp48_nsg_m <- check_bimodal(df_1week_rp48_nsg_m$cell_num)
res_1week_rp116_bl6_m <- check_bimodal(df_1week_rp116_bl6_m$cell_num)
res_1week_rp116_nsg_m <- check_bimodal(df_1week_rp116_nsg_m$cell_num)
res_3weeks_rp48_bl6_m <- check_bimodal(df_3weeks_rp48_bl6_m$cell_num)
res_3weeks_rp48_nsg_m <- check_bimodal(df_3weeks_rp48_nsg_m$cell_num)
res_3weeks_rp116_bl6_m <- check_bimodal(df_3weeks_rp116_bl6_m$cell_num)
res_3weeks_rp116_nsg_m <- check_bimodal(df_3weeks_rp116_nsg_m$cell_num)

#Female
res_1week_rp48_bl6_f <- check_bimodal(df_1week_rp48_bl6_f$cell_num)
res_1week_rp48_nsg_f <- check_bimodal(df_1week_rp48_nsg_f$cell_num)
res_1week_rp116_bl6_f <- check_bimodal(df_1week_rp116_bl6_f$cell_num)
res_1week_rp116_nsg_f <- check_bimodal(df_1week_rp116_nsg_f$cell_num)
res_3weeks_rp48_bl6_f <- check_bimodal(df_3weeks_rp48_bl6_f$cell_num)
res_3weeks_rp48_nsg_f <- check_bimodal(df_3weeks_rp48_nsg_f$cell_num)
res_3weeks_rp116_bl6_f <- check_bimodal(df_3weeks_rp116_bl6_f$cell_num)
res_3weeks_rp116_nsg_f <- check_bimodal(df_3weeks_rp116_nsg_f$cell_num)

# Make the plots
# Male
den_p1_m <- plot_bimodal(df_3weeks_rp116_bl6_m$cell_num, res_3weeks_rp116_bl6_m, show_fit = FALSE, show_valleys = FALSE, return_data = TRUE)
den_p2_m <- plot_bimodal(df_3weeks_rp48_bl6_m$cell_num, res_3weeks_rp48_bl6_m, show_fit = FALSE, show_valleys = FALSE, return_data = TRUE)
den_p3_m <- plot_bimodal(df_3weeks_rp116_nsg_m$cell_num, res_3weeks_rp116_nsg_m, show_fit = FALSE, show_valleys = FALSE, return_data = TRUE)
den_p4_m <- plot_bimodal(df_3weeks_rp48_nsg_m$cell_num, res_3weeks_rp48_nsg_m, show_fit = FALSE, show_valleys = FALSE, return_data = TRUE)

# Female
den_p1_f <- plot_bimodal(df_3weeks_rp116_bl6_f$cell_num, res_3weeks_rp116_bl6_f, show_fit = FALSE, show_valleys = FALSE, return_data = TRUE)
den_p2_f <- plot_bimodal(df_3weeks_rp48_bl6_f$cell_num, res_3weeks_rp48_bl6_f, show_fit = FALSE, show_valleys = FALSE, return_data = TRUE)
den_p3_f <- plot_bimodal(df_3weeks_rp116_nsg_f$cell_num, res_3weeks_rp116_nsg_f, show_fit = FALSE, show_valleys = FALSE, return_data = TRUE)
den_p4_f <- plot_bimodal(df_3weeks_rp48_nsg_f$cell_num, res_3weeks_rp48_nsg_f, show_fit = FALSE, show_valleys = FALSE, return_data = TRUE)

to_plot <- rbind(
    den_p1_m %>% mutate(Cell_Line = "rp116", Mouse_Genotype = "BL6", Sex = "M"),
    den_p2_m %>% mutate(Cell_Line = "rp48", Mouse_Genotype = "BL6", Sex = "M"),
    den_p3_m %>% mutate(Cell_Line = "rp116", Mouse_Genotype = "NSG", Sex = "M"),
    den_p4_m %>% mutate(Cell_Line = "rp48", Mouse_Genotype = "NSG", Sex = "M"),
    den_p1_f %>% mutate(Cell_Line = "rp116", Mouse_Genotype = "BL6", Sex = "F"),
    den_p2_f %>% mutate(Cell_Line = "rp48", Mouse_Genotype = "BL6", Sex = "F"),
    den_p3_f %>% mutate(Cell_Line = "rp116", Mouse_Genotype = "NSG", Sex = "F"),
    den_p4_f %>% mutate(Cell_Line = "rp48", Mouse_Genotype = "NSG", Sex = "F")
)

rp116_sex_plot <- to_plot %>%
    filter(Cell_Line == "rp116") %>%
    mutate(
        group = paste(Mouse_Genotype, Sex)
    ) %>%
    ggplot(
        aes(
            x = x,
            fill = group
        )
    ) +
    geom_density(alpha = 0.3) +
    labs(
        title = "Density Plot of Cell Number for Crebbp in Liver - rp116",
        x = "Log2 (Cell Number)",
        y = "Density"
    ) +
    theme_minimal()

rp48_sex_plot <- to_plot %>%
    filter(Cell_Line == "rp48") %>%
    mutate(
        group = paste(Mouse_Genotype, Sex)
    ) %>%
    ggplot(
        aes(
            x = x,
            fill = group
        )
    ) +
    geom_density(alpha = 0.3) +
    labs(
        title = "Density Plot of Cell Number for Crebbp in Liver - rp48",
        x = "Log2 (Cell Number)",
        y = "Density"
    ) +
    theme_minimal()

rp116_sex_plot + 
    # BL6 RP116 Male
    geom_text(x = res_3weeks_rp116_bl6_m$peaks[1], y = res_3weeks_rp116_bl6_m$peak_heights[1], label = paste0(round(den_p1_m$comp1_perc[1], 2), "%"), vjust = -1) +
    geom_text(x = res_3weeks_rp116_bl6_m$peaks[2], y = res_3weeks_rp116_bl6_m$peak_heights[2], label = paste0(round(den_p1_m$comp2_perc[1], 2), "%"), vjust = -1) +
    # BL6 RP116 Female
    geom_text(x = res_3weeks_rp116_bl6_f$peaks[1], y = res_3weeks_rp116_bl6_f$peak_heights[1], label = paste0(round(den_p1_f$comp1_perc[1], 2), "%"), vjust = -1) +
    geom_text(x = res_3weeks_rp116_bl6_f$peaks[2], y = res_3weeks_rp116_bl6_f$peak_heights[2], label = paste0(round(den_p1_f$comp2_perc[1], 2), "%"), vjust = -1) +
    # NSG RP116 Male
    geom_text(x = res_3weeks_rp116_nsg_m$peaks[1], y = res_3weeks_rp116_nsg_m$peak_heights[1], label = paste0(round(den_p3_m$comp1_perc[1], 2), "%"), vjust = -1) +
    geom_text(x = res_3weeks_rp116_nsg_m$peaks[2], y = res_3weeks_rp116_nsg_m$peak_heights[2], label = paste0(round(den_p3_m$comp2_perc[1], 2), "%"), vjust = -1) +
    # NSG RP116 Female
    geom_text(x = res_3weeks_rp116_nsg_f$peaks[1], y = res_3weeks_rp116_nsg_f$peak_heights[1], label = paste0(round(den_p3_f$comp1_perc[1], 2), "%"), vjust = -1) +
    geom_text(x = res_3weeks_rp116_nsg_f$peaks[2], y = res_3weeks_rp116_nsg_f$peak_heights[2], label = paste0(round(den_p3_f$comp2_perc[1], 2), "%"), vjust = -1) +
    theme_minimal()
rp48_sex_plot + 
    # BL6 RP48 Male
    geom_text(x = res_3weeks_rp48_bl6_m$peaks[1], y = res_3weeks_rp48_bl6_m$peak_heights[1], label = paste0(round(den_p2_m$comp1_perc[1], 2), "%"), vjust = -1) +
    geom_text(x = res_3weeks_rp48_bl6_m$peaks[2], y = res_3weeks_rp48_bl6_m$peak_heights[2], label = paste0(round(den_p2_m$comp2_perc[1], 2), "%"), vjust = -1) +
    # BL6 RP48 Female
    geom_text(x = res_3weeks_rp48_bl6_f$peaks[1], y = res_3weeks_rp48_bl6_f$peak_heights[1], label = paste0(round(den_p2_f$comp1_perc[1], 2), "%"), vjust = -1) +
    geom_text(x = res_3weeks_rp48_bl6_f$peaks[2], y = res_3weeks_rp48_bl6_f$peak_heights[2], label = paste0(round(den_p2_f$comp2_perc[1], 2), "%"), vjust = -1) +
    # NSG RP48 Male
    geom_text(x = res_3weeks_rp48_nsg_m$peaks[1], y = res_3weeks_rp48_nsg_m$peak_heights[1], label = paste0(round(den_p4_m$comp1_perc[1], 2), "%"), vjust = -1) +
    geom_text(x = res_3weeks_rp48_nsg_m$peaks[2], y = res_3weeks_rp48_nsg_m$peak_heights[2], label = paste0(round(den_p4_m$comp2_perc[1], 2), "%"), vjust = -1) +
    # NSG RP48 Female
    geom_text(x = res_3weeks_rp48_nsg_f$peaks[1], y = res_3weeks_rp48_nsg_f$peak_heights[1], label = paste0(round(den_p4_f$comp1_perc[1], 2), "%"), vjust = -1) +
    geom_text(x = res_3weeks_rp48_nsg_f$peaks[2], y = res_3weeks_rp48_nsg_f$peak_heights[2], label = paste0(round(den_p4_f$comp2_perc[1], 2), "%"), vjust = -1) +
    theme_minimal()

# Do the same for sgSafe
# Male
safe_df_1week_rp48_bl6_m <- df %>% filter(target == "sgSafe", Tissue == "Liver", Cell_Line == "rp48", Time_Point == "1week", Mouse_Genotype == "BL6", Sex == "M")
safe_df_1week_rp48_nsg_m <- df %>% filter(target == "sgSafe", Tissue == "Liver", Cell_Line == "rp48", Time_Point == "1week", Mouse_Genotype == "NSG", Sex == "M")
safe_df_1week_rp116_bl6_m <- df %>% filter(target == "sgSafe", Tissue == "Liver", Cell_Line == "rp116", Time_Point == "1week", Mouse_Genotype == "BL6", Sex == "M")
safe_df_1week_rp116_nsg_m <- df %>% filter(target == "sgSafe", Tissue == "Liver", Cell_Line == "rp116", Time_Point == "1week", Mouse_Genotype == "NSG", Sex == "M")
safe_df_3weeks_rp48_bl6_m <- df %>% filter(target == "sgSafe", Tissue == "Liver", Cell_Line == "rp48", Time_Point == "3weeks", Mouse_Genotype == "BL6", Sex == "M")
safe_df_3weeks_rp48_nsg_m <- df %>% filter(target == "sgSafe", Tissue == "Liver", Cell_Line == "rp48", Time_Point == "3weeks", Mouse_Genotype == "NSG", Sex == "M")
safe_df_3weeks_rp116_bl6_m <- df %>% filter(target == "sgSafe", Tissue == "Liver", Cell_Line == "rp116", Time_Point == "3weeks", Mouse_Genotype == "BL6", Sex == "M")
safe_df_3weeks_rp116_nsg_m <- df %>% filter(target == "sgSafe", Tissue == "Liver", Cell_Line == "rp116", Time_Point == "3weeks", Mouse_Genotype == "NSG", Sex == "M")

# Female
safe_df_1week_rp48_bl6_f <- df %>% filter(target == "sgSafe", Tissue == "Liver", Cell_Line == "rp48", Time_Point == "1week", Mouse_Genotype == "BL6", Sex == "F")
safe_df_1week_rp48_nsg_f <- df %>% filter(target == "sgSafe", Tissue == "Liver", Cell_Line == "rp48", Time_Point == "1week", Mouse_Genotype == "NSG", Sex == "F")
safe_df_1week_rp116_bl6_f <- df %>% filter(target == "sgSafe", Tissue == "Liver", Cell_Line == "rp116", Time_Point == "1week", Mouse_Genotype == "BL6", Sex == "F")
safe_df_1week_rp116_nsg_f <- df %>% filter(target == "sgSafe", Tissue == "Liver", Cell_Line == "rp116", Time_Point == "1week", Mouse_Genotype == "NSG", Sex == "F")
safe_df_3weeks_rp48_bl6_f <- df %>% filter(target == "sgSafe", Tissue == "Liver", Cell_Line == "rp48", Time_Point == "3weeks", Mouse_Genotype == "BL6", Sex == "F")
safe_df_3weeks_rp48_nsg_f <- df %>% filter(target == "sgSafe", Tissue == "Liver", Cell_Line == "rp48", Time_Point == "3weeks", Mouse_Genotype == "NSG", Sex == "F")
safe_df_3weeks_rp116_bl6_f <- df %>% filter(target == "sgSafe", Tissue == "Liver", Cell_Line == "rp116", Time_Point == "3weeks", Mouse_Genotype == "BL6", Sex == "F")
safe_df_3weeks_rp116_nsg_f <- df %>% filter(target == "sgSafe", Tissue == "Liver", Cell_Line == "rp116", Time_Point == "3weeks", Mouse_Genotype == "NSG", Sex == "F")

# Check for bimodality
# Male
safe_res_1week_rp48_bl6_m <- check_bimodal(safe_df_1week_rp48_bl6_m$cell_num)
safe_res_1week_rp48_nsg_m <- check_bimodal(safe_df_1week_rp48_nsg_m$cell_num)
safe_res_1week_rp116_bl6_m <- check_bimodal(safe_df_1week_rp116_bl6_m$cell_num)
safe_res_1week_rp116_nsg_m <- check_bimodal(safe_df_1week_rp116_nsg_m$cell_num)
safe_res_3weeks_rp48_bl6_m <- check_bimodal(safe_df_3weeks_rp48_bl6_m$cell_num)
safe_res_3weeks_rp48_nsg_m <- check_bimodal(safe_df_3weeks_rp48_nsg_m$cell_num)
safe_res_3weeks_rp116_bl6_m <- check_bimodal(safe_df_3weeks_rp116_bl6_m$cell_num)
safe_res_3weeks_rp116_nsg_m <- check_bimodal(safe_df_3weeks_rp116_nsg_m$cell_num)

# Female
safe_res_1week_rp48_bl6_f <- check_bimodal(safe_df_1week_rp48_bl6_f$cell_num)
safe_res_1week_rp48_nsg_f <- check_bimodal(safe_df_1week_rp48_nsg_f$cell_num)
safe_res_1week_rp116_bl6_f <- check_bimodal(safe_df_1week_rp116_bl6_f$cell_num)
safe_res_1week_rp116_nsg_f <- check_bimodal(safe_df_1week_rp116_nsg_f$cell_num)
safe_res_3weeks_rp48_bl6_f <- check_bimodal(safe_df_3weeks_rp48_bl6_f$cell_num)
safe_res_3weeks_rp48_nsg_f <- check_bimodal(safe_df_3weeks_rp48_nsg_f$cell_num)
safe_res_3weeks_rp116_bl6_f <- check_bimodal(safe_df_3weeks_rp116_bl6_f$cell_num)
safe_res_3weeks_rp116_nsg_f <- check_bimodal(safe_df_3weeks_rp116_nsg_f$cell_num)

# Male
safe_den_p1_m <- plot_bimodal(safe_df_3weeks_rp116_bl6_m$cell_num, safe_res_3weeks_rp116_bl6_m, show_fit = FALSE, show_valleys = FALSE, return_data = TRUE)
safe_den_p2_m <- plot_bimodal(safe_df_3weeks_rp48_bl6_m$cell_num, safe_res_3weeks_rp48_bl6_m, show_fit = FALSE, show_valleys = FALSE, return_data = TRUE)
safe_den_p3_m <- plot_bimodal(safe_df_3weeks_rp116_nsg_m$cell_num, safe_res_3weeks_rp116_nsg_m, show_fit = FALSE, show_valleys = FALSE, return_data = TRUE)
safe_den_p4_m <- plot_bimodal(safe_df_3weeks_rp48_nsg_m$cell_num, safe_res_3weeks_rp48_nsg_m, show_fit = FALSE, show_valleys = FALSE, return_data = TRUE)
#plot_bimodal(safe_df_3weeks_rp116_nsg_m$cell_num, safe_res_3weeks_rp116_nsg_m, show_fit = FALSE, show_valleys = FALSE, return_data = F)

# Female
safe_den_p1_f <- plot_bimodal(safe_df_3weeks_rp116_bl6_f$cell_num, safe_res_3weeks_rp116_bl6_f, show_fit = FALSE, show_valleys = FALSE, return_data = TRUE)
safe_den_p2_f <- plot_bimodal(safe_df_3weeks_rp48_bl6_f$cell_num, safe_res_3weeks_rp48_bl6_f, show_fit = FALSE, show_valleys = FALSE, return_data = TRUE)
safe_den_p3_f <- plot_bimodal(safe_df_3weeks_rp116_nsg_f$cell_num, safe_res_3weeks_rp116_nsg_f, show_fit = FALSE, show_valleys = FALSE, return_data = TRUE)
safe_den_p4_f <- plot_bimodal(safe_df_3weeks_rp48_nsg_f$cell_num, safe_res_3weeks_rp48_nsg_f, show_fit = FALSE, show_valleys = FALSE, return_data = TRUE)

# The Plot
to_plot <- rbind(
    safe_den_p1_m %>% mutate(Cell_Line = "rp116", Mouse_Genotype = "BL6", Sex = "M"),
    safe_den_p2_m %>% mutate(Cell_Line = "rp48", Mouse_Genotype = "BL6", Sex = "M"),
    safe_den_p3_m %>% mutate(Cell_Line = "rp116", Mouse_Genotype = "NSG", Sex = "M"),
    safe_den_p4_m %>% mutate(Cell_Line = "rp48", Mouse_Genotype = "NSG", Sex = "M"),
    safe_den_p1_f %>% mutate(Cell_Line = "rp116", Mouse_Genotype = "BL6", Sex = "F"),
    safe_den_p2_f %>% mutate(Cell_Line = "rp48", Mouse_Genotype = "BL6", Sex = "F"),
    safe_den_p3_f %>% mutate(Cell_Line = "rp116", Mouse_Genotype = "NSG", Sex = "F"),
    safe_den_p4_f %>% mutate(Cell_Line = "rp48", Mouse_Genotype = "NSG", Sex = "F")
)

safe_rp116_sex_plot <- to_plot %>%
    filter(Cell_Line == "rp116") %>%
    mutate(
        group = paste(Mouse_Genotype, Sex)
    ) %>%
    ggplot(
        aes(
            x = x,
            fill = group
        )
    ) +
    geom_density(alpha = 0.3) +
    labs(
        title = "Density Plot of Cell Number for sgSafe in Liver - rp116",
        x = "Log2 (Cell Number)",
        y = "Density"
    ) +
    theme_minimal()

safe_rp48_sex_plot <- to_plot %>%
    filter(Cell_Line == "rp48") %>%
    mutate(
        group = paste(Mouse_Genotype, Sex)
    ) %>%
    ggplot(
        aes(
            x = x,
            fill = group
        )
    ) +
    geom_density(alpha = 0.3) +
    labs(
        title = "Density Plot of Cell Number for sgSafe in Liver - rp48",
        x = "Log2 (Cell Number)",
        y = "Density"
    ) +
    theme_minimal()

safe_rp116_sex_plot +
    # BL6 RP116 Male
    geom_text(x = safe_res_3weeks_rp116_bl6_m$peaks[1], y = safe_res_3weeks_rp116_bl6_m$peak_heights[1], label = paste0(round(safe_den_p1_m$comp1_perc[1], 2), "%"), vjust = -1) +
    geom_text(x = safe_res_3weeks_rp116_bl6_m$peaks[2], y = safe_res_3weeks_rp116_bl6_m$peak_heights[2], label = paste0(round(safe_den_p1_m$comp2_perc[1], 2), "%"), vjust = -1) +
    # BL6 RP116 Female
    geom_text(x = safe_res_3weeks_rp116_bl6_f$peaks[1], y = safe_res_3weeks_rp116_bl6_f$peak_heights[1], label = paste0(round(safe_den_p1_f$comp1_perc[1], 2), "%"), vjust = -1) +
    geom_text(x = safe_res_3weeks_rp116_bl6_f$peaks[2], y = safe_res_3weeks_rp116_bl6_f$peak_heights[2], label = paste0(round(safe_den_p1_f$comp2_perc[1], 2), "%"), vjust = -1) +
    # NSG RP116 Male
    geom_text(x = safe_res_3weeks_rp116_nsg_m$peaks[1], y = safe_res_3weeks_rp116_nsg_m$peak_heights[1], label = paste0(round(safe_den_p3_m$comp1_perc[1], 2), "%"), vjust = -1) +
    geom_text(x = safe_res_3weeks_rp116_nsg_m$peaks[4], y = safe_res_3weeks_rp116_nsg_m$peak_heights[4], label = paste0(round(safe_den_p3_m$comp2_perc[1], 2), "%"), vjust = -1) +
    # NSG RP116 Female
    geom_text(x = safe_res_3weeks_rp116_nsg_f$peaks[1], y = safe_res_3weeks_rp116_nsg_f$peak_heights[1], label = paste0(round(safe_den_p3_f$comp1_perc[1], 2), "%"), vjust = -1) +
    geom_text(x = safe_res_3weeks_rp116_nsg_f$peaks[2], y = safe_res_3weeks_rp116_nsg_f$peak_heights[2], label = paste0(round(safe_den_p3_f$comp2_perc[1], 2), "%"), vjust = -1) +
    theme_minimal()

safe_rp48_sex_plot +
    # BL6 RP48 Male
    geom_text(x = safe_res_3weeks_rp48_bl6_m$peaks[1], y = safe_res_3weeks_rp48_bl6_m$peak_heights[1], label = paste0(round(safe_den_p2_m$comp1_perc[1], 2), "%"), vjust = -1) +
    geom_text(x = safe_res_3weeks_rp48_bl6_m$peaks[2], y = safe_res_3weeks_rp48_bl6_m$peak_heights[2], label = paste0(round(safe_den_p2_m$comp2_perc[1], 2), "%"), vjust = -1) +
    # BL6 RP48 Female
    geom_text(x = safe_res_3weeks_rp48_bl6_f$peaks[1], y = safe_res_3weeks_rp48_bl6_f$peak_heights[1], label = paste0(round(safe_den_p2_f$comp1_perc[1], 2), "%"), vjust = -1) +
    geom_text(x = safe_res_3weeks_rp48_bl6_f$peaks[2], y = safe_res_3weeks_rp48_bl6_f$peak_heights[2], label = paste0(round(safe_den_p2_f$comp2_perc[1], 2), "%"), vjust = -1) +
    # NSG RP48 Male
    geom_text(x = safe_res_3weeks_rp48_nsg_m$peaks[1], y = safe_res_3weeks_rp48_nsg_m$peak_heights[1], label = paste0(round(safe_den_p4_m$comp1_perc[1], 2), "%"), vjust = -1) +
    geom_text(x = safe_res_3weeks_rp48_nsg_m$peaks[2], y = safe_res_3weeks_rp48_nsg_m$peak_heights[2], label = paste0(round(safe_den_p4_m$comp2_perc[1], 2), "%"), vjust = -1) +
    # NSG RP48 Female
    geom_text(x = safe_res_3weeks_rp48_nsg_f$peaks[1], y = safe_res_3weeks_rp48_nsg_f$peak_heights[1], label = paste0(round(safe_den_p4_f$comp1_perc[1], 2), "%"), vjust = -1) +
    geom_text(x = safe_res_3weeks_rp48_nsg_f$peaks[2], y = safe_res_3weeks_rp48_nsg_f$peak_heights[2], label = paste0(round(safe_den_p4_f$comp2_perc[1], 2), "%"), vjust = -1) +
    theme_minimal()

# Plot for Rui
large_den_df <- rbind(
    den_p1_m %>% mutate(target = "Crebbp", Cell_Line = "rp116", Mouse_Genotype = "BL6", Sex = "M"),
    den_p2_m %>% mutate(target = "Crebbp", Cell_Line = "rp48", Mouse_Genotype = "BL6", Sex = "M"),
    den_p3_m %>% mutate(target = "Crebbp", Cell_Line = "rp116", Mouse_Genotype = "NSG", Sex = "M"),
    den_p4_m %>% mutate(target = "Crebbp", Cell_Line = "rp48", Mouse_Genotype = "NSG", Sex = "M"),
    den_p1_f %>% mutate(target = "Crebbp", Cell_Line = "rp116", Mouse_Genotype = "BL6", Sex = "F"),
    den_p2_f %>% mutate(target = "Crebbp", Cell_Line = "rp48", Mouse_Genotype = "BL6", Sex = "F"),
    den_p3_f %>% mutate(target = "Crebbp", Cell_Line = "rp116", Mouse_Genotype = "NSG", Sex = "F"),
    den_p4_f %>% mutate(target = "Crebbp", Cell_Line = "rp48", Mouse_Genotype = "NSG", Sex = "F"),
    safe_den_p1_m %>% mutate(target = "sgSafe", Cell_Line = "rp116", Mouse_Genotype = "BL6", Sex = "M"),
    safe_den_p2_m %>% mutate(target = "sgSafe", Cell_Line = "rp48", Mouse_Genotype = "BL6", Sex = "M"),
    safe_den_p3_m %>% mutate(target = "sgSafe", Cell_Line = "rp116", Mouse_Genotype = "NSG", Sex = "M"),
    safe_den_p4_m %>% mutate(target = "sgSafe", Cell_Line = "rp48", Mouse_Genotype = "NSG", Sex = "M"),
    safe_den_p1_f %>% mutate(target = "sgSafe", Cell_Line = "rp116", Mouse_Genotype = "BL6", Sex = "F"),
    safe_den_p2_f %>% mutate(target = "sgSafe", Cell_Line = "rp48", Mouse_Genotype = "BL6", Sex = "F"),
    safe_den_p3_f %>% mutate(target = "sgSafe", Cell_Line = "rp116", Mouse_Genotype = "NSG", Sex = "F"),
    safe_den_p4_f %>% mutate(target = "sgSafe", Cell_Line = "rp48", Mouse_Genotype = "NSG", Sex = "F")
)

peaks_df <- large_den_df %>%
    group_by(Mouse_Genotype, Sex, Cell_Line, target) %>%
    summarize(
        peak1 = mean(peak1),
        peak2 = mean(peak2),
        peak1_percentage = mean(comp1_perc),
        peak2_percentage = mean(comp2_perc),
        peak1_height = mean(peak1_height),
        peak2_height = mean(peak2_height)
    ) %>%
    mutate(
        group = paste(Mouse_Genotype, Sex, sep = " "),
        fill_group = case_when(
            target == "sgSafe" ~ "sgSafe",
            target == "Crebbp" & Mouse_Genotype == "BL6" ~ "Crebbp_BL6",
            target == "Crebbp" & Mouse_Genotype == "NSG" ~ "Crebbp_NSG"
        )
    ) %>%
    filter(Mouse_Genotype != "BL6")

den_plot_nsg <- large_den_df %>%
    mutate(
        group = paste(Mouse_Genotype, Sex, sep = " "),
        fill_group = case_when(
            target == "sgSafe" ~ "sgSafe",
            target == "Crebbp" & Mouse_Genotype == "BL6" ~ "Crebbp_BL6",
            target == "Crebbp" & Mouse_Genotype == "NSG" ~ "Crebbp_NSG"
        )
    ) %>%
    filter(Mouse_Genotype != "BL6") %>%
    ggplot(
        aes(
            x = x,
            fill = fill_group
        )
    ) +
    geom_density(alpha = 0.8) +
    geom_vline(
        data = peaks_df,
        aes(xintercept = peak1),
        linetype = "dashed",
        alpha = 0.3
    ) +
    geom_text_repel(
        data = peaks_df,
        aes(x = peak1, 
            y = peak1_height, 
            label = str_c(
                round(peak1_percentage, 2), "%\n",
                format(round(2^peak1), big.mark = ",")
            )
        ),
        vjust = -1
    ) +
    geom_vline(
        data = peaks_df,
        aes(xintercept = peak2),
        linetype = "dashed",
        alpha = 0.3
    ) +
    geom_text_repel(
        data = peaks_df,
        aes(x = peak2, 
            y = peak2_height, 
            label = str_c(
                round(peak2_percentage, 2), "%\n",
                format(round(2^peak2), big.mark = ",")
            )
        ),
        vjust = -1
    ) +
    facet_wrap(~ group + Cell_Line, nrow = 4, ncol = 2) +
    labs(
        title = "Density Plot of Cell Number for Crebbp and sgSafe in Liver",
        x = "Log2 (Cell Number)",
        y = "Density"
    ) +
    theme_bw() +
    ylim(NA, 0.25) +
        scale_fill_manual(
        name = "Time Point",
        values = c(
            "Crebbp_BL6" = "#0066CC",
            "Crebbp_NSG" = "#E60000",
            "sgSafe" = "#666666"
        ),
        labels = c(
            "Crebbp_BL6" = "Crebbp - BL6",
            "Crebbp_NSG" = "Crebbp - NSG",
            "sgSafe" = "sgSafe"
        )
    )
pdf("density_plot_liver_nsg_CrebbpVssgSafe.pdf", width = 10, height = 10)
print(den_plot_nsg)
dev.off()

peaks_df <- large_den_df %>%
    group_by(Mouse_Genotype, Sex, Cell_Line, target) %>%
    summarize(
        peak1 = mean(peak1),
        peak2 = mean(peak2),
        peak1_percentage = mean(comp1_perc),
        peak2_percentage = mean(comp2_perc),
        peak1_height = mean(peak1_height),
        peak2_height = mean(peak2_height)
    ) %>%
    mutate(
        peak2 = ifelse(peak2 == 0, NA, peak2),
        peak2_percentage = ifelse(is.na(peak2), NA, peak2_percentage),
        peak2_height = ifelse(is.na(peak2), NA, peak2_height),
        group = paste(Mouse_Genotype, Sex, sep = " "),
        fill_group = case_when(
            target == "sgSafe" ~ "sgSafe",
            target == "Crebbp" & Mouse_Genotype == "BL6" ~ "Crebbp_BL6",
            target == "Crebbp" & Mouse_Genotype == "NSG" ~ "Crebbp_NSG"
        )
    ) %>%
    filter(Mouse_Genotype != "NSG")

den_plot_bl6 <- large_den_df %>%
    mutate(
        group = paste(Mouse_Genotype, Sex, sep = " "),
        fill_group = case_when(
            target == "sgSafe" ~ "sgSafe",
            target == "Crebbp" & Mouse_Genotype == "BL6" ~ "Crebbp_BL6",
            target == "Crebbp" & Mouse_Genotype == "NSG" ~ "Crebbp_NSG"
        )
    ) %>%
    filter(Mouse_Genotype != "NSG") %>%
    ggplot(
        aes(
            x = x,
            fill = fill_group
        )
    ) +
    geom_density(alpha = 0.8) +
    geom_vline(
        data = peaks_df,
        aes(xintercept = peak1),
        linetype = "dashed",
        alpha = 0.3
    ) +
    geom_text_repel(
        data = peaks_df,
        aes(
            x = peak1,
            y = peak1_height,
            label = str_c(
                round(peak1_percentage, 2), "%\n",
                format(round(2^peak1), big.mark = ",")
            )
        ),
        vjust = -1
    ) +
    geom_vline(
        data = peaks_df,
        aes(xintercept = peak2),
        linetype = "dashed",
        alpha = 0.3
    ) +
    geom_text_repel(
        data = peaks_df,
        aes(
            x = peak2,
            y = peak2_height,
            label = str_c(
                round(peak2_percentage, 2), "%\n",
                format(round(2^peak2), big.mark = ",")
            )
        ),
        vjust = -1
    ) +
    facet_wrap(~ group + Cell_Line, nrow = 4, ncol = 2) +
    labs(
        title = "Density Plot of Cell Number for Crebbp and sgSafe in Liver",
        x = "Log2 (Cell Number)",
        y = "Density"
    ) +
    theme_bw() +
    ylim(NA, 0.25) +
    scale_fill_manual(
        name = "Time Point",
        values = c(
            "Crebbp_BL6" = "#0066CC",
            "Crebbp_NSG" = "#E60000",
            "sgSafe" = "#666666"
        ),
        labels = c(
            "Crebbp_BL6" = "Crebbp - BL6",
            "Crebbp_NSG" = "Crebbp - NSG",
            "sgSafe" = "sgSafe"
        )
    )
pdf("density_plot_liver_bl6_CrebbpVssgSafe.pdf", width = 10, height = 10)
print(den_plot_bl6)
dev.off()

# Crebbp Specific Analysis - For Gene, For Tissue, For Genotype, For Cell Line
# Notes: We have 2 different Cell Lines, but we are interested to see what environmental factors may affect the 5 given statistics
# Sex, Genotype
# For EVERY gene, in every Tissue, for BOTH Cell Lines
# We specify which Tissue we want, Liver first as a exploration
run_glm_models <- function(data, statistic, tissue = "Liver", formula_type = "Sex", time_point = "3weeks") {
    map_dfr(unique(data$target), function(sg) {       # This will do this for EVERY gene
        print(sg)
        results_list <- list()
        # Of course for each cell line we need to run the test twice
        for (cell_line in c("rp48", "rp116")) {
            print(cell_line)
            cell_data <- data %>%
                filter(
                    Time_Point == time_point,     # Time point is fixed, but we COULD run it for 1week in Liver just to see
                    Tissue == tissue,
                    Cell_Line == cell_line,
                    target == sg,
                    !is.na(!!sym(statistic)),
                    !!sym(statistic) != Inf,
                    !!sym(statistic) != -Inf
                )
            # Cell Data is now the dataframe specifically for that gene

            if (nrow(cell_data) > 0) {
                if (formula_type == "Sex") {
                    # If it's Sex, the variable needs to be swapped between Genotype...
                    for (genotype in c("BL6", "NSG")) {
                        geno_data <- cell_data %>%
                            filter(Mouse_Genotype == genotype)

                        if (nrow(geno_data) > 1 && length(unique(geno_data$Sex)) > 1) {
                            formula <- as.formula(paste(statistic, "~ Sex"))
                            results_list[[paste0(cell_line, "_", tolower(genotype))]] <-
                                lm(formula, data = geno_data) %>%
                                sjPlot::get_model_data(type = "est") %>%
                                mutate(
                                    Mouse_Genotype = genotype,
                                    Cell_Line = cell_line,
                                    Statistic = !!statistic
                                )
                        }
                    }
                } else if (formula_type == "Mouse_Genotype") {
                    # If it's Genotype, the variable needs to be swapped between Sex...
                    for (sex in c("F", "M")) {
                        sex_data <- cell_data %>%
                            filter(Sex == sex)
                        
                        if (nrow(sex_data) > 1 && length(unique(sex_data$Mouse_Genotype)) > 1) {
                            formula <- as.formula(paste(statistic, "~ Mouse_Genotype"))
                            results_list[[paste0(cell_line, "_", tolower(sex))]] <-
                                lm(formula, data = sex_data) %>%
                                sjPlot::get_model_data(type = "est") %>%
                                mutate(
                                    Sex = sex,
                                    Cell_Line = cell_line,
                                    Statistic = !!statistic
                                )
                        }
                    }
                }
            }
        }

        if (length(results_list) > 0) {
            results_df <- bind_rows(results_list)
            results_df$gene <- sg
            return(results_df)
        }
        return(NULL)
    })
}

# Run models for statistics
run_tissue_models <- function(final_df, tissue) {
    models_list <- list()

    # Statistics to analyze
    statistics <- c("Log2FC_MetBurden", "Log2FC_MetSeeding", "Log2FC_MetDormancy", "Log2FC_SizePercentile", "Log2FC_SizePercentile_50", "Log2FC_SizePercentile_20", "Log2FC_PeakMode")

    for (statistic in statistics) {
        # Run models for current statistic
        models_list[[paste0(statistic, "_Sex")]] <- run_glm_models(
            final_df,
            statistic = statistic,
            tissue = tissue,
            formula_type = "Sex"
        )

        models_list[[paste0(statistic, "_Genotype")]] <- run_glm_models(
            final_df,
            statistic = statistic,
            tissue = tissue,
            formula_type = "Mouse_Genotype"
        )

        # Optional: Run 1-week models
        models_list[[paste0(statistic, "_1week_Sex")]] <- run_glm_models(
            final_df,
            statistic = statistic,
            tissue = tissue,
            formula_type = "Sex",
            time_point = "1week"
        )

        models_list[[paste0(statistic, "_1week_Genotype")]] <- run_glm_models(
            final_df,
            statistic = statistic,
            tissue = tissue,
            formula_type = "Mouse_Genotype",
            time_point = "1week"
        )
    }

    return(models_list)
}

tissues <- unique(final_df$Tissue)
all_tissue_models <- lapply(tissues, function(tissue) {
    run_tissue_models(final_df, tissue)
})
names(all_tissue_models) <- tissues

# Create combined heatmap data for all tissues
all_tissue_heatmaps <- lapply(names(all_tissue_models), function(tissue) {
    print(tissue)
    models <- all_tissue_models[[tissue]]
    
    # Create heatmap data similar to existing code
    rbind(
        # MetBurden
        models[["Log2FC_MetBurden_Sex"]] %>%
            select(term, estimate, Mouse_Genotype, Cell_Line, gene, p.value) %>%
            mutate(
                Tissue = tissue,
                Statistic = "Log2FC_MetBurden",
                term = paste("Female vs Male"),
                Info = paste(Mouse_Genotype, Cell_Line, term, sep = " ")
            ) %>%
            select(-c(term, Mouse_Genotype, Cell_Line)) %>%
            pivot_wider(
                names_from = Info,
                values_from = c(estimate, p.value),
                values_fill = list(estimate = 0, p.value = NA)
            ) %>%
            left_join(
                models[["Log2FC_MetBurden_Genotype"]] %>%
                    select(term, estimate, Sex, Cell_Line, gene, p.value) %>%
                    mutate(
                        Statistic = "Log2FC_MetBurden",
                        term = paste("BL6 vs NSG"),
                        Sex = ifelse(Sex == "F", "Female", "Male"),
                        Info = paste(Sex, Cell_Line, term, sep = " ")
                    ) %>%
                    select(-c(term, Sex, Cell_Line)) %>%
                    pivot_wider(
                        names_from = Info,
                        values_from = c(estimate, p.value),
                        values_fill = list(estimate = 0, p.value = NA)
                    ),
                by = c("gene", "Statistic")
            ),
        
        # MetSeeding
        models[["Log2FC_MetSeeding_Sex"]] %>%
            select(term, estimate, Mouse_Genotype, Cell_Line, gene, p.value) %>%
            mutate(
                Tissue = tissue,
                Statistic = "Log2FC_MetSeeding",
                term = paste("Female vs Male"),
                Info = paste(Mouse_Genotype, Cell_Line, term, sep = " ")
            ) %>%
            select(-c(term, Mouse_Genotype, Cell_Line)) %>%
            pivot_wider(
                names_from = Info,
                values_from = c(estimate, p.value),
                values_fill = list(estimate = 0, p.value = NA)
            ) %>%
            left_join(
                models[["Log2FC_MetSeeding_Genotype"]] %>%
                    select(term, estimate, Sex, Cell_Line, gene, p.value) %>%
                    mutate(
                        Statistic = "Log2FC_MetSeeding",
                        term = paste("BL6 vs NSG"),
                        Sex = ifelse(Sex == "F", "Female", "Male"),
                        Info = paste(Sex, Cell_Line, term, sep = " ")
                    ) %>%
                    select(-c(term, Sex, Cell_Line)) %>%
                    pivot_wider(
                        names_from = Info,
                        values_from = c(estimate, p.value),
                        values_fill = list(estimate = 0, p.value = NA)
                    ),
                by = c("gene", "Statistic")
            ),

        # MetDormancy
        models[["Log2FC_MetDormancy_Sex"]] %>%
            select(term, estimate, Mouse_Genotype, Cell_Line, gene, p.value) %>%
            mutate(
                Tissue = tissue,
                Statistic = "Log2FC_MetDormancy",
                term = paste("Female vs Male"),
                Info = paste(Mouse_Genotype, Cell_Line, term, sep = " ")
            ) %>%
            select(-c(term, Mouse_Genotype, Cell_Line)) %>%
            pivot_wider(
                names_from = Info,
                values_from = c(estimate, p.value),
                values_fill = list(estimate = 0, p.value = NA)
            ) %>%
            left_join(
                models[["Log2FC_MetDormancy_Genotype"]] %>%
                    select(term, estimate, Sex, Cell_Line, gene, p.value) %>%
                    mutate(
                        Statistic = "Log2FC_MetDormancy",
                        term = paste("BL6 vs NSG"),
                        Sex = ifelse(Sex == "F", "Female", "Male"),
                        Info = paste(Sex, Cell_Line, term, sep = " ")
                    ) %>%
                    select(-c(term, Sex, Cell_Line)) %>%
                    pivot_wider(
                        names_from = Info,
                        values_from = c(estimate, p.value),
                        values_fill = list(estimate = 0, p.value = NA)
                    ),
                by = c("gene", "Statistic")
            ),

        # SizePercentile
        models[["Log2FC_SizePercentile_Sex"]] %>%
            select(term, estimate, Mouse_Genotype, Cell_Line, gene, p.value) %>%
            mutate(
                Tissue = tissue,
                Statistic = "Log2FC_SizePercentile",
                term = paste("Female vs Male"),
                Info = paste(Mouse_Genotype, Cell_Line, term, sep = " ")
            ) %>%
            select(-c(term, Mouse_Genotype, Cell_Line)) %>%
            pivot_wider(
                names_from = Info,
                values_from = c(estimate, p.value),
                values_fill = list(estimate = 0, p.value = NA)
            ) %>%
            left_join(
                models[["Log2FC_SizePercentile_Genotype"]] %>%
                    select(term, estimate, Sex, Cell_Line, gene, p.value) %>%
                    mutate(
                        Statistic = "Log2FC_SizePercentile",
                        term = paste("BL6 vs NSG"),
                        Sex = ifelse(Sex == "F", "Female", "Male"),
                        Info = paste(Sex, Cell_Line, term, sep = " ")
                    ) %>%
                    select(-c(term, Sex, Cell_Line)) %>%
                    pivot_wider(
                        names_from = Info,
                        values_from = c(estimate, p.value),
                        values_fill = list(estimate = 0, p.value = NA)
                    ),
                by = c("gene", "Statistic")
            ),
        
        # SizePercentile
        models[["Log2FC_SizePercentile_50_Sex"]] %>%
            select(term, estimate, Mouse_Genotype, Cell_Line, gene, p.value) %>%
            mutate(
                Tissue = tissue,
                Statistic = "Log2FC_SizePercentile_50",
                term = paste("Female vs Male"),
                Info = paste(Mouse_Genotype, Cell_Line, term, sep = " ")
            ) %>%
            select(-c(term, Mouse_Genotype, Cell_Line)) %>%
            pivot_wider(
                names_from = Info,
                values_from = c(estimate, p.value),
                values_fill = list(estimate = 0, p.value = NA)
            ) %>%
            left_join(
                models[["Log2FC_SizePercentile_50_Genotype"]] %>%
                    select(term, estimate, Sex, Cell_Line, gene, p.value) %>%
                    mutate(
                        Statistic = "Log2FC_SizePercentile_50",
                        term = paste("BL6 vs NSG"),
                        Sex = ifelse(Sex == "F", "Female", "Male"),
                        Info = paste(Sex, Cell_Line, term, sep = " ")
                    ) %>%
                    select(-c(term, Sex, Cell_Line)) %>%
                    pivot_wider(
                        names_from = Info,
                        values_from = c(estimate, p.value),
                        values_fill = list(estimate = 0, p.value = NA)
                    ),
                by = c("gene", "Statistic")
            ),

        # SizePercentile
        models[["Log2FC_SizePercentile_20_Sex"]] %>%
            select(term, estimate, Mouse_Genotype, Cell_Line, gene, p.value) %>%
            mutate(
                Tissue = tissue,
                Statistic = "Log2FC_SizePercentile_20",
                term = paste("Female vs Male"),
                Info = paste(Mouse_Genotype, Cell_Line, term, sep = " ")
            ) %>%
            select(-c(term, Mouse_Genotype, Cell_Line)) %>%
            pivot_wider(
                names_from = Info,
                values_from = c(estimate, p.value),
                values_fill = list(estimate = 0, p.value = NA)
            ) %>%
            left_join(
                models[["Log2FC_SizePercentile_20_Genotype"]] %>%
                    select(term, estimate, Sex, Cell_Line, gene, p.value) %>%
                    mutate(
                        Statistic = "Log2FC_SizePercentile_20",
                        term = paste("BL6 vs NSG"),
                        Sex = ifelse(Sex == "F", "Female", "Male"),
                        Info = paste(Sex, Cell_Line, term, sep = " ")
                    ) %>%
                    select(-c(term, Sex, Cell_Line)) %>%
                    pivot_wider(
                        names_from = Info,
                        values_from = c(estimate, p.value),
                        values_fill = list(estimate = 0, p.value = NA)
                    ),
                by = c("gene", "Statistic")
            )

        # PeakMode
        # models[["Log2FC_PeakMode_Sex"]] %>%
        #     select(term, estimate, Mouse_Genotype, Cell_Line, gene, p.value) %>%
        #     mutate(
        #         Tissue = tissue,
        #         Statistic = "Log2FC_PeakMode",
        #         term = paste("Female vs Male"),
        #         Info = paste(Mouse_Genotype, Cell_Line, term, sep = " ")
        #     ) %>%
        #     select(-c(term, Mouse_Genotype, Cell_Line)) %>%
        #     pivot_wider(
        #         names_from = Info,
        #         values_from = c(estimate, p.value),
        #         values_fill = 0
        #     ) %>%
        #     left_join(
        #         models[["Log2FC_PeakMode_Genotype"]] %>%
        #             select(term, estimate, Sex, Cell_Line, gene, p.value) %>%
        #             mutate(
        #                 Statistic = "Log2FC_PeakMode",
        #                 term = paste("BL6 vs NSG"),
        #                 Sex = ifelse(Sex == "F", "Female", "Male"),
        #                 Info = paste(Sex, Cell_Line, term, sep = " ")
        #             ) %>%
        #             select(-c(term, Sex, Cell_Line)) %>%
        #             pivot_wider(
        #                 names_from = Info,
        #                 values_from = c(estimate, p.value),
        #                 values_fill = 0
        #             ),
        #         by = c("gene", "Statistic")
        #     )
    )
})
# Combine all tissues
combined_heatmap_data <- bind_rows(all_tissue_heatmaps)

genes_with_enough_data <- df %>%
    mutate(
        target = gsub("mmu-sg[1-3]", "", sgID)
    ) %>%
    group_by(target, Mouse_ID) %>%
    summarize(
        n = n()
    ) %>% 
    ungroup() %>%
    group_by(target) %>%
    summarize(
        n = mean(n)
    ) %>% 
    filter(n > 100) %>%
    pull(target)

#fwrite(combined_heatmap_data, "combined_heatmap_data.csv")
combined_heatmap_data <- fread("combined_heatmap_data.csv")

statistic_heatmap_to_plot <- combined_heatmap_data %>%
    filter(Tissue == "Liver") %>%
    select(-Tissue) %>%
    #filter(Statistic != "Log2FC_PeakMode" & Statistic != "Log2FC_SizePercentile") %>%
    pivot_longer(
        cols = -c(gene, Statistic),
        names_to = c(".value", "comparison"),
        names_pattern = "(estimate|p.value)_(.*)",
        values_to = c("estimate", "p.value")
    ) %>%
    mutate(
        target = gsub("mmu-sg[1-3]", "", gene)
    ) %>%
    filter(target %in% genes_with_enough_data) %>%
    mutate(
        geneGroup = factor(case_when(
            gene == "sgSafe" ~ "Control",
            grepl("Nfib", gene) ~ "Nfib",
            gene %in% unique(gsub("mmu-sg[1-3]", "", sgOnco)) ~ "Oncogene",
            gene %in% unique(gsub("mmu-sg[1-3]", "", sgTSG)) ~ "TSG",
            gene %in% unique(gsub("mmu-sg[1-3]", "", sgBoth)) ~ "Both",
            gene %in% unique(gsub("mmu-sg[1-3]", "", sgOther)) ~ "Other",
            #TRUE ~ NA
        ), levels = c("Control", "Nfib", "TSG", "Oncogene", "Both", "Other")),
        # geneGroup = factor(case_when(
        #     gene %in% sgSafe ~ "Control",
        #     gene %in% sgNfib ~ "Nfib",
        #     gene %in% sgOnco ~ "Oncogene",
        #     gene %in% sgTSG ~ "TSG",
        #     gene %in% sgBoth ~ "Both",
        #     gene %in% sgOther ~ "Other",
        #     #TRUE ~ NA
        # ), levels = c("Control", "Nfib", "TSG", "Oncogene", "Both", "Other")),
        celllineGroup = factor(case_when(
            grepl("rp116", comparison) ~ "rp116",
            grepl("rp48", comparison) ~ "rp48",
            TRUE ~ NA
        ), levels = c("rp116", "rp48")),
        statisticGroup = factor(case_when(
            grepl("MetBurden", Statistic) ~ "MetBurden",
            grepl("MetSeeding", Statistic) ~ "MetSeeding",
            grepl("MetDormancy", Statistic) ~ "MetDormancy",
            grepl("SizePercentile_20", Statistic) ~ "SizePercentile_20",
            grepl("SizePercentile_50", Statistic) ~ "SizePercentile_50",
            grepl("SizePercentile", Statistic) ~ "SizePercentile",
            grepl("PeakMode", Statistic) ~ "PeakMode",
            TRUE ~ NA
        ), levels = c("MetBurden", "MetSeeding", "MetDormancy", "SizePercentile", "SizePercentile_50", "SizePercentile_20", "PeakMode")),
        comparison = case_when(
            grepl("rp116", comparison) ~ gsub("rp116 ", "", comparison),
            grepl("rp48", comparison) ~ gsub("rp48 ", "", comparison)
        ),
        pFDR = p.adjust(p.value, method = "BH")
    )
statistic_heatmap_plot <- statistic_heatmap_to_plot %>%
    ggplot(
        aes(
            x = comparison,
            y = gene,
            fill = estimate
        )
    ) +
    geom_tile(color = "black") +
    geom_point(
        data = subset(statistic_heatmap_to_plot, pFDR <= 0.05),
        aes(shape = "pFDR <= 0.05"),
        color = "black",
        size = 1
    ) +
    scale_fill_gradient2(
        low = "blue",
        mid = "white",
        high = "red",
        midpoint = 0
    ) +
    scale_y_discrete(expand = c(0, 0)) +
    scale_x_discrete(expand = c(0, 0)) +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 8),
        legend.text = element_text(face = "bold"),
        plot.background = element_blank(),
        panel.border = element_blank(),
        strip.text.y.left = element_text(angle = 0)
    ) +
    labs(
        x = "Comparison",
        y = "Gene",
        fill = "Coeff."
    ) +
    labs(title = "", x = "", y = "", z = "") +
    facet_grid(cols = vars(celllineGroup, statisticGroup), rows = vars(geneGroup), space = "free", scales = "free", switch = "both") +
    guides(
        shape = guide_legend(order = 2)
    ) +
    scale_shape_manual(name = "", values = c("pFDR <= 0.05" = 8))
pdf("statistic_heatmap_plot_liver.pdf", width = 18, height = 10)
print(statistic_heatmap_plot)
dev.off()