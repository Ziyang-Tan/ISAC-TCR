library(dplyr)
library(readr)
library(ggplot2)
library(ggrepel)

file_dir <- "data/DE_resample_results"
file_names <- list.files(file_dir, pattern = "*.csv")
file_list <- lapply(file_names, function(file) {
    read_csv(file.path(file_dir, file), show_col_types = FALSE) %>% mutate(file_name = file)
})

# gene_list <- lapply(file_list, function(file) file %>% select(gene, file_name, p_val_adj))
gene_list <- lapply(file_list, function(file) file %>% select(gene, file_name, rank_metric))
common_genes <- lapply(gene_list, function(x) x$gene) %>% Reduce(f = intersect)

common_gene_list <- lapply(gene_list, function(x) {
    x %>%
        filter(
            gene %in% common_genes
        ) %>%
        mutate(
            rank = as.numeric(rownames(.)),
            # significant = if_else(p_val_adj < 1e-2, TRUE, FALSE)
            show = if_else(abs(rank_metric) > 0.3, TRUE, FALSE)
        )
})

theme_slopegraph <- list(
    # move the x axis labels up top
    scale_x_discrete(position = "top", expand = c(1.3, 0.3)),
    theme_bw(),
    # Format tweaks
    # Remove the legend
    # theme(legend.position = "none"),
    guides(alpha = "none"),
    # Remove the panel background
    theme(panel.background = element_blank()),
    # Remove the panel border
    theme(panel.border = element_blank()),
    # Remove just about everything from the y axis
    theme(axis.title.y = element_blank()),
    theme(axis.text.y = element_blank()),
    theme(panel.grid.major.y = element_blank()),
    theme(panel.grid.minor.y = element_blank()),
    # Remove a few things from the x axis and increase font size
    theme(axis.title.x = element_blank()),
    theme(panel.grid.major.x = element_blank()),
    theme(axis.text.x.top = element_text(size = 12)),
    # Remove x & y tick marks
    theme(axis.ticks = element_blank()),
    # Format title & subtitle
    theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5)),
    theme(plot.subtitle = element_text(hjust = 0.5))
)

# test 
name1 <- common_gene_list[[2]]$file_name %>% unique()
name2 <- common_gene_list[[6]]$file_name %>% unique()
pair_name <- paste(
    name1,
    name2,
    sep = " vs "
)
df <- rbind(common_gene_list[[2]], common_gene_list[[6]]) %>%
    group_by(gene) %>%
    mutate(
        alpha_show = any(show),
        pos = sum(rank_metric > 0),
        neg = sum(rank_metric < 0),
        color_group = case_when(
            pos == 2 & all(show) ~ "both higher in expanded",
            neg == 2 & all(show) ~ "both higher in non-expanded",
            pos == 1 & neg == 1  & all(show) ~ "different in adults and isac",
            TRUE ~ "Others"
        )
    ) %>%
    ungroup()
ggplot(df, aes(x = file_name, y = -rank, group = gene)) + # reverse rank so positive logFC is on top
    geom_line(aes(alpha = alpha_show, color = color_group), linewidth = 1) +
    geom_text_repel(
        data = . %>% dplyr::filter(file_name == name1 & show),
        aes(label = gene, color = color_group),
        hjust = "right",
        fontface = "bold",
        size = 3,
        nudge_x = -.1,
        direction = "y",
        force = 0.3
    ) +
    geom_text_repel(
        data = . %>% dplyr::filter(file_name == name2 & show),
        aes(label = gene, color = color_group),
        hjust = "left",
        fontface = "bold",
        size = 3,
        nudge_x = .1,
        direction = "y",
        force = 0.3
    ) +
    scale_colour_manual(values = c("#ff8800", "#365690", "red", "grey")) +
    scale_alpha_discrete(range = c(0, 1)) +
    theme_slopegraph +
    ggtitle(pair_name)


for (i in c(1:3)) {
    for (j in c(4:8)) {
        name1 <- common_gene_list[[i]]$file_name %>% unique()
        name2 <- common_gene_list[[j]]$file_name %>% unique()
        pair_name <- paste(
            name1,
            name2,
            sep = " vs "
        )
        df <- rbind(common_gene_list[[i]], common_gene_list[[j]]) %>%
            group_by(gene) %>%
            mutate(
                alpha_show = any(show),
                pos = sum(rank_metric > 0),
                neg = sum(rank_metric < 0),
                color_group = case_when(
                    pos == 2 & all(show) ~ "both higher in expanded",
                    neg == 2 & all(show) ~ "both higher in non-expanded",
                    pos == 1 & neg == 1  & all(show) ~ "different in adults and isac",
                    TRUE ~ "Others"
                )
            ) %>%
            ungroup()
        ggplot(df, aes(x = file_name, y = -rank, group = gene)) + # reverse rank so positive logFC is on top
            geom_line(aes(alpha = alpha_show, color = color_group), linewidth = 1) +
            geom_text_repel(
                data = . %>% dplyr::filter(file_name == name1 & show),
                aes(label = gene, color = color_group),
                hjust = "right",
                fontface = "bold",
                size = 3,
                nudge_x = -.1,
                direction = "y",
                force = 0.3
            ) +
            geom_text_repel(
                data = . %>% dplyr::filter(file_name == name2 & show),
                aes(label = gene, color = color_group),
                hjust = "left",
                fontface = "bold",
                size = 3,
                nudge_x = .1,
                direction = "y",
                force = 0.3
            ) +
            scale_colour_manual(values = c(
                "both higher in expanded" = "#ff8800",
                "both higher in non-expanded" = "#365690",
                "different in adults and isac" = "#b85dbd",
                "Others" = "grey"
            )) +
            scale_alpha_discrete(range = c(0, 1)) +
            theme_slopegraph +
            ggtitle(pair_name)
        ggsave(paste0("figures/DE/slopegraph_resample/", pair_name, ".pdf"), width = 12, height = 10)
    }
}
