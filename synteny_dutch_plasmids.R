# Libraries ----

library(tidyverse)
library(gggenomes)
library(ggtext)
library(ggrepel)

# Plotting function ----

synteny_plot <- function(gbff, links, mobtyper, order) {
  theme <- theme_gggenomes_clean(base_size = 8) +
    theme()

  gbff <- dplyr::bind_rows(gbff, mobtyper) |>
    dplyr::mutate(
      gene_type = dplyr::case_when(
        !is.na(gene_type) ~ gene_type,
        str_detect(dbxref, "GO:0046677") ~ "response to antibiotic",
        str_detect(str_to_lower(product), "efflux") ~ "efflux",
        str_detect(str_to_lower(product), "transpos") ~ "transposable element",
        str_detect(product, "IS") ~ "transposable element",
        .default = "Other"
      ),
      seq_id = factor(seq_id, levels = order, ordered = TRUE)
    )

  plot <- gggenomes::gggenomes(gbff) |>
    gggenomes::add_links(links) +
    geom_link_curved(
      fill = "#C8C8C8",
      alpha = .5,
      linetype = 0
    ) +
    geom_seq(color = "#C8C8C8") +
    geom_gene(
      data = genes(.gene_types = "CDS"),
      aes(fill = factor(
        gene_type,
        levels = c(
          "response to antibiotic",
          "efflux",
          "transposable element",
          "replicon",
          "mate-pair-formation",
          "relaxase"
        ),
        ordered = TRUE
      )),
      position = "strand",
      linetype = 0,
      size = 3,
      show.legend = TRUE
    ) +
    scale_fill_manual(
      name = "Putative gene function",
      values = c(
        "response to antibiotic" = "#ca005d",
        "efflux" = "#ffb612",
        "transposable element" = "#6abda4",
        "replicon" = "#007bc7",
        "mate-pair-formation" = "#552c6f",
        "relaxase" = "#e17000"
      ),
      labels = function(breaks) {
        breaks[is.na(breaks)] <- "other"
        breaks
      },
      na.value = "#C8C8C8",
      drop = FALSE,
      guide = guide_legend(ncol = 1)
    ) +
    geom_text_repel(
      data = genes(
        .gene_types = "CDS",
        gene_type %in% c("response to antibiotic", "efflux")
      ),
      aes(
        label = name,
        x = x + ((xend - x) / 2),
        y = ifelse(strand == "+", y + 0.05, y - 0.05)
      ),
      max.overlaps = Inf,
      min.segment.length = 0,
      size = 7 / .pt,
      nudge_x = .5,
      nudge_y = 0.05,
      segment.curvature = -0.1,
      segment.ncp = 3,
      segment.angle = 20,
      segment.alpha = .5
    ) +
    geom_seq_label(size = 8 / .pt, nudge_y = -0.25) +
    scale_x_bp(suffix = "bp") +
    scale_y_discrete(expand = c(.05, .05)) +
    theme

  return(plot)
}

# Generate plot ----

gbff <- list.files(
  "synteny/cfr_rivm_mrsa/rotated",
  "*.gbff",
  full.names = TRUE
) |>
  map_dfr(read_gbk)

links <- read.csv("synteny/cfr_rivm_mrsa/blast_homology.csv")
mobtyper <- read.csv("synteny/cfr_rivm_mrsa/mob_markers.csv")
order <- readLines("synteny/cfr_rivm_mrsa/synteny_order.txt")

plot <- synteny_plot(
  gbff, links, mobtyper, order
)
ggsave(
  "synteny/rivm_mrsa_cfr_plasmids.pdf",
  plot = plot, width = 180, unit = "mm"
)