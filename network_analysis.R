# libraries ----

library(tidyverse)
library(tidygraph)
library(igraph)
library(ggnetwork)
library(ggforce)
library(ggtext)
library(patchwork)
library(ggplot2)
library(ggrepel)
library(mclust)

# read metadata ----

rivm_metadata <- read.delim(
  "network_analysis/rivm_network_metadata.csv",
  sep = ";"
)
public_metadata <- read.delim(
  "network_analysis/public_network_metadata.csv",
  sep = ";"
)

rivm_metadata <- rivm_metadata |>
  mutate(source = "RIVM")

public_metadata <- public_metadata |>
  mutate(source = "Public")

metadata <- rivm_metadata |>
  bind_rows(public_metadata) |>
  mutate(
    across(everything(), ~ na_if(.x, "")),
    Genus = ifelse(is.na(order), "Other|Unknown", genus),
    Country = ifelse(is.na(region), "Other|Unknown", country)
  ) |>
  group_by(order) |>
  mutate(n_order = n()) |>
  group_by(Genus) |>
  mutate(n_genus = n()) |>
  ungroup() |>
  arrange(desc(n_order), desc(n_genus), Genus) |>
  mutate(
    Genus = fct_inorder(Genus, ordered = TRUE),
    Genus = fct_relevel(Genus, "Other|Unknown", after = Inf)
  ) |>
  group_by(region) |>
  mutate(n_region = n()) |>
  group_by(Country) |>
  mutate(n_country = n()) |>
  ungroup() |>
  arrange(desc(n_region), desc(n_country), Country) |>
  mutate(
    Country = fct_inorder(Country, ordered = TRUE),
    Country = fct_relevel(Country, "Other|Unknown", after = Inf)
  )

# parse pling results ----

pling_subcom <- read.delim(
  "pling/dcj_thresh_4_graph/objects/typing.tsv",
  sep = "\t"
)

pling_containment <- read.delim(
  "pling/containment/all_pairs_containment_distance.tsv",
  sep = "\t"
)

pling_dcj <- read.delim(
  "pling/all_plasmids_distances.tsv",
  sep = "\t"
)

pling_subcom <- pling_subcom |>
  pull(type, plasmid)

network <- pling_containment |>
  left_join(
    pling_dcj,
    by = join_by(plasmid_1, plasmid_2),
    suffix = c("_containment", "_dcj")
  ) |>
  filter(
    distance_containment <= 0.5,
    distance_dcj <= 4,
    pling_subcom[plasmid_1] == pling_subcom[plasmid_2]
  )

pling_graph <- tbl_graph(
  nodes = metadata,
  edges = network,
  node_key = "accession",
  directed = FALSE
) |>
  activate(nodes) |>
  mutate(
    group = factor(group_components()),
    group = fct_lump_min(group, min = 2, other_level = "Singleton")
  )

# plot pling networks ----
plot_network <- function(graph_layout, category, colors) {
  theme <- ggnetwork::theme_blank(base_size = 7) +
    theme(
      legend.margin = margin(0),
      legend.text = element_text(size = 7),
      legend.title = element_text(size = 8),
      legend.key.size = unit(8, "pt")
    )
  group_label <- graph_layout |>
    filter(group != "Singleton") |>
    dplyr::group_by(group) |>
    dplyr::summarise(
      x = mean(x),
      y = mean(y)
    )
  plot <- ggplot(
    graph_layout,
    aes(
      x = x,
      y = y,
      xend = xend,
      yend = yend
    )
  ) +
    ggforce::geom_mark_hull(
      aes(
        group = group,
        filter = group != "Singleton",
        label = NA,
        description = NA
      ),
      show.legend = FALSE,
      color = "black",
      concavity = 10,
      expand = unit(2, "mm"),
      radius = unit(2, "mm")
    ) +
    ggnetwork::geom_edges(
      color = "#C6C6C6",
      curvature = 0.125
    ) +
    ggnetwork::geom_nodes(
      aes(
        shape = source,
        color = .data[[category]]
      ),
      size = 2
    ) +
    scale_color_manual(
      values = colors
    ) +
    scale_shape_manual(
      values = c(
        "Public" = 19,
        "RIVM" = 18
      )
    ) +
    geom_label(
      data = group_label,
      aes(
        x = x,
        y = y,
        label = group
      ),
      size = 5 / .pt,
      fontface = "bold",
      alpha = .5,
      inherit.aes = FALSE
    ) +
    theme
  return(plot)
}

country_colors <- c(
  "China" = "#005a32",
  "South Korea" = "#238b45",
  "India" = "#41ab5d",
  "Japan" = "#74c476",
  "Philippines" = "#a1d99b",
  "Taiwan" = "#c7e9c0",
  "Viet Nam" = "#edf8e9",
  "The Netherlands" = "#ca005d",
  "Italy" = "#8c2d04",
  "Germany" = "#d94801",
  "Ireland" = "#f16913",
  "Denmark" = "#fd8d3c",
  "France" = "#fdae6b",
  "Greece" = "#fdd0a2",
  "Poland" = "#feedde",
  "USA" = "#3182bd",
  "Canada" = "#9ecae1",
  "Brazil" = "#756bb1",
  "Chile" = "#bcbddc",
  "Australia" = "#ffb612"
)

genus_colors <- c(
  "Staphylococcus" = "#a50f15",
  "Mammaliicoccus" = "#de2d26",
  "Bacillus" = "#fb6a4a",
  "Paenibacillus" = "#fee5d9",
  "Lysinibacillus" = "#fcae91",
  "Escherichia" = "#005a32",
  "Proteus" = "#238b45",
  "Klebsiella" = "#41ab5d",
  "Citrobacter" = "#74c476",
  "Leclercia" = "#a1d99b",
  "Providencia" = "#c7e9c0",
  "Salmonella" = "#edf8e9",
  "Enterococcus" = "#756bb1",
  "Vagococcus" = "#bcbddc",
  "Streptococcus" = "#efedf5",
  "Pseudoalteromonas" = "#ffb612",
  "Campylobacter" = "#3182bd",
  "Other|Unknown" = "#C6C6C6"
)

set.seed(3721)
pling_layout <- fortify(
  pling_graph,
  layout = layout_components(pling_graph, layout_with_kk)
)

pling_genus <- plot_network(
  graph_layout = pling_layout,
  category = "Genus",
  colors = genus_colors
) +
  labs(title = "B. Pling subcommunities, by genus")
ggsave(
  "network_analysis/pling_genus.pdf",
  plot = pling_genus, width = 180, unit = "mm"
)

pling_country <- plot_network(
  graph_layout = pling_layout,
  category = "Country",
  colors = country_colors
) +
  labs(title = "D. Pling subcommunities, by country")
ggsave(
  "network_analysis/pling_country.pdf",
  plot = pling_country, width = 180, unit = "mm"
)

# parse mash results ----

mash <- read.delim(
  "mash/mash_edges.tsv",
  sep = "\t",
  header = FALSE,
  col.names = c("seq1", "seq2", "distance", "p_val", "shared_hashes")
) |>
  mutate(
    seq1 = str_extract(seq1, "([^/]*).fasta$", group = 1),
    seq2 = str_extract(seq2, "([^/]*).fasta$", group = 1)
  )

mash_nearest_neighbour <- mash |>
  group_by(seq1) |>
  arrange(distance) |>
  slice_head(n = 1) |>
  select(seq1, seq2, distance) |>
  rename(
    accession = seq1,
    nearest_neighbour = seq2,
    mash_distance = distance
  )

mash <- mash |>
  filter(distance <= 0.05)

mash_graph <- tbl_graph(
  nodes = metadata,
  edges = mash,
  node_key = "accession",
  directed = FALSE
) |>
  activate(nodes) |>
  mutate(
    group = factor(group_components()),
    group = fct_lump_min(group, min = 2, other_level = "Singleton")
  )

set.seed(3721)
mash_layout <- fortify(
  mash_graph,
  layout = layout_components(mash_graph, layout_with_kk)
)

mash_genus <- plot_network(
  graph_layout = mash_layout,
  category = "Genus",
  colors = genus_colors
) +
  labs(title = "A. Mash clusters, by genus")
ggsave(
  "network_analysis/mash_genus.pdf",
  plot = mash_genus, width = 180, unit = "mm"
)

mash_country <- plot_network(
  graph_layout = mash_layout,
  category = "Country",
  colors = country_colors
) +
  labs(title = "C. Mash clusters, by country")
ggsave(
  "network_analysis/mash_country.pdf",
  plot = mash_country, width = 180, unit = "mm"
)

design <- "
AAAEBBB
AAAEBBB
CCCEDDD
CCCEDDD
"

multi_panel <- mash_genus +
  pling_genus +
  mash_country +
  pling_country +
  guide_area() +
  plot_layout(design = design, guides = "collect")
ggsave(
  "network_analysis/network_multipanel.pdf",
  plot = multi_panel, width = 180, unit = "mm"
)

node_table <- pling_graph |>
  activate(nodes) |>
  as_tibble() |>
  select(
    accession, order, genus, Genus, region, country, Country, source, group
  ) |>
  left_join(
    mash_graph |>
      activate(nodes) |>
      as_tibble() |>
      select(accession, group),
    by = "accession",
    suffix = c("_pling", "_mash")
  ) |>
  left_join(
    mash_nearest_neighbour,
    by = "accession"
  )

write.csv(
  node_table,
  file = "network_analysis/network_node_table.csv",
  row.names = FALSE
  )

bind_cols(
  node_table |>
    summarise(singleton_incl_ari = adjustedRandIndex(group_pling, group_mash)),
  node_table |>
    filter(group_pling != "Singleton" & group_mash != "Singleton") |>
    summarise(singleton_rm_ari = adjustedRandIndex(group_pling, group_mash)),
  node_table |>
    mutate(
      group_pling = as.character(group_pling),
      group_pling = if_else(
        group_pling == "Singleton",
        str_glue("Singleton_{row_number()}"),
        group_pling
      ),
      group_mash = as.character(group_mash),
      group_mash = if_else(
        group_mash == "Singleton",
        str_glue("Singleton_{row_number() * 10}"),
        group_mash
      )
    ) |>
    summarise(singleton_idx_ari = adjustedRandIndex(group_pling, group_mash))
) |>
  write.csv(
    file = "network_analysis/ARI.csv",
    row.names = FALSE
)
