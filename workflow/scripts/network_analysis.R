# libraries ----

library(tidyverse)
library(countrycode)
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

rivm_metadata <- read.csv(snakemake@input[["rivm_metadata"]]) |>
  mutate(sequence_source = "RIVM")

public_metadata <- read.csv(snakemake@input[["public_metadata"]]) |>
  mutate(sequence_source = "Public")

metadata <- rivm_metadata |>
  bind_rows(public_metadata) |>
  mutate(
    country = replace_values(
      country, c("not applicable", "", "Arctic Ocean") ~ NA
    ),
    continent = fct_infreq(
      countrycode(country, "country.name", "continent"),
      ordered = TRUE
    )
  ) |>
  group_by(continent) |>
  mutate(continent_count = n()) |>
  group_by(country) |>
  mutate(country_count = n()) |>
  group_by(family) |>
  mutate(family_count = n()) |>
  group_by(genus) |>
  mutate(genus_count = n()) |>
  ungroup() |>
  arrange(desc(continent_count), desc(country_count)) |>
  mutate(
    continent = fct_inorder(continent, ordered = TRUE),
    country = fct_inorder(country, ordered = TRUE)
  ) |>
  arrange(desc(family_count), desc(genus_count)) |>
  mutate(
    family = fct_inorder(family, ordered = TRUE),
    family = fct_lump_n(family, 5),
    genus = ifelse(family == "Other", NA, genus),
    genus = fct_inorder(genus, ordered = TRUE)
  )

# parse pling results ----

pling_subcommunity <- read.delim(
  file.path(
    snakemake@input[["pling"]], "dcj_thresh_4_graph/objects/typing.tsv"
  ),
  sep = "\t"
) |>
  pull(type, plasmid)

pling_containment <- read.delim(
  file.path(
    snakemake@input[["pling"]], "containment/all_pairs_containment_distance.tsv"
  ),
  sep = "\t"
)

pling_dcj <- read.delim(
  file.path(
    snakemake@input[["pling"]], "all_plasmids_distances.tsv"
  ),
  sep = "\t"
)

network <- pling_containment |>
  left_join(
    pling_dcj,
    by = join_by(plasmid_1, plasmid_2),
    suffix = c("_containment", "_dcj")
  ) |>
  filter(
    distance_containment <= snakemake@params[["pling_containment_cutoff"]],
    distance_dcj <= snakemake@params[["pling_dcj_cutoff"]],
    pling_subcommunity[plasmid_1] == pling_subcommunity[plasmid_2]
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
        shape = sequence_source,
        fill = .data[[category]]
      ),
      color = "#C6C6C6",
      size = 2,
      linewidth = .125
    ) +
    scale_fill_manual(
      values = colors,
      labels = function(breaks) {
        breaks[is.na(breaks)] <- "other|unknown"
        breaks
      },
      na.value = "#C6C6C6",
      guide = guide_legend(ncol = 1, override.aes = list(shape = 21))
    ) +
    scale_shape_manual(
      values = c(
        "Public" = 21,
        "RIVM" = 23
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

rivm_colors <- c(
  "#007bc7",
  "#ca005d",
  "#ffb612",
  "#552c6f",
  "#6abda4",
  "#e17000",
  "#777b00",
  "#94710a",
  "#81b6d0",
  "#673327",
  "#01689B",
  "#275837",
  "#d883b9",
  "#a90061"
)

generate_palette <- function(df, category, group_levels, level) {
  group_number <- which(group_levels == level)
  base_color <- rivm_colors[group_number]
  levels <- levels(forcats::fct_drop(df[[category]]))
  n <- length(levels)
  hues <- colorRampPalette(
    c(base_color, "#FFFFFF"),
    interpolate = "linear"
  )(n + 1)
  hues <- hues[1:n]
  names(hues) <- levels
  hues
}

genus_colors <- metadata |>
  group_by(family) |>
  group_map(
    ~ generate_palette(.x, "genus", levels(metadata$family), .y$family)
  ) |>
  flatten()

country_colors <- metadata |>
  group_by(continent, .drop = FALSE) |>
  group_map(
    ~ generate_palette(.x, "country", levels(metadata$continent), .y$continent)
  ) |>
  flatten()

set.seed(3721)
pling_layout <- fortify(
  pling_graph,
  layout = layout_components(pling_graph, layout_with_kk)
)

pling_genus <- plot_network(
  graph_layout = pling_layout,
  category = "genus",
  colors = genus_colors
) +
  labs(title = "B. Pling subcommunities, by genus")

pling_country <- plot_network(
  graph_layout = pling_layout,
  category = "country",
  colors = country_colors
) +
  labs(title = "D. Pling subcommunities, by country")

# parse mash results ----

mash <- read.delim(
  snakemake@input[["mash"]],
  sep = "\t",
  header = FALSE,
  col.names = c("seq1", "seq2", "distance", "p_val", "shared_hashes")
) |>
  mutate(
    seq1 = str_extract(seq1, "([^/]*).fasta$", group = 1),
    seq2 = str_extract(seq2, "([^/]*).fasta$", group = 1)
  )

mash_l <- mash |>
  select(seq1, seq2, distance) |>
  rename(
    accession = seq1,
    nearest_neighbour = seq2,
    mash_distance = distance
  )

mash_r <-
  mash |>
  select(seq1, seq2, distance) |>
  rename(
    accession = seq2,
    nearest_neighbour = seq1,
    mash_distance = distance
  )

mash_nearest_neighbour <- bind_rows(mash_l, mash_r) |>
  arrange(mash_distance) |>
  group_by(accession) |>
  group_modify(~ head(.x, n = 1L))

mash <- mash |>
  filter(distance <= snakemake@params[["mash_cluster_cutoff"]])

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
  category = "genus",
  colors = genus_colors
) +
  labs(title = "A. Mash clusters, by genus")

mash_country <- plot_network(
  graph_layout = mash_layout,
  category = "country",
  colors = country_colors
) +
  labs(title = "C. Mash clusters, by country")

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
  snakemake@output[["figure"]],
  plot = multi_panel, width = 180, unit = "mm"
)

node_table <- pling_graph |>
  activate(nodes) |>
  as_tibble() |>
  select(
    accession, group
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
  file = snakemake@output[["node_table"]],
  row.names = FALSE
)

bind_cols(
  node_table |>
    summarise(ari_incl_singletons = adjustedRandIndex(group_pling, group_mash)),
  node_table |>
    filter(group_pling != "Singleton" & group_mash != "Singleton") |>
    summarise(ari_excl_singletons = adjustedRandIndex(group_pling, group_mash)),
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
    summarise(
      ari_numbered_singletons = adjustedRandIndex(group_pling, group_mash)
    )
) |>
  write.csv(
    file = snakemake@output[["ari"]],
    row.names = FALSE
  )
