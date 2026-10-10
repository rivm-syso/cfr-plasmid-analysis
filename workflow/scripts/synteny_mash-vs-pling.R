# Libraries ----

library(tidyverse)
library(gggenomes)
library(ggtext)
library(ggrepel)
library(patchwork)

# logging ----

log <- file(snakemake@log[[1]], open = "wt")
sink(log)

# Plotting function ----

synteny_plot <- function(gbff, links, mobtyper) {
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
      seq_id = factor(seq_id, ordered = TRUE)
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
    geom_seq_label(size = 8 / .pt) +
    scale_x_bp(suffix = "bp") +
    scale_y_discrete(expand = c(.05, .05)) +
    theme

  return(plot)
}

# file parsing ----

parse_biomarkers <- function(file) {
  sample_id <- stringr::str_remove(basename(file), fixed(".tab"))
  df <- tryCatch(
    {
      read.delim(file, sep = "\t")
    },
    error = function(cond) {
      print(cond)
      return(tibble())
    }
  )
  if (nrow(df) != 0) {
    df |>
      mutate(
        seq_id = sample_id[[1]][1],
        strand = case_when(
          sstrand == "plus" ~ "+",
          sstrand == "minus" ~ "-"
        ),
        type = "CDS"
      ) |>
      rename(
        start = sstart,
        end = send,
        gene_type = biomarker
      ) |>
      select(
        seq_id, start, end, strand, gene_type, type
      )
  } else {
    df
  }
}

parse_links <- function(file) {
  df <- read.csv(file) |>
    rename(
      seq_id = qseqid,
      length = qlen,
      start = qstart,
      end = qend,
      seq_id2 = sseqid,
      length2 = slen,
      start2 = sstart,
      end2 = send,
      strand2 = sstrand
    ) |>
    mutate(
      strand = "+",
      strand2 = case_when(
        strand2 == "plus" ~ "+",
        strand2 == "minus" ~ "-"
      )
    )
  df
}

# Panel A ----

gbff <- snakemake@input[["genbank_a"]] |>
  map_dfr(read_gbk)

links <- snakemake@input[["links_a"]] |>
  map_dfr(parse_links)

mobtyper <- snakemake@input[["biomarkers_a"]] |>
  map_dfr(parse_biomarkers)

panel_a <- synteny_plot(
  gbff, links, mobtyper
)

# Panel B ----

gbff <- snakemake@input[["genbank_b"]] |>
  map_dfr(read_gbk)

links <- snakemake@input[["links_b"]] |>
  map_dfr(parse_links)

mobtyper <- snakemake@input[["biomarkers_b"]] |>
  map_dfr(parse_biomarkers)

panel_b <- synteny_plot(
  gbff, links, mobtyper
)


# Panel C ----

gbff <- snakemake@input[["genbank_c"]] |>
  map_dfr(read_gbk)

links <- snakemake@input[["links_c"]] |>
  map_dfr(parse_links)

mobtyper <- snakemake@input[["biomarkers_c"]] |>
  map_dfr(parse_biomarkers)

panel_c <- synteny_plot(
  gbff, links, mobtyper
)

# Multipanel ----

design <- "
AAAAAADD
AAAAAADD
BBBBBBDD
BBBBBBDD
CCCCCCCC
CCCCCCCC
"
multipanel <- panel_a + panel_b + panel_c + guide_area() +
  plot_layout(guides = "collect", design = design) +
  plot_annotation(tag_levels = "A")

ggsave(snakemake@output[[1]], plot = multipanel, width = 180, unit = "mm")
