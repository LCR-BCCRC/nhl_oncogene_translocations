prepare_circos <- function(breakpoints,
                           genes,
                           region_order = NULL,
                           sample_ids = NULL,
                           maf_id_col = "Tumor_Sample_Barcode",
                           breakpoint_id_col = "sample_id",
                           maf = NULL,
                           partner_upstream_buffer = 1e5,
                           partner_downstream_buffer = 1e5,
                           target_upstream_buffer = 0,
                           target_downstream_buffer = 5e5) {
  circos <- list()

  # Get data for circos links
  circos$links <- as.data.table(breakpoints)

  circos$target <- unique(circos$links$target)
  circos$chrom_target <- unique(circos$links$chrom_target)

  # Store genes in the circos object for use by plot_circos
  circos$genes <- data.table(genes)

  # Obtain the regions that encompass translocations
  # and define their correspondence to the chromosomes
  circos$correspondance <-
    circos$genes %>%
    filter(region != circos$target) %>%
    group_by(region, chrom) %>%
    summarise(
      start = min(start) - partner_upstream_buffer,
      end = max(end) + partner_downstream_buffer
    ) %>%
    ungroup() %>%
    add_row(
      chrom = circos$chrom_target,
      start = min(circos$genes[circos$genes$region == circos$target, ]$start) - target_upstream_buffer,
      end = max(circos$genes[circos$genes$region == circos$target, ]$end) + target_downstream_buffer,
      region = circos$target
    ) %>%
    mutate(
      start2 = start,
      end2 = end
    ) %>%
    select(chrom, start, end, region, start2, end2) %>%
    arrange(as.numeric(str_remove(chrom, "chr")), start) %>%
    as.data.table()

  # Obtain the chromosome regions for the outer ring
  circos$regions <- circos$correspondance$chrom

  # Handle custom region_order specification
  if (!is.null(region_order)) {
    circos$regions <- circos$regions[
      order(match(circos$regions, region_order))
    ]
    circos$correspondance <- circos$correspondance[
      order(match(circos$correspondance$chrom, region_order))
    ]
  }

  # Get data for circos links
  # Ensure they don't go outside of specified plotting regions
  region_dt <- circos$correspondance[, 1:3]
  setkey(region_dt, chrom, start, end)
  breakpoints <- as.data.table(breakpoints)
  breakpoints_1 <- foverlaps(
    breakpoints,
    region_dt,
    by.x = c("chrom_target", "start_target", "end_target")
  ) %>%
    filter(!is.na(start)) %>%
    select(-start, -end)

  breakpoints_2 <- foverlaps(
    breakpoints_1,
    region_dt,
    by.x = c("chrom_partner", "start_partner", "end_partner")
  ) %>%
    filter(!is.na(start)) %>%
    filter(partner %in% circos$correspondance$region) %>%
    select(-start, -end)

  circos$links <- data.table(breakpoints_2)

  if (!is.null(sample_ids)) {
    circos$links <- circos$links[
      str_detect(
        circos$links[[breakpoint_id_col]],
        paste0(sample_ids, collapse = "|")
      ),
    ]
  }


  # Rearrange the gene bed file to use region instead of chromosome
  circos$genes <- circos$genes[, c("region", "start", "end", "name", "colour")]

  # Subset maf to plotting regions
  if (!is.null(maf)) {
    maf <- maf %>%
      select(
        chrom = Chromosome,
        start = Start_Position,
        end = End_Position,
        Mutation_Overlap_WRCY,
        all_of(maf_id_col)
      ) %>%
      data.table()
    region_dt <- circos$correspondance[, 1:4]
    setkey(region_dt, chrom, start, end)
    subset_maf <- foverlaps(maf, region_dt)
    subset_maf <- subset_maf %>%
      filter(!is.na(start)) %>%
      select(
        region,
        start = i.start,
        end = i.end,
        Mutation_Overlap_WRCY,
        all_of(maf_id_col)
      )

    if (!is.null(sample_ids)) {
      subset_maf <- subset_maf[subset_maf[[maf_id_col]] %in% sample_ids, ]
    }
    circos$rainfall <- list()
    circos$rainfall$other <- filter(subset_maf, Mutation_Overlap_WRCY == "FALSE")
    circos$rainfall$WRCY <- filter(subset_maf, Mutation_Overlap_WRCY != "FALSE")
  }

  circos
}
