process_regions <- function(
    regions = NULL,
    regions_df = NULL,
    region_padding = 0,
    skip_regions = NULL,
    only_regions = NULL,
    projection = "grch37"
) {
    # Use default ashm region table if no regions are provided
    if (is.null(regions)) {
        if (is.null(regions_df)) {
            message("Using default GAMBLR aSHM regions. ")
            if (projection == "grch37") {
                regions_df <- grch37_ashm_regions %>%
                    dplyr::mutate(chr_name = str_remove(chr_name, "chr"))
            } else {
                regions_df <- hg38_ashm_regions
            }
            # Fix column names
            regions_df <- regions_df %>%
                dplyr::rename_with(
                    ~ str_remove(.x, "^hg.*_")
                ) %>%
                dplyr::rename(chrom = chr_name) 
            if (!is.null(skip_regions)) {
                # drop user-specified regions
                regions_df = regions_df %>%
                    dplyr::filter(!gene %in% skip_regions)
            }
            if (!is.null(only_regions)) {
                # keep only user-specified regions
                regions_df = regions_df %>%
                    dplyr::filter(gene %in% only_regions)
                
            }
            regions_df = dplyr::rename(regions_df, name = gene)
        }

        required_cols <- c("chrom", "start", "end", "name")
        if (min(required_cols %in% colnames(regions_df)) == 0) {
            stop("Provided regions_df lacks required column names. Ensure columns chrom, start, end, and name are present. ")
        }

        # gene column is required for later joins
        if (!"gene" %in% colnames(regions_df)) {
            regions_df <- mutate(regions_df, gene = name)
        }
    } else {
        # Convert character vector of regions to df
        regions_df <- bind_rows(lapply(regions, function(x) {
            chunks <- GAMBLR:::region_to_chunks(x)
            df <- data.frame(
                chromosome = chunks$chromosome,
                start = as.numeric(chunks$start),
                end = as.numeric(chunks$end)
            )
        }))
        if (!is.null(names(regions))) {
            regions_df$name = names(regions)
            regions_df$gene = names(regions)
        } else {
            message("WARNING: Regions provided as an unnamed character vector. It is strongly recommended to provide a named vector with a unique name per region, or a bed file-formatted data frame. ")
            regions_df$name = regions_df$chromosome
            regions_df$gene = regions_df$chromosome
        }
    }

    # Collapse regions with duplicate names
    if (length(unique(regions_df$name)) < length(regions_df$name)) {
        message("Warning: Multiple regions in the provided data frame have the same name. Merging these entries based on min(start) and max(end) per name value. ")
        regions_df <- regions_df %>%
            group_by(name) %>%
            mutate(
                start = min(start),
                end = max(end)
            ) %>%
            ungroup() %>%
            distinct()
    }

    regions = unlist(apply(
        regions_df,
        1,
        function(x) {
            # add specified padding around each region
            paste0(x[1], ":", as.numeric(x[2]) - region_padding, "-", as.numeric(x[3]) + region_padding)
        }
    ))
    names(regions) <- regions_df$name

    return(
        list(
            regions = regions,
            regions_df = regions_df
        )
    )
}
