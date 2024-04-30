calc_mut_freq_sliding_windows <- function(
    this_region,
    chromosome,
    start_pos,
    end_pos,
    these_samples_metadata,
    these_sample_ids,
    maf,
    seq_type = "genome",
    slide_by = 100,
    window_size = 1000,
    plot_type = "none",
    sortByColumns = "pathology",
    classification_column = "lymphgen",
    custom_colours,
    return_format = "long-simple",
    min_count_per_bin = 0,
    return_count = TRUE,
    drop_unmutated = FALSE,
    from_indexed_flatfile = TRUE,
    mode = "slms-3") {
    # Create objects to describe region both as string and individual objects
    try(if (missing(this_region) & missing(chromosome)) {
        stop("No region information provided. Please provide a region as a string in the chrom:start-end format, or as individual arguments. ")
    })

    if (missing(this_region)) {
        this_region <- paste0(
            chromosome, ":", start_pos, "-",
            end_pos
        )
    } else {
        chunks <- GAMBLR:::region_to_chunks(this_region)
        chromosome <- chunks$chromosome
        start_pos <- as.numeric(chunks$start)
        end_pos <- as.numeric(chunks$end)
    }

    # Check for provided metadata, else use GAMBL metadata
    if (missing(these_samples_metadata)) {
        metadata <- get_gambl_metadata()
    } else {
        metadata <- these_samples_metadata
    }
    # Ensure metadata is subset to specified sample IDs
    if (!missing(these_sample_ids)) {
        metadata <- dplyr::filter(
            metadata,
            sample_id %in% these_sample_ids
        )
    } else {
        these_sample_ids <- metadata$sample_id
    }

    # Subset the metadata to specified columns
    metadata_tidy <- data.frame(
        sample_id = metadata$sample_id,
        sort_column = metadata[[sortByColumns]],
        classification_column = metadata[[classification_column]]
    ) %>%
        mutate(region_start = start_pos)
    metadata <- metadata_tidy

    # Check region size and compare to max region size
    # Is this really needed?
    max_region <- 5e+06

    region_size <- end_pos - start_pos
    if (region_size < max_region) {
        message(paste(
            "processing bins of size", window_size,
            "across", region_size, "bp region"
        ))
    } else {
        message(paste("CAUTION!\n", region_size, "exceeds maximum size recommended by this function."))
    }

    # Split region into windows
    windows <- data.frame(
        start = seq(start_pos, end_pos, by = slide_by)
    ) %>%
        mutate(end = start + window_size - 1)
    windows.dt <- as.data.table(windows)

    # Obtain SSM coordinates from GAMBL if no maf was provided
    if (missing(maf)) {
        region_ssm <- GAMBLR::get_ssm_by_region(
            region = this_region,
            streamlined = FALSE,
            seq_type = seq_type,
            from_indexed_flatfile = from_indexed_flatfile,
            mode = mode
        ) %>%
            dplyr::select(
                start = Start_Position,
                sample_id = Tumor_Sample_Barcode
            ) %>%
            mutate(mutated = 1) %>%
            filter(sample_id %in% these_sample_ids)
    } else {
        #  Subset provided maf to specified region
        maf.dt <- data.table(maf)
        region_bed <- data.table(
            "Chromosome" = str_split(this_region, ":")[[1]][1],
            "Start_Position" = as.numeric(str_split(str_split(this_region, ":")[[1]][2], "-")[[1]][1]),
            "End_Position" = as.numeric(str_split(str_split(this_region, ":")[[1]][2], "-")[[1]][2])
        )
        setkey(region_bed)
        region_ssm <- foverlaps(maf.dt, region_bed) %>%
            dplyr::filter(!is.na(Start_Position)) %>%
            dplyr::select(
                start = i.Start_Position,
                sample_id = Tumor_Sample_Barcode
            ) %>%
            dplyr::mutate(mutated = 1)
    }

    # Ensure SSM data are in data.table format
    region.dt <- region_ssm %>%
        dplyr::mutate(
            start = as.numeric(as.character(start)),
            end = start + 1, end = as.numeric(as.character(end))
        ) %>%
        dplyr::relocate(start, .before = end) %>%
        as.data.table()

    # Overlap SSM data with windows table
    setkey(windows.dt, start, end)
    setkey(region.dt, start, end)

    # Count mutations per window
    windows_tallied <- foverlaps(windows.dt, region.dt) %>%
        dplyr::filter(!is.na(start)) %>%
        dplyr::rename(
            window_start = "i.start",
            mutation_position = "start"
        ) %>%
        dplyr::select(-i.end, -end, -mutation_position) %>%
        as.data.frame() %>%
        group_by(
            sample_id,
            window_start
        ) %>%
        tally() %>%
        arrange(sample_id) %>%
        as.data.frame()

    # Join with metadata
    windows_tallied_md <- windows_tallied %>%
        full_join(metadata) %>%
        full_join(windows.dt[, "start"], by = c("window_start" = "start")) %>%
        pivot_wider(
            names_from = window_start,
            values_from = n,
            values_fill = 0
        ) %>%
        select(-`NA`) %>%
        pivot_longer(
            -c(sample_id, sort_column, classification_column, region_start),
            names_to = "window_start",
            values_to = "n"
        ) %>%
        distinct()

    # Remove unmutated windows if requested
    if (drop_unmutated | min_count_per_bin > 0) {
        windows_tallied_md <- windows_tallied_md %>%
            filter(n >= min_count_per_bin)
    }

    # Convert sample_id to a factor and ensure levels are consistent across dfs
    metadata$sort_column <- factor(metadata$sort_column)
    metadata$classification_column <- factor(metadata$classification_column)
    metadata$sample_id <- factor(
        metadata$sample_id,
        levels = metadata[order(
            rev(metadata$sort_column),
            rev(metadata$classification_column)
        ), ]$sample_id
    )

    windows_tallied_md$sample_id <- factor(
        windows_tallied_md$sample_id,
        levels = levels(metadata$sample_id)
    )
    windows_tallied_md$sort_column <- factor(
        windows_tallied_md$sort_column,
        levels = levels(metadata$sort_column)
    )
    windows_tallied_md$classification_column <- factor(
        windows_tallied_md$classification_column,
        levels = levels(metadata$classification_column)
    )

    # Specify colours
    if (missing(custom_colours)) {
        message("using GAMBLR colours")
        class_cols <- get_gambl_colours()[
            names(get_gambl_colours()) %in% levels(metadata$classification_column)
        ]
        class_cols <- class_cols[levels(metadata$classification_column)]
        print(class_cols)
        sort_cols <- get_gambl_colours()[
            names(get_gambl_colours()) %in% levels(metadata$sort_column)
        ]
        sort_cols <- sort_cols[levels(metadata$sort_column)]
    } else {
        class_cols <- custom_colours[
            levels(metadata$classification_column)
        ]
        class_cols <- class_cols[levels(metadata$classification_column)]
        message(class_cols)
        sort_cols <- custom_colours[
            levels(metadata$sort_column)
        ]
        sort_cols <- sort_cols[levels(metadata$sort_column)]
    }

    tile_width <- (end_pos - start_pos) / 100

    if (plot_type %in% c("points", "point")) {
        # Dot plot of binary mutated/unmutated per window
        p <- windows_tallied_md %>%
            filter(!is.na(sample_id)) %>%
            ggplot2::ggplot() +
            geom_tile(
                aes(
                    x = as.numeric(region_start - tile_width),
                    y = sample_id,
                    fill = sort_column,
                    width = tile_width
                )
            ) +
            geom_point(aes(
                x = ifelse(n > 0, as.numeric(window_start), NA),
                y = sample_id,
                colour = classification_column
            )) +
            theme(axis.text = element_text(size = 4)) +
            theme(axis.text.y = element_blank()) +
            xlim(start_pos - (2 * tile_width), end_pos) +
            xlab("Position") +
            ylab("Sample")

        if (length(class_cols) == length(levels(metadata$classification_column))) {
            p <- p +
                scale_colour_manual(
                    values = c(class_cols), name = paste0(classification_column)
                )
        } else {
            message("Insufficient colours in custom scale for classification_column. Using ggplot2 default colours. ")
        }

        if (length(sort_cols) == length(levels(metadata$sort_column))) {
            p <- p +
                scale_fill_manual(
                    values = c(sort_cols), name = sortByColumns
                )
        } else {
            message("Insufficient colours in custom scale for sortByColumns. Using ggplot2 default colours. ")
        }
    } else if (plot_type %in% c("tile", "tiled")) {
        # Heatmap-style plot with gradient to represent number of mutations per window
        p <- windows_tallied_md %>%
            filter(!is.na(sample_id)) %>%
            ggplot2::ggplot() +
            geom_point(
                aes(
                    x = as.numeric(region_start - 2 * tile_width),
                    y = sample_id,
                    colour = sort_column
                ),
                shape = "square",
                size = 4
            ) +
            geom_point(
                aes(
                    x = as.numeric(region_start - tile_width),
                    y = sample_id,
                    colour = classification_column
                ),
                shape = "square",
                size = 4
            ) +
            geom_tile(
                aes(
                    x = as.numeric(window_start),
                    y = sample_id,
                    fill = as.numeric(n)
                )
            ) +
            scale_fill_gradient2(
                low = "white", mid = "orange", high = "red",
                limits = c(0, max(windows_tallied_md$n, na.rm = TRUE)),
                midpoint = max(windows_tallied_md$n, na.rm = TRUE) / 2,
                na.value = NA
            ) +
            cowplot::theme_cowplot() +
            theme(axis.text.y = element_blank()) +
            xlim(start_pos - (3 * tile_width), end_pos) +
            xlab("Position") +
            ylab("Sample")

        if (
            length(class_cols) == length(levels(metadata$classification_column)) &
                length(sort_cols) == length(levels(metadata$sort_column))
        ) {
            p <- p +
                scale_colour_manual(
                    values = c(class_cols, sort_cols), name = paste0(classification_column, "/", sortByColumns)
                )
        } else {
            message("Insufficient colours in custom scale for classification_column and/or sortByColumns. Using ggplot2 default colours. ")
        }
    }
    # Store the outputs in a list
    output <- list()
    if (plot_type != "none") {
        output$plot <- p
    }
    # Create requested data output format
    if (return_count) {
        # Return table of mutation counts per bin
        windows_tallied_final <- mutate(
            windows_tallied_md,
            bin = paste0(window_start, "_", chromosome)
        ) %>%
            mutate(mutated = n)
    } else {
        # Return table of binary mutated/unmutated status per bin
        windows_tallied_final <- mutate(
            windows_tallied_md,
            bin = paste0(window_start, "_", chromosome)
        ) %>%
            mutate(mutated = ifelse(n > 0, 1, 0))
    }

    if (return_format == "long") {
        output$data <- windows_tallied_final
    } else if (return_format == "long-simple") {
        windows_simple <- dplyr::select(
            windows_tallied_final,
            sample_id,
            bin,
            mutated
        )
        output$data <- windows_simple
    } else {
        widened <- windows_tallied_final %>%
            pivot_wider(
                names_from = bin,
                values_from = mutated,
                values_fill = 0
            )
        output$data <- widened
    }

    return(output)
}
