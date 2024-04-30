
calc_mut_freq_sliding_windows <- function(
        this_region,
        chromosome,
        start_pos,
        end_pos,
        these_samples_metadata,
        these_sample_ids,
        maf = NULL,
        projection = "grch37",
        slide_by = 100,
        window_size = 1000,
        return_format = "long-simple",
        min_count_per_bin = 0,
        return_count = TRUE,
        drop_unmutated = FALSE,
        from_indexed_flatfile = TRUE,
        mode = "slms-3"
    ) {
    # Create objects to describe region both as string and individual objects
    try(if (missing(this_region) & missing(chromosome)) {
        stop("No region information provided. Please provide a region as a string in the chrom:start-end format, or as individual arguments. ")
    })

    if((drop_unmutated | min_count_per_bin > 0) & return_format == "wide"){
        message("To return a wide table, all samples and windows must be kept. Ignoring drop_unmutated and min_count_per_bin flags. ")
    }

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
        metadata <- get_gambl_metadata(
            seq_type_filter = c("capture", "genome"),
            tissue_status_filter = "tumour"
        )
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

    if(
        (str_detect(chromosome, "chr") & projection == "grch37") |
            (!str_detect(chromosome, "chr") & projection == "hg38")
    ) {
        stop("chr prefixing status of region and specified projection don't match. ")
    }

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
        chrom = chromosome,
        window_start = seq(start_pos, end_pos, by = slide_by)
    ) %>%
        mutate(window_end = window_start + window_size - 1) %>%
        select(chrom, window_start, window_end)

    # Option to return full region count instead of sliding window
    if(window_size == 0){
        windows <- data.frame(
            chrom = chromosome,
            window_start = start_pos,
            window_end = end_pos
        )
    }

    # Obtain SSM coordinates from GAMBL if no maf was provided
    if (is.null(maf)) {
        try(
            if(!"seq_type" %in% colnames(metadata))
            stop("seq_type must be present in metdata for compatibility with get_ssm_by_sample")
        )
        message("Using GAMBLR::get_ssm_by_region...")
        region_ssm <- list()
        for (st in unique(metadata$seq_type)){
            this_seq_type <- GAMBLR::get_ssm_by_region(
                region = this_region,
                projection = projection,
                streamlined = FALSE,
                seq_type = st,
                from_indexed_flatfile = TRUE,
                mode = "slms-3"
            ) %>%
                dplyr::mutate(end = Start_Position + 1) %>%
                dplyr::select(
                    chrom = Chromosome,
                    start = Start_Position,
                    end,
                    sample_id = Tumor_Sample_Barcode
                ) %>%
                mutate(mutated = 1, seq_type = st) %>%
                filter(sample_id %in% these_sample_ids)

            region_ssm[[st]] <- data.frame(metadata) %>%
                select(sample_id, seq_type) %>%
                filter(seq_type == st) %>%
                left_join(this_seq_type) %>%
                filter(!is.na(mutated)) %>%
                select(-seq_type)
        }
        region_ssm <- bind_rows(region_ssm)
    } else {
        #  Subset provided maf to specified region
        message("Using provided maf...")
        maf.dt <- data.table(maf)
        region_bed <- data.table(
            "Chromosome" = as.character(chromosome),
            "Start_Position" = as.numeric(start_pos),
            "End_Position" = as.numeric(end_pos)
        )
        setkey(region_bed)
        region_ssm <- foverlaps(maf.dt, region_bed) %>%
            dplyr::filter(!is.na(Start_Position)) %>%
            mutate(end = i.Start_Position - 1) %>%
            dplyr::select(
                chrom = Chromosome,
                start = i.Start_Position,
                end,
                sample_id = Tumor_Sample_Barcode
            ) %>%
            dplyr::mutate(mutated = 1)

        region_ssm <- data.frame(metadata) %>%
            select(sample_id) %>%
            left_join(region_ssm) %>%
            filter(!is.na(mutated))
    }

    # Check if the region is empty.
    # If yes return NULL so that running this function with lapply will allow bind_rows to run on the output.
    if(nrow(region_ssm) == 0 & (drop_unmutated | min_count_per_bin > 0)){
        message(paste0("No mutations found in region ", this_region, " for this sample set. "))
        return(NULL)
    }

    # Count mutations per window
    windows_tallied <- inner_join(
        windows,
        region_ssm,
        by = "chrom"
    ) %>%
        dplyr::filter(
            start >= window_start,
            start <= window_end
        ) %>%
        group_by(
            sample_id,
            window_start
        ) %>%
        tally() %>%
        ungroup() %>%
        full_join(select(metadata, sample_id)) %>%
        arrange(sample_id) %>%
        full_join(select(windows, window_start)) %>%
        distinct() %>%
        pivot_wider(
            names_from = window_start,
            values_from = n,
            values_fill = 0
        ) %>%
        select(-matches("^NA$")) %>%
        pivot_longer(
            -c(sample_id),
            names_to = "window_start",
            values_to = "n"
        ) %>%
        distinct() %>%
        filter(!is.na(sample_id))

    # Remove unmutated windows if requested
    if (drop_unmutated | min_count_per_bin > 0) {
        windows_tallied <- windows_tallied %>%
            filter(n >= min_count_per_bin)
        if(drop_unmutated & min_count_per_bin == 0){
            windows_tallied %>%
                filter(n > 0)
        }
    }

    # Create requested data output format
    if (return_count) {
        # Return table of mutation counts per bin
        windows_tallied_final <- mutate(
            windows_tallied,
            bin = paste0(chromosome, "_", window_start)
        ) %>%
            mutate(mutation_count = n) %>%
            select(
                sample_id,
                bin,
                mutation_count
            )
    } else {
        # Return table of binary mutated/unmutated status per bin
        windows_tallied_final <- mutate(
            windows_tallied,
            bin = paste0(chromosome, "_", window_start)
        ) %>%
            mutate(mutated = ifelse(n > 0, 1, 0)) %>%
            select(
                sample_id,
                bin,
                mutated
            )
    }

    if (return_format == "wide") {
        widened <- windows_tallied_final %>%
            pivot_wider(
                names_from = bin,
                values_from = matches("mutat"),
                values_fill = 0
            )
        return(widened)
    } else {
        return(windows_tallied_final)
    }
}
