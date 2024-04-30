heatmap_mutation_frequency_bin = function(
    regions = NULL,
    regions_df = NULL,
    these_samples_metadata,
    these_sample_ids,
    maf,
    mut_freq_matrix,
    seq_type = "genome",
    projection = "grch37",
    region_padding = 1000,
    drop_unmutated = FALSE,
    metadataColumns = c("pathology"),
    sortByColumns = c("pathology"),
    expressionColumns = NULL,
    orientation = "sample_rows",
    skip_regions,
    only_regions,
    customColours = NULL,
    backgroundColour = "grey90",
    slide_by = 100,
    window_size = 500,
    min_count_per_bin = 0,
    min_bin_recurrence = 5,
    min_mut_tumour = 0,
    region_fontsize = 8,
    cluster_rows_heatmap = FALSE,
    cluster_cols_heatmap = FALSE,
    show_gene_colours = FALSE,
    label_regions_by = "gene",
    legend_row = 3,
    legend_col = 3,
    legend_direction = "horizontal",
    legendFontSize = 10,
    legend_side = "bottom",
    return_heatmap_obj = FALSE,
    from_indexed_flatfile = TRUE,
    mode = "slms-3"
){

    # Get region specifications
    if(missing(skip_regions)){skip_regions = NULL}
    if(missing(only_regions)){only_regions = NULL}
    regions <- process_regions(
        regions = regions,
        regions_df = regions_df,
        region_padding = region_padding,
        skip_regions,
        only_regions
    )
    regions_df <- regions$regions_df
    regions <- regions$regions


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


    # Ensure all requested metadata columns are present in the metadata
    allMetaCols <- unique(c(metadataColumns, sortByColumns, expressionColumns))

    if(!min(allMetaCols %in% colnames(metadata))){
        stop("Not all requested columns are present in the metadata table. ")
    }

    if(!missing(mut_freq_matrix)){
        samples_in_colnames = max(colnames(mut_freq_matrix) %in% these_sample_ids)
        samples_in_rownames = max(rownames(mut_freq_matrix) %in% these_sample_ids)
        if(samples_in_colnames == 1){
            all_matrix <- mut_freq_matrix
            matrix_samps <- colnames(mut_freq_matrix)
        } else if(samples_in_rownames == 1){
            all_matrix <- t(mut_freq_matrix)
            rownames(all_matrix) <- colnames(mut_freq_matrix)
            colnames(all_matrix) <- rownames(mut_freq_matrix)
        } else {
            stop("Error: sample IDs are not present in either rownames or colnames of provided matrix. ")
        }
    } else {
        # Obtain sliding window mutation frequencies for all regions
        if(missing(maf)){maf = NULL}
        all_wide <- get_mutation_frequency_bins(
            maf = maf,
            regions_df = regions_df,
            these_samples_metadata = metadata,
            projection = projection,
            slide_by = slide_by,
            window_size = window_size,
            drop_unmutated = drop_unmutated,
            return_format = "matrix",
            from_indexed_flatfile = from_indexed_flatfile,
            mode = mode
        )

        # Convert to a matrix with samples in colnames and bins in rownames
        all_matrix <- data.frame(t(select(all_wide, -sample_id)))
        colnames(all_matrix) <- all_wide$sample_id
        rownames(all_matrix) <- colnames(all_wide)[2:length(colnames(all_wide))]

    }


    # Normalize the expression columns
    if(!is.null(expressionColumns)){
        these_samples_metadata = these_samples_metadata %>%
        mutate(across(all_of(expressionColumns), ~ trim_scale_expression(.x)))
    }

    # Subset metadata to specified display columns

    meta_show = metadata %>%
        select(sample_id, all_of(allMetaCols)) %>%
        drop_na() %>%
        mutate(across(all_of(allMetaCols), ~ factor(.x))) %>%
        arrange(across(all_of(sortByColumns))) %>%
        dplyr::filter(sample_id %in% colnames(all_matrix)) %>%
        column_to_rownames(var = "sample_id")

    message(paste("starting with", length(colnames(all_matrix)), "samples"))
    samples_show <- colnames(all_matrix)[which(colSums(all_matrix)>= min_mut_tumour)]
    message(paste("returning matrix with", length(samples_show), "samples"))
    if(length(samples_show) < 2){
        stop("Insufficient samples remaining after filtering. This function requries at least two samples. ")
    }
    meta_show <- meta_show[rownames(meta_show) %in% samples_show, , drop = FALSE]
    matrix_show <- all_matrix[which(rowSums(all_matrix) > min_bin_recurrence), rownames(meta_show)]
    identical(rownames(meta_show), colnames(matrix_show))

    # Set heatmap colour function
    bin_col_fun = colorRamp2(
        c(0, 3, 6, 9),
        c(backgroundColour, "orange", "red", "purple")
    )

    # Handle custom colours
    annoColumns <- unique(c(metadataColumns, sortByColumns))
    # Get columns with no custom colours specified
    needsColour <- annoColumns[!annoColumns %in% names(customColours)]

    gamblColours <- NULL
    if(length(needsColour) > 0){
        gamblColours <- lapply(needsColour, function(x){
            colours <- get_gambl_colours()[levels(meta_show[[x]])]
            colours <- colours[unique(names(colours))][!is.na(names(colours))]
        })
        names(gamblColours) <- needsColour
    }
    annoColoursTmp <- append(gamblColours, customColours)

    # Check that there are enough colour values for each annotation
    annoColours <- lapply(annoColumns, function(x){
        if(length(levels(meta_show[[x]])[
            !levels(meta_show[[x]]) %in% names(annoColoursTmp[[x]])
        ] > 1)) {
            message(paste(
                "Warning: Insufficient values available for annotation ",
                x,
                "- using default colours. "
            ))
            colours <- get_gambl_colours("blood")[1:length(levels(meta_show[[x]]))]
            names(colours) <- levels(meta_show[[x]])
        } else {
            return(annoColoursTmp[[x]])
        }
        return(colours)
    })
    names(annoColours) <- annoColumns

    # Add colour functions for expression columns
    if(!is.null(expressionColumns)){
        col_fun = circlize::colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))
        exprColours <- lapply(expressionColumns, function(x){
        return(col_fun)
        })
        names(exprColours) <- expressionColumns
        annoColours <- append(annoColours, exprColours)
    }

    # assign bins back to regions for better annotation
    assign_bins_to_region = function(bin_names, rdf, label_by = "name"){
        bin_df = data.frame(bin_name = bin_names) %>%
        separate(bin_name, into = c("chrom", "start")) %>%
        mutate(start = as.integer(start)) %>%
        mutate(end = start + 1) %>%
        mutate(bin_name = bin_names)

        regions.dt = rdf %>%
        mutate(
            start = start - window_size,
            end = end + window_size
        ) %>%
        as.data.table()
        setkey(regions.dt, chrom, start, end)

        bin.dt = as.data.table(bin_df)
        setkey(bin.dt, chrom, start, end)
        bin_overlapped = foverlaps(bin.dt, regions.dt) %>%
        as.data.frame() %>%
        dplyr::rename(label = !!sym(label_by)) %>%
        arrange(label, start) %>%
        select(bin_name, label) %>% 
        distinct(bin_name, .keep_all = TRUE) %>%
        column_to_rownames(var = "bin_name") %>%
        drop_na() 

        return(bin_overlapped)
    }

    bin_annot = assign_bins_to_region(
        bin_names = rownames(matrix_show),
        rdf = regions_df,
        label_by = label_regions_by
    )
    heatmap_legend_param = list(
        title = "Mutation count",
        at = c(0, 2, 4, 6, 8, 10),
        nrow = legend_row,
        ncol = legend_col,
        legend_direction = legend_direction,
        labels_gp = gpar(fontsize = legendFontSize)
    )

    annotation_legend_param = list(
        nrow = legend_row,
        ncol = legend_col,
        direction = legend_direction,
        labels_gp = gpar(fontsize = legendFontSize)
    )

    if(orientation == "sample_rows"){

        to_show_t = t(matrix_show[rownames(bin_annot), ])

        row_annot = HeatmapAnnotation(
        df = meta_show,
        show_legend = T,
        which = 'row',
        col = annoColours,
        annotation_legend_param = annotation_legend_param
        )
        if(show_gene_colours){
        col_annot = HeatmapAnnotation(
            df = bin_annot,
            show_legend = F,
            which = 'col',
            annotation_legend_param = annotation_legend_param
        )
        }else{
        col_annot = HeatmapAnnotation(value = anno_empty(border = FALSE))
        }
        ht <- Heatmap(
        to_show_t[rownames(meta_show), rownames(bin_annot)],
                cluster_columns = cluster_cols_heatmap,
                cluster_rows = cluster_rows_heatmap,
                col = bin_col_fun,
                bottom_annotation = col_annot,
                left_annotation = row_annot,
                show_row_names = F,
                show_column_names = F,
                column_split = factor(bin_annot$label),
                column_title_gp = gpar(fontsize = region_fontsize),
                column_title_rot = 0,
                row_title_gp = gpar(fontsize = 10),
                heatmap_legend_param = heatmap_legend_param
        )

    }else{
        col_annot = HeatmapAnnotation(
        df = meta_show,
        show_legend = T,
        which = 'col',
        col = annoColours,
        annotation_legend_param = annotation_legend_param
        )
        if(show_gene_colours){
        row_annot = HeatmapAnnotation(
            df = bin_annot,
            show_legend = F,
            which = 'row',
            annotation_legend_param = annotation_legend_param
        )
        }else{
        row_annot = rowAnnotation(value = anno_empty(border = FALSE))
        }
        ht <- Heatmap(
        as.matrix(matrix_show)[rownames(bin_annot),rownames(meta_show)],
        show_heatmap_legend = T,
        cluster_columns = cluster_rows_heatmap,
        cluster_rows = cluster_cols_heatmap,
        col = bin_col_fun,
        bottom_annotation = col_annot,
        left_annotation = row_annot,
        show_row_names = F,
        show_column_names = F,
        row_split = factor(bin_annot$label),
        row_title_gp = gpar(fontsize = region_fontsize),
        row_title_rot = 0,
        column_title_gp = gpar(fontsize = 8),
        heatmap_legend_param = heatmap_legend_param)
    }

    if(return_heatmap_obj){
        return(ht)
    } else {
        draw(ht, heatmap_legend_side = legend_side, annotation_legend_side = legend_side)
    }

}
