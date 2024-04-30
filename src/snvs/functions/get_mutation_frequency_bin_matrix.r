get_mutation_frequency_bins = function(
  regions = NULL,
  regions_df = NULL,
  these_samples_metadata,
  these_sample_ids,
  maf,
  projection = "grch37",
  region_padding = 1000,
  drop_unmutated = FALSE,
  skip_regions = NULL,
  only_regions = NULL,
  slide_by = 100,
  window_size = 500,
  return_format = "matrix",
  from_indexed_flatfile = TRUE,
  mode = "slms-3"
){

  regions <- process_regions(
    regions = regions,
    regions_df = regions_df,
    region_padding = region_padding,
    skip_regions = skip_regions,
    only_regions = only_regions
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

  # Obtain sliding window mutation frequencies for all regions
  # Base R lapply crashes with the IGH region on a maf with only 1000 samples.
  # mclapply is not affected by this issue. It's also faster and more efficient.
  if(missing(maf)){maf = NULL}
  dfs = parallel::mclapply(names(regions), function(x){
    df <- calc_mut_freq_sliding_windows(
        this_region = regions[x],
        these_samples_metadata = metadata,
        maf = maf,
        projection = projection,
        drop_unmutated = drop_unmutated,
        slide_by = slide_by,
        window_size = window_size,
        min_count_per_bin = 0,
        return_count = TRUE,
        from_indexed_flatfile = from_indexed_flatfile,
        mode = mode
    ) %>%
    mutate(name = x)
    return(df)
  })

  all <- dplyr::bind_rows(dfs) %>%
    dplyr::distinct(bin, sample_id, .keep_all = TRUE)

  # If none of the samples are mutated, return the mutation frequency df and exit.
  if(max(all$mutation_count) == 0){
    message("No mutations found in specified regions for specified samples. Exiting. ")
    return(all)
  }

  if(return_format == "matrix"){
    # Convert mutation frequency table to a matrix
    all_wide <- all %>%
      select(sample_id, mutation_count, bin) %>%
      pivot_wider(
          names_from = bin,
          values_from = mutation_count,
          values_fill = 0
      )
      return(all_wide)
    } else {
      return(all)
    }
}
