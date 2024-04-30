plot_circos <- function(circos,
                        prop_chroms,
                        prop_regions = prop_chroms,
                        angle = 155,
                        link_colours,
                        nest_colours,
                        colour_links_by,
                        type_links_by,
                        title = NULL,
                        legend = FALSE) {
  circos.clear()

  # For outer ring:
  #   - Plot ideogram for MYC and IG chromosomes
  outer_ring_fun <- function() {
    circos.par(gap.degree = 5, start.degree = angle)

    circos.initializeWithIdeogram(
      plotType = c("ideogram", "labels"),
      species = "hg19", chromosome.index = unique(circos$regions),
      sector.width = prop_chroms
    )

    title(title)
  }

  # For inner ring:
  #   - Plot genes in MYC and IG regions
  #   - Plot rainfall plot, highlighting mutations in AID motif (in red)
  #   - Plot IG-MYC translocations as links
  inner_ring_fun <- function() {
    circos.par(track.height = 0.2, gap.degree = 3, start.degree = angle)

    circos.genomicInitialize(circos$correspondance[, 4:6],
      plotType = NULL,
      sector.width = prop_regions, 
      tickLabelsStartFromZero = FALSE
    )
    
    circos.genomicTrackPlotRegion(
      circos$correspondance[, 4:6],
      ylim = c(0, 1),
      track.height = 0.08,
      bg.border = NA,
      panel.fun = function(region, value, ...) {
        seq_by = 1e5
        if(max(region$end) - min(region$start) > 6e5){
          seq_by = 2e5
        }  
        if(max(region$end) - min(region$start) > 1e6){
          seq_by = 5e5
        }
        axis_start = (round(min(region$start) / 1e6, digits = 1) * 1e6) - seq_by
        axis_values = seq(axis_start, max(region$end) + seq_by, by = seq_by)
        axis_labels = paste0(axis_values / 1e6, "MB")
        circos.genomicAxis(
          tickLabelsStartFromZero = FALSE,
          major.at = axis_values,
          labels = axis_labels,
          direction = "outside",
          h = "bottom"
        )
      }
    )

    circos.genomicTrackPlotRegion(circos$genes,
      ylim = c(0, 1),
      track.height = 0.2,
      panel.fun = function(region, value, ...) {
        all_tx = unique(value$name)
        for(i in seq_along(all_tx)) {
            l = value$name == all_tx[i]
            # for each transcript
            current_tx_start = min(region[l, 1])
            current_tx_end = max(region[l, 2]) 
            label_pos = mean(c(current_tx_start, current_tx_end))
            circos.lines(c(current_tx_start, current_tx_end), 
                c(0.2, 0.2), col = "black")
            circos.text(
              x = label_pos, y = 0.7, 
              labels = all_tx[i], cex = 0.8,
              niceFacing = TRUE, facing = "bending"
            )
            circos.genomicRect(region, ytop = 0.35, 
                ybottom = 0.05, col = value$colour, border = "black")
        }
      }
    )
    
    if ("rainfall" %in% names(circos)) {
      circos.genomicRainfall(
          circos$rainfall,
          pch = 16, cex = 0.1, ylim = c(0, 4.5),
          col = c("black", "firebrick"), 
          track.height = 0.2
      )
    }


    linetypes <- c("1", "2")
    names(linetypes) <- c("TRUE", "FALSE")

    circos.genomicLink(select(circos$links, target, start_target, end_target),
      select(circos$links, partner, start_partner, end_partner),
      col = link_colours[circos$links[[colour_links_by]]],
      lwd = 0.5
      # lty = circos$links[[type_links_by]]
    )
  }

  circos.nested(outer_ring_fun, inner_ring_fun, circos$correspondance,
    connection_col = nest_colours[circos$correspondance[[4]]],
    adjust_start_degree = FALSE,
    connection_height = 0.08
  )

  if (legend) {
    lgd <- Legend(
      labels = names(link_colours),
      type = "lines",
      legend_gp = gpar(
        col = link_colours,
        lex = 2
      )
    )

    draw(lgd,
      x = unit(dev.size("in")[1] * 0.05, "in"),
      y = unit(dev.size("in")[2] * 0.95, "in"),
      just = c("left", "top")
    )
  }
}
