## ggmanh_vec: Manhattan plots as vector graphics

## Description:
## Function to make ggplot object of Manhattan plot suitable for exporting as vector graphic.

## Usage:
## ggmanh_vec(df, chr = "CHR", pos = "POS", log_pv = "LOG10_P")

## Arguments:
## df:            a data.frame with columns for chromosome, SNP position and p-value (-log10-transformed)
## chr:           name of chromosome column in input data.frame
## pos:           name of SNP position column in input data.frame
## log_pv:        a string with name of -log10(p-value) column in input data.frame
## hw:            height-to-width ratio of plot
## buffer_dist:   distance parameter passed to st_buffer(); determines point size
## gap:           gap between chromsomes in Manhattan plot in physical position units
## main_cols:     vector with colour codes for odd and even chromosomes, resp.
## sig_col:       colour code for line of genome-wide significance
## sig_type:      linetype for line of genome-wide significance
## title:         plot title
## title_sz:      font size of plot title
## axis_title_sz: font size of axis tick labels
## ncores:        number of cores passed to registerDoParallel() to be used by doParallel

## Value:
## A ggplot object containing the Manhattan plot.



## Load dependencies:
library(doParallel)
library(dplyr)
library(foreach)
library(ggplot2)
library(sf)  # requires GDAL to be available on system




ggmanh_vec <- function(df,
                       chr = "CHR",
                       pos = "POS",
                       log_pv = "LOG10_P",
                       hw = 2.35,
                       buffer_dist = 1.2,
                       gap = 2e7,
                       main_cols = c("#00468B", "#0095AF"),
                       sig_col = "black",
                       sig_type = "solid",
                       title = "",
                       title_sz = 8,
                       axis_title_sz = 7,
                       axis_text_sz = 6,
                       ncores = 1) {
    
    ## Make the first SNP of each chromosome have position 0
    df <- df %>%
        group_by(.data[[chr]]) %>%
        mutate(POS_0 = .data[[pos]] - .data[[pos]][1])

    df <- df %>%
        ## Compute chromosome size
        group_by(.data[[chr]]) %>%
        summarise(chr_len = max(POS_0)) %>%
        ## Calculate cumulative position of each chromosome
        mutate(gap_add = seq(0, gap * (n() - 1), gap)) %>%
        mutate(tot = cumsum(as.numeric(chr_len)) - chr_len + gap_add) %>%
        select(-chr_len) %>%  
        ## Add this info to the initial dataset
        left_join(df, ., by = .data[[chr]]) %>%
        ## Add a cumulative position of each SNP
        arrange(.data[[chr]], POS_0) %>%
        mutate(POScum = POS_0 + tot) %>%
        select(.data[[chr]], POScum, .data[[log_pv]])

    ## Rescale SNP positions so plot has desired H/W ratio
    maxpos <- max(df$POScum)
    maxp <- max(df[, log_pv])
    df <- df %>%
        mutate(POScum = POScum / maxpos * maxp * hw)

    ## Get chromosome center positions for x-axis
    axis_pos <- df %>%
        group_by(.data[[chr]]) %>%
        summarize(centre = (max(POScum) + min(POScum)) / 2)
    

    registerDoParallel(ncores)
    ## Odd/even chr must end up in separate polygons so can have different colours
    odd_chr <- list()
    P <- list()
    st <- list(); en <- list()
    odd_chr_seg <- list()

    odd_out <- foreach (o=seq(1, 21, 2)) %dopar% {
        P[[o]] <- nrow(df %>% filter(CHR == o))

        if (P[[o]] > 10000) {
            if ((P[[o]] %% 10000) != 0) {
                ## Seq of intervals into which to slice data for faster processing
                st[[o]] <- seq(1, P[[o]], 10000)
                en[[o]] <- c(seq(10000, P[[o]], 10000), P[[o]])
            } else {
                st[[o]] <- seq(1, P[[o]], 10000)
                en[[o]] <- seq(10000, P[[o]], 10000)
            }

            odd_chr[[o]] <- df %>%
                filter(CHR == o) %>%
                slice(st[[o]][1]:en[[o]][1]) %>%
                sf::st_as_sf(coords = c(2, 3)) %>%
                sf::st_buffer(dist = buffer_dist) %>%
                sf::st_union()
            if (length(st[[o]]) > 1) {
                for (s in 2:length(st[[o]])) {
                    odd_chr_seg[[o]] <- df %>%
                        filter(CHR == o) %>%
                        slice(st[[o]][s]:en[[o]][s]) %>%
                        sf::st_as_sf(coords = c(2, 3)) %>%
                        sf::st_buffer(dist = buffer_dist) %>%
                        sf::st_union()
                    odd_chr[[o]] <- sf::st_union(odd_chr[[o]], odd_chr_seg[[o]])
                }
            }
        } else {
            odd_chr[[o]] <- df %>%
                filter(CHR %in% o) %>%
                sf::st_as_sf(coords = c(2, 3)) %>%
                sf::st_buffer(dist = buffer_dist) %>%
                sf::st_union()
        }
        return(odd_chr[[o]])
    } 
 
    for (o in 1:length(odd_out)) {
        if (o == 1) {
            df_com_odd <- odd_out[[o]]
        } else {
            df_com_odd <- sf::st_union(df_com_odd, odd_out[[o]])
        }
    }
    df_pol_odd <- sf::st_cast(df_com_odd, "POLYGON")


    even_chr <- list()
    even_chr_seg <- list()

    even_out <- foreach (e=seq(2, 22, 2)) %dopar% {
        P[[e]] <- nrow(df %>% filter(CHR == e))

        if (P[[e]] > 10000) {
            if ((P[[e]] %% 10000) != 0) {
                st[[e]] <- seq(1, P[[e]], 10000)
                en[[e]] <- c(seq(10000, P[[e]], 10000), P[[e]])
            } else {
                st[[e]] <- seq(1, P[[e]], 10000)
                en[[e]] <- seq(10000, P[[e]], 10000)
            }

            even_chr[[e]] <- df %>%
                filter(CHR == e) %>%
                slice(st[[e]][1]:en[[e]][1]) %>%
                sf::st_as_sf(coords = c(2, 3)) %>%
                sf::st_buffer(dist = buffer_dist) %>%
                sf::st_union()
            if (length(st[[e]]) > 1) {
                for (s in 2:length(st[[e]])) {
                    even_chr_seg[[e]] <- df %>%
                        filter(CHR == e) %>%
                        slice(st[[e]][s]:en[[e]][s]) %>%
                        sf::st_as_sf(coords = c(2, 3)) %>%
                        sf::st_buffer(dist = buffer_dist) %>%
                        sf::st_union()
                    even_chr[[e]] <- sf::st_union(even_chr[[e]], even_chr_seg[[e]])
                }
            }
        } else {
            even_chr[[e]] <- df %>%
                filter(CHR %in% e) %>%
                sf::st_as_sf(coords = c(2, 3)) %>%
                sf::st_buffer(dist = buffer_dist) %>%
                sf::st_union()
        }
        return(even_chr[[e]])
    } 
 
    for (e in 1:length(even_out)) {
        if (e == 1) {
            df_com_even <- even_out[[e]]
        } else {
            df_com_even <- sf::st_union(df_com_even, even_out[[e]])
        }
    }
    df_pol_even <- sf::st_cast(df_com_even, "POLYGON")

    
    manh <- ggplot() +
        geom_sf(data = df_pol_odd,  fill = main_cols[1], colour = NA) +
        geom_sf(data = df_pol_even, fill = main_cols[2], colour = NA) +
        scale_x_continuous(label = axis_pos$CHR, breaks = axis_pos$centre,
                           expand = expansion(mult = c(0.01, 0.01))) +  # Remove some white space space at the edges
        scale_y_continuous(expand = expansion(mult = c(0.001, 0.05))) +
        labs(x = "Chromosome", y = expression("-log"[10]*"(p-value)")) +
        ## Use geom_line instead of geom_vline because the latter produces warnings:
        ## "Warning message: range backtransformation not implemented in this coord; results may be wrong."
        ## See https://github.com/tidyverse/ggplot2/issues/2820
        geom_line(data = data.frame(x = c(0 - 0.01 * df$POScum[nrow(df)],
                                          df$POScum[nrow(df)] + 0.01 * df$POScum[nrow(df)]),
                                    y = -log10(5e-8)),
            aes(x, y), size = 0.4, linetype = sig_type, colour = sig_col) +
        coord_sf(xlim = c(0, df$POScum[nrow(df)])) +
        theme_bw()

    if (nchar(title) > 0) manh <- manh + ggtitle(title)  # Add title if available

    manh <- manh + 
        theme(plot.title = element_text(size = 8),
              axis.title = element_text(size = 7),
              axis.text  = element_text(size = 6),
              legend.position = "none",
              panel.grid.minor = element_blank(),
              panel.grid.major = element_blank())
}
