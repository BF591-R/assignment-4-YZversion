library('tidyverse')
library('RColorBrewer')

#' Read the expression data "csv" file as a dataframe, not tibble
#'
#' @param filename (str): the path of the file to read
#' @param delimiter (str): generalize the function so it can read in data with
#'   your choice of delimiter
#'
#' @return A dataframe containing the example intensity data with rows as probes
#'   and columns as samples
#' @export
#'
#' @examples
read_data <- function(intensity_data, delimiter) {
    intensity_df <- utils::read.table(
        file = intensity_data,
        sep = delimiter,
        header = TRUE,
        row.names = 1,
        check.names = FALSE,
        stringsAsFactors = FALSE
    )

    return(as.data.frame(intensity_df))
}

#' Define a function to calculate the proportion of variance explained by each PC
#'
#' @param pca_results (obj): the results returned by `prcomp()`
#'
#' @return A vector containing the values of the variance explained by each PC
#' @export
#'
#' @examples
calculate_variance_explained <- function(pca_results) {
    stopifnot(!is.null(pca_results$sdev))

    variances <- pca_results$sdev^2
    variance_explained <- variances / sum(variances)

    return(variance_explained)
}

#' Define a function that takes in the variance values and the PCA results to
#' make a tibble with PC names, variance explained by each PC, and the
#' cumulative sum of variance explained. These columns should be named 
#' "principal_components", "variance_explained", and "cumulative", respectively.
#' 
#'
#'
#' @param pca_ve (vector): the vector generated in the previous function with
#'   the variance explained values
#' @param pca_results (object): the results returned by `prcomp()`
#' @return A tibble that contains the names of the PCs, the individual variance
#'   explained, and the cumulative variance explained with names described above
#' @export
#' @examples 
make_variance_tibble <- function(pca_ve, pca_results) {
    pc_names <- colnames(pca_results$x)
    if (is.null(pc_names) || length(pc_names) < length(pca_ve)) {
        pc_names <- paste0("PC", seq_along(pca_ve))
    } else {
        pc_names <- pc_names[seq_along(pca_ve)]
    }

    variance_tibble <- tibble(
        principal_components = factor(pc_names, levels = pc_names),
        variance_explained = pca_ve,
        cumulative = cumsum(pca_ve)
    )

    return(variance_tibble)
}



#' Define a function to create a biplot of PC1 vs. PC2 labeled by
#' SixSubTypesClassification
#'
#' @param metadata (str): The path to the proj_metadata.csv file
#' @param pca_results (obj): The results returned by `prcomp()`
#'
#' @return A ggplot consisting of a scatter plot of PC1 vs PC2 labeled by
#'   SixSubTypesClassification found in the metadata
#' @export
#'
#' @examples
make_biplot <- function(metadata, pca_results) {
    metadata_tbl <- readr::read_csv(metadata, show_col_types = FALSE)

    scores <- as.data.frame(pca_results$x)
    scores$geo_accession <- rownames(scores)

    merged_tbl <- dplyr::left_join(
        scores,
        metadata_tbl %>% dplyr::select(geo_accession, SixSubtypesClassification),
        by = "geo_accession"
    ) %>%
        dplyr::filter(!is.na(SixSubtypesClassification))

    pc_cols <- intersect(c("PC1", "PC2"), colnames(scores))
    if (length(pc_cols) < 2) {
        pc_cols <- colnames(scores)[seq_len(min(2, ncol(scores)))]
    }

    merged_tbl$pc_x <- merged_tbl[[pc_cols[1]]]
    merged_tbl$pc_y <- merged_tbl[[pc_cols[2]]]

    biplot <- ggplot(
        merged_tbl,
        aes(x = pc_x, y = pc_y, color = SixSubtypesClassification)
    ) +
        geom_point(size = 3, alpha = 0.85) +
        labs(
            x = pc_cols[1],
            y = pc_cols[2],
            color = "Subtype",
            title = "PC1 vs PC2 Biplot"
        ) +
        theme_minimal()

    return(biplot)
}

#' Define a function to return a list of probeids filtered by signifiance
#'
#' @param diff_exp_tibble (tibble): A tibble containing the differential expression results
#' @param fdr_threshold (float): an appropriate FDR threshold, we will use a
#'   value of .01. This is the column "padj" in the tibble.
#'
#' @return A list with the names of the probeids passing the fdr_threshold
#' @export
#'
#' @examples
list_significant_probes <- function(diff_exp_tibble, fdr_threshold) {
    significant_ids <- diff_exp_tibble %>%
        dplyr::filter(!is.na(padj), padj < fdr_threshold) %>%
        dplyr::pull(probeid)

    return(significant_ids)
}

#' Define a function that uses the list of significant probeids to return a
#' matrix with the intensity values for only those probeids.
#' @param intensity (dataframe): The dataframe of intensity data generated in
#'   part 1
#' @param sig_ids_list (list/vector): The list of differentially expressed
#'   probes generated in part 6
#'
#' @return A `matrix()` of the probe intensities for probes in the list of
#'   significant probes by FDR determined in the previous function.
#'
#' @export
#'
#' @examples
return_de_intensity <- function(intensity, sig_ids_list) {
    intensity_df <- as.data.frame(intensity)
    available_ids <- sig_ids_list[sig_ids_list %in% rownames(intensity_df)]

    de_subset <- intensity_df[available_ids, , drop = FALSE]

    return(as.matrix(de_subset))
}

#' Define a function that takes the intensity values for significant probes and
#' creates a color-blind friendly heatmap
#'
#' @param de_intensity (matrix): The matrix of intensity values for significant
#'   differentially expressed probes returned in part 7
#' @param num_colors (int): The number of colors in a specificed RColorBrewer
#'   palette
#' @param palette (str): The name of the chosen RColorBrewer palette
#'
#' @return A heatmap displaying the intensity values for the differentially
#'   expressed probes
#' @export
#'
#' @examples
plot_heatmap <- function(de_intensity, num_colors, palette) {
    intensity_matrix <- as.matrix(de_intensity)

    if (!palette %in% rownames(RColorBrewer::brewer.pal.info)) {
        stop("Palette not found in RColorBrewer.")
    }

    max_colors <- RColorBrewer::brewer.pal.info[palette, "maxcolors"]
    base_colors <- RColorBrewer::brewer.pal(max_colors, palette)
    palette_fn <- grDevices::colorRampPalette(base_colors)

    heatmap(
        intensity_matrix,
        col = palette_fn(num_colors),
        scale = "row"
    )
}

