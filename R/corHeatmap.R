#' plot heatmap with normalized expression
#' @export corHeatmap
#'
#' @return NULL
#'
#' @import dplyr
#' @import RColorBrewer
#' @import ComplexHeatmap
#' @import circlize
#' @import grid

corHeatmap <- function(cor_table, circ_matrix, gene_matrix, posNum = 50, negNum = 50, needLog = TRUE) {

  # gene filter
  pos_cor <- cor_table %>%
    dplyr::filter(cor > 0) %>%
    dplyr::arrange(desc(cor)) %>%
    head(posNum)

  neg_cor <- cor_table %>%
    dplyr::filter(cor < 0) %>%
    dplyr::arrange(desc(cor)) %>%
    tail(negNum)
  # gene matrix
  gn_data <- gene_matrix[c(pos_cor$target, neg_cor$target), ]
  if (needLog)
    gn_data <- log2(gn_data + 0.01)
  gn_scale <- scale(t(as.matrix(gn_data)))
  # circ matrix
  circ_data <- circ_matrix
  CIRC_NAME <- rownames(circ_data)
  if (needLog)
    circ_data <- log2(circ_data + 0.01)
  circ_scale <- scale(t(as.matrix(circ_data)))
  circ_order <- circ_scale[, 1]
  names(circ_order) <- rownames(circ_scale)
  circ_order <- sort(circ_order)
  circ_df    <- data.frame(circ_order)
  colnames(circ_df) <- CIRC_NAME

  # heatmap -------------------------
  # color
  # `annotation_label` only for `ComplexHeatmap > v 2.3.3`
  CIRC_LIM <- 0.5*(floor(max(abs(circ_order)) / 0.5) + 1)
  col_fun <- circlize::colorRamp2(c(- CIRC_LIM, 0, CIRC_LIM),
                                  c("#377EB8", "white", "#E41A1C"))
  circ_col <- list(col_fun)
  names(circ_col) <- CIRC_NAME
  heat_circ <- ComplexHeatmap::HeatmapAnnotation(df = circ_df,
                                                 col = circ_col,
                                                 annotation_legend_param = list(direction = "horizontal"),
                                                 annotation_name_side = "right",
                                                 annotation_name_gp = gpar(fontsize = 7),
                                                 border = TRUE)

  heat_all <- ComplexHeatmap::Heatmap(t(gn_scale)[, rownames(circ_df)],
                                      name = "gene",
                                      col = col_fun,
                                      top_annotation = heat_circ,
                                      column_order = rownames(circ_df),
                                      show_column_names = FALSE,
                                      show_row_names = TRUE,
                                      row_names_gp = gpar(fontsize = 4),
                                      heatmap_legend_param = list(direction = "horizontal")
                                      # row_title = NULL,
                                      # row_split = c(rep("A", posNum), rep("B", negNum)),
                                      # column_dend_height = unit(2, "cm"),row_dend_width = unit(2, "cm")
  )
  draw(heat_all,
       merge_legend = TRUE,
       heatmap_legend_side = "bottom",
       annotation_legend_side = "bottom")

}

