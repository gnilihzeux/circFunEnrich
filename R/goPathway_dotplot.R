#' biological pathway dotplot
#' @export pathDotplot
#'
#' @return NULL
#'
#' @import ggplot2
#' @import stringr
#'
pathDotplot <- function(x, Num = 10){

  df_full <- head(x, Num)
  gene_ratio <- as.numeric(sapply(strsplit(df_full$geneRatio, "/"), "[", 1)) /
    as.numeric(sapply(strsplit(df_full$geneRatio, "/"), "[", 2))
  df_dt   <- data.frame(geneRatio = gene_ratio,
                        Description= df_full$Description,
                        Count      = df_full$overlapGeneCount,
                        logP       = - log10(df_full$pvalue),
                        stringsAsFactors = FALSE)
  df_dt   <- df_dt[order(df_dt$Count, df_dt$logP), ]
  df_dt$Description <- factor(df_dt$Description, levels = df_dt$Description)

  plt <- ggplot(data = df_dt,
                aes(x = Description, y = geneRatio, color = logP, size = Count)) +
    geom_point() +
    theme_classic() +
    scale_color_gradient(low = "red", high = "blue",
                         guide = guide_colorbar(order =1)) +
    labs(x = "", y = "Gene Ratio" ) +
    theme(axis.text.x = element_text(size=8),
          legend.key.size = unit(0.5, "cm"),
          legend.text = element_text(size = 5)) +
    scale_size(range = c(1, 5)) +
    coord_flip() +
    scale_x_discrete(labels = function(x) str_wrap(x, width = 48))

  # @issue legend
  plt$labels$colour <- "- logP"
  plt
}

