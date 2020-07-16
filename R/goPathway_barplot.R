#' biological pathway barplot
#' @export pathBarplot
#'
#' @return NULL
#'
#' @import ggplot2
#' @import stringr

pathBarplot <- function(x, Num = 10){

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

  go_plot <- ggplot(data = df_dt, aes(x = Description, y = Count, fill = logP)) +
    geom_bar(stat = "identity") +
    theme_classic() +
    scale_fill_gradient(low = "red", high = "blue") +
    labs(x = "", y="Count")+
    theme(axis.text.y = element_text(size= 8),
          axis.title.x = element_text(size = 8),
          legend.key.size = unit(0.5, "cm"),
          legend.text = element_text( size = 5)) +
    coord_flip() +
    scale_x_discrete(labels = function(x) str_wrap(x, width = 48))
  go_plot$labels$fill <- "- logP"
  go_plot
}

