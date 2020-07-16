load("data/MSigDB_partial.rda")
#' GSEA main function
#' @export gseaPerform
#'
#' @return NULL
#' @import fgsea
#' @import dplyr

gseaPerform <- function(m_list, rank_vector, pValue = 0.05, FDR = 1) {
  # fgsea
  fg_res <- fgsea::fgsea(pathways = m_list,
                         stats    = rank_vector,
                         minSize  = 5,
                         maxSize  = 500,
                         nperm    = 10000)
  fg_sig <- fg_res[(fg_res$pval <= pValue) & (fg_res$padj <= FDR), ]

  # collpase
  fg_up  <- fgsea::collapsePathways(fg_sig[fg_sig$NES > 0, ],
                                    pathways = m_list,
                                    stats    = rank_vector,
                                    nperm    = 500,
                                    pval.threshold = 0.05)
  fg_dn  <- fgsea::collapsePathways(fg_sig[fg_sig$NES < 0, ],
                                    pathways = m_list,
                                    stats    = rank_vector,
                                    nperm    = 500,
                                    pval.threshold = 0.05)
  fg_tbl <- fg_res[fg_res$pathway %in% c(fg_up$mainPathways, fg_dn$mainPathways), ] %>%
    dplyr::arrange(dplyr::desc(NES))
  fg_tbl$Enrichment <- ifelse(fg_tbl$NES > 0, "Activated", "Suppressed")

  fg_tbl
}

#' pathway name modification
#' @export pathNameModify
#'
#' @return NULL
#'

pathNameModify <- function(fg_tbl) {
  paths <- as.character(fg_tbl$pathway)
  fg_tbl$pathway <- paste0(sub("_.*", "", paths),
                           rep(" | ", length(paths)),
                           tolower(gsub("_", " ", sub(".*?_", "", paths))))
  fg_tbl
}

lollipop_plot <- function(fg_tbl) {
  fg_tbl <- pathNameModify(fg_tbl)
  library(ggplot2)
  ac_num <- min(10, sum(fg_tbl$Enrichment == "Activated"))
  sp_num <- min(10, sum(fg_tbl$Enrichment == "Suppressed"))
  fg_plot <- rbind(head(fg_tbl, ac_num),
                   tail(fg_tbl, sp_num))
  fg_plot$pathway <- factor(fg_plot$pathway, levels = fg_plot$pathway)
  ggplot(fg_plot, aes(x = reorder(pathway, NES), y = NES, label = round(NES, 2))) +
    geom_segment(aes(reorder(pathway, NES), xend = pathway, y = 0, yend = NES),
                 color = "gray75", size = 2) +
    geom_point(aes(color = Enrichment), size = 7) +
    scale_color_manual(values = c("Suppressed" = "#377EB8",
                                  "Activated" = "#E41A1C") ) +
    coord_flip() +
    geom_hline(yintercept = 0, linetype = 2, color = "gray75") +
    labs(x="GeneSets", y="Normalized Enrichment Score") +
    theme_classic() +
    theme(legend.position = "top",
          legend.title = element_blank(),
          legend.justification = "center",
          axis.title = element_text(size = 15),
          axis.text.y = element_text(size=10,
                                     colour= c(rep("#377EB8", sp_num),
                                               rep("#E41A1C", ac_num)))) +
    geom_text(color="white", size=2) +
    scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 48))
}

#' gsea plot
#' @export gseaPlot
#'
#' @return NULL
#'

gseaPlot <- function(m_list, rank_vector, fg_tbl, rowNum = 1) {
  if (rowNum <= nrow(fg_tbl)) {
    pathways <- fg_tbl$pathway
    fg_tbl   <- pathNameModify(fg_tbl)
    path_table <- data.frame(new = fg_tbl$pathway,
                             raw = pathways,
                             stringsAsFactors = FALSE)

    plotEnrichment(pathway = m_list[[path_table[rowNum, "raw"]]],
                   stats   = rank_vector) +
      labs(title = path_table[rowNum, "new"])
  } else {
    message("Please input true Row-Number in the GSEA result table.")
  }
}

