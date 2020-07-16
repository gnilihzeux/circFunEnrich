#' shinyAppServer
#' @export shinyAppServer
#'
#' gene_expression
#'
#' @return NULL
#' @import shiny
#' @import data.table
#' @import bsplus
#' @import openxlsx
#' @import clusterProfiler

circbank <- readRDS("data/circbank.rds")

shinyAppServer <- function(input, output, session) {

  # ----------------------------------------
  # circBank Annotation                     \
  # -----------------------------------------
  # circRNA coordinates input
  circbank_filter <- reactive({
    if (input$circ_coord == "")
      return(NULL)
    itms <- unique(unlist(strsplit(unlist(strsplit(input$circ_coord, "\n")), " ")))
    index <- circbank$position %in% itms
    if (!any(index)) {
      return("<strong>Sorry, no items matched.<br>Please check your format.</strong>")
    } else {
      return(merge(circbank, data.frame(position = itms), all.y = TRUE))
    }
  })


  ct_value <- reactiveValues(data = NULL)

  observeEvent(input$coord_button, {
    ct_value$data <- circbank_filter()
  })

  # render ---------------------------
  output$circbk_txt <- renderText({
    if (!is.character(ct_value$data))
      return(NULL)
    ct_value$data
  })
  output$circbk_tbl <- DT::renderDataTable({
    if (!is.character(ct_value$data))
      DT::datatable(ct_value$data,
                    extensions = "Buttons",
                    options = list(paging      = TRUE,
                                   searching   = TRUE,
                                   fixedColums = TRUE,
                                   autoWidth   = TRUE,
                                   ordering    = TRUE)
      )
  })
  # download
  output$circbk_download <- renderUI({
    req(ct_value$data)
    downloadButton('circbk_download_item',
                   label = "Download Output")
  })
  output$circbk_download_item <- downloadHandler(
    filename = function(){
      "circBank_annotation_table.csv"
    },
    content = function(file) {
      write.csv(ct_value$data, file, row.names = FALSE)
    }
  )

  # ----------------------------------------
  # CFEA                                    \
  # -----------------------------------------
  # main ---------------------
  cfea_react <- reactiveValues()
  observeEvent(input$cfea_button,
               {
                 input$circrna_matrix
                 input$gene_matrix
                 # circrna matrix
                 circrna_matrix <- data.frame(fread(input$circrna_matrix$datapath))
                 exp_m1      <- circrna_matrix[- 1]
                 rownames(exp_m1) <- circrna_matrix[, 1]
                 cfea_react$exp_m1 <- exp_m1

                 # gene matrix
                 gene_matrix <- data.frame(fread(input$gene_matrix$datapath))
                 exp_m2      <- gene_matrix[- 1]
                 rownames(exp_m2) <- gene_matrix[, 1]
                 cfea_react$exp_m2 <- exp_m2
                 # cor
                 withProgress(message = 'correlation calculating ...', value = 30, {
                   cfea_react$coexp_reslt <- corMatrix2LongTable(pearsonTestMatrix(t(cfea_react$exp_m1),
                                                                                   t(cfea_react$exp_m2)),
                                                                 symmetric = FALSE)
                 })
                 cfea_react$coexp_result <- cfea_react$coexp_reslt %>%
                   dplyr::arrange(desc(cor))
                 # overall panel ---------------
                 index <- (cfea_react$coexp_result$cor >= input$coexp_cor) & (cfea_react$coexp_result$pvalue <= input$coexp_p) & (cfea_react$coexp_result$padj <= input$coexp_padj)
                 if (any(index)) {
                   cfea_react$coexp_filter <- cfea_react$coexp_result[index, ]
                 } else {
                   cfea_react$sorry <- "<strong>Sorry, No Results.Please set lower thresholds."
                 }
                 # CFEA
                 withProgress(message = 'function enrichment process ...', value = 90, {
                   cfea_react$cfea_enrich <- goPathwayEnrichment(cfea_react$coexp_filter$target,
                                                                 species = "human",
                                                                 pvalueCutoff = 1,
                                                                 qvalueCutoff = 1,
                                                                 minOverlapNum = 3)
                 }
                 )

               })
  # items
  observeEvent(input$cfea_items,
               {
                 cfea_react$cfea_which_item <- switch(input$cfea_items,
                                                      go_BP = cfea_react$cfea_enrich$format_bp,
                                                      go_CC = cfea_react$cfea_enrich$format_cc,
                                                      go_MF = cfea_react$cfea_enrich$format_mf,
                                                      KEGG  = cfea_react$cfea_enrich$format_kegg,
                                                      Reactome = cfea_react$cfea_enrich$format_reactome
                 )
               })

  # ui ----------------------
  # format
  output$cfea_format <- downloadHandler(
    filename = function(){
      "file_format.csv"
    },
    content = function(file) {
      file_format <- read.csv(system.file("extdata",
                                          "gene_expression.csv",
                                          package = "circFunEnrich.beta1.1"))
      write.csv(file_format, file, row.names = FALSE)
    }
  )
  # coexpression
  output$coexp_tbl <- DT::renderDataTable({
    DT::datatable(cfea_react$coexp_filter,
                  extensions = "Buttons",
                  options = list(paging      = TRUE,
                                 lengthChange = TRUE,
                                 pageLength  = 6,
                                 digits      = 3,
                                 searching   = TRUE,
                                 fixedColums = TRUE,
                                 autoWidth   = TRUE,
                                 ordering    = TRUE
                  )
    )
  }
  )

  output$coexp_download <- renderUI({
    req(cfea_react$coexp_filter)
    downloadButton('coexp_download_item',
                   label = "Download Output")
  })
  output$coexp_download_item <- downloadHandler(
    filename = function(){
      "coexpressiton_correlation.csv"
    },
    content = function(file) {
      write.csv(cfea_react$coexp_filter, file, row.names = FALSE)
    }
  )
  # cfea
  output$coexp_cfea_items <- renderUI({
    req(cfea_react$cfea_enrich)
    selectInput("cfea_items",
                label = "Canical Function Sets",
                choices = c("GO: Biological Pathway" = "go_BP",
                            "GO: Cellular Component" = "go_CC",
                            "GO: Molecular Function" = "go_MF",
                            "KEGG" = "KEGG",
                            "Reactome" = "Reactome"),
                width = "250px"
    )
  })

  output$cfea_tbl <- renderUI({
    req(is.data.frame(cfea_react$cfea_which_item))
    DT::dataTableOutput("cfea_tbl_item")
  })
  output$cfea_tbl_item <- DT::renderDataTable({
    req(is.data.frame(cfea_react$cfea_which_item))
    DT::datatable(cfea_react$cfea_which_item,
                  extensions = "Buttons",
                  options = list(paging      = TRUE,
                                 scrollX     = TRUE,
                                 lengthChange = TRUE,
                                 pageLength  = 5,
                                 digits      = 3,
                                 searching   = TRUE,
                                 fixedColums = TRUE,
                                 autoWidth   = TRUE,
                                 ordering    = TRUE
                  )
    )
  })
  output$cfea_download <- renderUI({
    req(is.data.frame(cfea_react$cfea_which_item))
    downloadButton('cfea_download_item',
                   label = "Download Output")
  })
  output$cfea_download_item <- downloadHandler(
    filename = function(){
      "CFEA_results.csv"
    },
    content = function(file) {
      req(is.data.frame(cfea_react$cfea_which_item))
      write.csv(cfea_react$cfea_which_item, file, row.names = FALSE)
    }
  )
  output$coexp_sorry <- renderText({
    req(is.character(cfea_react$cfea_which_item))
    cfea_react$sorry
  })
  # dotplot
  output$cfea_width <- renderUI({
    req(is.data.frame(cfea_react$cfea_which_item))
    numericInput("cfea_width_num",
                 label = "Width(px)",
                 value = 850,
                 min   = 400,
                 max   = 1500,
                 step = 1,
                 width = "150px")
  })
  output$cfea_height <- renderUI({
    req(is.data.frame(cfea_react$cfea_which_item))
    numericInput("cfea_height_num",
                 value = 450,
                 min   = 400,
                 max   = 1500,
                 step = 1,
                 label = "Height(px)",
                 width = "150px")
  })
  output$cfea_plot <- renderUI({
    req(is.data.frame(cfea_react$cfea_which_item))
    numericInput("cfea_plot_num",
                 value = 10,
                 min   = 1,
                 max   = nrow(cfea_react$cfea_which_item),
                 step = 1,
                 label = "No.of TOP",
                 width = "150px")
  })

  output$cfea_dotplot <- renderUI({
    req(is.data.frame(cfea_react$cfea_which_item))
    plotOutput("dotplot",
               width = input$cfea_width_num,
               height = input$cfea_height_num)
  })
  output$dotplot <- renderPlot({
    pathDotplot(cfea_react$cfea_which_item, input$cfea_plot_num)
  }, width = "auto", height = "auto", res = 150)

  output$cfea_barplot <- renderUI({
    req(is.data.frame(cfea_react$cfea_which_item))
    plotOutput("barplot",
               width = input$cfea_width_num,
               height = input$cfea_height_num)
  })
  output$barplot <- renderPlot({
    pathBarplot(cfea_react$cfea_which_item, input$cfea_plot_num)
  }, width = "auto", height = "auto", res = 150)

  # ----------------------------------------
  # GSEA                                    \
  # -----------------------------------------

  # main -------------
  gsea_react <- reactiveValues()

  observeEvent(input$gsea_button,
               {
                 input$circrna_matrix_gsea
                 input$gene_matrix_gsea
                 # circrna matrix
                 circrna_matrix <- data.frame(fread(input$circrna_matrix_gsea$datapath))
                 exp_m1      <- circrna_matrix[- 1]
                 rownames(exp_m1) <- circrna_matrix[, 1]
                 exp_m1 <- exp_m1[1, ]
                 gsea_react$exp_m1 <- exp_m1

                 # gene matrix
                 gene_matrix <- data.frame(fread(input$gene_matrix_gsea$datapath))
                 exp_m2      <- gene_matrix[- 1]
                 rownames(exp_m2) <- gene_matrix[, 1]
                 gsea_react$exp_m2 <- exp_m2
                 # cor
                 gsea_react$coexp_reslt <- corMatrix2LongTable(pearsonTestMatrix(t(gsea_react$exp_m1),
                                                                                 t(gsea_react$exp_m2)),
                                                               symmetric = FALSE)
                 gsea_react$coexp_result <- gsea_react$coexp_reslt %>%
                   dplyr::arrange(desc(cor))

                 # GSEA
                 coexp_rslt <- gsea_react$coexp_result %>%
                   mutate(rank = base::rank(- cor, ties.method = "random"))
                 rank_vector <- coexp_rslt$rank
                 names(rank_vector) <- coexp_rslt$target
                 gsea_react$rank_vector <- rank_vector
                 gsea_react$msigdb_list <- rec_list
               })
  observeEvent(input$gsea_p,
               {
                 gsea_react$msigdb_list
                 gsea_react$rank_vector
                 input$gsea_padj
                 gsea_react$gsea_table <- gseaPerform(gsea_react$msigdb_list,
                                                      gsea_react$rank_vector,
                                                      input$gsea_p,
                                                      input$gsea_padj)
                 # change name
                 if (is.data.frame(gsea_react$gsea_table)) {
                   gsea_react$gsea_result <- pathNameModify(gsea_react$gsea_table)
                 } else {
                   gsea_react$sorry <- "Threshold may be too strigency."
                 }
               }
  )

  # ui --------------------
  # format download
  output$gsea_format <- downloadHandler(
    filename = function(){
      "file_format.csv"
    },
    content = function(file) {
      #file_format <- read.csv("data/gene_expression.csv")
      file_format <- read.csv(system.file("extdata",
                                          "gene_expression.csv",
                                          package = "circFunEnrich.beta1.1"))
      write.csv(file_format, file, row.names = FALSE)
    }
  )
  # correlation table
  output$icor_tbl <- renderUI({
    req(is.data.frame(gsea_react$coexp_result))
    DT::dataTableOutput("icor_tbl_item")
  })
  output$icor_tbl_item <- DT::renderDataTable({
    req(is.data.frame(gsea_react$coexp_result))
    DT::datatable(gsea_react$coexp_result,
                  extensions = "Buttons",
                  options = list(paging      = TRUE,
                                 scrollX     = TRUE,
                                 lengthChange = TRUE,
                                 pageLength  = 5,
                                 digits      = 3,
                                 searching   = TRUE,
                                 fixedColums = TRUE,
                                 autoWidth   = TRUE,
                                 ordering    = TRUE
                  )
    )
  })
  output$icor_download <- renderUI({
    req(is.data.frame(gsea_react$coexp_result))
    downloadButton('icor_download_item',
                   label = "Download Output")
  })
  output$icor_download_item <- downloadHandler(
    filename = function(){
      "coexpressiton_correlation.csv"
    },
    content = function(file) {
      req(is.data.frame(gsea_react$coexp_result))
      write.csv(gsea_react$coexp_result, file, row.names = FALSE)
    }
  )
  # heatmap
  output$gsea_heatmap <- renderUI({
    req(is.data.frame(gsea_react$coexp_result))
    plotOutput("cor_heatmap",
               width = "960px",
               height = "825px")
  })
  output$cor_heatmap <- renderPlot({
    corHeatmap(gsea_react$coexp_result,
               gsea_react$exp_m1,
               gsea_react$exp_m2)
  }, res = 150)
  # gsea table
  output$gsea_descript <- renderUI({
    req(gsea_react$rank_vector)
    htmlOutput("gsea_descript_id")
  })
  output$gsea_descript_id <- renderText({
    "<br/>
     <br/>
      <p><strong>GSEA Functional Sets</strong> integrates with
        <i><b>H:hallmark gene sets & C2:pathways without CGP</b></i>
        derived from
        <a href='https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp'>MSigDB</a>.
      </p>
     <p>GSEA cutoffs</p>
    "
  })

  output$gsea_p_desc <- renderUI({
    textOutput("p_desc")
  })
  output$p_desc <- renderText({
    req(gsea_react$rank_vector)
    "p-value <"
  })
  output$gsea_padj_desc <- renderUI({
    textOutput("padj_desc")
  })
  output$padj_desc <- renderText({
    req(gsea_react$rank_vector)
    "p-adj <"
  })

  output$gsea_p_pos <- renderUI({
    req(gsea_react$rank_vector)
    numericInput("gsea_p",
                 label = "",
                 value = 0.05,
                 min   = 0,
                 max   = 1,
                 width = "100px"
    )
  })
  output$gsea_padj_pos <- renderUI({
    req(gsea_react$rank_vector)
    numericInput("gsea_padj",
                 label = "",
                 value = 1,
                 min   = 0,
                 max   = 1,
                 width = "100px"
    )
  })
  output$gsea_tbl <- renderUI({
    req(is.data.frame(gsea_react$gsea_result))
    DT::dataTableOutput("gsea_tbl_item")
  })
  output$gsea_tbl_item <- DT::renderDataTable({
    DT::datatable(gsea_react$gsea_result,
                  extensions = "Buttons",
                  options = list(paging      = TRUE,
                                 scrollX     = TRUE,
                                 lengthChange = TRUE,
                                 pageLength  = 3,
                                 searching   = TRUE,
                                 fixedColums = TRUE,
                                 autoWidth   = TRUE,
                                 ordering    = TRUE,
                                 digits      = 3
                  )
    )
  }
  )
  output$gsea_sorry <- renderText({
    gsea_react$sorry
  })
  # download
  output$gsea_download <- renderUI({
    req(gsea_react$gsea_result)
    downloadButton('gsea_download_item',
                   label = "Download Output")
  })
  output$gsea_download_item <- downloadHandler(
    filename = function(){
      "gsea_result_table.xlsx"
    },
    content = function(file) {
      write.xlsx(gsea_react$gsea_result, file)
    }
  )

  # lolliplot
  output$gsea_lolliplot <- renderUI({
    req(gsea_react$gsea_result)
    plotOutput("lolliplot",
               width = '1050px',
               height = '800px')
  })
  output$lolliplot <- renderPlot({
    #' note: not gsea_react$gsea_result
    lollipop_plot(gsea_react$gsea_table)
  }, width = "auto", height = "auto", res = 150)

  # gseaplot
  output$gseaplot_num <- renderUI({
    req(gsea_react$gsea_result)
    numericInput("gseaplot_rownum",
                 label = "Row number of pathway to plot",
                 value = 1,
                 min   = 1,
                 max   = nrow(gsea_react$gsea_result),
                 step  = 1)
  })

  output$gsea_gseaplot <- renderUI({
    req(gsea_react$gsea_result)
    plotOutput("gseaplot",
               width = '850px',
               height = '750px')
  })
  output$gseaplot <- renderPlot({
    gseaPlot(gsea_react$msigdb_list,
             gsea_react$rank_vector,
             gsea_react$gsea_table,
             rowNum = input$gseaplot_rownum)
  }, res = 125)

}
