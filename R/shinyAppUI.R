#' shinyAppServer
#' @export shinyAppUI
#'
#' @return NULL
#' @import shiny
#' @import shinythemes

shinyAppUI <- navbarPage(
  # title = img(src="logo.png", height = "80px"),
  title = "circFunEnrich",
  id = "navBar",

  #theme = "readable.css",
  theme = shinythemes::shinytheme("cerulean"),
  collapsible = TRUE,
  inverse = FALSE,
  windowTitle = "human circular RNA tookits for functions & Pearson Correlation",
  position = "fixed-top",
  header = tags$style(
    ".navbar-right {
      float: right !important;
    }",
    ".navbar .navbar-header {
      position: absolute;
      left:18%
    }",
    ".navbar .navbar-nav {
      position: absolute;
      left:30%
    }",
    "body {padding-top: 180px;}"),
  footer = includeHTML("./www/footer.html"),

  tabPanel("Home", value = "home",

           includeHTML("./www/home.html")

  ),

  # ----------------------------------------
  # circBank Annotation                     \
  # -----------------------------------------
  tabPanel("circBank Annotation", value = "circ_anno",

           fluidPage(title = "circBank Annotation",
                     shiny::HTML("
                                      <p style='font-size:18px;text-align:center'>
                                        <strong>circBank Annotation</strong>
                                        could annotate circular RNA genomic coordinates
                                        with <a href:http://www.circbank.cn/>circBank</a>
                                        annotation.
                                      </p>
                                                "),

                     fluidRow(
                       column(2, offset = 2,
                              textAreaInput("circ_coord", label = "",
                                            width = "300px",
                                            height = "200px")
                       ),
                       column(4, offset = 2,
                              br(),
                              helpText(
                                "NOTE"
                              ),
                              helpText(
                                "Please input genomic coordinates, i.e. chr11:33307958-33361080"),
                              helpText(
                                "Multiple entries should be one-line-one-entry."),

                              actionButton("coord_button",
                                           label = "Submmit",
                                           icon  = icon("play-circle"))

                       )
                     ),
                     br(),
                     br(),

                     fluidRow(
                       column(4, offset = 4,

                              htmlOutput("circbk_txt")
                       )
                     ),
                     fluidRow(
                       column(12, offset = 0,

                              DT::dataTableOutput("circbk_tbl"))
                     ),
                     fluidRow(
                       column(2, offset = 8,
                              uiOutput('circbk_download')
                       )
                     )
           )
  ),

  # ----------------------------------------
  # CFEA                                    \
  # -----------------------------------------
  tabPanel("CFEA",
           fluidRow(
             column(8, offset = 2,
                    includeMarkdown("www/CFEA_main.md")
             )
           ),

           # input ----------------------
           fluidRow(br(),


                    column(8, offset = 2,

                           shiny::HTML("
                                          <div style='color:#28A1E6'>
                                            <p><strong>INPUT files:</strong></p>
                                            <p style='color:black'>&emsp;&emsp;Regular Expression Matrix,
                                            Rows are genes(mRNA/lncRNA/circRNA/...),</p>
                                            <p style='color:black'>&emsp;&emsp;Colums are samples,
                                            and Values are <strong>normalized expression values</strong>.</p>
                                          </div>
                                          "),
                           shiny::HTML("
                                          <div style='color:#28A1E6'>
                                            <p><strong>Allowed Input Format</strong></p>
                                            <p style='color:black'>&emsp;&emsp;txt / csv</p>
                                          </div>
                                          "),
                           shiny::HTML("
                                          <div style='color:red'>
                                            <p><strong>NOTE!</strong></p>
                                            <p style='color:black'>&emsp;&emsp;The first column is the gene name by default.</p>
                                          </div>
                                          ")

                    )
           ),
           fluidRow(
             column(2, offset = 2,
                    downloadButton("cfea_format",
                                   label = "File Format",
                                   icon  = "file-csv")
             ),
             column(8, offset = 2,
                    hr()
             )
           ),

           fluidRow(
             br(),

             column(3, offset = 2,

                    fileInput("circrna_matrix", label = "circRNA Expression Matrix",
                              multiple = FALSE,
                              buttonLabel = "browser...",
                              accept = c(
                                "text/txt",
                                "text/tab-separated-values,text/plain",
                                ".txt",
                                "text/csv",
                                "text/comma-separated-values,text/plain",
                                ".csv"),
                              placeholder = "No file selected",
                              width = "100%")
             )

           ),
           fluidRow(
             br(),

             column(3, offset = 2,

                    fileInput("gene_matrix", label = "Gene Expression Matrix",
                              multiple = FALSE,
                              buttonLabel = "browser...",
                              accept = c(
                                "text/txt",
                                "text/tab-separated-values,text/plain",
                                ".txt",
                                "text/csv",
                                "text/comma-separated-values,text/plain",
                                ".csv"),
                              placeholder = "No file selected",
                              width = "100%")
             )
           ),
           fluidRow(
             br(),
             column(6, offset = 2,
                    helpText("Make sure that the columns of expression matrix are consistent bewteen circRNA and gene.")
             )
           ),

           # result filter
           fluidRow(

             br(),
             column(4, offset = 2,
                    p(strong("Co-expression Filter"))
             ),
             column(4, offset = 0,
                    helpText("Correlation Coefficient must > 0."))
           ),

           fluidRow(
             column(3, offset = 2, style='padding-top:30px;padding-right:10px',
                    p("correlation coefficient >")
             ),
             column(1, offset = 0, style = 'padding-left:0px;margin-left:-4em',
                    numericInput("coexp_cor",
                                 label = "",
                                 value = 0.9,
                                 min   = 0,
                                 max   = 1,
                                 width = "100px"
                    )
             )
           ),

           fluidRow(
             column(3, offset = 2, style='padding-top:30px;padding-right:10px',
                    p("p-value <")
             ),
             column(1, offset = 0, style = 'padding-left:0px;margin-left:-4em',
                    numericInput("coexp_p",
                                 label = "",
                                 value = 0.05,
                                 min   = 0,
                                 max   = 1,
                                 width = "100px"
                    )
             )
           ),

           fluidRow(
             column(3, offset = 2, style='padding-top:30px;padding-right:10px',
                    p("p-adj <")
             ),
             column(1, offset = 0, style = 'padding-left:0px;margin-left:-4em',
                    numericInput("coexp_padj",
                                 label = "",
                                 value = 1,
                                 min   = 0,
                                 max   = 1,
                                 width = "100px"
                    )
             )
           ),

           # submit
           fluidRow(
             column(4, offset = 2,
                    actionButton("cfea_button",
                                 label = "Submmit",
                                 icon  = icon("play-circle"))
             )
           ),

           fluidRow(
             column(12, offset = 0,

                    DT::dataTableOutput("coexp_tbl"))
           ),

           fluidRow(
             column(7, offset = 3,
                    br(),
                    htmlOutput("coexp_sorry")
             )
           ),

           fluidRow(
             column(2, offset = 8,
                    uiOutput('coexp_download')
             )
           ),
           # GO / KEGG / Reactome
           fluidRow(
             column(4, offset = 0,
                    uiOutput('coexp_cfea_items')
             )
           ),
           fluidRow(
             column(12, offset = 0,
                    uiOutput('cfea_tbl'))
           ),
           fluidRow(
             column(2, offset = 8,
                    uiOutput('cfea_download'))
           ),
           # dotplot
           fluidRow(
             column(2, offset = 2,
                    br(),
                    uiOutput('cfea_width')),
             column(2, offset = 0,
                    br(),
                    uiOutput('cfea_height')),
             column(2, offset = 0,
                    br(),
                    uiOutput('cfea_plot'))

           ),

           fluidRow(
             column(8, offset = 0,
                    br(),
                    br(),

                    uiOutput('cfea_dotplot')
             )
           ),
           fluidRow(
             column(8, offset = 0,
                    br(),
                    br(),

                    uiOutput('cfea_barplot')
             )
           )


  ),

  # ----------------------------------------
  # GSEA                                    \
  # -----------------------------------------
  tabPanel("GSEA",

           fluidRow(
             column(8, offset = 2,
                    includeMarkdown("www/GSEA_main.md")
             )
           ),

           # input ----------------------
           fluidRow(br(),


                    column(8, offset = 2,

                           shiny::HTML("
                                          <div style='color:#28A1E6'>
                                            <p><strong>INPUT files:</strong></p>
                                            <p style='color:black'>&emsp;&emsp;Regular Expression Matrix,
                                            Rows are genes(mRNA/lncRNA/circRNA/...),</p>
                                            <p style='color:black'>&emsp;&emsp;Colums are samples,
                                            and Values are <strong>normalized expression values</strong>.</p>
                                          </div>
                                          "),
                           shiny::HTML("
                                          <div style='color:#28A1E6'>
                                            <p><strong>Allowed Input Format</strong></p>
                                            <p style='color:black'>&emsp;&emsp;txt / csv</p>
                                          </div>
                                          "),
                           shiny::HTML("
                                          <div style='color:red'>
                                            <p><strong>NOTE!</strong></p>
                                            <p style='color:black'>&emsp;&emsp;The first column is the gene name by default.</p>
                                          </div>
                                          ")

                    )
           ),
           fluidRow(
             column(2, offset = 2,
                    downloadButton("gsea_format",
                                   label = "File Format",
                                   icon  = "file-csv")
             ),
             column(8, offset = 2,
                    hr()
             )
           ),

           fluidRow(
             br(),

             column(3, offset = 2,

                    fileInput("circrna_matrix_gsea", label = "circRNA Expression Matrix",
                              multiple = FALSE,
                              buttonLabel = "browser...",
                              accept = c(
                                "text/txt",
                                "text/tab-separated-values,text/plain",
                                ".txt",
                                "text/csv",
                                "text/comma-separated-values,text/plain",
                                ".csv"),
                              placeholder = "No file selected",
                              width = "100%")
             )

           ),
           fluidRow(
             br(),

             column(3, offset = 2,

                    fileInput("gene_matrix_gsea", label = "Gene Expression Matrix",
                              multiple = FALSE,
                              buttonLabel = "browser...",
                              accept = c(
                                "text/txt",
                                "text/tab-separated-values,text/plain",
                                ".txt",
                                "text/csv",
                                "text/comma-separated-values,text/plain",
                                ".csv"),
                              placeholder = "No file selected",
                              width = "100%")
             )
           ),
           fluidRow(
             br(),
             column(6, offset = 2,
                    helpText("Make sure that the columns of expression matrix are consistent bewteen circRNA and gene."),
                    helpText("If your circRNA matrix has multiple rows, only the first will been as input.")
             )
           ),
           # submit button ------------------------
           fluidRow(
             column(4, offset = 2,
                    br(),
                    br(),
                    actionButton("gsea_button",
                                 label = "Submmit",
                                 icon  = icon("play-circle"))
             )
           ),
           # sorry ------------------
           fluidRow(
             column(7, offset = 3,
                    br(),
                    htmlOutput("gsea_sorry")
             )
           ),

           # correlation block ---------
           fluidRow(
             column(12, offset = 0,

                    uiOutput("icor_tbl"))
           ),
           fluidRow(
             column(2, offset = 8,
                    uiOutput('icor_download'))
           ),
           fluidRow(
             column(8, offset = 2,
                    uiOutput('gsea_heatmap')
             )
           ),

           # gsea block --------------
           fluidRow(
             column(8, offset = 2,
                    htmlOutput('gsea_descript')
             ),
             column(3, offset = 2, style='padding-top:30px;padding-right:10px',
                    uiOutput('gsea_p_desc')
             ),
             column(1, offset = 0, style = 'padding-left:0px;margin-left:-4em',
                    uiOutput('gsea_p_pos')
             )
           ),
           fluidRow(
             column(3, offset = 2, style='padding-top:30px;padding-right:10px',
                    uiOutput('gsea_padj_desc')
             ),
             column(1, offset = 0, style = 'padding-left:0px;margin-left:-4em',
                    uiOutput('gsea_padj_pos')
             )
           ),

           fluidRow(
             column(12, offset = 0,

                    uiOutput("gsea_tbl"))
           ),

           fluidRow(

             column(2, offset = 8,
                    br(),
                    br(),
                    uiOutput('gsea_download')
             )

           ),

           # plot block ------------
           fluidRow(
             column(8, offset = 0,
                    br(),
                    br(),
                    uiOutput('gsea_lolliplot')
             )
           ),
           fluidRow(
             column(10, offset = 1,
                    br(),
                    uiOutput('gseaplot_num'))
           ),
           fluidRow(
             column(10, offset = 1,
                    br(),
                    uiOutput('gsea_gseaplot'))
           )
  )
)
