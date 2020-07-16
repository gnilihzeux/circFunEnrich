#' enrichment score
#' @export enrichScore
#'
#' @return NULL
#'
#' @import clusterProfiler
#' @import ReactomePA


enrichScore <- function(enrichReslt){
  # enrichment score = overlapGeneCount*bgGeneNum / (diffGeneNum*termGeneNum)
  #
  # Args:
  #   enrichReslt: enrichGO's or enrichKEGG's result, class-enrichReslt
  #
  # Returns:
  #   enrichment score, class-numeric
  overlapGeneCount <- as.numeric(sapply(strsplit(enrichReslt$GeneRatio, "/"), "[", 1))
  diffGeneNum      <- as.numeric(sapply(strsplit(enrichReslt$GeneRatio, "/"), "[", 2))
  bgGeneNum        <- as.numeric(sapply(strsplit(enrichReslt$BgRatio, "/"), "[", 2))
  termGeneNum      <- as.numeric(sapply(strsplit(enrichReslt$BgRatio, "/"), "[", 1))

  overlapGeneCount*bgGeneNum / (diffGeneNum*termGeneNum)
}

#' enrichment analysis for GO
#' @export goEn
#'
#' @return NULL
#'
#'
goEn <- function(gene,
                 Ont          = "BP",
                 universeList = NULL,
                 orgdb        = NULL,
                 keyType      = "SYMBOL",
                 minGeneNum   = 5,
                 maxGeneNum   = 500,
                 pvalueCutoff = 1,
                 qvalueCutoff = 1,
                 resultPvalueCutoff = 0.05,
                 minOverlapNum = 3){
  entrez_tbl <- gene
  if (keyType == "ENTREZID") {
    gene <- entrez_tbl$ENTREZID
  } else if (keyType == "SYMBOL") {
    gene <- entrez_tbl$SYMBOL
  }

  enrich_go <- enrichGO(gene          = gene,
                        OrgDb         = orgdb,
                        keyType       = keyType, # keytype or keyType
                        ont           = Ont,
                        pvalueCutoff  = pvalueCutoff,
                        pAdjustMethod = "BH",
                        qvalueCutoff = qvalueCutoff,
                        minGSSize     = minGeneNum,
                        maxGSSize     = maxGeneNum,
                        readable      = FALSE
  )@result
  type <- switch(Ont,
                 BP = "biological process",
                 CC = "cellular component",
                 MF = "molecular function")
  format_go <- data.frame(databaseID    = enrich_go$ID,
                          Description   = enrich_go$Description,
                          type          = type,
                          geneRatio     = enrich_go$GeneRatio,
                          bgRatio       = enrich_go$BgRatio,
                          pvalue        = enrich_go$pvalue,
                          padj          = enrich_go$p.adjust,
                          qvalue        = enrich_go$qvalue,
                          enrichScore   = enrichScore(enrich_go),
                          overlapGeneList = enrich_go$geneID,
                          overlapGeneCount = enrich_go$Count,
                          stringsAsFactors = FALSE
  )

  # transform
  if (keyType == "ENTREZID") {
    format_go$overlapGeneList <- sapply(strsplit(format_go$overlapGeneList, "/"),
                                        function(x)
                                          paste(entrez_tbl$SYMBOL[match(x, entrez_tbl$ENTREZID)], collapse = "/")
    )
  }

  format_go   = format_go[order(format_go$pvalue, - format_go$overlapGeneCount), ]
  format_go   = format_go[format_go$pvalue < resultPvalueCutoff, ]
  fb          = format_go[format_go$overlapGeneCount >= minOverlapNum, ]
  if (nrow(fb) >= 3)
    format_go <- fb
  format_go
}

#' enrichment analysis for KEGG
#' @export keggPath
#'
#' @return NULL
#'

keggPath <- function(gene,
                     universeList = NULL,
                     species      = "hsa",
                     keyType      = "SYMBOL",
                     minGeneNum   = 5,
                     maxGeneNum   = 500,
                     pvalueCutoff = 1,
                     qvalueCutoff = 1,
                     resultPvalueCutoff = 0.05,
                     minOverlapNum = 3,
                     use_internal_data = FALSE
){
  entrez_tbl <- gene
  #' if (species == "hsa" && use_internal_data) {
  #'   #' @note `quitely` must been add
  #'   require(KEGG.db, quietly = TRUE)
  #'   use_internal_data = TRUE
  #' }
  #'
  if (keyType == "ENTREZID") {
    gene <- entrez_tbl$ENTREZID
  } else if (keyType == "SYMBOL") {
    gene <- entrez_tbl$SYMBOL
  }

  enrich_kegg <- enrichKEGG(gene          = gene,
                            organism      = species,
                            keyType       = "kegg", # keytype or keyType
                            pvalueCutoff  = pvalueCutoff,
                            pAdjustMethod = "BH",
                            minGSSize     = minGeneNum,
                            maxGSSize     = maxGeneNum,
                            qvalueCutoff = qvalueCutoff,
                            use_internal_data = use_internal_data
  )@result
  format_kegg <- data.frame(databaseID          = enrich_kegg$ID,
                            Description = enrich_kegg$Description,
                            geneRatio     = enrich_kegg$GeneRatio,
                            bgRatio       = enrich_kegg$BgRatio,
                            pvalue        = enrich_kegg$pvalue,
                            padj          = enrich_kegg$p.adjust,
                            qvalue        = enrich_kegg$qvalue,
                            enrichScore   = enrichScore(enrich_kegg),
                            overlapGeneList = enrich_kegg$geneID,
                            overlapGeneCount = enrich_kegg$Count,
                            stringsAsFactors = FALSE
  )

  # transform
  if (keyType == "ENTREZID") {
    format_kegg$overlapGeneList <- sapply(strsplit(format_kegg$overlapGeneList, "/"),
                                          function(x)
                                            paste(entrez_tbl$SYMBOL[match(x, entrez_tbl$ENTREZID)], collapse = "/")
    )
  }

  format_kegg = format_kegg[order(format_kegg$pvalue, - format_kegg$overlapGeneCount), ]
  format_kegg = format_kegg[format_kegg$pvalue < resultPvalueCutoff, ]
  fk          = format_kegg[format_kegg$overlapGeneCount >= minOverlapNum, ]
  if (nrow(fk) >= 3)
    format_kegg <- fk
  format_kegg
}

#' enrichment analysis for reactome
#' @export reacPath
#'
#' @return NULL
#'

reacPath <- function(gene,
                     universeList = NULL,
                     species      = "human",
                     keyType      = "SYMBOL",
                     minGeneNum   = 5,
                     maxGeneNum   = 500,
                     pvalueCutoff = 1,
                     qvalueCutoff = 1,
                     resultPvalueCutoff = 0.05,
                     minOverlapNum = 3
){
  entrez_tbl <- gene
  if (keyType == "ENTREZID") {
    gene <- entrez_tbl$ENTREZID
  } else if (keyType == "SYMBOL") {
    gene <- entrez_tbl$SYMBOL
  }

  enrich_reactome <- enrichPathway(gene         = gene,
                                   organism     = species,
                                   pAdjustMethod = "BH",
                                   pvalueCutoff = pvalueCutoff,
                                   qvalueCutoff = qvalueCutoff,
                                   minGSSize    = minGeneNum,
                                   maxGSSize    = maxGeneNum,
                                   readable     = FALSE)@result
  format_reactome <- data.frame(databaseID        = enrich_reactome$ID,
                                Description   = enrich_reactome$Description,
                                geneRatio     = enrich_reactome$GeneRatio,
                                bgRatio       = enrich_reactome$BgRatio,
                                pvalue        = enrich_reactome$pvalue,
                                padj          = enrich_reactome$p.adjust,
                                qvalue        = enrich_reactome$qvalue,
                                enrichScore   = enrichScore(enrich_reactome),
                                overlapGeneList = enrich_reactome$geneID,
                                overlapGeneCount = enrich_reactome$Count,
                                stringsAsFactors = FALSE)
  # transform
  if (keyType == "ENTREZID") {
    format_reactome$overlapGeneList <- sapply(strsplit(format_reactome$overlapGeneList, "/"),
                                              function(x)
                                                paste(entrez_tbl$SYMBOL[match(x, entrez_tbl$ENTREZID)], collapse = "/")
    )
  }

  format_reactome = format_reactome[order(format_reactome$pvalue, - format_reactome$overlapGeneCount), ]
  format_reactome = format_reactome[format_reactome$pvalue < resultPvalueCutoff, ]
  fr          = format_reactome[format_reactome$overlapGeneCount >= minOverlapNum, ]
  if (nrow(fr) >= 3)
    format_reactome <- fr
  format_reactome
}

#' the main function
#' @export goPathwayEnrichment
#'
#' @return NULL
#'

goPathwayEnrichment <- function(degList,
                                universeList = NULL,
                                species      = "hsa",
                                keyType      = "SYMBOL",
                                species_db   = NULL,
                                minGeneNum   = 5,
                                maxGeneNum   = 500,
                                pvalueCutoff = 1,
                                qvalueCutoff = 1,
                                resultPvalueCutoff = 0.05,
                                minOverlapNum = 3
){
  # GO and KEGG pathway enrichment analysis with zebrafisher's test based on local database
  #
  # Args:
  #   degList     : differential expression gene list
  #                 class-character
  #   universeList: background gene list
  #                 class-character
  #   species     : human, mouse or rat, correspongding to 'hsa', 'mmu' or 'rno'
  #                 class-character
  #   minGeneNum  : go/pathway item contains at least 5 genes
  #                 class-numeric
  #   maxGeneNum  : go/pathway item contains at most 500 genes
  #                 class-numeric
  #
  # Returns:
  #   4 tables with columns: go/pahID, go/pathDescription, goType, geneRatio, bgRatio, pvalue, padj, overlapGeneList, overlapGeneCount
  #                 class-list
  #           1. go:bp class-data.frame
  #           2. go:cc class-data.frame
  #           3. go:mf class-data.frame
  #           4. kegg  class-data.frame
  # @issue package version
  #   there are some differences between 'v3.4' and 'v3.6' of 'clusterProfiler':
  #           1. Bioconductor version
  #                Bioc v3.4 --> v3.4
  #                Bioc v3.6 --> v3.6
  #           2. R version
  #                R v3.3.x --> v3.4
  #                R >= v3.4.2 --> v3.6
  #           3. enrichGO
  #                keytype --> v3.4
  #                keyType --> v3.6
  #  @strategy
  #    omit the argue name 'keytype/keyType' but input parameters in order

  # header --------------------------------------
    suppressMessages(require(clusterProfiler))
    suppressMessages(require(ReactomePA))

  options(stringsAsFactors = FALSE)
  options(digits = 7)

  # input ---------------------------------------
  # species
  species <- switch(species,
                    human = "hsa",
                    hsa   = "hsa",
                    mouse = "mmu",
                    mmu   = "mmu",
                    rat   = "rno",
                    rno   = "rno",
                    dre   = "dre",
                    zebrafish  = "dre",
                    aalb  = "aalb",
                    aedes_albopictus = "aalb"
  )
  SPECIES <- switch(species,
                    human = "human",
                    hsa   = "human",
                    mouse = "mouse",
                    mmu   = "mouse",
                    rat   = "rat",
                    rno   = "rat",
                    dre   = "zebrafish",
                    zebrafish  = "zebrafish",
                    aalb  = "aedes_albopictus",
                    aedes_albopictus = "aedes_albopictus"
  )
  if (is.null(species_db)) {
    species_db <- switch(species,
                         hsa = "org.Hs.eg.db",
                         mmu = "org.Mm.eg.db",
                         rno = "org.Rn.eg.db",
                         dre = "org.Dr.eg.db"
    )
    require(species_db, character.only = TRUE)
  }


  # id transform
  entrez_tbl <- bitr(degList, fromType = keyType, toType = c("ENTREZID", "SYMBOL"), OrgDb = species_db)

  format_bp <- tryCatch(goEn(gene          = entrez_tbl,
                             Ont           = "BP",
                             keyType       = "ENTREZID", # keytype or keyType
                             orgdb         = species_db,
                             pvalueCutoff  = pvalueCutoff,
                             qvalueCutoff  = qvalueCutoff,
                             resultPvalueCutoff = resultPvalueCutoff,
                             minGeneNum     = minGeneNum,
                             maxGeneNum     = maxGeneNum,
                             minOverlapNum = minOverlapNum),
                        error = function(e)e)
  if (inherits(format_bp, "simpleError")) {
    if(sum(grep("with no slots", format_bp$message) != 0))
      format_bp <- "--> No gene can be mapped...\n"
  }

  format_cc <- tryCatch(goEn(gene          = entrez_tbl,
                             Ont           = "CC",
                             keyType       = "ENTREZID", # keytype or keyType
                             orgdb         = species_db,
                             pvalueCutoff  = pvalueCutoff,
                             qvalueCutoff  = qvalueCutoff,
                             resultPvalueCutoff = resultPvalueCutoff,
                             minGeneNum     = minGeneNum,
                             maxGeneNum     = maxGeneNum,
                             minOverlapNum = minOverlapNum),
                        error = function(e)e)
  if (inherits(format_cc, "simpleError")) {
    if(sum(grep("with no slots", format_cc$message) != 0))
      format_cc <- "--> No gene can be mapped...\n"
  }


  format_mf <- tryCatch(goEn(gene          = entrez_tbl,
                             Ont           = "MF",
                             keyType       = "ENTREZID", # keytype or keyType
                             orgdb         = species_db,
                             pvalueCutoff  = pvalueCutoff,
                             qvalueCutoff  = qvalueCutoff,
                             resultPvalueCutoff = resultPvalueCutoff,
                             minGeneNum     = minGeneNum,
                             maxGeneNum     = maxGeneNum,
                             minOverlapNum = minOverlapNum),
                        error = function(e)e)
  if (inherits(format_mf, "simpleError")) {
    if(sum(grep("with no slots", format_mf$message) != 0))
      format_mf <- "--> No gene can be mapped...\n"
  }

  format_kegg <- tryCatch(keggPath(gene          = entrez_tbl,
                                   species       = species,
                                   keyType       = "ENTREZID", # keytype or keyType
                                   pvalueCutoff  = pvalueCutoff,
                                   resultPvalueCutoff = resultPvalueCutoff,
                                   minGeneNum     = minGeneNum,
                                   maxGeneNum     = maxGeneNum,
                                   qvalueCutoff  = qvalueCutoff,
                                   minOverlapNum = 1,
                                   use_internal_data = FALSE),
                          error = function(e)e)
  if (inherits(format_kegg, "simpleError")) {
    if(sum(grep("with no slots", format_kegg$message) != 0))
      format_kegg <- "--> No gene can be mapped...\n"
  }


  format_reactome <- tryCatch(reacPath(gene         = entrez_tbl,
                                       keyType      = "ENTREZID",
                                       species      = SPECIES,
                                       pvalueCutoff = pvalueCutoff,
                                       qvalueCutoff = qvalueCutoff,
                                       resultPvalueCutoff = resultPvalueCutoff,
                                       minGeneNum    = minGeneNum,
                                       maxGeneNum    = maxGeneNum,
                                       minOverlapNum = minOverlapNum),
                              message = function(e)e)
  if (inherits(format_reactome, "message"))
    format_reactome <- format_reactome$message

  list(format_bp       = format_bp,
       format_cc       = format_cc,
       format_mf       = format_mf,
       format_kegg     = format_kegg,
       format_reactome = format_reactome)
}
