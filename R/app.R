library(shiny)
library(lpSolve)
library(shinyWidgets)
library(splatter)
library(SingleCellExperiment)
library(DropletUtils)
library(Seurat)
library(scater)
library(igraph)
library(ggplot2)
library(ggnetwork)
library(plotly)
library(DescTools)
library(dplyr)
library(Matrix)

#' The length of a string
#'
#' A function to load TREK barcodes into memory.
#'
#' @return Loads the TREK barcodes into memory.
#' 
#' @export
#' @examples
#' message('correct usage') #load_TREK()
load_TREK <- function() {
  master <- read.table("~/master_sc_variants.tsv", sep="\t", header=T)

  get_combined_names <- function(tab){
    paste(tab[,1], tab[,2], tab[,3], tab[,4], tab[,5], tab[,6])
  }
  combined_names <- get_combined_names(master)
  cells <- paste0(master[,8], "-1")

  master2 <- data.frame(var=combined_names, cell=cells, n=master[,9])

  meb <- master2 %>% group_by(var) %>% summarize(amount = n())

  meb2 <- meb[meb$amount > 500,]
  targs <- meb2 %>% pull(var)

  master2 <- master2[master2$var %in% targs,]
  master2['n'] <- as.numeric(match(master2[,3], unique(master2[,3])))
  master2['rown'] <- as.numeric(match(master2[,1], unique(master2[,1])))
  master2['cown'] <- as.numeric(match(master2[,2], unique(master2[,2])))

  master2 <- master2[order(master2[,5], master2[,4]),]

  matrixify <- function(mat, expression_matrix){
    message('pass1')
    dgc <- Matrix(ncol = ncol(expression_matrix), nrow=length(targs), sparse=TRUE, data=0)
    colnames(dgc) <- colnames(expression_matrix)
    rownames(dgc) <- targs
    dgc <- r_to_Cpp_to_R(dgc, as.numeric(mat[,4] - 1), as.numeric(mat[,5] - 1), as.numeric(mat[,3]), ncol = ncol(dgc))
    dgc
  }

  res <- matrixify(master2, expression_matrix)
  res
}

run <- function(dat, weights) {
  browser()
  f.obj <- weights
  f.con <- dat
  len1 <- nrow(dat) - 2 * length(weights)
  len2 <- length(weights)
  f.dir <- c(rep('>=', len1 + len2), rep('<=', len2))
  f.rhs <- c(rep(1, len1), rep(0, len2), rep(1, len2))

  solution <- lp ("min", f.obj, f.con, f.dir, f.rhs)
  res <- solution$solution
  npick <- runif(len2)
  pick <- ifelse(res >= npick, 1, 0)
  pick
}

run1 <- function(dat, weights) {
  f.obj <- weights
  len1 <- ncol(dat) - 2 * length(weights)
  len2 <- length(weights)
  f.con <- dat[1:len1, 1:len2]
  f.dir <- c(rep('>=', len1))
  f.rhs <- c(rep(1, len1))

  solution <- lp ("min", f.obj, f.con, f.dir, f.rhs)
  res <- solution$solution
  npick <- runif(len2)
  pick <- ifelse(res >= npick, 1, 0)
  pick
}

iterate <- function(dat, weights, times) {
  len1 <- nrow(dat) - 2 * length(weights)
  len2 <- length(weights)
  runs <- lapply(1:times, function(i) {
    run(dat, weights)
  })
  runs <- do.call(rbind, runs)
  runs <- colSums(runs)
  res <- ifelse(runs > 0, 1, 0)
  dat100 <- dat[1:len1, 1:len2]
  message("Selected ", sum(res), " antigens.")
  message("Total Covered: ", sum(rowSums(dat100[,res==1])>0))
  message("Total Weight: ", sum(weights * res))
  list(which(res>0), which(rowSums(dat100[,res==1]) > 0))
}

run_lp_sim <- function(ncells, nantigens, k, times) {
  dat <- sim_antigens(ncells, nantigens, k)
  dat <- do.call(rbind, dat)
  weights <- runif(nantigens)
  iterate(dat, weights, times)
}

if(FALSE) {
ncells <- 1000
nantigens <- 100
k <- 30

dat <- sim_antigens(ncells, nantigens, k)
dat <- do.call(rbind, dat)

dat <- as.matrix(mat)
dat <- t(dat)
dat[dat > 0] <- 1
tab2 <- lapply(1:ncol(dat), function(j) {
  res <- rep(0, ncol(dat))
  res[j] <- 1
  res
})

tab2 <- do.call(rbind, tab2)

dat <- rbind(dat, tab2)
dat <- rbind(dat, tab2)

weights <- rep(1, nrow(mmat))
ncells <- ncol(mat)

iter <- 10

res <- lapply (1:iter, function(i) {

  coverage <- c()
  antigens <- c()
  times <- 0

  while(length(coverage) != ncells) {
    message("Start coverage: ", length(coverage))
    li <- iterate(dat, weights, 1)
    antigens <- c(antigens, as.character(li[[1]]))
    antigens <- antigens[!duplicated(antigens)]
    coverage <- c(coverage, as.character(li[[2]]))
    coverage <- coverage[!duplicated(coverage)]
    times <- times + 1
    message("Coverage: ", length(coverage))
    message("Antigens: ", length(antigens))
    message("Times: ", times)
  }
  list(coverage, antigens, times)
})
}

run_sim <- function(expression_matrix, nantigens, max_antigens, ncells, dropout, libloc){

  #expression_matrix <- Seurat::Read10X(data.dir = expression_matrix)
  rn <- rownames(expression_matrix)

  gs <- sample(rn, nantigens) #c('CACNA1S', 'CAPN3', 'CENPF', 'CXCR2',	'DPEP1', 'FABP6', 'FAM216A', 'FBXL2', 'GABBR1', 'GREB1', 'IVD', 'JADE3', 'LACTB2', 'LRP8', 'MAGEA11', 'MAT2A', 'MMP11', 'MOK',	'MPP3',	'MVK', 'NPTX2', 'NUP43', 'OVGP1',	'PDXDC1',	'PTGER4',	'RRAD',	'SLC7A5',	'STK31',	'TDRD9',	'TFDP',	'TRIP13',	'TUBB4A',	'UBE4A',	'WFS1',	'ZNF200')
  ww <- which(rownames(expression_matrix) %in% gs)
  ee <- expression_matrix[ww,]

  #ncells <- 100
  params <- newSplatParams(group.prob = c(0.5, 0.5), nGenes = length(gs), batchCells=ncells, dropout.type='experiment', dropout.mid=dropout, lib.loc = libloc)
  sim.groups <- splatSimulate(params, method="groups",
                              verbose = FALSE)

  sim.groups1 <- sim.groups[,1:(ncells/2)] #sim.groups[,which(col.group == "Group1")]
  sim.groups2 <- sim.groups[,((ncells/2) + 1):(ncells)] #which(col.group == "Group2")

  tc <- assays(sim.groups)['TrueCounts'][[1]]
  rownames(tc) <- gs
  ass <- assay(sim.groups)
  ass[ass < 1] <- 0
  col.group <- colData(sim.groups)$Group


  dd <- as.data.frame(ee)
  dat <- dd

  ww <- which(rownames(tc) %in% gs)
  ee <- tc[ww,]

  dat <- data.frame(ee)

  cols <- sapply(dat, is.logical)
  dat[,cols] <- lapply(dat[,cols], as.numeric)
  lll <- as.list(dat)

  group_dat <- colData(sim.groups)
  vv <- vapply(lll[seq_along(lll)], function(x){
    nnn <- as.character(as.numeric(x > 4))
    paste(nnn, collapse="")
  }, character(1))

  gc <- tc
  colnames(gc) <- group_dat$Group


  gg <- vv
  names(gg) <- group_dat$Group

  df <- data.frame(group = names(gg), profile = unname(gg))

  dfList <- split(df, df$profile)

  dfList= dfList[order(sapply(dfList,nrow),decreasing = T)]

  groupo <- vapply(dfList, function(x) {
    names(which.max(table(x$group)))
  }, character(1))

  df_len <- vapply(dfList, function(x) {nrow(x)}, numeric(1))
  df_ooo <- order(df_len, decreasing=T)
  df_ooo <- order(df_len, decreasing=T)

  ttt <- table(vv)


  num <- max_antigens
  ooo <- ttt[order(ttt, decreasing=T)]
  num <- min(num, length(ooo[!is.na(ooo)]))
  ooo <- ooo[1:num]

  uooo <- unname(ooo)
  rrr <- rownames(ooo)

  dister <- lapply(rownames(ooo), function(x){
    res <- lapply(rownames(ooo), function(y){
      if(compare(x, y) %in% c(0, 1))
        data.frame(x, y, text = x)
    })
    do.call(rbind, res)
  })

  dister <- lapply(seq_along(dister), function(x) {
    temp <- dister[[x]]
    temp$group <- groupo[x]
    temp
  })


  lengo <- lapply(dister, function(x) {nrow(x)})
  lengo <- unlist(lengo)

  roro <- rep(rrr, times = lengo)
  uouo <- rep(uooo, times = lengo)

  distero <- do.call(rbind, dister)
  distero <- distero[!duplicated(distero),]
  distero$size <- uouo

  disto <- unname(unlist(dister))


  net <- graph.data.frame(distero, directed = T)
  df_net <- list()
  set.seed(1)
  dfneto <- ggnetwork(net)
  dfneto <- dfneto[!is.na(dfneto$size),]
  df_net$dfnet <- dfneto
  df_net$gs <- gs
  df_net
}

run_real <- function(expression_matrix, trek, max_antigens, targets, find_minimum_coverage_set, remove0){

  #expression_matrix <- Seurat::Read10X(data.dir = expression_matrix)
  rn <- rownames(expression_matrix)

  expression_matrix <- trek

  browser()

  #gs <- vap2
  #gs <- c(unique(master2[,1]))
#  message(gs)
  gs <- rownames(trek)
  #gs <- c('CACNA1S', 'CAPN3', 'CENPF', 'CXCR2',	'DPEP1', 'FABP6', 'FAM216A', 'FBXL2', 'GABBR1', 'GREB1', 'IVD', 'JADE3', 'LACTB2', 'LRP8', 'MAGEA11', 'MAT2A', 'MMP11', 'MOK',	'MPP3',	'MVK', 'NPTX2', 'NUP43', 'OVGP1',	'PDXDC1',	'PTGER4',	'RRAD',	'SLC7A5',	'STK31',	'TDRD9',	'TFDP',	'TRIP13',	'TUBB4A',	'UBE4A',	'WFS1',	'ZNF200')
  targets <- which(gs %in% targets)
  ww <- which(rownames(expression_matrix) %in% gs)
  gs <- gs[gs %in% rownames(expression_matrix)]
  ee <- expression_matrix[ww,]

  dd <- as.data.frame(ee)
  dat <- dd

  cols <- sapply(dat, is.logical)
  dat[,cols] <- lapply(dat[,cols], as.numeric)
  lll <- as.list(dat)

  vv <- vapply(lll[seq_along(lll)], function(x){
    nnn <- as.character(as.numeric(x > 1))
    paste(nnn, collapse="")
  }, character(1))

  ttt <- table(vv)

  num <- max_antigens
  ooo <- ttt[order(ttt, decreasing=T)]
  num <- min(num, length(ooo[!is.na(ooo)]))
  ooo <- ooo[1:num]

  uooo <- unname(ooo)
  rrr <- rownames(ooo)

  dister <- lapply(rownames(ooo), function(x){
    res <- lapply(rownames(ooo), function(y){
      if(compare(x, y) %in% c(0, 1))
        data.frame(x, y, text = x)
    })
    do.call(rbind, res)
  })

  lengo <- lapply(dister, function(x) {nrow(x)})
  lengo <- unlist(lengo)

  roro <- rep(rrr, times = lengo)
  print(lengo)
  uouo <- rep(uooo, times = lengo)

  distero <- do.call(rbind, dister)
  distero <- distero[!duplicated(distero),]
  distero$size <- uouo

  if(remove0){
    distero <- distero[!distero$x==paste0(as.character(rep(0, length(gs))), collapse=""),]
    distero <- distero[!distero$y==paste0(as.character(rep(0, length(gs))), collapse=""),]
  }

  disto <- unname(unlist(dister))

  values <- vapply(1:nrow(distero), function(x){
    ret <- "Not targetted"
    for(y in targets){
      if(substr(distero[x,1], y, y) == "1")
        if(ret == "Not targetted")
          ret <- gs[y]
    }
    ret
  }, character(1))

  distero$targets <- values

  net <- graph.data.frame(distero, directed = T)
  df_net <- list()
  set.seed(1)
  dfneto <- ggnetwork(net)
  dfneto <- dfneto[!is.na(dfneto$size),]

  pepe <- paste0(dfneto$targets, " (", 100 * round(dfneto$size/ncol(expression_matrix), 4) ,"%)")

  colos <- lapply(seq_along(dfneto$name), function(x) {
    i <- dfneto$name[x]
    ul <- as.logical(as.numeric(unlist(strsplit(i, split=""))))
    wh <- which(ul)
    data.frame(cluster=rep(x, length(wh)), element=wh)
  })



  if(find_minimum_coverage_set) {
  colos <- lapply(seq_along(dfneto$name), function(x) {
    i <- dfneto$name[x]
    ul <- unlist(strsplit(i, split=""))
    ul <- as.logical(as.numeric(ul))
    ul <- gs[ul]
    ul <- which(gs %in% ul)
    df <- data.frame(n = ul, sets=rep(paste0('cluster', x), length(ul)))
    df
  })

  #colos <- do.call(cbind, list(data.frame(antigens = seq_along(gs)), colos))
  colos <- do.call(rbind, colos)
  #colos <- colos[-1]

  items.mat <- table(colos$sets, colos$n)
  cos <- colnames(items.mat)
  gss <- gs[as.numeric(cos)]
  dimnames(items.mat) = list( items=rep(paste0('cluster', seq_len(nrow(items.mat)))), antigens=cos)

  #colos <- data.matrix(colos)

  f.obj <-  rep(1,ncol(items.mat))
  f.dir <- rep(">=", nrow(items.mat))
  f.rhs <- rep(1, nrow(items.mat))

  sol <- lp ("min", f.obj, items.mat, f.dir, f.rhs)$solution
  covering <- gss[as.logical(sol)]

  covering_w <- which(gs %in% covering)

  values2 <- vapply(1:nrow(dfneto), function(x){
    ret <- "Not targetted"
    for(y in covering_w){
      if(substr(dfneto[x,3], y, y) == "1")
        if(ret == "Not targetted")
          ret <- gs[y]
    }
    ret
  }, character(1))

  dfneto$targets <- values2

  }

  dfneto2 <- dfneto[!duplicated(dfneto$name),]
  dfneto2 <- aggregate(dfneto2$size, by=list(Category=dfneto2$targets), FUN=sum)
  colnames(dfneto2) <- c("targets", "Neoantigen_gene_targets")
  total_cells <- sum(dfneto2$Neoantigen_gene_targets)
  dfneto2$Neoantigen_gene_targets <- paste0(dfneto2$targets, " (", 100 * round(dfneto2$Neoantigen_gene_targets/total_cells, 4) ,"%)")
  dfneto <- left_join(dfneto, dfneto2, by="targets")

  df_net$dfnet <- dfneto

#  colos <- lapply(seq_along(df_net$dfnet$name), function(x) {
#    i <- df_net$dfnet$name[x]
#    ul <- unlist(strsplit(i, split=""))
#    ul <- as.numeric(ul)
#    df <- data.frame(ul)
#    names(df) <- paste0('cluster', x)
#    rownames(df) <- seq_along(gs)
#    df
#  })




  df_net$covering <- covering
  df_net$gs <- gs
  df_net
}

if(FALSE){

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Tumor Neoantigen Profiler"),

    # Sidebar with a slider input for number of bins
    sidebarLayout(
        sidebarPanel(
            checkboxInput("simulated", label = "use simulated data", value = TRUE),
            textInput("expression_matrix", label = h4("Expression matrix file"), value='/Users/vantwisk/filtered_feature_bc_matrix'),
            textInput("TREK output", label = h4("TREK output file"), value='/Users/vantwisk/master_sc_variants.tsv'),
            checkboxInput("find_minimum_coverage_set", label = "Find Minimum Coverage Set", value = TRUE),
            conditionalPanel(condition = "input.find_minimum_coverage_set == false",
            pickerInput(
              inputId = "targets",
              label = "Top Genes to target",
              choices = unique(master2[,1]),
              options = list(
                `actions-box` = TRUE,
                size = 10,
                `selected-text-format` = "count > 3"
              )),
              multiple = TRUE
            ),
            numericInput("max_antigens", label = h4("Max # antigens profiles"), value=30),
            conditionalPanel( condition = "input.simulated == true",
                              numericInput("nantigens", label = h4("# antigens"), value=30),
                              numericInput("min_cell_app", label = h4("Minimum occurence in cells"), value=30),
                              numericInput("ncells", label = h4("# simulated cells"), value=100),
                              numericInput("dropout", label = h4("Dropout"), value=2),
                              numericInput("libloc", label = h4("Library Location"), value=11)),
            checkboxInput("remove0", label = "Remove Cells Without Antigen Targets", value = FALSE),
            checkboxInput("logscale", label = "Log Scale Size of graph", value = FALSE)
        ),

        # Show a plot of the generated distribution
        mainPanel(
           plotlyOutput("graphPlot"),
           tableOutput("antigenTable")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {

  getExpressionMatrix <- reactive({
    Seurat::Read10X(data.dir = input$expression_matrix)
  })

  getTREKVarients <- reactive({
    load_TREK()
  })

  datasetInput <- reactive({
    #expression_matrix <- getExpressionMatrix()
    if(!input$simulated) {
      res <- run_real(getExpressionMatrix(), getTREKVarients(), input$max_antigens, input$targets, input$find_minimum_coverage_set, input$remove0)
    } else {
      res <- run_sim(expression_matrix, input$nantigens, input$max_antigens, input$ncells, input$dropout, input$libloc)
    }
    res
  })

#  output$dfnet <- reactive({
#    datasetInput()
#  })

  output$graphPlot <- renderPlotly({
    res <-datasetInput()
    df_net <- res$dfnet
    Percent_of_cells_targetted <- df_net$size
    is_log <- ""
    if(input$logscale){
      Percent_of_cells_targetted <- log(Percent_of_cells_targetted)
      is_log <- "Log "
    }
    plot3 <- ggplot(df_net, aes(x = x, y = y, xend = xend, yend = yend)) +
      geom_edges(size=0.4, alpha = 0.25) +
      geom_nodes(aes(color = Neoantigen_gene_targets, size = Percent_of_cells_targetted, text = text), show.legend=T) +
      theme_blank() +
      labs(color = "", size="AML antigen gene (% of cells)") +
      ggtitle(paste0("Single-Cell Transcriptomic Network of ", input$max_antigens, " Most Common Antigen Profiles"))

    plot3 %>% ggplotly(tooltip="text")
  })

  output$antigenTable <- renderTable(
    {
      res <- datasetInput()
      covering <- res$covering
      X <- data.frame(covering_neoantigens = covering)

      return(X)
    }
  )

  outputOptions(output, "graphPlot", suspendWhenHidden = FALSE)
}

#    output$graphPlot <- renderPlotly({

#      if(!input$simulated) {
#        res <- run_real(expression_matrix, input$max_antigens, input$targets)
#      } else {
#        res <- run_sim(expression_matrix, input$nantigens, input$max_antigens, input$ncells, input$dropout, input$libloc)
#      }

#      df_net <- res$dfnet
#      df_net <- df_net[!is.na(df_net$text),]

#        plot3 <- ggplot(df_net, aes(x = x, y = y, xend = xend, yend = yend)) +
#          geom_edges(size=0.4, alpha = 0.25) +
#          geom_nodes(aes(color = targets, size = size, text = text)) +
#          theme_blank() +
#          ggtitle(paste0("Simulated Single-Cell Transcriptomic Network of ", num, " Most Common Antigen Profiles"))

#        plot3 %>% ggplotly(tooltip="text")
#    })

#    output$graphplot <- renderPlotly({

#    })
#}

# Run the application
shinyApp(ui = ui, server = server)
}
