library(monocle)
library(scater)
library(dplyr)
library(stringr)
library(biomaRt)
library(M3Drop)
library(ggthemes)
library(ggbeeswarm)
library(DT)

ui = fluidPage(title = "monocle analyses",
               tabsetPanel(
                 tabPanel( title = "Sample Selector",
                           fluidRow(column(1,actionButton("deleteRows", "Delete Rows"))),
                           fluidRow(column(3,dataTableOutput("table1"))),
                           fluidRow(column(1,actionButton("recalculate","Calculate")))
                 ),
                 tabPanel( title = "M3Drop",
                           p("Identifies significantly highly variable genes"),
                           a("See M3Drop" , href = "https://www.bioconductor.org/packages/release/bioc/vignettes/M3Drop/inst/doc/M3Drop_Vignette.pdf" ),
                           plotOutput("m3drop")
                 ),
                 tabPanel( title = "Diff genes",
                           tableOutput("diff_genes")
                 ),
                 tabPanel( title = "Pseudotime course",
                           plotOutput("pseudotime")
                 ),
                 tabPanel(title = "fetal cells",
                          plotOutput("fetal_cells")
                 )
               )
)

# read qc scater object
reads.qc<-readRDS('data/reads.monocle.rds')
reads.qc <- reads.qc[rowData(reads.qc)$use, colData(reads.qc)$use]
endog_genes <- !rowData(reads.qc)$is_feature_control
reads.qc <- reads.qc[endog_genes, ]
# create read matrix
reads.mat<-assay(reads.qc)
#list object created by M3Drop
clean_list<- M3DropCleanData(reads.mat,suppress.plot = T)
# create empty vector to receive genes from M3Drop 
DE_genes <- vector()
# data frame of sample names
df <-data.frame(colData(reads.qc)[,c("dev_stage","expt","gate")]) #d)[,c("dev_stage","expt","gate")]
df$samples<-rownames(df)
rownames(df)<-1:nrow(df)
# temporary dataframe for monocle data
monocle_df <- data.frame(type = vector(), pseudotime = vector(), State = vector, stringsAsFactors = F)

server = function(input, output) {
values <- reactiveValues(df_data = df,                # reactive vector of sample names
                         mat_data = clean_list$data,  # reactive matrix for M3Drop 
                         genes = DE_genes,            # reactive differentially expressed gene list
                         reads_object = reads.qc,     # reactive scater object reads.qc 
                         monocle_time = monocle_df    # data frame for monocle states
  )  
  
##On press delete rows button
observeEvent(input$deleteRows, {
    if (!is.null(input$table1_rows_selected)){
      values$df_data<-values$df_data[-as.numeric(input$table1_rows_selected),]
    }
  })
  
  output$table1 <- renderDataTable({
    values$df_data
  })
  # press recalculate button
  observeEvent(input$recalculate, {
    output$m3drop <- renderPlot( {
      # on pressing recalculate run m3dop clean data on selected samples
      temp_list<- M3DropCleanData(values$mat_data[,values$df_data$samples],suppress.plot = T)
      temp_mat<-temp_list$data
      values$genes<-M3DropDifferentialExpression(temp_mat,
                                                mt_method = "fdr",
                                                mt_threshold = 0.01,suppress.plot = F)
    })
    
    output$diff_genes <- renderTable({
      values$genes
    })
    

    output$pseudotime <- renderPlot ({
      d <- values$reads_object[which(rownames(values$reads_object) %in% values$genes$Gene), as.character(values$df_data$samples)]
      colnames(d) <- 1:ncol(d)
      geneNames <- rownames(d)
      rownames(d) <- 1:nrow(d)
      pd <- data.frame(colData(d)[,c("dev_stage","expt","gate")])
     # pd <- data.frame(type = colData(d)$type)
      pd <- new("AnnotatedDataFrame", data=pd)
      fd <- data.frame(gene_short_name = geneNames)
      fd <- new("AnnotatedDataFrame", data=fd)
      dCellData <- newCellDataSet(assay(d), phenoData = pd, featureData = fd, expressionFamily = negbinomial())
      dCellData <- setOrderingFilter(dCellData, which(geneNames %in% values$genes$Gene))
      dCellData <- estimateSizeFactors(dCellData)
      dCellDataSet <- reduceDimension(dCellData, pseudo_expr = 1)
      dCellDataSet <- orderCells(dCellDataSet, reverse = FALSE)
      
      temp_df <-
        data.frame(
          dev_stage = phenoData(dCellDataSet)$dev_stage,
          pseudotime = phenoData(dCellDataSet)$Pseudotime,
          State = phenoData(dCellDataSet)$State,
          gate = phenoData(dCellDataSet)$gate,
          expt = phenoData(dCellDataSet)$expt
        )
      # 
      # temp_df<-data.frame(type=phenoData(dCellDataSet)$type,
      #                     pseudotime=phenoData(dCellDataSet)$Pseudotime,
      #                     State=phenoData(dCellDataSet)$State
      #                     )
      rownames(temp_df) <- 1:ncol(d)
      values$monocle_time<-temp_df
      plot_cell_trajectory(dCellDataSet)
    })

    output$fetal_cells <- renderPlot ({
              temp_reads.qc <- values$reads_object[,as.character(values$df_data$samples)]
              temp_reads.qc$pseudotime_monocle <- values$monocle_time$pseudotime
              temp_reads.qc$pseudotime_state <- values$monocle_time$State
              colData(temp_reads.qc)
              ggplot(as.data.frame(colData(temp_reads.qc)),
                     aes(x = values$monocle_time$pseudotime, y = colData(temp_reads.qc)$dev_stage, colour = factor(pseudotime_state),shape = expt)) +
                     geom_quasirandom(groupOnX = FALSE) +
                     scale_color_tableau() +
                     #scale_color_manual( values=c("green4",  "tomato" ,"dodgerblue")) +
                     theme_classic() +
                     xlab("monocle pseudotime") + ylab("stage") +
                     ggtitle("Cells ordered by monocle pseudotime")
    })

  })
}

shinyApp(ui = ui, server = server)