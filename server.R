# an App that uses ISB-CGC uploaded databases
# This App will use PanCan(TCGA), BRCA (breast cancer) cohort
# App will display cross cohort between TCGA, and user data
# All data must be in google cloud platform, BigQuery

library(bigrquery)
library(data.table)

server <- function(input, output) {
  
  
  # OUTPUTS -----------------------------------------------------------------
  
  
  output$sqlresults = renderPlot({
    
    mypick2 = input$myPicker2
    mypick3 = input$myGene
    # App Auto runs, this is to catch on first empty run
    if (all(c(is.null(mypick2), nchar(mypick3)==0))) {
      output$graphtext = renderText("No Criteria selected")
    } 
    else if (length(mypick2) < 2) {
      output$graphtext = renderText("Please select at least two cohort subtypes")
    }
    else if (nchar(mypick3) < 1){
      output$graphtext = renderText("Please enter a gene")
    }
    # Everything else is for code after selection
    else {
      output$graphtext = renderText("")
      globdata()
    }
  })
  

  # Input Data --------------------------------------------------------------
  
  
  globdata = reactive({
    project = "<<Your Google Cloud Platform Project>>"
    
    # SQL Query for TCGA Datasets ---------------------------------------------
    mylist = input$myPicker2
    addthis <- ""
    for (i in 1:length(mylist)) {
      x = switch(mylist[i], 
                 "Basal" = " 'BRCA.Basal' ",
                 "HER2" = " 'BRCA.Her2' ", 
                 "LuminalA" = "'BRCA.LumA' ",
                 "LuminalB" = " 'BRCA.LumB' ",
                 "Normal" = " 'BRCA.Normal' "
      )
      addthis <- paste0(addthis, x)
      if (i != length(mylist)) {
        addthis <- paste0(addthis, ", ")
      }
    }
    
    # This SQL Query searches 2 ISB-CGC TCGA uploaded tables
    # Both tables are needed to get expression data values
    # and user specified cohorts
    sql1 = paste0("SELECT `pancancer-atlas.Individual_Manuscript_Tables.Pan_Immune_Feature_Matrix_mmc2`.string_field_0, 
                  `pancancer-atlas.Individual_Manuscript_Tables.Pan_Immune_Feature_Matrix_mmc2`.string_field_3, 
                  `pancancer-atlas.Filtered.EBpp_AdjustPANCAN_IlluminaHiSeq_RNASeqV2_genExp_filtered`.normalized_count 
                  FROM (`pancancer-atlas.Filtered.EBpp_AdjustPANCAN_IlluminaHiSeq_RNASeqV2_genExp_filtered` INNER JOIN 
                  `pancancer-atlas.Individual_Manuscript_Tables.Pan_Immune_Feature_Matrix_mmc2` ON 
                  `pancancer-atlas.Filtered.EBpp_AdjustPANCAN_IlluminaHiSeq_RNASeqV2_genExp_filtered`.ParticipantBarcode = 
                  `pancancer-atlas.Individual_Manuscript_Tables.Pan_Immune_Feature_Matrix_mmc2`.string_field_0) WHERE 
                  `pancancer-atlas.Filtered.EBpp_AdjustPANCAN_IlluminaHiSeq_RNASeqV2_genExp_filtered`.Symbol = '", 
                  input$myGene, "' AND `pancancer-atlas.Individual_Manuscript_Tables.Pan_Immune_Feature_Matrix_mmc2`.string_field_1 = '", 
                  input$myPicker1, "' AND `pancancer-atlas.Individual_Manuscript_Tables.Pan_Immune_Feature_Matrix_mmc2`.string_field_3 
                  IN (", addthis, ")") 
    
    tb1 = bq_project_query(project, sql1)
    # This turns the Query into a Dataframe
    df1 = data.frame(bq_table_download(tb1))
    n = length(unique(df1$string_field_3))
    
    # Normalize values same way as ISPY Data
    df1$normalized_count = log2(df1$normalized_count) + 1
    tcgaval = data.frame(df1)
    colnames(tcgaval) = c("ParticipantBarcode", "Study", "Gene_Expression_Values")
    # Remove leading 'BRCA' from Subtype
    tcgaval$Study = sapply(strsplit(tcgaval$Study,"\\."), `[`, 2)
    
    
    #  SQL Query for ISPY Data ------------------------------------------------
    
    # This code was originally built for ISPY specific data
    # User will have to pull their own expression and subtype data tables, 
    # in order to filter with TCGA data above
    # SQl query will pull expression data, subtype, and filter on input$myGene
    sql3 = paste0("SELECT `<<subtype data table>>`.PID, `<<subtype data table>>`.Call, 
    `<<expression data table>>`.VALUE 
    FROM (`<<expression data table>>` INNER JOIN `<<subtype data table>>` 
    ON `<<expression data table>>`.PID=`<<subtype data table>>`.PID) 
                  WHERE `<<expression data table>>`.GENE = '", input$myGene, "'")
    tb3 = bq_project_query(project, sql3)
    mydata3 = data.frame(bq_table_download(tb3))
    
    # Human readable on left
    # Actual Columns name for subtype on right, user table dependent
    myfilter = input$myPicker2
    for (i in 1:length(myfilter)) {
      x = switch(myfilter[i], 
                 "Basal" = "Basal",
                 "HER2" = "Her2", 
                 "LuminalA" = "LumA",
                 "LuminalB" = "LumB",
                 "Normal" = "Normal"
      )
      myfilter[i] = x
    }
    
    # Index of which subtypes were selected
    ipam50 = which(mydata3$Call %in% myfilter)
    mdgb = mydata3[ipam50,]
    colnames(mdgb) = c("PID", "Study", "Gene_Expression_Values")
    
    # Prep Values, remove infinites or NULLs ----------------------------------
    
    mdgb$Gene_Expression_Values = as.vector(mdgb$Gene_Expression_Values)
    mdgb$Gene_Expression_Values = as.numeric(mdgb$Gene_Expression_Values)
    tcgaval$Gene_Expression_Values = as.numeric(tcgaval$Gene_Expression_Values)
    
    myi = which(is.infinite((mdgb$Gene_Expression_Values)))
    if (length(myi) > 0) {
      mdgb = mdgb[-myi,]
    }
    myi = which(is.infinite(tcgaval$Gene_Expression_Values))
    if (length(myi) > 0) {
      tcgaval = tcgaval[-myi,]
    }
    
    
    # BoxPlots ----------------------------------------------------------------
    
    
    par(mfrow=c(2,1), pin = c(3.5,2.8), mar = c(2.9,4,2.9,4))
    
    bp = boxplot(Gene_Expression_Values~Study, data = mdgb, horizontal = TRUE, main = paste0("ISPY2 Distribution of Gene Expression Values [",input$myGene, "]"), col=(c("gray","deepskyblue3")), las = 1)
    nbGroup <- length(unique(mdgb$Study))
    text( 
      y=c(1:nbGroup)-0.3, 
      x=bp$stats[nrow(bp$stats),] + .50, 
      paste("n = ",table(mdgb$Study),sep="")  
    )
    
    bp = boxplot(Gene_Expression_Values~Study, data = tcgaval, horizontal = TRUE, main = paste0("Pan-Cancer Atlas Distribution of Gene Expression Values [",input$myGene, "]"), col=(c("gray","deepskyblue3")), las = 1)
    nbGroup <- length(unique(tcgaval$Study))
    text( 
      y=c(1:nbGroup)-0.3, 
      x=bp$stats[nrow(bp$stats),] + .03, 
      paste("n = ",table(tcgaval$Study),sep="")  
    )
    
    
    # P Values Section --------------------------------------------------------
    
    # This is a formula to show how many combinations (rows) of 2, for cohorts selected
    myrows = factorial(n) / (2*factorial(n - 2))
    ptab = data.frame(matrix(NA, nrow = myrows, ncol = 2))
    colnames(ptab) = c("ISPY2", "TCGA")
    mypicks = sort(unique(mdgb$Study))
    
    # Create P values test table
    count = 1
    for (i in 1:(n-1)) {
      for (j in (i+1):n) {
        rownames(ptab)[count] = paste0(mypicks[i],"-", mypicks[j])
        pval_ISPY = t.test(as.numeric(mdgb[mdgb$Study %in% c(mypicks[i],mypicks[j]),"Gene_Expression_Values"])~mdgb[mdgb$Study %in% c(mypicks[i],mypicks[j]),"Study"])
        pval_TCGA = t.test(as.numeric(tcgaval[tcgaval$Study %in% c(mypicks[i],mypicks[j]),"Gene_Expression_Values"])~tcgaval[tcgaval$Study %in% c(mypicks[i],mypicks[j]),"Study"])
        ptab[count,1] = pval_ISPY$p.value
        ptab[count,2] = pval_TCGA$p.value
        count = count + 1
      }
    }
    
    ptab$Study = rownames(ptab)
    ptab = ptab[c(3,1,2)]
    output$tabletext = renderText(paste0("T-Test, P-Values:"))
    output$pvalues = renderTable(ptab, digits = 4)
    
    
  } # end globdata
  ) # end globdata
} # end server

