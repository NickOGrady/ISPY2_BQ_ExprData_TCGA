# an App that uses ISB-CGC uploaded databases
# This App will use PanCan(TCGA), BRCA (breast cancer) cohort
# App will display cross cohort between TCGA, and user data
# All data must be in google cloud platform, BigQuery

library(shiny)
library(shinyWidgets)

# Breast Cancer subtypes
mysubtype = c("Normal", "Basal", "HER2", "LuminalA", "LuminalB")

# Interface page
ui <- fluidPage(
  
  titlePanel(HTML(
    paste(
      h3("ISPY2 PRoBE:"),h4("Pan-Cancer Atlas Search")))),
  
  sidebarLayout(
    sidebarPanel(
      
      # Option here to add more cohorts studies, LUSC, PAAD, KIRP, LUAD etc 
      pickerInput(
        inputId = "myPicker1", 
        label = "Select TCGA Cohort", 
        choices = "BRCA",
        selected = NULL,
        options = list(
          `actions-box` = TRUE, 
          size = 10,
          `selected-text-format` = "count > 5"
        )
      ),
      pickerInput(
        inputId = "myPicker2", 
        label = "Select Cohort Subtype", 
        choices = mysubtype, 
        options = list(
          `actions-box` = TRUE, 
          size = 10,
          `selected-text-format` = "count > 5"
        ), multiple = TRUE), 
      
      textInput(inputId = "myGene", label = 'Enter Gene of Interest', value = "", width = NULL,
                placeholder = NULL),
      HTML(paste("Enter single gene, e.g. 'CD274', or 'TP53'", "\n")),
      submitButton("Submit"),
      textOutput("tabletext"),
      tableOutput(outputId = "pvalues")
    ),
    
    mainPanel( HTML(paste(h4("Results:"))),
               textOutput("graphtext"),
               plotOutput(outputId = "sqlresults"))
    
  ) # End Sidebar layout
) # End UI fluid page
