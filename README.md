# ISPY2_BigQuery_ExprData_TCGA

This is an RShiny application that displays a blox plot and T-Test P-values for the BRCA cohort, various subtypes (Basal, Luminal A, Luminal B), and on a gene of interest.  Data will be displayed for both user data and publicaly available data from the TCGA. The application pulls data from BigQuery.  User selects cohort, subtypes, and specific gene of interest. Use of application is highlighted in "PRoBE the Cloud Toolkit: Finding the Best Biomarkers of Drug Response within a Breast Cancer Clinical Trial."

User will need a Google Cloud Platform project, and relevant data uploaded to BigQuery.  Application authorizes user persmissions through email.  User will need to enter specific project variables inside << >> notation.  
