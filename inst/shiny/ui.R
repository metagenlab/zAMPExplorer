
# In ui.R
shinyUI(
  shinydashboard::dashboardPage(
    shinydashboard::dashboardHeader(
      title = "zAMPExplorer",
      tags$li(class = "dropdown",
              tags$a(href = "javascript:void(0);",
                     icon("book"),
                     "Documentation",
                     title = "Documentation",
                     style = "color: white; font-size: 18px; padding: 10px;",
                     onclick = "window.open('https://github.com/metagenlab/zAMExplorer', '_blank');")
      )
    ),
    shinydashboard::dashboardSidebar(
      shinydashboard::sidebarMenu(
        menuItem("Upload Data", tabName = "upload", icon = icon("upload")),
        menuItem("Phyloseq components", tabName = "phyloseqcomponents", icon = icon("project-diagram")),
        menuItem("Reads QC", tabName = "readsQC", icon = icon("chart-bar")),
        menuItem("Taxa Overview", tabName = "taxaoverview", icon = icon("eye")),
        menuItem("Compositional Barplot", tabName = "barplot", icon = icon("chart-area")),
        menuItem("Heatmap", tabName = "heatmap", icon = icon("th-large")),
        menuItem("Alpha Diversity", tabName = "alphaDiversity", icon = icon("dna")),
        menuItem("Beta Diversity", tabName = "betaDiversity", icon = icon("retweet")),
        menuItem("Differential Abundance", tabName = "differentialAbundance", icon = icon("chart-pie")),
        menuItem("Community Typing (DMM)", tabName = "communityTyping", icon = icon("users")),
        menuItem("RDA Plot", tabName = "rdaPlot", icon = icon("chart-line"))
      )
    ),
    
    shinydashboard::dashboardBody(
      tags$head(tags$style(HTML('
        .dropdown .btn {
          background-color: #3498db;
          color: white;
          border: none;
          margin: 10px;
          padding: 10px 20px;
          border-radius: 5px;
          font-size: 18px;
        }
      '))),
      tabItems(
        tabItem(tabName = "upload",
                fluidRow(
                  box(
                    title = "Upload Phyloseq Object", width = 12, status = "primary",
                    fileInput("physeqFile", "Choose Phyloseq RDS File", accept = c(".rds")),
                    helpText("Upload the phyloseq object for analysis.")
                  )
                )
        ),
        tabItem(tabName = "phyloseqcomponents",
                fluidRow(
                  box(title = "Phyloseq Components", width = 12, status = "primary",
                      actionButton("show_metadata", "Metadata"),
                      actionButton("show_combined_table", "Abundance/Taxonomy Table"),
                      actionButton("show_summary_statistics", "Summary Statistics")
                  )
                ),
                uiOutput("dynamic_tables"),
                uiOutput("phylogenetic_tree_controls")
        ),
        tabItem(
          tabName = "readsQC",
          fluidRow(
            box(
              title = "Reads Distribution", width = 12, status = "primary",
              actionButton("show_reads_distribution_samples", "Reads Distribution Across Samples"),
              actionButton("show_reads_distribution_groups", "Reads Distribution Across Groups"),
              actionButton("show_number_of_reads", "Number of Reads per Sample"),
              actionButton("show_rarefaction_curves", "Rarefaction Curves"),
              actionButton(
                "explanation_button_reads_qc", 
                label = NULL, 
                icon = icon("info-circle"),
                style = "color: #31708f; border: none; background: none; font-size: 30px; padding: 10px; width: 50px; height: 50px;"
              )
            )
          ),
          # Dynamically render the content based on user selection
          fluidRow(
            uiOutput("reads_qc_content")
          )
        ),
        tabItem(tabName = "taxaoverview",
                fluidRow(
                  box(title = "Taxa Overview", width = 12, status = "primary",
                      actionButton("show_taxa_prevalence_samples", "Taxa Prevalence Across Samples"),
                      actionButton("explain_taxa_prevalence_samples", label = NULL, icon = icon("info-circle"),
                                   class = "btn-info", style = "color: #31708f; border: none; background: none; font-size: 30px; padding: 10px; width: 50px; height: 50px;"),
                      actionButton("show_taxa_prevalence_groups", "Taxa Prevalence Across Groups"),
                      actionButton("explain_taxa_prevalence_groups", label = NULL, icon = icon("info-circle"),
                                   class = "btn-info", style = "color: #31708f; border: none; background: none; font-size: 30px; padding: 10px; width: 50px; height: 50px;"),
                      actionButton("show_dominant_taxa_samples", "Dominant Taxa Per Sample"),
                      actionButton("explain_dominant_taxa_samples", label = NULL, icon = icon("info-circle"),
                                   class = "btn-info", style = "color: #31708f; border: none; background: none; font-size: 30px; padding: 10px; width: 50px; height: 50px;"),
                      actionButton("show_prevalence_abundance_plot", "Prevalence vs Abundance Plot"),
                      actionButton("explain_prevalence_abundance_plot", label = NULL, icon = icon("info-circle"),
                                   class = "btn-info", style = "color: #31708f; border: none; background: none; font-size: 30px; padding: 10px; width: 50px; height: 50px;"),
                      actionButton("show_upset_plot", "Upset Plot"),
                      actionButton("explain_upset_plot", label = NULL, icon = icon("info-circle"),
                                   class = "btn-info", style = "color: #31708f; border: none; background: none; font-size: 30px; padding: 10px; width: 50px; height: 50px;"),
                      actionButton("show_venn_diagram_plot", "Venn Diagram"),
                      actionButton("explain_venn_diagram_plot", label = NULL, icon = icon("info-circle"),
                                   class = "btn-info", style = "color: #31708f; border: none; background: none; font-size: 30px; padding: 10px; width: 50px; height: 50px;"),
                  )
                ),
                fluidRow(
                  uiOutput("taxa_overview_content")
                )
        ),

        # Tab for compositional bar plot
        tabItem(tabName = "barplot",
                fluidRow(
                  box(title = "Compositional Barplot",
                      actionButton("show_compositional_explanation", label = NULL, icon = icon("info-circle"),
                                   style = "color: #31708f; border: none; background: none; font-size: 30px; padding: 10px; width: 50px; height: 50px;"),
                      width = 12, status = "primary",
                      selectInput("hueRank", "Select Primary Rank", choices = NULL),
                      selectInput("shadeRank", "Select Lower Rank", choices = NULL),
                      selectInput("facetBy", "Facet by", choices = NULL),
                      sliderInput("nHues", "Number of Primary Taxa", min = 1, max = 10, value = 3),
                      sliderInput("nShades", "Number of Lower Taxa in the main rank", min = 1, max = 10, value = 4)
                  )
                ),
                fluidRow(
                  box(title = "Compositional Barplot", width = 12, status = "primary",
                      plotlyOutput("compositionalBarplot", height = "1000px", width = "100%"),
                      selectInput("filetype", "Select file type", choices = c("png", "pdf", "jpeg", "svg")),
                      numericInput("plot_width", "Plot Width (in pixels)", value = 1200, min = 400),
                      numericInput("plot_height", "Plot Height (in pixels)", value = 800, min = 400),
                      downloadButton("download_compositionalBarplot", "Download Plot")
                  )
                )
        ),

        # Tab for relative abundance heatmap
        tabItem(tabName = "heatmap",
                fluidRow(
                  box(title = "Relative Abundance Heatmap",
                      actionButton("show_heatmap_explanation", label = NULL, icon = icon("info-circle"),
                                   style = "color: #31708f; border: none; background: none; font-size: 30px; padding: 10px; width: 50px; height: 50px;"),
                      width = 12, status = "primary",

                      selectInput("normalizationMethod", "Select Normalization Method", choices = c("compositional", "Z", "log10", "log10p", "hellinger", "identity", "clr")),
                      selectInput("heatmapRank", "Select Taxonomic Rank",
                                  choices = c("ASV", "Phylum", "Family", "Genus", "Species"),
                                  selected = "ASV"),
                      selectInput("annotationColumn1", "Select First Annotation Column", choices = NULL),
                      selectInput("annotationColumn2", "Select Second Annotation Column", choices = NULL),
                      checkboxInput("clusterRows", "Cluster Rows", value = TRUE),
                      checkboxInput("clusterColumns", "Cluster Columns", value = TRUE),
                      numericInput("topTaxa", "Number of Top Taxa", value = 10, min = 1),
                      actionButton("plotHeatmap", "Plot Heatmap"),
                      InteractiveComplexHeatmapOutput("relativeAbundanceHeatmap")
                  )
                )
        ),
        # Tab for alpha diversity
        tabItem(tabName = "alphaDiversity",
                fluidRow(
                  box(title = "Alpha Diversity Analysis", width = 12, status = "primary",
                      actionButton("show_alpha_explanation", label = NULL, icon = icon("info-circle"),
                                   style = "color: #31708f; border: none; background: none; font-size: 30px; padding: 10px; width: 50px; height: 50px;"),
                      selectInput("alphaMetric", "Select Diversity Metric", choices = NULL, multiple = TRUE),
                      selectInput("alphaGroupingColumn", "Grouping Column", choices = NULL),
                      actionButton("plotAlpha", "Plot Alpha Diversity"),
                      actionButton("showStats", "Show Stats"),
                      plotlyOutput("alphaDiversityPlot", height = "900px"), #, width = "100%"
                      DTOutput("statsTable"),
                      verbatimTextOutput("alphaDiversityNote"),
                      selectInput("alpha_filetype", "Select file type", choices = c("html", "png", "pdf", "svg")),
                      numericInput("plot_width", "Plot Width (in pixels)", value = 8, min = 12),
                      numericInput("plot_height", "Plot Height (in pixels)", value = 5, min = 9),
                      downloadButton("download_alphaDiversityPlot", "Download Plot"),
                      downloadButton("downloadStats", "Download Stats")
                  )
                )
        ),
        
        # Tab for beta diversity
        tabItem(
          tabName = "betaDiversity",
          fluidRow(
            box(
              title = "Beta Diversity Analysis", 
              width = 12, 
              status = "primary",
              actionButton("show_beta_explanation", label = NULL, icon = icon("info-circle"),
                           style = "color: #31708f; border: none; background: none; font-size: 30px; padding: 10px; width: 50px; height: 50px;"),
              
              # Create tabset panel for Plotting and Statistics tabs
              tabsetPanel(
                id = "betaDiversityTabs",
                
                # Plotting Tab
                tabPanel(
                  title = "Plotting",
                  selectInput("betaMethod", "Select Method", choices = c("PCoA", "PCA", "NMDS")),
                  selectInput("distanceMetric", "Select Distance Metric", choices = c("robust.aitchison", "bray", "jaccard")),
                  selectInput("normalizationMethod", "Select Normalization Method", choices = c("identity", "compositional", "binary", "clr", "hellinger", "log10", 'log10p', 'Z')),
                  selectInput("taxRank", "Select Taxonomic Rank", choices = NULL),
                  selectInput("groupingColumnBeta_plotting", "Grouping Column (Plotting)", choices = NULL),
                  #selectInput("groupingColumnBeta", "Grouping Column", choices = NULL),
                  selectInput("shapeColumnBeta", "Shape Column", choices = NULL),
                  
                  # Plot buttons
                  actionButton("plotPCoA", "Plot PCoA"),
                  actionButton("plotPCA", "Plot PCA"),
                  actionButton("plotNMDS", "Plot NMDS"),
                  
                  # Plot output
                  plotOutput("betaDiversityPlot", height = "800px", width = "100%"),
                  verbatimTextOutput("betaDiversityNote"),  # For displaying any notes or errors
                  
                  # Plot download options
                  selectInput("beta_filetype", "Select file type", choices = c("html", "png", "pdf", "svg")),
                  numericInput("plot_width", "Plot Width (in inches)", value = 8, min = 12),
                  numericInput("plot_height", "Plot Height (in inches)", value = 5, min = 9),
                  downloadButton("download_PCoAPlot", "Download PCoA Plot"),
                  downloadButton("download_PCAPlot", "Download PCA Plot"),
                  downloadButton("download_NMDSPlot", "Download NMDS Plot")
                ),
        
        # Statistics Tab
        tabPanel(
          title = "Statistics",
          
          # Variables selection for PERMANOVA
          selectInput("groupingColumnBeta_statistics", "Grouping Column (Statistics)", choices = NULL),
          # selectInput("groupingColumnBeta", "Grouping Column", choices = NULL),  # Choices will be updated dynamically based on metadata columns
          # Additional settings for PERMANOVA
          selectInput("distanceMetric", "Distance Metric", choices = c("bray", "jaccard")),
          selectInput("strata", "Strata (for stratified analysis)", choices = NULL, selected = "None"),
          
          # Additional settings
          selectInput("significance_by", "Significance By", choices = c("term", "margin", "onedf")),
          selectInput("taxRank", "Select Taxonomic Rank", choices = c("Phylum","Family", "Genus", "Species"), selected = "Genus"),
          
          # Statistics buttons
          actionButton("calculatePERMANOVA", "Calculate PERMANOVA"),
          actionButton("calculateDispersionTest", "Calculate Dispersion Test"),
          actionButton("plotCentroids", "Plot Centroids"),
          
          # Dynamic output for PERMANOVA Test results
          uiOutput("PERMANOVAResultsUI"),
          
          # Dynamic output for Dispersion Test results
          uiOutput("dispersionTestResultsUI"),
          
          # Dynamic output for Centroid Plot
          uiOutput("centroidPlotUI")
        )
              )
            )
          )
        ),
        tabItem(
          tabName = "differentialAbundance",
          fluidRow(
            box(
              title = "Differential Abundance Analysis",
              width = 12,
              status = "primary",
              actionButton("show_maaslin2_explanation", label = NULL, icon = icon("info-circle"),
                           style = "color: #31708f; border: none; background: none; font-size: 30px; padding: 10px; width: 50px; height: 50px;"),
              tabsetPanel(
                id = "daTabs",
                
                # Maaslin2 Tab
                tabPanel(
                  title = "Maaslin2",
                  fluidRow(
                    # First column: Output folder and core inputs
                    column(
                      6, # Specifies the width of the column in a 12-column grid layout
                      shinyDirButton("select_folder", "Select Output Folder", "Please select a folder", icon = icon("folder")), # Folder selection button
                      verbatimTextOutput("selected_folder"),  # Displays selected folder path
                      textInput("Maaslin2_output", "Output Directory", placeholder = "Output folder path will appear here"),
                      selectizeInput("fixed_variables", "Fixed Effects", choices = NULL, multiple = TRUE), # Fixed effects input
                      selectizeInput("taxonomy_level", "Select Taxonomy Level", choices = c("Phylum", "Family", "Genus", "Species"), selected = "Genus"),
                      selectizeInput("random_variables", "Random Effects", choices = NULL, multiple = TRUE), # Random effects input
                      fluidRow(
                        column(
                          6,
                          selectInput("reference_variable", "Reference Variable", choices = NULL) # New dropdown for reference variable
                        ),
                        column(
                          6,
                          selectInput("reference_level", "Reference Level", choices = NULL) # New dropdown for reference level
                        )
                      ),
                      numericInput("min_abundance", "Minimum Abundance (Relative)", value = 0.0001, min = 0, step = 0.0001), # Minimum abundance
                      numericInput("min_prevalence", "Minimum Prevalence (Proportion)", value = 0.1, min = 0, max = 1, step = 0.01), # Minimum prevalence
                      numericInput("max_significance", "Maximum Significance (q-value)", value = 0.05, min = 0, max = 1, step = 0.01) # Maximum significance
                    ),
                    # Second column: Maaslin2 advanced settings
                    column(
                      6,
                      selectInput(
                        "normalization", "Normalization Method",
                        choices = c("NONE", "CSS", "TSS", "CLR", "TMM"), selected = "CSS"
                      ),
                      selectInput(
                        "transform", "Transformation",
                        choices = c("NONE", "LOG", "LOGIT", "AST"), selected = "NONE"
                      ),
                      selectInput(
                        "analysis_method", "Analysis Method",
                        choices = c("NEGBIN", "ZINB", "LM", "CPLM", "ZICP"), selected = "NEGBIN"
                      ),
                      selectInput(
                        "correction", "Correction Method",
                        choices = c("BH", "BY", "holm", "hochberg"), selected = "BH"
                      ),
                      checkboxInput("standardize", "Standardize Continuous Metadata", value = TRUE), # Standardization option
                      checkboxInput("plot_heatmap", "Plot Heatmap", value = TRUE), # Heatmap plotting
                      numericInput("heatmap_first_n", "Top N Features in Heatmap", value = 50, min = 1), # Top features in heatmap
                      checkboxInput("plot_scatter", "Plot Scatter Plots", value = TRUE), # Scatter plot option
                      numericInput("cores", "Number of Cores", value = 1, min = 1) # Number of cores
                    )
                  ),
                  actionButton("run_maaslin2", "Run Maaslin2"), # Button to run Maaslin2
                  DT::dataTableOutput("maaslin2_results") # Display results in a DataTable
                )
              ) # Close tabsetPanel
            ) # Close box
          ) # Close fluidRow
        ), # Close tabItem for differentialAbundance
        
                
        # Tab for community typing
        tabItem(tabName = "communityTyping",
                fluidRow(
                  box(title = "Dirichlet Multinomial Mixture (DMM) Model",
                      actionButton("show_community_typing_explanation", label = NULL, icon = icon("info-circle"),
                                   style = "color: #31708f; border: none; background: none; font-size: 30px; padding: 10px; width: 50px; height: 50px;"),
                      status = "primary", width = 12,
                      selectInput("taxRankDMM", "Select Taxonomic Rank", choices = NULL),
                      numericInput("detectionThreshold", "Detection Threshold (%)", value = 0.1, min = 0.01, max = 100, step = 0.01),
                      numericInput("prevalenceThreshold", "Prevalence Threshold (%)", value = 50, min = 0.01, max = 100, step = 0.01),
                      numericInput("numComponents", "Number of Dirichlet Components", value = 2, min = 1),
                      actionButton("runDMM", "Run DMM Analysis")
                  )
                ),

                fluidRow(
                  # Add buttons for each result
                  box(width = 12,
                      actionButton("showModelFit", "Show Model Fit Plot"),
                      actionButton("showMixtureParams", "Show Mixture Parameters"),
                      actionButton("showSampleAssignments", "Show Sample Assignments"),
                      actionButton("showDriverPlots", "Show Driver Plots")
                  )
                ),

                # Create space for each result to be displayed conditionally
                fluidRow(
                  conditionalPanel(
                    condition = "input.showModelFit == 1",
                    box(title = "Model Fit", width = 12, plotOutput("modelFitPlot"))
                  ),
                  conditionalPanel(
                    condition = "input.showMixtureParams == 1",
                    box(title = "Mixture Parameters", width = 12, verbatimTextOutput("mixtureParams"))
                  ),
                  conditionalPanel(
                    condition = "input.showSampleAssignments == 1",
                    box(title = "Sample Assignments", width = 12, verbatimTextOutput("sampleAssignments"))
                  ),
                  conditionalPanel(
                    condition = "input.showDriverPlots == 1",
                    box(title = "Driver Plots", width = 12, uiOutput("driverPlotsUI"))
                  )
                ),

                # Add download buttons
                fluidRow(
                  box(title = "Download Results", width = 12, status = "primary",
                      downloadButton("downloadModelFitPlot", "Download Model Fit Plot"),
                      downloadButton("downloadMixtureParams", "Download Mixture Parameters"),
                      DTOutput("sampleAssignments"),
                      downloadButton("downloadSampleAssignments", "Download Sample Assignments"),
                      downloadButton("downloadDriverPlots", "Download Driver Plots"),
                      downloadButton("downloadUpdatedPhyseq", "Download Updated Phyloseq Object"),
                      downloadButton("downloadMetadata", "Download Metadata")
                  )
                )
        ),

        tabItem(tabName = "rdaPlot",
                fluidRow(
                  box(
                    title = "RDA Plot Settings",
                    width = 12,
                    status = "primary",
                    actionButton(
                      "show_rda_explanation",
                      label = NULL,
                      icon = icon("info-circle"),
                      style = "color: #31708f; border: none; background: none; font-size: 30px; padding: 10px; width: 50px; height: 50px;"),
                    selectInput("taxRankRDA", "Select Taxonomic Rank", choices = NULL),
                    selectInput("normMethod", "Select Normalization Method", choices = c("compositional", "hellinger", "log10", "log10p",
                                                                                         "identity", "clr")),
                    selectInput("RDAdistanceMetric", "Distance Metric", choices = c("bray", "jaccard")),
                    selectInput("constraints", "Select Constraints", choices = NULL, multiple = TRUE),
                    selectInput("conditions", "Select Conditions", choices = NULL, multiple = TRUE, selected = NULL),
                    selectInput("colorBy", "Color By", choices = NULL),
                    selectInput("shapeBy", "Select Shape Column", choices = NULL),
                    actionButton("plotRDA", "Plot RDA"),
                    plotOutput("rdaPlot", height = "800px", width = "100%"),  # Use plotOutput for ggplot
                    selectInput("rda_filetype", "Select file type", choices = c("pdf", "png", "svg", "jpeg")),
                    downloadButton("download_rdaPlot", "Download Plot")
                  )
                ),
                fluidRow(
                  box(
                    title = "CAP Analysis Statistics",
                    width = 12,
                    status = "primary",
                    tabsetPanel(
                      tabPanel("Summary", DT::dataTableOutput("inertiaTable"), downloadButton("downloadInertia", "Download Summary")),
                      tabPanel("Eigenvalues", DT::dataTableOutput("eigenTable"), downloadButton("downloadEigen", "Download Eigenvalues")),
                      tabPanel("Adjusted R²", DT::dataTableOutput("r2Table"), downloadButton("downloadR2", "Download Adjusted R²")),
                      tabPanel("Regression Coefficients", DT::dataTableOutput("coeffTable"), downloadButton("downloadCoeff", "Download Coefficients")),
                      tabPanel("Site Scores", DT::dataTableOutput("siteScoresTable"), downloadButton("downloadSiteScores", "Download Site Scores")),
                      tabPanel("Species Scores", DT::dataTableOutput("speciesScoresTable"), downloadButton("downloadSpeciesScores", "Download Species Scores")),
                      tabPanel("ANOVA Overall", DT::dataTableOutput("anovaOverallTable"), downloadButton("downloadAnovaOverall", "Download Overall ANOVA")),
                      tabPanel("ANOVA Terms", DT::dataTableOutput("anovaTermsTable"), downloadButton("downloadAnovaTerms", "Download ANOVA Terms")),
                      tabPanel("ANOVA Axes", DT::dataTableOutput("anovaAxesTable"), downloadButton("downloadAnovaAxes", "Download ANOVA Axes"))
                    )
                  )
                )
        )
      ) # End of tabItems
    ) # End of dashboardBody
  ) # End of dashboardPage
) # End of shinyUI

