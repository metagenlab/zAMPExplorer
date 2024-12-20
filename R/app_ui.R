#' Application User Interface
#'
#' Defines the user interface for the Shiny application.
#'
#' @param request Internal parameter for `{shiny}`. DO NOT REMOVE.
#' @return A Shiny UI definition created using `{shinydashboard}`.
#' @import shiny
#' @importFrom shinydashboard dashboardPage
#' @noRd

app_ui <- function(request) {
  shinydashboard::dashboardPage(
    shinydashboard::dashboardHeader(
      title = "zAMPExplorer",
      tags$li(
        class = "dropdown",
        tags$a(
          href = "javascript:void(0);",
          icon("book"),
          "Documentation",
          title = "Documentation",
          style = "color: white; font-size: 18px; padding: 10px;",
          onclick = "window.open('https://github.com/metagenlab/zAMPExplorer', '_blank');"
        )
      )
    ),
    shinydashboard::dashboardSidebar(
      shinydashboard::sidebarMenu(
        shinydashboard::menuItem("Upload Data", tabName = "upload", icon = icon("upload")),
        shinydashboard::menuItem("Phyloseq components", tabName = "phyloseqcomponents", icon = icon("project-diagram")),
        shinydashboard::menuItem("Reads QC", tabName = "readsQC", icon = icon("chart-bar")),
        shinydashboard::menuItem("Taxa Overview", tabName = "taxaoverview", icon = icon("eye")),
        shinydashboard::menuItem("Compositional Barplot", tabName = "barplot", icon = icon("chart-area")),
        shinydashboard::menuItem("Alpha Diversity", tabName = "alphaDiversity", icon = icon("dna")),
        shinydashboard::menuItem("Beta Diversity", tabName = "betaDiversity", icon = icon("retweet")),
        shinydashboard::menuItem("Community Typing (DMM)", tabName = "communityTyping", icon = icon("users")),
        shinydashboard::menuItem("RDA Plot", tabName = "rdaPlot", icon = icon("chart-line"))
      )
    ),
    shinydashboard::dashboardBody(
      tags$head(tags$style(HTML("
        .dropdown .btn {
          background-color: #3498db;
          color: white;
          border: none;
          margin: 10px;
          padding: 10px 20px;
          border-radius: 5px;
          font-size: 18px;
        }
      "))),
      shinydashboard::tabItems(
        shinydashboard::tabItem(
          tabName = "upload",
          fluidRow(
            shinydashboard::box(
              title = "Upload Phyloseq Object", width = 12, status = "primary",
              fileInput("physeqFile", "Choose Phyloseq RDS File", accept = c(".rds")),
              helpText("Upload the phyloseq object for analysis.")
            )
          )
        ),
        shinydashboard::tabItem(
          tabName = "phyloseqcomponents",
          fluidRow(
            shinydashboard::box(
              title = "Phyloseq Components", width = 12, status = "primary",
              actionButton("show_metadata", "Metadata"),
              actionButton("show_combined_table", "Abundance/Taxonomy Table"),
              actionButton("show_summary_statistics", "Summary Statistics")
            )
          ),
          uiOutput("dynamic_tables"),
          uiOutput("phylogenetic_tree_controls")
        ),
        shinydashboard::tabItem(
          tabName = "readsQC",
          fluidRow(
            shinydashboard::box(
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
        shinydashboard::tabItem(
          tabName = "taxaoverview",
          fluidRow(
            shinydashboard::box(
              title = "Taxa Overview", width = 12, status = "primary",
              actionButton("show_taxa_prevalence_samples", "Taxa Prevalence Across Samples"),
              actionButton("explain_taxa_prevalence_samples",
                label = NULL, icon = icon("info-circle"),
                class = "btn-info", style = "color: #31708f; border: none; background: none; font-size: 30px; padding: 10px; width: 50px; height: 50px;"
              ),
              actionButton("show_taxa_prevalence_groups", "Taxa Prevalence Across Groups"),
              actionButton("explain_taxa_prevalence_groups",
                label = NULL, icon = icon("info-circle"),
                class = "btn-info", style = "color: #31708f; border: none; background: none; font-size: 30px; padding: 10px; width: 50px; height: 50px;"
              ),
              actionButton("show_prevalence_abundance_plot", "Prevalence vs Abundance Plot"),
              actionButton("explain_prevalence_abundance_plot",
                label = NULL, icon = icon("info-circle"),
                class = "btn-info", style = "color: #31708f; border: none; background: none; font-size: 30px; padding: 10px; width: 50px; height: 50px;"
              ),
              actionButton("show_upset_plot", "Upset Plot"),
              actionButton("explain_upset_plot",
                label = NULL, icon = icon("info-circle"),
                class = "btn-info", style = "color: #31708f; border: none; background: none; font-size: 30px; padding: 10px; width: 50px; height: 50px;"
              ),
              actionButton("show_venn_diagram_plot", "Venn Diagram"),
              actionButton("explain_venn_diagram_plot",
                label = NULL, icon = icon("info-circle"),
                class = "btn-info", style = "color: #31708f; border: none; background: none; font-size: 30px; padding: 10px; width: 50px; height: 50px;"
              ),
            )
          ),
          fluidRow(
            uiOutput("taxa_overview_content")
          )
        ),

        # Tab for compositional bar plot
        shinydashboard::tabItem(
          tabName = "barplot",
          fluidRow(
            shinydashboard::box(
              title = "Compositional Barplot",
              actionButton("show_compositional_explanation",
                label = NULL, icon = icon("info-circle"),
                style = "color: #31708f; border: none; background: none; font-size: 30px; padding: 10px; width: 50px; height: 50px;"
              ),
              width = 12, status = "primary",
              selectInput("TaxRankcomp", "Select Rank", choices = NULL),
              numericInput("detectionThreshold", "Detection Threshold (%)", value = 0.1, min = 0.01, max = 100, step = 0.01),
              numericInput("prevalenceThreshold", "Prevalence Threshold (%)", value = 50, min = 0.01, max = 100, step = 0.01),
            )
          ),
          fluidRow(
            shinydashboard::box(
              title = "Compositional Barplot", width = 12, status = "primary",
              plotlyOutput("compositionalBarplot", height = "1000px", width = "100%"),
              selectInput("filetype", "Select file type", choices = c("png", "pdf", "jpeg", "svg")),
              numericInput("plot_width", "Plot Width (in pixels)", value = 1200, min = 400),
              numericInput("plot_height", "Plot Height (in pixels)", value = 800, min = 400),
              downloadButton("download_compositionalBarplot", "Download Plot")
            )
          )
        ),
        # Tab for alpha diversity
        shinydashboard::tabItem(
          tabName = "alphaDiversity",
          fluidRow(
            shinydashboard::box(
              title = "Alpha Diversity Analysis", width = 12, status = "primary",
              actionButton("show_alpha_explanation",
                label = NULL, icon = icon("info-circle"),
                style = "color: #31708f; border: none; background: none; font-size: 30px; padding: 10px; width: 50px; height: 50px;"
              ),
              selectInput("alphaMetric", "Select Diversity Metric", choices = NULL, multiple = TRUE),
              selectInput("alphaGroupingColumn", "Grouping Column", choices = NULL),
              actionButton("plotAlpha", "Plot Alpha Diversity"),
              actionButton("showStats", "Show Stats"),
              plotlyOutput("alphaDiversityPlot", height = "900px"), # , width = "100%"
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
        shinydashboard::tabItem(
          tabName = "betaDiversity",
          fluidRow(
            shinydashboard::box(
              title = "Beta Diversity Analysis",
              width = 12,
              status = "primary",
              actionButton("show_beta_explanation",
                label = NULL, icon = icon("info-circle"),
                style = "color: #31708f; border: none; background: none; font-size: 30px; padding: 10px; width: 50px; height: 50px;"
              ),

              # Create tabset panel for Plotting and Statistics tabs
              tabsetPanel(
                id = "betaDiversityTabs",

                # Plotting Tab
                tabPanel(
                  title = "Plotting",
                  selectInput("betaMethod", "Select Method", choices = c("PCoA", "PCA", "NMDS")),
                  selectInput("distanceMetric", "Select Distance Metric", choices = c("aitchison", "bray", "jaccard")),
                  selectInput("normalizationMethod", "Select Normalization Method", choices = c("identity", "compositional", "binary", "clr", "hellinger", "log10", "log10p", "Z")),
                  selectInput("taxRank", "Select Taxonomic Rank", choices = NULL),
                  selectInput("groupingColumnBeta_plotting", "Grouping Column (Plotting)", choices = NULL),
                  # selectInput("groupingColumnBeta", "Grouping Column", choices = NULL),
                  selectInput("shapeColumnBeta", "Shape Column", choices = NULL),

                  # Plot buttons
                  actionButton("plotPCoA", "Plot PCoA"),
                  actionButton("plotPCA", "Plot PCA"),
                  actionButton("plotNMDS", "Plot NMDS"),

                  # Plot output
                  plotOutput("betaDiversityPlot", height = "800px", width = "100%"),
                  verbatimTextOutput("betaDiversityNote"), # For displaying any notes or errors

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
                  selectInput("taxRank", "Select Taxonomic Rank", choices = c("Phylum", "Family", "Genus", "Species"), selected = "Genus"),

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

        # Tab for community typing
        shinydashboard::tabItem(
          tabName = "communityTyping",
          fluidRow(
            shinydashboard::box(
              title = "Dirichlet Multinomial Mixture (DMM) Model",
              actionButton("show_community_typing_explanation",
                label = NULL, icon = icon("info-circle"),
                style = "color: #31708f; border: none; background: none; font-size: 30px; padding: 10px; width: 50px; height: 50px;"
              ),
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
            shinydashboard::box(
              width = 12,
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
              shinydashboard::box(title = "Model Fit", width = 12, plotOutput("modelFitPlot"))
            ),
            conditionalPanel(
              condition = "input.showMixtureParams == 1",
              shinydashboard::box(title = "Mixture Parameters", width = 12, verbatimTextOutput("mixtureParams"))
            ),
            conditionalPanel(
              condition = "input.showSampleAssignments == 1",
              shinydashboard::box(title = "Sample Assignments", width = 12, verbatimTextOutput("sampleAssignments"))
            ),
            conditionalPanel(
              condition = "input.showDriverPlots == 1",
              shinydashboard::box(title = "Driver Plots", width = 12, uiOutput("driverPlotsUI"))
            )
          ),

          # Add download buttons
          fluidRow(
            shinydashboard::box(
              title = "Download Results", width = 12, status = "primary",
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
        shinydashboard::tabItem(
          tabName = "rdaPlot",
          fluidRow(
            shinydashboard::box(
              title = "RDA Plot Settings",
              width = 12,
              status = "primary",
              actionButton(
                "show_rda_explanation",
                label = NULL,
                icon = icon("info-circle"),
                style = "color: #31708f; border: none; background: none; font-size: 30px; padding: 10px; width: 50px; height: 50px;"
              ),
              selectInput("taxRankRDA", "Select Taxonomic Rank", choices = NULL),
              selectInput("normMethod", "Select Normalization Method", choices = c(
                "compositional", "hellinger", "log10", "log10p",
                "identity", "clr"
              )),
              selectInput("RDAdistanceMetric", "Distance Metric", choices = c("bray", "jaccard")),
              selectInput("constraints", "Select Constraints", choices = NULL, multiple = TRUE),
              selectInput("conditions", "Select Conditions", choices = NULL, multiple = TRUE, selected = NULL),
              selectInput("colorBy", "Color By", choices = NULL),
              selectInput("shapeBy", "Select Shape Column", choices = NULL),
              actionButton("plotRDA", "Plot RDA"),
              plotOutput("rdaPlot", height = "800px", width = "100%"), # Use plotOutput for ggplot
              selectInput("rda_filetype", "Select file type", choices = c("pdf", "png", "svg", "jpeg")),
              downloadButton("download_rdaPlot", "Download Plot")
            )
          ),
          fluidRow(
            shinydashboard::box(
              title = "CAP Analysis Statistics",
              width = 12,
              status = "primary",
              tabsetPanel(
                tabPanel("Summary", DTOutput("inertiaTable"), downloadButton("downloadInertia", "Download Summary")),
                tabPanel("Eigenvalues", DTOutput("eigenTable"), downloadButton("downloadEigen", "Download Eigenvalues")),
                tabPanel("Adjusted R^2", DTOutput("r2Table"), downloadButton("downloadR2", "Download Adjusted R^2")),
                tabPanel("Regression Coefficients", DTOutput("coeffTable"), downloadButton("downloadCoeff", "Download Coefficients")),
                tabPanel("Site Scores", DTOutput("siteScoresTable"), downloadButton("downloadSiteScores", "Download Site Scores")),
                tabPanel("Species Scores", DTOutput("speciesScoresTable"), downloadButton("downloadSpeciesScores", "Download Species Scores")),
                tabPanel("ANOVA Overall", DTOutput("anovaOverallTable"), downloadButton("downloadAnovaOverall", "Download Overall ANOVA")),
                tabPanel("ANOVA Terms", DTOutput("anovaTermsTable"), downloadButton("downloadAnovaTerms", "Download ANOVA Terms")),
                tabPanel("ANOVA Axes", DTOutput("anovaAxesTable"), downloadButton("downloadAnovaAxes", "Download ANOVA Axes"))
              )
            )
          )
        )
      ) # End of tabItems
    ) # End of dashboardBody
  ) # End of dashboardPage
} # End of app_ui function
