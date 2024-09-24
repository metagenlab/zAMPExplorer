

# In ui.R



library(shiny)
library(shinydashboard)
library(shinyWidgets)

library(DT)
library(ape)
library(mice)
library(writexl)
library(phyloseq)
library(microViz)
library(microbiome)
library(ggplot2)
library(cowplot)
library(circlize)
library(merTools)
library(reshape2)
library(dplyr)
library(tidyr)
library(vegan)
library(gridExtra)
library(ggpubr)
library(mia)
library(eulerr)
library(miaViz)
library(MicEco)
library(plotly)
library(metagMisc)
library(reticulate)
library(webshot2)
library(ggplotify)
library(ComplexUpset)
library(htmlwidgets)
library(colourpicker)
library(RColorBrewer)
library(radiant.data)
library(VennDiagram)
library(ComplexHeatmap)
library(microbiomeMarker)
library(microbiomeutilities)
library(DirichletMultinomial)
library(InteractiveComplexHeatmap)
library(TreeSummarizedExperiment)







shinyUI(
  dashboardPage(
    dashboardHeader(
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
    dashboardSidebar(
      sidebarMenu(
        menuItem("Upload Data", tabName = "upload", icon = icon("upload")),
        menuItem("Phyloseq components", tabName = "phyloseqcomponents", icon = icon("list-alt")),
        menuItem("Reads QC", tabName = "readsQC", icon = icon("list-alt")),
        menuItem("Taxa Overview", tabName = "taxaoverview", icon = icon("list-alt")),
        menuItem("Compositional Barplot", tabName = "barplot", icon = icon("chart-bar")),
        menuItem("Heatmap", tabName = "heatmap", icon = icon("th")),
        #menuItem("Venn Diagram", tabName = "vennDiagram", icon = icon("project-diagram")),
        menuItem("Alpha Diversity", tabName = "alphaDiversity", icon = icon("bar-chart")),
        menuItem("Beta Diversity", tabName = "betaDiversity", icon = icon("project-diagram")),
        menuItem("Community Typing (DMM)", tabName = "communityTyping", icon = icon("users")),
        menuItem("RDA Plot", tabName = "rdaPlot", icon = icon("chart-line"))
      )
    ),
    dashboardBody(
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
        # Tab for uploading the phyloseq object
        tabItem(tabName = "upload",
                fluidRow(
                  box(title = "Upload Phyloseq Object", width = 12, status = "primary",
                      fileInput("physeqFile", "Choose RDS File", accept = c(".rds"))
                  )
                )

        ),
        # Tab for phyloseq components
        tabItem(tabName = "phyloseqcomponents",
                fluidRow(
                  box(title = "Phyloseq Components", width = 12, status = "primary",
                      actionButton("show_metadata", "Metadata"),
                      actionButton("show_combined_table", "Abundance/Taxonomy Table"),
                      actionButton("show_phylogenetic_tree", "Phylogenetic Tree"),
                      actionButton("show_summary_statistics", "Summary Statistics"),
                  )
                ),
                uiOutput("dynamic_tables")  # Render the selected table or plot dynamically
        ),

        # Tab for reads QC
        tabItem(tabName = "readsQC",
                fluidRow(
                  box(title = "Reads Distribution", width = 12, status = "primary",
                      actionButton("show_reads_distribution_samples", "Read Distribution Across Samples"),
                      actionButton("show_reads_distribution_groups", "Read Distribution Across Groups"),
                      actionButton("show_number_of_reads", "Number of Reads per Sample"),
                      actionButton("explanation_button_reads_qc", label = NULL, icon = icon("info-circle"),
                                   style = "color: #fff; background-color: #007bff; border-color: #007bff;")
                  )
                ),
                fluidRow(
                  uiOutput("reads_qc_content")
                )
        ),

        # Tab for uploading taxa statistics
        tabItem(tabName = "taxaoverview",
                fluidRow(
                  box(title = "Taxa Overview", width = 12, status = "primary",
                      actionButton("show_taxa_prevalence_samples", "Taxa Prevalence Across Samples"),
                      actionButton("explain_taxa_prevalence_samples", label = NULL, icon = icon("info-circle"), class = "btn-info",
                                   style = "color: #fff; background-color: #007bff; border-color: #007bff;"),

                      actionButton("show_taxa_prevalence_groups", "Taxa Prevalence Across Groups"),
                      actionButton("explain_taxa_prevalence_groups", label = NULL, icon = icon("info-circle"), class = "btn-info",
                      style = "color: #fff; background-color: #007bff; border-color: #007bff;"),

                      actionButton("show_dominant_taxa_groups", "Dominant Taxa Across Groups"),
                      actionButton("explain_dominant_taxa_groups", label = NULL, icon = icon("info-circle"), class = "btn-info",
                                   style = "color: #fff; background-color: #007bff; border-color: #007bff;"),

                      actionButton("show_dominant_taxa_samples", "Dominant Taxa Per Sample"),
                      actionButton("explain_dominant_taxa_samples", label = NULL, icon = icon("info-circle"), class = "btn-info",
                      style = "color: #fff; background-color: #007bff; border-color: #007bff;")
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
                                  choices = c("ASV", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
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
        tabItem(tabName = "betaDiversity",
                fluidRow(
                  box(title = "Beta Diversity Analysis", width = 12, status = "primary",
                      actionButton("show_beta_explanation", label = NULL, icon = icon("info-circle"),
                                   style = "color: #31708f; border: none; background: none; font-size: 30px; padding: 10px; width: 50px; height: 50px;"),
                      selectInput("betaMethod", "Select Method", choices = c("PCoA", "PCA", "NMDS")),
                      selectInput("distanceMetric", "Select Distance Metric", choices = c("robust.aitchison", "bray", "jaccard")),
                      selectInput("normalizationMethod", "Select Normalization Method", choices = c("identity", "compositional", "binary", "clr", "hellinger", "log10", 'log10p', 'Z')),
                      selectInput("taxRank", "Select Taxonomic Rank", choices = NULL),
                      selectInput("groupingColumnBeta", "Color Column", choices = NULL),
                      selectInput("shapeColumnBeta", "Shape Column", choices = NULL),
                      actionButton("plotPCoA", "Plot PCoA"),
                      actionButton("plotPCA", "Plot PCA"),
                      actionButton("plotNMDS", "Plot NMDS"),
                      actionButton("Statistics", "Statistics"),
                      plotOutput("betaDiversityPlot", height = "800px", width = "100%"),
                      verbatimTextOutput("betaDiversityNote"),  # Add this line for the note
                      selectInput("beta_filetype", "Select file type", choices = c("html", "png", "pdf", "svg")),
                      numericInput("plot_width", "Plot Width (in inches)", value = 8, min = 12),
                      numericInput("plot_height", "Plot Height (in inches)", value = 5, min = 9),
                      downloadButton("download_PCoAPlot", "Download PCoA Plot"),
                      downloadButton("download_PCAPlot", "Download PCA Plot"),
                      downloadButton("download_NMDSPlot", "Download NMDS Plot")
                  )
                ),
                fluidRow(
                  box(title = "PERMANOVA Results", width = 12, status = "primary",
                      DT::dataTableOutput("permanovaResults"),
                      downloadButton("download_permanovaResults", "Download PERMANOVA Results")
                  )
                ),
                fluidRow(
                  box(title = "Dispersion Test Results", width = 12, status = "primary",
                      DT::dataTableOutput("dispersionTestResults"),
                      downloadButton("download_dispersionTestResults", "Download Dispersion Test Results")
                  )
                ),
                fluidRow(
                  box(title = "Distances to Centroid", width = 12, status = "primary",
                      plotOutput("centroidBoxplot")
                  )
                )
        ),
    #   fluidRow(
    #     box(
    #       title = "Beta Diversity Statistics",
    #       width = 12,
    #       status = "primary",
    #       DT::dataTableOutput("permanovaTable"),  # Output for PERMANOVA results
    #       DT::dataTableOutput("betadisperTable"),  # Output for betadisper results
    #       plotOutput("dispersionBoxplot"),  # Output for the boxplot
    #       downloadButton("download_permanovaTable", "Download PERMANOVA Table"),
    #       downloadButton("download_betadisperTable", "Download Dispersion Table")
    #     )
    #   )
    # ),
    #
    #
    #


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

        # tabItem(
        #   tabName = "rdaPlot",
        #
        #   # Settings for RDA and dbRDA
        #   fluidRow(
        #     box(
        #       title = "RDA/dbRDA Plot Settings", width = 12, status = "primary",
        #       actionButton("show_rda_explanation", label = NULL, icon = icon("info-circle"),
        #                    style = "color: #31708f; border: none; background: none; font-size: 30px; padding: 10px; width: 50px; height: 50px;"
        #       ),
        #       selectInput("taxRankRDA", "Select Taxonomic Rank", choices = NULL),
        #       selectInput("normMethod", "Select Normalization Method", choices = c("hellinger", "log10", "normalize", "relabundance",
        #                                                                            "standardize", "total", "z")),
        #       selectInput("distanceMetric", "Select Distance Metric", choices = c("bray", "jaccard")),
        #       selectInput("constraints", "Select Constraints", choices = NULL, multiple = TRUE),
        #       selectInput("conditions", "Select Conditions", choices = NULL, multiple = TRUE, selected = NULL),
        #       selectInput("colorBy", "Color By", choices = NULL),
        #       selectInput("shapeBy", "Select Shape Column", choices = NULL)
        #     )
        #   ),
        #
        #   # Buttons to generate plots
        #   fluidRow(
        #     column(
        #       width = 6,
        #       actionButton("plotRDA", "Plot RDA", class = "btn btn-primary btn-block")
        #     ),
        #     column(
        #       width = 6,
        #       actionButton("plotDbRDA", "Plot dbRDA", class = "btn btn-primary btn-block")
        #     )
        #   ),
        #
        #   # RDA Plot, Tables, and Download
        #   fluidRow(
        #     conditionalPanel(
        #       condition = "input.plotRDA > 0",
        #       box(
        #         title = "RDA Plot",
        #         plotlyOutput("rdaPlot", height = "600px", width = "100%")
        #       ),
        #       fluidRow(
        #         column(
        #           width = 6,
        #           selectInput("rda_filetype", "Select file type", choices = c("pdf", "png", "svg", "jpeg", "html"))
        #         ),
        #         column(
        #           width = 6,
        #           downloadButton("download_rdaPlot", "Download Plot")
        #         )
        #       ),
        #       fluidRow(
        #         box(
        #           title = "RDA Statistics",
        #           width = 12,
        #           status = "primary",
        #           DT::dataTableOutput("rdaInfoTable"),
        #           downloadButton("download_rdaInfoTable", "Download RDA Statistics Table"),
        #           uiOutput("downloadPhyseqUI")
        #         )
        #       )
        #     )
        #   ),

        #   # dbRDA Plot, Tables, and Download
        #   fluidRow(
        #     conditionalPanel(
        #       condition = "input.plotDbRDA > 0",
        #       box(
        #         title = "dbRDA Plot",
        #         plotOutput("dbRdaPlot", height = "600px", width = "100%")
        #       ),
        #       fluidRow(
        #         column(
        #           width = 6,
        #           selectInput("dbrda_filetype", "Select file type", choices = c("pdf", "png", "svg", "jpeg"))
        #         ),
        #         column(
        #           width = 6,
        #           downloadButton("download_dbRdaPlot", "Download dbRDA Plot")
        #         )
        #       ),
        #       fluidRow(
        #         box(
        #           title = "dbRDA Statistics",
        #           width = 12,
        #           status = "primary",
        #           DT::dataTableOutput("dbRdaInfoTable"),
        #           downloadButton("download_dbRdaInfoTable", "Download dbRDA Statistics Table"),
        #           DT::dataTableOutput("termsSignificanceTable"),
        #           downloadButton("download_termsSignificanceTable", "Download Terms Significance Table")
        #         )
        #       ),
        #       fluidRow(
        #         box(
        #           title = "Species Scores",
        #           width = 12,
        #           status = "primary",
        #           DT::dataTableOutput("speciesScoresTable"),
        #           downloadButton("download_speciesScoresTable", "Download Species Scores")
        #         )
        #       )
        #     )
        #   )
        # )
        #works

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
                      style = "color: #31708f; border: none; background: none; font-size: 30px; padding: 10px; width: 50px; height: 50px;"
                    ),
                    selectInput("taxRankRDA", "Select Taxonomic Rank", choices = NULL),
                    selectInput("normMethod", "Select Normalization Method", choices = c("hellinger", "log10", "normalize","relabundance",
                      "standardize", "total", "z")),
                    #selectInput("distanceMetric", "Select Distance Metric", choices = c("bray", "jaccard")),
                    selectInput("constraints", "Select Constraints", choices = NULL, multiple = TRUE),
                    selectInput("conditions", "Select Conditions", choices = NULL, multiple = TRUE, selected = NULL),
                    selectInput("colorBy", "Color By", choices = NULL),
                    selectInput("shapeBy", "Select Shape Column", choices = NULL),
                    actionButton("plotRDA", "Plot RDA"),
                    plotlyOutput("rdaPlot", height = "800px", width = "100%"),  # Use plotlyOutput
                    selectInput("rda_filetype", "Select file type", choices = c("pdf", "png", "svg", "jpeg", "html")),  # Allow HTML for plotly
                    downloadButton("download_rdaPlot", "Download Plot")
                  )
                ),
                #c("alr", "chi.square", "clr", "hellinger", "log", "log10", "normalize", "pa", "range", "rank", "rclr", "relabundance", "rrank",
                #"standardize", "total", "z")),
                fluidRow(
                  box(
                    title = "PERMANOVA Statistics",
                    width = 12,
                    status = "primary",
                    DT::dataTableOutput("rdaInfoTable"),  # This will render the PERMANOVA table
                    downloadButton("download_rdaInfoTable", "Download Table"),
                    uiOutput("downloadPhyseqUI")  # Download button for the updated phyloseq object with imputing NA values
                  )
                )
        )


        # # Tab for redundancy analysis
        # tabItem(tabName = "rdaPlot",
        #         fluidRow(
        #           box(title = "RDA Plot Settings", width = 12, status = "primary",
        #               actionButton("show_rda_explanation", label = NULL, icon = icon("info-circle"),
        #                            style = "color: #31708f; border: none; background: none; font-size: 30px; padding: 10px; width: 50px; height: 50px;"),
        #               selectInput("taxRankRDA", "Select Taxonomic Rank", choices = NULL),
        #               selectInput("normMethod", "Select Normalization Method", choices = c("identity", "compositional", "clr", "log10")),
        #               selectInput("constraints", "Select Constraints", choices = NULL, multiple = TRUE),
        #               selectInput("conditions", "Select Conditions", choices = NULL, multiple = TRUE, selected = NULL),
        #               selectInput("colorBy", "Color By", choices = NULL),
        #               selectInput("shapeBy", "Select Shape Column", choices = NULL),
        #               actionButton("plotRDA", "Plot RDA"),
        #               plotOutput("rdaPlot", height = "800px", width = "100%"),  # Changed from plotlyOutput to plotOutput
        #               selectInput("rda_filetype", "Select file type", choices = c("pdf", "png", "svg", "jpeg")),  # Updated file type choices
        #               downloadButton("download_rdaPlot", "Download Plot")
        #           )
        #         ),
        #         fluidRow(
        #           box(
        #             title = "PERMANOVA Statistics",
        #             width = 12,
        #             status = "primary",
        #             DT::dataTableOutput("rdaInfoTable"),  # This will render the PERMANOVA table
        #             downloadButton("download_rdaInfoTable", "Download Table"),
        #             uiOutput("downloadPhyseqUI")  # Download button for the updated phyloseq object with imputing NA values
        #           )
        #         )
        # )





      # # Tab for redundancy analysis
      #   tabItem(tabName = "rdaPlot",
      #           fluidRow(
      #             box(title = "RDA Plot Settings", width = 12, status = "primary",
      #                 selectInput("taxRankRDA", "Select Taxonomic Rank", choices = NULL),
      #                 selectInput("normMethod", "Select Normalization Method", choices = c("clr", "log10", "identity", "compositional")),
      #                 selectInput("distanceMetric", "Select Distance Metric", choices = c("bray", "jaccard")),
      #                 selectInput("explanatoryvariables", "Select explanatory variables (RDA Formula)", choices = NULL, multiple = TRUE),
      #                 selectInput("colorBy", "Color By", choices = NULL),
      #                 selectInput("shapeColumn", "Select Shape Column", choices = NULL),
      #                 actionButton("plotRDA", "Plot RDA"),
      #                 plotlyOutput("rdaPlot"),
      #                 selectInput("rda_filetype", "Select file type", choices = c("html", "pdf", "svg")),
      #                 downloadButton("download_rdaPlot", "Download Plot")
      #             )
      #           ),
      #           fluidRow(
      #             box(title = "PERMANOVA Statistics", width = 12, status = "primary",
      #                 DT::dataTableOutput("rdaInfoTable"),  # This will render the PERMANOVA table
      #                 downloadButton("download_rdaInfoTable", "Download Table")
      #             )
      #           )
      #   )
      )  # Close tabItems
    )  # Close dashboardBody
  )  # Close dashboardPage
)  # Close shinyUI
