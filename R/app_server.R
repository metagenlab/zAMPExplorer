#' Application Server Logic
#'
#' Defines the server-side logic for the Shiny application, including reactive expressions,
#' outputs, and event observers.
#'
#' @param input,output,session Default Shiny server parameters.
#' @return None. This function is called for its side effects.
#' @import shiny
#' @importFrom ggplot2 ggplot aes geom_bar geom_histogram geom_boxplot theme
#' @importFrom utils combn write.csv write.table
#' @importFrom DT datatable DTOutput
#' @importFrom plotly plotlyOutput renderPlotly
#' @importFrom InteractiveComplexHeatmap InteractiveComplexHeatmapOutput
#' @importFrom MicrobiotaProcess ggrarecurve
#' @importFrom RColorBrewer brewer.pal
#' @importFrom phyloseq distance otu_table sample_data rank_names sample_variables taxa_sums nsamples
#' @importFrom phyloseq tax_table prune_samples sample_sums sample_names sample_data
#' @importFrom dplyr filter mutate select
#' @importFrom tidyr pivot_longer pivot_wider
#' @importFrom vegan capscale adonis vegdist
#' @importFrom writexl write_xlsx
#' @importFrom microbiome transform aggregate_taxa core core_members
#' @importFrom grDevices colorRampPalette dev.off pdf png svg
#' @importFrom graphics box layout legend lines
#' @importFrom stats anova as.formula coef fitted median p.adjust quantile setNames update wilcox.test
#' @importFrom plotly config subplot
#' @importFrom shinyFiles shinyFilesButton shinyDirButton
#' @importFrom shinyWidgets pickerInput
#' @importFrom magrittr %>%
#' @importFrom reshape2 melt

#' @noRd



app_server <- function(input, output, session) {

  # Reactive value to store the uploaded phyloseq object
  physeq <- reactive({
    req(input$physeqFile)
    pseq <- readRDS(input$physeqFile$datapath)

    # Preprocess the phyloseq object
    pseq <- phyloseq::prune_samples(sample_sums(pseq) > 0, pseq)
    pseq <- phyloseq::prune_taxa(taxa_sums(pseq) > 0, pseq)

    # Drop existing diversity metrics from sample data
    drop <- c("Observed", "Chao1", "se.chao1", "ACE", "se.ACE", "Shannon", "Simpson", "InvSimpson", "Fisher", "Observed_min_1")
    phyloseq::sample_data(pseq) <- phyloseq::sample_data(pseq)[, !(names(phyloseq::sample_data(pseq)) %in% drop)]

    return(pseq)  # Return the processed phyloseq object
  })

  # Debug output to track which file is being used
  observe({
    if (!is.null(physeq())) {
      cat("Phyloseq object has been uploaded and processed.\n")
    } else {
      cat("No phyloseq object uploaded yet.\n")
    }
  })

  if (!requireNamespace("microViz", quietly = TRUE)) {
    stop("The package 'microViz' is required for this functionality. Please install it using install.packages('microViz').")
  }



  # Update downstream components based on the uploaded phyloseq object
  observeEvent(physeq(), {
    req(physeq())  # Ensure physeq is available before proceeding

    showNotification("Phyloseq object successfully uploaded and processed.", type = "message")

    # Update downstream components such as selectInputs
    updateSelectInput(session, "phylum", choices = unique(phyloseq::tax_table(physeq())[, "Phylum"]))
    updateSelectInput(session,"label_tips", "Label Tips By (taxonomy):", choices = rank_names(physeq()))
    updateSelectInput(session,"color_by", "Choose Color By (metadata):", choices = sample_variables(physeq()))
    updateSelectInput(session,"shape_by", "Choose Shape By (taxonomy):", choices = rank_names(physeq()))
    updateSelectInput(session, "color_by", choices = sample_variables(physeq()))
    updateSelectInput(session, "taxaLevel", choices = rank_names(physeq()))
    updateSelectInput(session, "sampleGroup", choices = sample_variables(physeq()))
    updateSelectInput(session, "hueRank", choices = rank_names(physeq()))
    updateSelectInput(session, "shadeRank", choices = rank_names(physeq()))
    updateSelectInput(session, "facetBy", choices = c(sample_variables(physeq())))
    updateSelectInput(session, "heatmapRank", choices = rank_names(physeq()))
    updateSelectInput(session, "annotationColumn1", choices = sample_variables(physeq()))
    updateSelectInput(session, "annotationColumn2", choices = sample_variables(physeq()))
    updateSelectInput(session, "alphaMetric", choices = colnames(microbiome::alpha(physeq(), index = "all")))
    updateSelectInput(session, "alphaGroupingColumn", choices = sample_variables(physeq()))
    updateSelectInput(session, "groupingColumn", choices = sample_variables(physeq()))
    updateSelectInput(session, "groupingColumnBeta", choices = sample_variables(physeq()))
    updateSelectInput(session, "shapeColumnBeta", choices = sample_variables(physeq()))
    updateSelectInput(session, "taxRank", choices = rank_names(physeq()))
    updateSelectInput(session, "dmmTaxonomicRank", choices = c("ASV", rank_names(physeq())))
    updateSelectInput(session, "rdaTaxonomicRank", choices = c("ASV", rank_names(physeq())))
    updateSelectInput(session, "shapeColumn", choices = sample_variables(physeq()))
    updateSelectInput(session, "explanatoryvariables", choices = sample_variables(physeq()))
    updateSelectInput(session, "colorBy", choices = sample_variables(physeq()))
  })




  #generate_download_handler Function
  generate_download_handler <- function(render_plot, filename_prefix, vwidth = 1200, vheight = 800, delay = 0.1) {
    downloadHandler(
      filename = function() {
        paste0(filename_prefix, Sys.Date(), ".", input$filetype)
      },
      content = function(file) {
        # Render the plotly object
        p <- render_plot()

        # Save the plotly plot to an HTML file first
        tempfile <- tempfile(fileext = ".html")
        htmlwidgets::saveWidget(as_widget(p), tempfile)

        # Use webshot to convert the HTML to the selected format
        webshot::webshot(tempfile, file = file, vwidth = vwidth, vheight = vheight, delay = delay)
      }
    )
  }


  ######### phyloseq components


  #value to track which button is pressed
  output_view <- reactiveVal("")

  # Observers for button clicks
  observeEvent(input$show_metadata, {
    output_view("metadata")
  })

  observeEvent(input$show_count_table, {
    output_view("count/Abundance_table")
  })

  observeEvent(input$show_summary_statistics, {
    output_view("summary_statistics")
  })


  observeEvent(input$show_combined_table, {
    output_view("combined_table")
  })

  # Show the selected data section

  output$dynamic_tables <- renderUI({
    req(output_view())  # Ensure the button was clicked
    req(physeq())       # Ensure physeq() is initialized

    if (output_view() == "metadata") {
      fluidRow(
        box(title = "Metadata Structure", width = 12, status = "primary",
            DT::dataTableOutput("metadata_structure")
        )
      )
    } else if (output_view() == "combined_table") {
      fluidRow(
        box(title = "Abundance/Taxonomy Table", width = 12, status = "primary",
            DT::dataTableOutput("combined_table"),
            selectInput("filetype", "Choose file type:", choices = c("csv", "xlsx", "tsv")),
            downloadButton("download_combined_table", "Download Abundance/Taxonomy Table")
        )
      )
    } else if (output_view() == "summary_statistics") {
      fluidRow(
        box(title = "Summary Statistics", width = 12, status = "primary",
            DT::dataTableOutput("summary_statistics")
        )
      )
    } else {
      # Placeholder to prevent errors
      fluidRow(
        box(title = "No Data", width = 12, status = "warning",
            p("Please upload a valid phyloseq object and select an option.")
        )
      )
    }
  })

  # output$dynamic_tables <- renderUI({
  #   if (output_view() == "metadata") {
  #     fluidRow(
  #       box(title = "Metadata Structure", width = 12, status = "primary",
  #           DT::dataTableOutput("metadata_structure")
  #       )
  #     )
  #   } else if (output_view() == "combined_table") {
  #     fluidRow(
  #       box(title = "Abundance/Taxonomy Table", width = 12, status = "primary",
  #           # Display combined table
  #           DT::dataTableOutput("combined_table"),
  #           # File format selection for download (under the table)
  #           selectInput("filetype", "Choose file type:", choices = c("csv", "xlsx", "tsv")),
  #           # Download button for the combined table
  #           downloadButton("download_combined_table", "Download Abundance/Taxonomy Table")
  #       )
  #     )
  #   } else if (output_view() == "summary_statistics") {
  #     fluidRow(
  #       box(title = "Summary Statistics", width = 12, status = "primary",
  #           DT::dataTableOutput("summary_statistics")
  #       )
  #     )
  #
  #   }
  #
  # })



  # Render Metadata Structure when button is clicked
  output$metadata_structure <- DT::renderDT({
    req(physeq())  # Ensure physeq is loaded
    req(output_view() == "metadata")  # Ensure correct tab is active

    as.data.frame(phyloseq::sample_data(physeq()))
  }, options = list(pageLength = 10, scrollX = TRUE))


  # Store the combined table in a reactive expression to avoid duplication
  combined_table <- reactive({
    # Extract the OTU/ASV count table
    otu_table <- as.data.frame(otu_table(physeq()))

    # Extract the taxonomy table
    taxonomy_table <- as.data.frame(phyloseq::tax_table(physeq()))

    # Ensure that the OTU table and taxonomy table align by rows (taxa/ASVs)
    if (nrow(taxonomy_table) != nrow(otu_table)) {
      stop("Number of rows in the taxonomy table does not match the OTU table")
    }

    # Combine the OTU/ASV count table with the taxonomy table
    combined <- merge(otu_table, taxonomy_table, by = 'row.names', all = TRUE)

    # Rename the Row.names column to ASVs
    colnames(combined)[1] <- "ASVs"

    return(combined)
  })

  # Render Combined Abundance/Taxonomy Table
  output$combined_table <- DT::renderDT({
    req(physeq())  # Ensure physeq is initialized
    req(output_view() == "combined_table")  # Ensure the correct button is active

    validate(
      need(!is.null(phyloseq::otu_table(physeq())), "Error: OTU table is missing.")
    )

    # Render the combined table
    DT::datatable(combined_table(), options = list(pageLength = 10, scrollX = TRUE), rownames = FALSE)
  })

  # Generate download handler for the Abundance/Taxonomy table
  output$download_combined_table <- downloadHandler(
    filename = function() {
      paste0("abundance_taxonomy_table_", Sys.Date(), ".", input$filetype)
    },
    content = function(file) {
      combined <- combined_table()  # Use the precomputed combined table

      # Save the file based on the selected filetype
      if (input$filetype == "csv") {
        write.csv(combined, file, row.names = FALSE)
      } else if (input$filetype == "xlsx") {
        writexl::write_xlsx(combined, path = file)
      } else if (input$filetype == "tsv") {
        write.table(combined, file, sep = "\t", row.names = FALSE, quote = FALSE)
      }
    }
  )

  # Render Summary Statistics
  output$summary_statistics <- DT::renderDT({
    req(physeq())
    validate(
      need(!is.null(physeq()), "Error: No phyloseq object found."),
      need(!is.null(phyloseq::tax_table(physeq())), "Error: Taxonomy table is missing.")
    )

    num_samples <- nsamples(physeq())
    num_asvs <- ntaxa(physeq())
    num_phylum <- length(unique(phyloseq::tax_table(physeq())[,"Phylum"]))
    num_family <- length(unique(phyloseq::tax_table(physeq())[,"Family"]))
    num_genus <- length(unique(phyloseq::tax_table(physeq())[,"Genus"]))
    num_species <- length(unique(phyloseq::tax_table(physeq())[,"Species"]))

    summary_df <- data.frame(
      "Metric" = c("Number of Samples", "Number of ASVs", "Number of Phyla", "Number of Families", "Number of Genera", "Number of Species"),
      "Count" = c(num_samples, num_asvs, num_phylum, num_family, num_genus, num_species)
    )

    DT::datatable(summary_df, options = list(pageLength = 10, scrollX = TRUE), rownames = FALSE)
  })




  # output$summary_statistics <- DT::renderDT({
  #   req(output_view() == "summary_statistics")
  #   num_samples <- nsamples(physeq())
  #   num_asvs <- ntaxa(physeq())
  #   num_phylum <- length(unique(phyloseq::tax_table(physeq())[,"Phylum"]))
  #   num_family <- length(unique(phyloseq::tax_table(physeq())[,"Family"]))
  #   num_genus <- length(unique(phyloseq::tax_table(physeq())[,"Genus"]))
  #   num_species <- length(unique(phyloseq::tax_table(physeq())[,"Species"]))
  #
  #   summary_df <- data.frame(
  #     "Metric" = c("Number of Samples", "Number of ASVs", "Number of Phyla", "Number of Families", "Number of Genera", "Number of Species"),
  #     "Count" = c(num_samples, num_asvs, num_phylum, num_family, num_genus, num_species)
  #   )
  #
  #   DT::datatable(summary_df, options = list(pageLength = 10, scrollX = TRUE), rownames = FALSE)
  # })





  ######### Reads QC

  # Explanation popup for Reads QC
  observeEvent(input$explanation_button_reads_qc, {
    showModal(modalDialog(
      title = "Reads QC Section Overview",
      p("This section provides quality control analysis for sequencing reads."),
      p("You can analyze the distribution of reads across samples or groups, view the number of reads per sample, or plot rarefaction curves to assess sequencing depth sufficiency."),
      easyClose = TRUE,
      footer = modalButton("Close")
    ))
  })

  # Track which button is pressed in Reads QC tab
  output_reads_qc_view <- reactiveVal("")

  observeEvent(input$show_reads_distribution_samples, {
    output_reads_qc_view("reads_distribution_samples")
  })

  observeEvent(input$show_reads_distribution_groups, {
    output_reads_qc_view("reads_distribution_groups")
  })

  observeEvent(input$show_number_of_reads, {
    output_reads_qc_view("number_of_reads")
  })

  observeEvent(input$show_rarefaction_curves, {
    output_reads_qc_view("rarefaction_curves")
  })



  # Dynamically render UI content based on the button clicked
  output$reads_qc_content <- renderUI({
    if (output_reads_qc_view() == "reads_distribution_groups") {
      fluidRow(
        column(
          width = 12,
          selectInput(
            "group_column",
            "Select Grouping Column",
            choices = sample_variables(physeq())
          ),
          div(
            style = "padding: 20px; border: 1px solid #ddd; background-color: #f9f9f9;",
            plotlyOutput("reads_distribution_boxplot_groups", height = "600px", width = "90%")
          )
        )
      )
    } else if (output_reads_qc_view() == "reads_distribution_samples") {
      fluidRow(
        column(
          width = 12,
          plotlyOutput("reads_distribution_plot_samples"),
          selectInput("plot_filetype_samples", "Select File Type", choices = c("png", "pdf", "svg")),
          numericInput("plot_width_samples", "Plot Width (inches)", value = 8, min = 4),
          numericInput("plot_height_samples", "Plot Height (inches)", value = 6, min = 4),
          downloadButton("download_plot_samples", "Download Plot")
        )
      )
    } else if (output_reads_qc_view() == "number_of_reads") {
      fluidRow(
        column(
          width = 12,
          DT::dataTableOutput("number_of_reads_table"),
          numericInput("reads_threshold", "Reads Threshold for Sample Filtration", value = 1000, min = 0),
          selectInput("sample_name_filter", "Filter by Sample Name(s)", choices = sample_names(physeq()), multiple = TRUE),
          actionButton("filter_samples", "Filter Samples"),
          DT::dataTableOutput("filtered_samples_table"),
          downloadButton("download_filtered_physeq", "Download Filtered Phyloseq")
        )
      )
    } else if (output_reads_qc_view() == "rarefaction_curves") {
      fluidRow(
        column(
          width = 12,
          selectInput("group_column", "Select Grouping Column", choices = sample_variables(physeq())),
          plotlyOutput("rarefaction_curves_plot", height = "600px", width = "100%")
        )
      )
    } else {
      NULL  # Render nothing if no button is clicked
    }
  })


  # Define common download options
  downloadOptions <- function(prefix) {
    list(
      selectInput(paste0("plot_filetype_", prefix), "Select File Type", choices = c("png", "pdf", "svg")),
      numericInput(paste0("plot_width_", prefix), "Plot Width (inches)", value = 8, min = 4),
      numericInput(paste0("plot_height_", prefix), "Plot Height (inches)", value = 6, min = 4),
      downloadButton(paste0("download_", prefix), "Download Plot")
    )
  }


  # Render Reads Distribution Plot Across Samples
  output$reads_distribution_plot_samples <- renderPlotly({
    req(physeq())
    reads_per_sample <- sample_sums(physeq())
    p <- ggplot2::ggplot(data.frame(Reads = reads_per_sample), aes(x = Reads)) +
      ggplot2::geom_histogram(binwidth = 10000, fill = "blue", color = "black", alpha = 0.7) +
      labs(title = "Reads Distribution Across Samples", x = "Number of Reads", y = "Frequency") +
      theme_minimal()
    ggplotly(p)
  })

  # Download handler for Reads Distribution Across Samples plot
  output$download_plot_samples <- downloadHandler(
    filename = function() {
      paste0("reads_distribution_samples_", Sys.Date(), ".", input$plot_filetype_samples)
    },
    content = function(file) {
      p <- ggplot2::ggplot(data.frame(Reads = sample_sums(physeq())), aes(x = Reads)) +
        ggplot2::geom_histogram(binwidth = 10000, fill = "blue", color = "black", alpha = 0.7) +
        labs(title = "Reads Distribution Across Samples", x = "Number of Reads", y = "Frequency") +
        theme_minimal()

      # Use ggsave with the specified device
      ggsave(filename = file, plot = p, device = input$plot_filetype_samples, width = input$plot_width_samples, height = input$plot_height_samples)
    }
  )



  # Plot and download handlers for Reads Distribution Across Groups with Boxplot
  # Render Reads Distribution Plot Across Groups with Boxplot
  output$reads_distribution_boxplot_groups <- renderPlotly({
    req(physeq(), input$group_column)  # Ensure physeq and group_column are available

    validate(
      need(input$group_column != "", "Please select a grouping column to generate the plot.")
    )

    reads_df <- data.frame(
      Sample = sample_names(physeq()),
      Group = phyloseq::sample_data(physeq())[[input$group_column]],
      Reads = sample_sums(physeq())
    )

    p <- ggplot2::ggplot(reads_df, aes(x = Group, y = Reads, fill = Group)) +
      geom_boxplot(alpha = 0.7, outlier.color = "red") +
      geom_jitter(aes(text = paste("Sample:", Sample, "<br>Reads:", Reads)),
                  color = "black", size = 1.5, width = 0.2, alpha = 0.6) +
      labs(
        title = "Reads Distribution per Sample Group",
        x = "Group",
        y = "Number of Reads"
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.margin = margin(20, 20, 20, 20)  # Add margins around the plot
      )

    ggplotly(p, tooltip = "text")
  })




  # Download handler for Reads Distribution Box Plot Across Groups
  output$download_boxplot_groups <- downloadHandler(
    filename = function() { paste0("reads_distribution_boxplot_", Sys.Date(), ".", input$plot_filetype_groups) },
    content = function(file) {
      # Regenerate the ggplot2 object
      req(physeq(), input$group_column)
      reads_df <- data.frame(Sample = sample_names(physeq()), Group = phyloseq::sample_data(physeq())[[input$group_column]], Reads = sample_sums(physeq()))

      p <- ggplot2::ggplot(reads_df, aes(x = Group, y = Reads, fill = Group)) +
        geom_boxplot(alpha = 0.7, outlier.color = "red") +
        geom_jitter(aes(text = paste("Sample:", Sample, "<br>Reads:", Reads)), color = "black", size = 1.5, width = 0.2, alpha = 0.6) +
        labs(title = "Reads Distribution per Sample Group", x = "Group", y = "Number of Reads") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))

      # Use ggsave with the specified device and dimensions
      ggsave(filename = file, plot = p, device = input$plot_filetype_groups, width = input$plot_width_groups, height = input$plot_height_groups)
    }
  )


  # Reactive value to store the filtered phyloseq object, initially NULL
  filtered_pseq <- reactiveVal(NULL)

  # Define a proxy function to use the filtered or original phyloseq object
  current_physeq <- reactive({
    if (!is.null(filtered_pseq())) {
      filtered_pseq()  # Return filtered version if available
    } else {
      physeq()  # Otherwise, return the original
    }
  })

  # Reactive data table for Number of Reads per Sample
  output$number_of_reads_table <- DT::renderDT({
    req(current_physeq())  # Use current_physeq() instead of physeq()
    reads_df <- data.frame(Sample = sample_names(current_physeq()), Reads = sample_sums(current_physeq()))
    DT::datatable(reads_df, options = list(pageLength = 10, scrollX = TRUE))
  })

  # Sample filtering based on reads threshold and sample names
  observeEvent(input$filter_samples, {
    req(input$reads_threshold, physeq())  # Ensure threshold and original physeq are available

    # Apply filtering on the phyloseq object
    filtered <- physeq() %>%
      phyloseq::prune_samples(sample_sums(.) >= input$reads_threshold, .) %>%
      phyloseq::prune_taxa(taxa_sums(.) > 0, .)

    if (!is.null(input$sample_name_filter) && length(input$sample_name_filter) > 0) {
      filtered <- phyloseq::prune_samples(!sample_names(filtered) %in% input$sample_name_filter, filtered)
      filtered <- phyloseq::prune_taxa(taxa_sums(filtered) > 0, filtered)
    }

    # Update the reactive variable to hold the filtered phyloseq object
    filtered_pseq(filtered)  # Now, current_physeq() will reference this filtered object

    # Display filtered samples in a data table
    output$filtered_samples_table <- DT::renderDT({
      filtered_samples_df <- data.frame(Sample = sample_names(current_physeq()), Reads = sample_sums(current_physeq()))
      DT::datatable(filtered_samples_df, options = list(pageLength = 10, scrollX = TRUE))
    })

    # Download filtered phyloseq object
    output$download_filtered_physeq <- downloadHandler(
      filename = function() { paste("filtered_phyloseq_", Sys.Date(), ".rds", sep = "") },
      content = function(file) { saveRDS(current_physeq(), file) }
    )
  })



  # Now, replace `physeq()` with `current_physeq()` in the rest of the app components
  # Example for Reads per Sample Plot
  output$reads_per_sample_plot <- renderPlotly({
    req(current_physeq(), input$generate_reads_plot)  # Use current_physeq instead of physeq

    reads_per_sample <- sample_sums(current_physeq())
    samples_df <- data.frame(Sample = sample_names(current_physeq()), Reads = reads_per_sample)

    p <- ggplot2::ggplot(samples_df, aes(x = Sample, y = Reads)) +
      geom_bar(stat = "identity", fill = "steelblue") +
      labs(title = "Total Number of Reads per Sample", x = "Sample", y = "Total Reads") +
      theme_minimal() + theme(axis.text.x = element_text(angle = 90, hjust = 1))

    ggplotly(p)
  })

  # Download handler for Reads per Sample plot using `current_physeq()`
  output$download_reads_plot <- downloadHandler(
    filename = function() { paste0("reads_per_sample_plot_", Sys.Date(), ".", input$plot_filetype_reads) },
    content = function(file) {
      reads_per_sample <- sample_sums(current_physeq())
      samples_df <- data.frame(Sample = sample_names(current_physeq()), Reads = reads_per_sample)

      p <- ggplot2::ggplot(samples_df, aes(x = Sample, y = Reads)) +
        geom_bar(stat = "identity", fill = "steelblue") +
        labs(title = "Total Number of Reads per Sample", x = "Sample", y = "Total Reads") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1))

      ggsave(filename = file, plot = p, device = input$plot_filetype_reads, width = input$plot_width_reads, height = input$plot_height_reads)
    }
  )



  # # Reactive expression for rarefaction curve data
  rarefaction_data <- reactive({
    req(current_physeq())  # Ensure phyloseq object is loaded

    # Validate that a grouping column has been selected
    validate(
      need(input$group_column != "", "Please select a grouping column to generate the plot.")
    )

    # Filter samples with sufficient reads
    filtered_physeq <- phyloseq::prune_samples(sample_sums(current_physeq()) >= 100, current_physeq())  # Keep samples with at least 100 reads

    # Check if filtering removed all samples
    if (nsamples(filtered_physeq) == 0) {
      showNotification("No samples have sufficient reads for rarefaction. Please check your data.", type = "error")
      return(NULL)
    }

    # Adjust chunk size based on minimum read count
    min_count <- min(sample_sums(filtered_physeq))
    chunk_size <- ifelse(min_count < 99, min_count, 99)

    # Generate rarefaction curve data
    set.seed(1024)
    tryCatch({
      get_rarecurve(obj = filtered_physeq, chunks = chunk_size)
    }, error = function(e) {
      showNotification("Error generating rarefaction curves. Please check data and chunk size.", type = "error")
      NULL
    })
  })


  if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
    showNotification("The required package 'RColorBrewer' is not installed. Please install it to use this feature.", type = "error")
    return()
  }


  # Reactive expression to create the rarefaction plot
  rarefaction_plot <- reactive({
    req(current_physeq(), input$group_column)  # Ensure phyloseq object and group column are selected

    # Generate the rarefaction data
    rareres <- rarefaction_data()  # Use the reactive rarefaction data

    # Create the ggplot2 rarefaction plot
    prare2 <- ggrarecurve(
      obj = rareres,
      factorNames = input$group_column,
      shadow = FALSE,
      indexNames = c("Observe", "Chao1", "ACE")
    ) +
      scale_color_manual(values = colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(20))+
      ggplot2::theme_bw() +
      theme(
        axis.text = element_text(size = 8),
        panel.grid = element_blank(),
        strip.background = element_rect(colour = NA, fill = "grey"),
        strip.text.x = element_text(face = "bold")
      )

    prare2  # Return the ggplot object
  })

  # Render the rarefaction plot using Plotly for interactivity
  output$rarefaction_curves_plot <- renderPlotly({
    ggplotly(rarefaction_plot())  # Convert the ggplot object to a Plotly object
  })

  # Download handler for the rarefaction plot
  output$download_rarefaction_curves <- downloadHandler(
    filename = function() {
      paste0("rarefaction_curves_", Sys.Date(), ".", input$plot_filetype_groups)
    },
    content = function(file) {
      # Use ggsave to save the plot generated by the reactive expression
      ggsave(
        filename = file,
        plot = rarefaction_plot(),  # Use the reactive plot directly
        device = input$plot_filetype_groups,
        width = input$plot_width_groups,
        height = input$plot_height_groups
      )
    }
  )





  ######### Taxa overview

  # Reactive value to track which button is pressed in Taxa Overview tab
  output_taxa_overview_view <- reactiveVal("")

  # Observe button clicks and update reactive value accordingly
  observeEvent(input$show_taxa_prevalence_samples, {
    output_taxa_overview_view("taxa_prevalence_samples")
  })

  observeEvent(input$show_taxa_prevalence_groups, {
    output_taxa_overview_view("taxa_prevalence_groups")
  })

  observeEvent(input$show_dominant_taxa_samples, {
    output_taxa_overview_view("dominant_taxa_samples")
  })

  # Observe button clicks and update reactive value accordingly
  observeEvent(input$show_prevalence_abundance_plot, {
    output_taxa_overview_view("prevalence_abundance_plot")
  })


  observeEvent(input$show_upset_plot, {
    output_taxa_overview_view("upset_plot")
  })

  observeEvent(input$show_core_microbiome_plot, {
    output_taxa_overview_view("core_microbiome_plot")
  })
  #
  #   observeEvent(input$show_sunburst_plot, {
  #     output_taxa_overview_view("sunburst_plot")
  #   })




  # Dynamically render UI content based on the button clicked
  output$taxa_overview_content <- renderUI({
    if (output_taxa_overview_view() == "taxa_prevalence_samples") {
      fluidRow(
        box(title = "Taxa Prevalence Across Samples", width = 12, status = "primary",
            selectInput("taxa_level_samples", "Select Taxonomic Level", choices = rank_names(physeq())),
            DT::dataTableOutput("taxa_prevalence_table_samples"),
            downloadButton("download_taxa_prevalence_samples", "Download Taxa Prevalence Data")
        )
      )
    } else if (output_taxa_overview_view() == "taxa_prevalence_groups") {
      fluidRow(
        box(title = "Taxa Prevalence Across Groups", width = 12, status = "primary",
            selectInput("taxa_level_groups", "Select Taxonomic Level", choices = rank_names(physeq())),
            selectInput("group_column", "Select Grouping Column", choices = sample_variables(physeq())),
            DT::dataTableOutput("taxa_prevalence_table_groups"),
            downloadButton("download_taxa_prevalence_groups", "Download Taxa Prevalence Data")
        )
      )

    } else if (output_taxa_overview_view() == "dominant_taxa_samples") {
      fluidRow(
        box(title = "Dominant Taxa Per Sample", width = 12, status = "primary",
            selectInput("dominant_taxa_level_samples", "Select Taxonomic Level", choices = rank_names(physeq())),
            selectInput("group_column", "Select Grouping Column", choices = sample_variables(physeq())),
            numericInput("Abundance_threshold", "Set Abundance Threshold (0 to 1)", value = 0.001, min = 0, max = 1, step = 0.001),
            numericInput("min_prevalence", "Minimum Prevalence (%)", value = 20, min = 0, max = 100, step = 5),
            numericInput("n_max_dominant", "Set Maximum Number of Dominant Taxa", value = 10, min = 3),
            actionButton("calculate_dominant_taxa", "Calculate Dominant Taxa"),
            DT::dataTableOutput("dominant_taxa_table_groups"),
            plotOutput("dominant_taxa_barplot",  height = "600px", width = "100%"),
            selectInput("plot_filetype_groups", "Select File Type", choices = c("png", "pdf", "svg")),
            numericInput("plot_width_groups", "Plot Width (inches)", value = 8, min = 4),
            numericInput("plot_height_groups", "Plot Height (inches)", value = 6, min = 4),
            downloadButton("download_dominant_taxa_plot", "Download Dominant Taxa Plot")
        )
      )
    } else if (output_taxa_overview_view() == "prevalence_abundance_plot") {
      fluidRow(
        box(title = "Prevalence vs Abundance Plot", width = 12, status = "primary",
            plotlyOutput("prevalence_abundance_plot", height = "600px", width = "100%"),
            selectInput("plot_filetype_groups", "Select File Type", choices = c("png", "pdf", "svg")),
            numericInput("plot_width_groups", "Plot Width (inches)", value = 8, min = 4),
            numericInput("plot_height_groups", "Plot Height (inches)", value = 6, min = 4),
            downloadButton("download_prevalence_abundance_plot", "Download Prevalence vs Abundance Plot")
        )
      )
    } else if (output_taxa_overview_view() == "upset_plot") {
      fluidRow(
        box(title = "Upset Plot", width = 12, status = "primary",
            # Add inputs for detection and prevalence threshold
            numericInput("detection_threshold_upset", "Detection Threshold", value = 0.001, min = 0, max = 1, step = 0.001),
            numericInput("prevalence_threshold_upset", "Prevalence Threshold", value = 0.1, min = 0, max = 1, step = 0.01),
            selectInput("group_column", "Select Grouping Column", choices = sample_variables(physeq())),
            plotOutput("upset_plot", height = "600px", width = "100%"),
            selectInput("plot_filetype_upset", "Select File Type", choices = c("png", "pdf", "svg")),
            numericInput("plot_width_upset", "Plot Width (inches)", value = 8, min = 4),
            numericInput("plot_height_upset", "Plot Height (inches)", value = 6, min = 4),
            downloadButton("download_upset_plot", "Download Upset Plot")
        )
      )
    } else if (output_taxa_overview_view() == "core_microbiome_plot") {
      fluidRow(
        box(title = "Core Microbiome Plot", width = 12, status = "primary",
            plotlyOutput("core_microbiome_plot", height = "600px", width = "100%"),
            downloadButton("download_core_microbiome_plot", "Download Core Microbiome Plot")
        )
      )
    }
  })


  # Server-side code to handle explanation modals
  observeEvent(input$explain_taxa_prevalence_samples, {
    showModal(modalDialog(
      title = "Taxa Prevalence Across Samples - Explanation",
      p("This section allows you to explore the prevalence of taxa across all samples at different taxonomic levels."),
      easyClose = TRUE,
      footer = modalButton("Close")
    ))
  })

  observeEvent(input$explain_taxa_prevalence_groups, {
    showModal(modalDialog(
      title = "Taxa Prevalence Across Groups - Explanation",
      p("This section allows you to compare the prevalence of taxa across different groups at different taxonomic levels."),
      easyClose = TRUE,
      footer = modalButton("Close")
    ))
  })


  observeEvent(input$explain_dominant_taxa_samples, {
    showModal(modalDialog(
      title = "Dominant Taxa Per Sample - Explanation",
      p("This section allows you to determine the dominant taxa within each sample at different taxonomic levels."),
      easyClose = TRUE,
      footer = modalButton("Close")
    ))
  })



  # Taxa Prevalence Across samples
  # Reactive value to store the prevalence table for samples
  prevalence_table_samples <- reactiveVal(NULL)

  # Calculate taxa prevalence across samples
  output$taxa_prevalence_table_samples <- DT::renderDT({
    req(input$taxa_level_samples, current_physeq())
    taxa_level <- input$taxa_level_samples

    # Aggregate OTU table at the selected taxonomic level
    physeq_agg <- aggregate_taxa(current_physeq(), taxa_level)

    prevalence_df <- apply(as.data.frame(otu_table(physeq_agg)), 1, function(x) sum(x > 0) / nsamples(physeq_agg) * 100)
    prevalence_df <- round(prevalence_df, 3)  # Round to 3 decimal places

    absolute_abundance <- rowSums(otu_table(physeq_agg))

    prevalence_table <- data.frame(Species = phyloseq::tax_table(physeq_agg)[, taxa_level], Prevalence = prevalence_df, Abundance = absolute_abundance)
    prevalence_table_samples(prevalence_table)  # Store the table in reactive value

    # Render the table without row names
    DT::datatable(prevalence_table, options = list(pageLength = 10, scrollX = TRUE), rownames = FALSE)
  })

  # Download handler for taxa prevalence across samples as TSV
  output$download_taxa_prevalence_samples <- downloadHandler(
    filename = function() {
      paste("taxa_prevalence_samples_", Sys.Date(), ".tsv", sep = "")
    },
    content = function(file) {
      write.table(prevalence_table_samples(), file, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
    }
  )



  # Taxa Prevalence Across Groups
  # Reactive value to store the prevalence table
  prevalence_table_reactive <- reactiveVal(NULL)


  # Render the taxa prevalence table across groups
  output$taxa_prevalence_table_groups <- DT::renderDT({
    req(input$taxa_level_groups, input$group_column, current_physeq())
    taxa_level <- input$taxa_level_groups
    group_column <- input$group_column

    # Aggregate OTU table at the selected taxonomic level
    physeq_agg <- aggregate_taxa(current_physeq(), taxa_level)

    # Get the OTU table and the grouping variable
    otu_table_agg <- as.data.frame(otu_table(physeq_agg))
    groups <- as.factor(phyloseq::sample_data(physeq_agg)[[group_column]])

    # Calculate prevalence within each group with rounding to two decimal places
    prevalence_list <- lapply(levels(groups), function(group) {
      group_samples <- which(groups == group)
      apply(otu_table_agg[, group_samples], 1, function(x) round(sum(x > 0) / length(group_samples) * 100, 2))
    })

    # Combine the prevalence results into a single data frame
    prevalence_df <- do.call(cbind, prevalence_list)
    colnames(prevalence_df) <- levels(groups)

    # Create the prevalence table
    prevalence_table <- data.frame(Species = phyloseq::tax_table(physeq_agg)[, taxa_level], prevalence_df)
    prevalence_table_reactive(prevalence_table)  # Store the table in reactive value

    # Render the table without row names
    DT::datatable(prevalence_table, options = list(pageLength = 10, scrollX = TRUE), rownames = FALSE)
  })

  # Download handler for taxa prevalence across groups as TSV
  output$download_taxa_prevalence_groups <- downloadHandler(
    filename = function() {
      paste("taxa_prevalence_groups_", Sys.Date(), ".tsv", sep = "")
    },
    content = function(file) {
      write.table(prevalence_table_reactive(), file, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
    }
  )



  # Dominant taxa in each sample
  # Reactive plot object for dominant taxa (globally defined)
  dominant_taxa_barplot_reactive <- reactiveVal()

  # Calculate dominant taxa when the button is clicked
  observeEvent(input$calculate_dominant_taxa, {
    req(current_physeq(), input$dominant_taxa_level_samples, input$Abundance_threshold, input$n_max_dominant, input$group_column, input$min_prevalence)

    taxa_level <- input$dominant_taxa_level_samples
    threshold <- input$Abundance_threshold   # minimum proportion at which to consider a sample dominated by a taxon
    n_max <- input$n_max_dominant            # maximum number of taxa that can be listed as dominant taxa
    min_prevalence <- input$min_prevalence / 100  # Convert percentage to a proportion

    # Calculate dominant taxa plot and store it in the reactive value
    tryCatch({
      physeq_dominant <- current_physeq() %>%
        tax_fix() %>%
        phyloseq_validate() %>%
        ps_calc_dominant(
          rank = taxa_level,
          threshold = threshold,
          n_max = n_max,
          other = "Other",
          none = "Not dominated",
          var = paste("dominant", taxa_level, sep = "_")
        )

      # Create the dominant taxa plot
      p <- microViz::comp_barplot(physeq_dominant, tax_level = taxa_level, facet_by = input$group_column,
                                  label = paste("dominant", taxa_level, sep = "_"), n_taxa = n_max) +
        coord_flip() +
        labs(title = paste("Dominant Taxa Composition in each Sample"))

      # Store the plot in the reactive value
      dominant_taxa_barplot_reactive(p)
    }, error = function(e) {
      showNotification("Error in generating dominant taxa barplot: Check your inputs.", type = "error")
      dominant_taxa_barplot_reactive(NULL)
    })
  })

  # Render the dominant taxa bar plot
  output$dominant_taxa_barplot <- renderPlot({
    req(dominant_taxa_barplot_reactive())
    dominant_taxa_barplot_reactive()
  })

  # Download handler for dominant taxa bar plot
  output$download_dominant_taxa_plot <- downloadHandler(
    filename = function() {
      paste("dominant_taxa_plot_", Sys.Date(), ".png", sep = "")
    },
    content = function(file) {
      # Retrieve the plot stored in the reactive value
      p <- dominant_taxa_barplot_reactive()
      req(p)  # Ensure that the plot exists

      # Save the plot to the specified file
      ggsave(file, plot = p, device = "png", width = 8, height = 6)
    }
  )




  # Reactive expression to store the Prevalence vs Abundance plot
  prevalence_abundance_plot <- reactive({
    req(current_physeq())

    # Transform to relative abundance (compositional data)
    ps_rel <- microbiome::transform(current_physeq(), "compositional")

    # Aggregate data at the Genus level
    physeq_agg <- aggregate_taxa(ps_rel, "Genus")

    # Calculate prevalence and abundance
    prevalence_df <- apply(as.data.frame(otu_table(physeq_agg)), 1, function(x) sum(x > 0) / nsamples(physeq_agg) * 100)
    abundance_df <- rowSums(otu_table(physeq_agg))

    # Extract Genus and Phylum from tax_table
    taxa_df <- as.data.frame(phyloseq::tax_table(physeq_agg))

    # Create a data frame for plotting
    plot_df <- data.frame(
      Prevalence = prevalence_df,
      Abundance = abundance_df,
      Genus = taxa_df$Genus,       # Add Genus as a column
      Phylum = taxa_df$Phylum      # Add Phylum as a column
    )

    # Create the ggplot object
    p <- ggplot2::ggplot(plot_df, aes(x = Prevalence, y = Abundance, color = Phylum, label = Genus)) +
      ggplot2::geom_point(size = 3, alpha = 0.8) +
      scale_y_log10() +  # Log scale for abundance
      theme_minimal() +
      labs(title = "Prevalence vs Abundance", x = "Prevalence (%)", y = "Abundance") +
      theme(legend.position = "right")

    return(p)
  })

  # Render the plot
  output$prevalence_abundance_plot <- renderPlotly({
    ggplotly(prevalence_abundance_plot())  # Convert the ggplot object to an interactive plotly object
  })

  # Download handler for Prevalence vs Abundance Plot
  output$download_prevalence_abundance_plot <- downloadHandler(
    filename = function() {
      paste0("prevalence_vs_abundance_", Sys.Date(), ".", input$plot_filetype_groups)  # Dynamic filename based on date and file type
    },
    content = function(file) {
      # Reuse the ggplot object stored in the reactive expression
      ggsave(file, plot = prevalence_abundance_plot(),
             device = input$plot_filetype_groups,  # Dynamic file type
             width = input$plot_width_groups,      # Dynamic width
             height = input$plot_height_groups)    # Dynamic height
    }
  )


  # Add explanation for Prevalence vs Abundance Plot
  observeEvent(input$explain_prevalence_abundance_plot, {
    showModal(modalDialog(
      title = "Prevalence vs Abundance Plot - Explanation",
      p("This plot compares the prevalence of each Genus (percentage of samples in which the Genus is present) with its total abundance across all samples. Points are colored by Phylum."),
      easyClose = TRUE,
      footer = modalButton("Close")
    ))
  })



  ### upset_plot

  # Reusable function to generate the upset plot using get_upset
  generate_upset_plot <- function(physeq, detection_threshold, prevalence_threshold, group_var) {
    # Get the core taxa based on the thresholds
    physeq_core <- core(physeq, detection = detection_threshold, prevalence = prevalence_threshold)

    # Check if there are any taxa remaining
    if (ntaxa(physeq_core) == 0) {
      stop("No taxa meet the detection and prevalence thresholds.")
    }

    # Check for missing values in the grouping column
    if (any(is.na(phyloseq::sample_data(physeq_core)[[group_var]]))) {
      stop("The selected grouping column contains missing values. Please choose a different column or remove missing values.")
    }

    # Use get_upset to generate the binary matrix for upset plot
    upset_data <- get_upset(physeq_core, factorNames = group_var)

    # Create the upset plot using sample groupings
    upset_plot <- UpSetR::upset(
      upset_data,
      sets = unique(as.vector(phyloseq::sample_data(physeq_core)[[group_var]])),  # Grouping by sample data
      sets.bar.color = "#56B4E9",
      keep.order = TRUE,              # Keep the order of the sets
      order.by = "freq",              # Order by frequency of occurrence
      nsets = 10,                     # Number of sets to display
      nintersects = 20,               # Number of intersections to display
      main.bar.color = "blue",        # Customize bar color
      empty.intersections = "on"      # Show empty intersections
    )

    return(upset_plot)
  }


  # Render the Upset Plot using the reusable function
  output$upset_plot <- renderPlot({
    req(current_physeq(), input$detection_threshold_upset, input$prevalence_threshold_upset, input$group_column)

    # Use the reusable function to generate the plot
    upset_plot <- generate_upset_plot(current_physeq(), input$detection_threshold_upset, input$prevalence_threshold_upset, input$group_column)

    return(upset_plot)
  })



  # Download handler for Upset Plot using the reusable function
  output$download_upset_plot <- downloadHandler(
    filename = function() {
      paste0("upset_plot_", Sys.Date(), ".", input$plot_filetype_upset)
    },
    content = function(file) {
      req(current_physeq(), input$detection_threshold_upset, input$prevalence_threshold_upset, input$group_column)

      # Open the appropriate device based on the selected file type
      if (input$plot_filetype_upset == "png") {
        png(file, width = input$plot_width_upset, height = input$plot_height_upset, units = "in", res = 300)
      } else if (input$plot_filetype_upset == "pdf") {
        pdf(file, width = input$plot_width_upset, height = input$plot_height_upset)
      } else if (input$plot_filetype_upset == "svg") {
        svg(file, width = input$plot_width_upset, height = input$plot_height_upset)
      }

      # Generate the plot and print it to the device
      print(generate_upset_plot(current_physeq(), input$detection_threshold_upset, input$prevalence_threshold_upset, input$group_column))

      # Close the device
      dev.off()
    }
  )


  # Enhanced Explanation for Upset Plot
  observeEvent(input$explain_upset_plot, {
    showModal(modalDialog(
      title = "Upset Plot - Explanation",
      p("This Upset plot shows the intersections of core taxa across samples."),
      p("It helps to identify shared or unique taxa across multiple groups or conditions."),
      p("Core taxa are defined based on user-adjustable detection and prevalence thresholds."),
      tags$ul(
        tags$li("**Set Size (left panel)**: The horizontal bars on the left represent the total number of taxa found in each group."),
        tags$li("**Intersection Size (top panel)**: The vertical bars at the top represent the number of taxa shared between the selected groups or unique to a specific group. Each dot and line combination below the bars shows the groups that are being compared in that intersection. For instance, a bar with a single dot below it indicates taxa unique to one group, whereas bars with lines connecting multiple dots indicate shared taxa across multiple groups."),
        tags$li("**Interpretation**: You can use this plot to see which taxa are unique to certain conditions (groups) or common to multiple groups. The higher the intersection size, the more core taxa are shared among the connected groups.")
      ),
      p("Adjust the detection and prevalence thresholds to fine-tune the definition of core taxa based on your analysis needs."),
      p("Empty intersections in an UpSet plot can occur when certain combinations of groups (sets) have no overlapping/shared taxa. This means that, for the given threshold values of detection and prevalence, there are no taxa shared among those groups, so the intersection is effectively empty."),
      p("You can prevent empty intersections from being displayed in the UpSet plot by using the 'empty.intersections' = 'off' parameter in the server code when creating the plot."),
      easyClose = TRUE,
      footer = modalButton("Close")
    ))
  })



  # UI part for Venn Diagram with selectable grouping column
  observeEvent(input$show_venn_diagram_plot, {
    output$taxa_overview_content <- renderUI({
      fluidRow(
        box(title = "Venn Diagram - Unique and Shared Taxa", width = 12, status = "primary",
            selectInput("taxonomy_level", "Select Taxonomy Level", choices = rank_names(current_physeq())),
            selectInput("grouping_column_venn", "Select Grouping Column", choices = colnames(phyloseq::sample_data(current_physeq()))),
            numericInput("detection_threshold_venn", "Detection Threshold", value = 0.001, min = 0, max = 1, step = 0.001),
            numericInput("prevalence_threshold_venn", "Prevalence Threshold", value = 0.1, min = 0, max = 1, step = 0.01),
            plotOutput("venn_diagram_plot", height = "600px", width = "100%"),
            downloadButton("download_unique_taxa", "Download Unique Taxa (TSV)"),
            downloadButton("download_shared_pairwise_taxa", "Download Pairwise Shared Taxa (TSV)"),
            downloadButton("download_shared_all_taxa", "Download Shared Across All Groups (TSV)"),
            selectInput("plot_filetype_venn", "Select File Type", choices = c("png", "pdf", "svg")),
            numericInput("plot_width_venn", "Plot Width (inches)", value = 8, min = 4),
            numericInput("plot_height_venn", "Plot Height (inches)", value = 6, min = 4),
            downloadButton("download_venn_plot", "Download Venn Diagram")
        )
      )
    })
  })

  # Helper function to get the taxonomy name above the selected level
  get_higher_taxonomy <- function(physeq, taxonomy_level) {
    ranks <- rank_names(physeq)
    level_idx <- which(ranks == taxonomy_level)
    if (level_idx == 1) {
      return(taxonomy_level)
    }
    return(ranks[level_idx - 1])
  }

  # Reactive variable to store the grouped OTU data for Venn plot with enhanced taxonomy names
  group_otus_reactive <- reactive({
    req(current_physeq(), input$taxonomy_level, input$detection_threshold_venn, input$prevalence_threshold_venn, input$grouping_column_venn)

    # Validate grouping column existence and handle missing values
    if (!input$grouping_column_venn %in% colnames(phyloseq::sample_data(current_physeq()))) {
      stop("Selected grouping column does not exist in the sample data.")
    }
    if (any(is.na(phyloseq::sample_data(current_physeq())[[input$grouping_column_venn]]))) {
      stop("The selected grouping column contains missing values. Please choose a different column or clean the data.")
    }

    # Get the higher taxonomy level
    higher_taxonomy_level <- get_higher_taxonomy(current_physeq(), input$taxonomy_level)

    # Aggregate physeq data at the selected taxonomy level
    physeq_agg <- microbiome::aggregate_taxa(current_physeq(), input$taxonomy_level)
    tax_table_data <- phyloseq::tax_table(physeq_agg)

    # Get core microbiome based on thresholds
    physeq_core <- core(physeq_agg, detection = input$detection_threshold_venn, prevalence = input$prevalence_threshold_venn)

    # Extract OTU data and metadata
    otu_data <- as(otu_table(physeq_core), "matrix")
    meta_data <- data.frame(phyloseq::sample_data(physeq_core))
    sample_groups <- unique(meta_data[[input$grouping_column_venn]])

    # Create list of OTUs for each group with enhanced taxonomy names
    group_otus <- lapply(sample_groups, function(group) {
      group_samples <- rownames(meta_data[meta_data[[input$grouping_column_venn]] == group, ])
      otus_in_group <- rownames(otu_data[, group_samples][rowSums(otu_data[, group_samples]) > 0, ])

      # Append the higher taxonomy level to the OTU names
      otus_with_higher_tax <- sapply(otus_in_group, function(otu) {
        higher_tax <- as.character(tax_table_data[otu, higher_taxonomy_level])
        otu_name <- as.character(tax_table_data[otu, input$taxonomy_level])

        if (is.na(higher_tax)) {
          return(otu_name)
        } else {
          return(paste0(higher_tax, "_", otu_name))
        }
      })

      return(otus_with_higher_tax)
    })

    # Name the list elements with group names
    names(group_otus) <- sample_groups

    return(group_otus)
  })

  # Reactive variable to store the Venn plot
  venn_plot_reactive <- reactive({
    group_otus <- group_otus_reactive()

    venn_plot <- ggvenn::ggvenn(
      group_otus,
      fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF", "darkgreen", "orchid4"),
      stroke_size = 0.5, set_name_size = 4
    ) +
      ggtitle(paste("Venn Diagram - Unique and Shared", input$taxonomy_level))

    return(venn_plot)
  })

  # Render the Venn Diagram using the reactive variable
  output$venn_diagram_plot <- renderPlot({
    venn_plot_reactive()
  })

  # Download handler for the Venn Diagram using the same reactive variable
  output$download_venn_plot <- downloadHandler(
    filename = function() {
      paste0("venn_diagram_", Sys.Date(), ".", input$plot_filetype_venn)
    },
    content = function(file) {
      venn_plot <- venn_plot_reactive()
      ggsave(file, plot = venn_plot, device = input$plot_filetype_venn,
             width = input$plot_width_venn, height = input$plot_height_venn)
    }
  )

  # Download handler for unique taxa in each group
  output$download_unique_taxa <- downloadHandler(
    filename = function() {
      paste0("pure_unique_taxa_", Sys.Date(), ".tsv")
    },
    content = function(file) {
      group_otus <- group_otus_reactive()

      # Identify unique taxa for each group
      unique_taxa <- lapply(names(group_otus), function(group) {
        other_groups <- setdiff(names(group_otus), group)
        unique_to_group <- setdiff(group_otus[[group]], Reduce(union, group_otus[other_groups]))
        if (length(unique_to_group) > 0) {
          return(data.frame(Group = group, Taxa = unique_to_group))
        } else {
          return(NULL)
        }
      })

      # Combine results and write to file
      unique_taxa_df <- do.call(rbind, unique_taxa)
      if (is.null(unique_taxa_df) || nrow(unique_taxa_df) == 0) {
        write.table(data.frame(Group = character(), Taxa = character()), file, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
      } else {
        write.table(unique_taxa_df, file, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
      }
    }
  )

  # Download handler for pairwise shared taxa
  output$download_shared_pairwise_taxa <- downloadHandler(
    filename = function() {
      paste0("pairwise_shared_taxa_", Sys.Date(), ".tsv")
    },
    content = function(file) {
      group_otus <- group_otus_reactive()

      # Get combinations and find shared taxa for each pair
      group_combinations <- combn(names(group_otus), 2, simplify = FALSE)
      shared_taxa_list <- lapply(group_combinations, function(groups) {
        shared_otus <- intersect(group_otus[[groups[1]]], group_otus[[groups[2]]])
        if (length(shared_otus) > 0) {
          return(data.frame(Group1 = groups[1], Group2 = groups[2], Shared_Taxa = shared_otus))
        } else {
          return(NULL)
        }
      })

      # Write combined results to file
      shared_taxa_df <- do.call(rbind, shared_taxa_list)
      if (is.null(shared_taxa_df) || nrow(shared_taxa_df) == 0) {
        write.table(data.frame(Group1 = character(), Group2 = character(), Shared_Taxa = character()), file, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
      } else {
        write.table(shared_taxa_df, file, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
      }
    }
  )

  # Download handler for shared taxa across all groups
  output$download_shared_all_taxa <- downloadHandler(
    filename = function() {
      paste0("shared_all_taxa_", Sys.Date(), ".tsv")
    },
    content = function(file) {
      group_otus <- group_otus_reactive()

      # Identify shared taxa across all groups (intersection of all lists)
      shared_taxa <- Reduce(intersect, group_otus)

      # Convert to a data frame and write to TSV
      if (length(shared_taxa) == 0) {
        write.table(data.frame(Taxa = character()), file, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
      } else {
        shared_taxa_df <- data.frame(Shared_Taxa = shared_taxa)
        write.table(shared_taxa_df, file, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
      }
    }
  )




  observeEvent(input$explain_venn_diagram_plot, {
    showModal(modalDialog(
      title = "Venn Diagram Explanation",
      p("This Venn diagram displays the unique and shared taxa across multiple sample groups."),
      p("The analysis is performed at the taxonomy level selected by the user (e.g., Genus, Species). Additionally, core taxa are filtered based on the detection and prevalence thresholds provided."),

      p("Here's how to interpret and utilize the results:"),

      h4("Unique Taxa"),
      p("Unique taxa refer to those that are found exclusively in a particular sample group and not shared with any other group. You can download this data using the 'Download Unique Taxa' button. This is useful when trying to identify biomarkers or specific organisms that are distinct to a condition or sample type."),

      h4("Pairwise Shared Taxa"),
      p("This identifies the taxa shared between pairs of groups. Use the 'Download Pairwise Shared Taxa' button to obtain this data. It's valuable for understanding the common microbiota between two specific groups and can help in comparative studies."),

      h4("Shared Taxa Across All Groups"),
      p("These taxa are shared across all sample groups. Use the 'Download Shared Across All Groups' button to download this data. Identifying the core microbiota present in all groups can be helpful for finding broadly conserved organisms across different conditions or environments."),

      h4("Venn Diagram Visualization"),
      p("The Venn diagram provides a visual representation of the shared and unique taxa between the sample groups. Each section of the diagram corresponds to different sets of shared taxa between the groups. The overlap between two or more groups shows the taxa that are shared between them."),

      p("This analysis is particularly useful when comparing microbial communities from different environments, treatment groups, or conditions. By visualizing and downloading these insights, researchers can easily identify both common and unique taxa among their samples."),

      easyClose = TRUE,
      footer = modalButton("Close")
    ))
  })



  ########### compositional barplot

  compositional_barplot_reactive <- reactive({
    req(current_physeq(), input$hueRank, input$shadeRank, input$facetBy)

    hueRank <- input$hueRank
    shadeRank <- input$shadeRank

    pseq2 <- current_physeq() %>%
      microViz::tax_sort(by = sum, at = shadeRank) %>%
      microViz::tax_sort(by = sum, at = hueRank) %>%
      microViz::tax_agg(rank = shadeRank)

    nHues <- input$nHues
    nShades <- input$nShades

    hierarchicalPalInfo <- data.frame(
      hue = as.vector(tt_get(pseq2)[, hueRank]),
      shade = as.vector(tt_get(pseq2)[, shadeRank]),
      counts = phyloseq::taxa_sums(otu_get(pseq2))
    )

    hierarchicalPalInfo <- hierarchicalPalInfo %>%
      dplyr::mutate(
        hue = forcats::fct_other(
          f = hue, keep = unique(hue)[seq_len(nHues)],
          other_level = paste("Other", hueRank)
        ),
        nChrHue = nchar(as.character(hue)), padHue = max(nChrHue) - nChrHue
      ) %>%
      dplyr::group_by(hue) %>%
      dplyr::mutate(
        shade = forcats::fct_other(
          f = shade, keep = unique(shade)[seq_len(nShades)],
          other_level = "Other"
        )
      ) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(
        nChrShade = nchar(as.character(shade)), padShade = max(nChrShade) - nChrShade,
        Taxa = paste0(hue, ": ", strrep(" ", padHue), shade, strrep(" ", padShade))
      )

    unique_taxa <- unique(hierarchicalPalInfo$Taxa)
    palette_colors <- colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(length(unique_taxa))
    names(palette_colors) <- unique_taxa

    facet_by <- if (input$facetBy == "None") NULL else input$facetBy

    p <- pseq2 %>%
      ps_get() %>%
      tax_mutate("Taxa" = hierarchicalPalInfo$Taxa, .keep = "none") %>%
      microViz::comp_barplot(
        tax_level = "Taxa", n_taxa = length(palette_colors),
        facet_by = facet_by,
        tax_order = "asis", palette = palette_colors, bar_width = 0.975
      ) +
      coord_flip() +
      theme(legend.text = element_text(family = "mono"))

    ggplotly(p)
  })

  #render the plot UI
  output$compositionalBarplot <- renderPlotly({
    compositional_barplot_reactive()
  })

  # Implement download handler for the compositional barplot
  output$download_compositionalBarplot <- downloadHandler(
    filename = function() {
      paste0("compositional_barplot_", Sys.Date(), ".", input$filetype)
    },
    content = function(file) {
      # Save the plotly object to an HTML file first
      tempfile <- tempfile(fileext = ".html")
      p <- compositional_barplot_reactive()
      htmlwidgets::saveWidget(as_widget(p), tempfile)

      # Convert the HTML to the selected format using user-specified dimensions
      webshot::webshot(tempfile, file = file,
                       vwidth = input$plot_width,
                       vheight = input$plot_height,
                       delay = 0.2)
    }
  )

  # Observer for showing the modal with the explanation for the compositional barplot
  observeEvent(input$show_compositional_explanation, {
    showModal(modalDialog(
      title = "Compositional Barplot Explanation",
      HTML("<p>The Compositional Barplot visualizes the distribution of different taxa across samples, organized by selected taxonomic ranks.</p>
              <p><strong>Key Features:</strong></p>
              <ul>
                <li><strong>Hue Rank:</strong> Select the taxonomic level for the main color groups (e.g., Phylum).</li>
                <li><strong>Shade Rank:</strong> Choose a more detailed taxonomic level for sub-grouping within each color group (e.g., Genus).</li>
                <li><strong>Faceting:</strong> Split the barplot by a sample variable to compare groups (e.g., treatment groups).</li>
                <li><strong>Number of Hues/Shades:</strong> Control how many taxa are shown, allowing for focus on the most abundant groups.</li>
              </ul>
              <p>This plot helps in understanding the relative abundances of different taxa and comparing these distributions across samples or groups.</p>"),
      easyClose = TRUE,
      footer = modalButton("Close")
    ))
  })




  ########### heatmap

  observeEvent(input$plotHeatmap, {
    req(current_physeq(), input$normalizationMethod, input$heatmapRank, input$annotationColumn1, input$annotationColumn2, input$topTaxa)

    # Conditionally aggregate the phyloseq object based on the user's selection
    aggregated_physeq <- reactive({
      if (input$heatmapRank == "ASV") {
        # Only normalize the physeq object without aggregation
        microbiome::transform(current_physeq(), input$normalizationMethod)
      } else {
        # Both normalize and aggregate at the selected taxonomic level
        microViz::tax_transform(current_physeq(), input$normalizationMethod, rank = input$heatmapRank)
      }
    })

    # Ensure aggregated_physeq has data
    pseq <- aggregated_physeq()
    req(nrow(otu_table(pseq)) > 0)

    # Select top taxa
    top_taxa <- names(sort(phyloseq::taxa_sums(pseq), decreasing = TRUE))[1:input$topTaxa]
    psq_normalized_pruned <- phyloseq::prune_taxa(top_taxa, pseq)

    # Define annotation colors dynamically based on the number of unique values
    unique_vals1 <- unique(phyloseq::sample_data(psq_normalized_pruned)[[input$annotationColumn1]])
    cols1 <- distinct_palette(n = length(unique_vals1), add = NA)
    names(cols1) <- unique_vals1

    unique_vals2 <- unique(phyloseq::sample_data(psq_normalized_pruned)[[input$annotationColumn2]])
    cols2 <- distinct_palette(n = length(unique_vals2), add = NA)
    names(cols2) <- unique_vals2


    # Prepare sample annotations, with legend title
    sample_anno <- sampleAnnotation(
      State1 = anno_sample_cat(input$annotationColumn1, legend_title = input$annotationColumn1),
      col = list(State1 = cols1, State2 = cols2), border = FALSE,
      State2 = anno_sample_cat(input$annotationColumn2, col = cols2, legend_title = input$annotationColumn2),
      annotation_label = c(input$annotationColumn1, input$annotationColumn2)
    )

    # Check if ComplexHeatmap and InteractiveComplexHeatmap are available
    if (!requireNamespace("ComplexHeatmap", quietly = TRUE) || !requireNamespace("InteractiveComplexHeatmap", quietly = TRUE)) {
      showNotification("The required packages 'ComplexHeatmap' and 'InteractiveComplexHeatmap' are not installed. Please install them to use this feature.", type = "error")
      return()
    }

    if (requireNamespace("microViz", quietly = TRUE)) {
      # Use microViz-related functionality
    } else {
      stop("The 'microViz' package is required for this feature. Please install it.")
    }

    # Generate the heatmap using comp_heatmap
    heatmap_obj <- microViz::comp_heatmap(
      psq_normalized_pruned,
      taxa = top_taxa,
      tax_anno = taxAnnotation(Prev. = anno_tax_prev(bar_width = 0.3, size = grid::unit(1, "cm"))),
      sample_anno = sample_anno,
      sample_seriation = "OLO_ward",
      tax_seriation = "OLO_ward",
      colors = heat_palette(palette = "Rocket", rev = TRUE),
      cluster_rows = input$clusterRows,
      cluster_columns = input$clusterColumns,
      sample_names_show = TRUE
    )

    # Draw the heatmap for interactive use, with annotation legends
    ht1 <- ComplexHeatmap::draw(
      object = heatmap_obj,
      annotation_legend_list = attr(heatmap_obj, "AnnoLegends"), merge_legends = TRUE
    )

    # Render the interactive heatmap using InteractiveComplexHeatmap
    InteractiveComplexHeatmap::makeInteractiveComplexHeatmap(input, output, session, ht1,
                                  heatmap_id = "relativeAbundanceHeatmap")
  })


  # Server-side code to show the explanation
  observeEvent(input$show_heatmap_explanation, {
    showModal(modalDialog(
      title = "Heatmap Explanation",
      HTML("<p>The Heatmap is a visualization tool used to display the relative abundance of taxa across samples, with colors representing the abundance levels.</p>
              <p><strong>Key Features:</strong></p>
              <ul>
                <li><strong>Normalization Method:</strong> The data can be normalized using different methods to adjust for varying sample sizes or sequencing depths.</li>
                <li><strong>Taxonomic Level:</strong> You can aggregate the data at different taxonomic levels (e.g., Genus, Family) for visualization.</li>
                <li><strong>Top Taxa Selection:</strong> The heatmap can focus on the most abundant taxa, helping you to identify key players in your data.</li>
                <li><strong>Annotations:</strong> You can annotate the heatmap with sample metadata, such as groupings or experimental conditions.</li>
                <li><strong>Interactivity:</strong> The zoomed heatmap on the right allows for detailed inspection of specific regions of interest.</li>
              </ul>
              <p>This visualization helps in identifying patterns and differences in microbial communities across samples.</p>"),
      easyClose = TRUE,
      footer = modalButton("Close")
    ))
  })




  ###### Alpha diversity

  # Reactive expression to store the ggplot2 and plotly plots
  plotly_plots <- reactiveVal(list())
  ggplot_plots <- reactiveVal(list())

  # Observing for available metrics and updating inputs
  observe({
    req(current_physeq())
    pseq <- current_physeq()

    # Calculate diversity metrics
    tab <- microbiome::alpha(pseq, index = "all")

    # Ensure that metrics are available
    available_metrics <- colnames(tab)
    if (length(available_metrics) == 0) {
      showNotification("No diversity metrics available. Please check your data.", type = "error")
      return()
    }

    # Add diversity metrics to sample data
    meta <- as(phyloseq::sample_data(pseq), "data.frame")
    for (metric in available_metrics) {
      meta[[metric]] <- tab[[metric]]
    }
    phyloseq::sample_data(pseq) <- phyloseq::sample_data(meta)  # Update sample data in phyloseq

    # Update selectInput choices
    updateSelectInput(session, "alphaMetric", choices = available_metrics, selected = available_metrics[1])
    updateSelectInput(session, "alphaGroupingColumn", choices = sample_variables(pseq))

    # Store updated phyloseq in a reactive value if needed
    filtered_pseq(pseq)
  })

  # Generate Alpha Diversity Plot
  observeEvent(input$plotAlpha, {
    req(current_physeq(), input$alphaMetric, input$alphaGroupingColumn)
    pseq <- current_physeq()

    meta <- as(phyloseq::sample_data(pseq), "data.frame")
    plotly_lst <- list()
    ggplot_lst <- list()

    clusters <- levels(factor(meta[[input$alphaGroupingColumn]]))
    comparisons <- combn(seq_along(clusters), 2, simplify = FALSE, FUN = function(i) clusters[i])

    custom_colors <- c("#0072B2", "#009E73", "#F0E442", "#56B4E9", "#D55E00", "#CC79A7",
                       "indianred", "darkblue", "seagreen", "skyblue", "lightsteelblue2",
                       "plum2", "gray", "indianred", "blue")

    for (metric in input$alphaMetric) {
      if (metric %in% colnames(meta)) {  # Check if metric exists in meta
        p <- ggplot2::ggplot(meta, ggplot2::aes_string(x = input$alphaGroupingColumn, y = metric, color = input$alphaGroupingColumn)) +
          geom_boxplot(outlier.shape = NA, width = 0.4, alpha = 0.75) +
          geom_jitter(height = 0, width = 0.2, alpha = 0.5, size = 2) +
          scale_color_manual(values = custom_colors) +
          theme_minimal(base_size = 15) +
          labs(y = metric, x = "Groups") +
          theme(
            legend.position = "right",
            plot.title = element_text(hjust = 0.5),
            axis.title.x = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1),
            plot.margin = margin(t = 40, r = 40, b = 40, l = 60, unit = "pt")
          )

        # Store ggplot and plotly versions
        ggplot_lst[[metric]] <- p
        plotly_lst[[metric]] <- ggplotly(p)
      }
    }

    # Update the reactive values
    ggplot_plots(ggplot_lst)
    plotly_plots(plotly_lst)
  })

  # Render the Alpha Diversity Plot
  output$alphaDiversityPlot <- renderPlotly({
    plots <- plotly_plots()
    req(length(plots) > 0)

    if (length(plots) == 1) {
      plots[[1]] %>%
        plotly::layout(dragmode = "zoom", margin = list(l = 120, r = 40, t = 40, b = 40), showlegend = TRUE) %>%
        plotly::config(displayModeBar = TRUE, scrollZoom = TRUE)
    } else {
      plotly::subplot(plots, nrows = ceiling(length(plots) / 2), shareX = TRUE, shareY = FALSE, titleX = TRUE, titleY = TRUE) %>%
        plotly::layout(dragmode = "zoom", margin = list(l = 120, r = 40, t = 40, b = 40), showlegend = TRUE) %>%
        plotly::config(displayModeBar = TRUE, scrollZoom = TRUE)
    }
  })

  # Download handler for the alpha diversity plot
  output$download_alphaDiversityPlot <- downloadHandler(
    filename = function() {
      paste0("alpha_diversity_", Sys.Date(), ".", input$alpha_filetype)
    },
    content = function(file) {
      if (input$alpha_filetype == "html") {
        # Save Plotly as HTML
        plots <- plotly_plots()  # Ensure this function returns valid Plotly objects
        p <- do.call(subplot, c(plots, list(nrows = ceiling(length(plots) / 2), shareX = TRUE, shareY = FALSE)))
        htmlwidgets::saveWidget(as_widget(p), file)
      } else {
        # Save ggplot objects as other formats
        plots <- ggplot_plots()  # Ensure this function returns valid ggplot objects

        # Check if `plots` is a list of ggplot objects
        if (length(plots) > 1) {
          p <- marrangeGrob(plots, ncol = 2, nrow = ceiling(length(plots) / 2))
        } else {
          p <- plots[[1]]  # Single plot case
        }

        # Save the plot with ggsave
        ggsave(
          filename = file,
          plot = p,
          width = input$plot_width,
          height = input$plot_height,
          units = "in",
          device = input$alpha_filetype,  # Ensure this matches supported devices
          limitsize = FALSE  # Avoid issues with large plots
        )
      }
    }
  )


  # Reactive expression to store the Wilcoxon test results
  alpha_stats <- reactiveVal(data.frame())

  # Calculate Wilcoxon test statistics when the button is clicked
  observeEvent(input$showStats, {
    print("showStats button clicked")  # Check if button click is detected
    req(current_physeq(), input$alphaMetric, input$alphaGroupingColumn)

    pseq <- current_physeq()
    meta <- as(phyloseq::sample_data(pseq), "data.frame")

    # Define clusters and pairwise comparisons
    clusters <- levels(factor(meta[[input$alphaGroupingColumn]]))
    comparisons <- combn(seq_along(clusters), 2, simplify = FALSE, FUN = function(i) clusters[i])

    stats_list <- list()

    # Calculate statistics for each selected metric
    for (metric in input$alphaMetric) {
      # Initialize p-values with NA in case there are issues with the comparisons
      p_values <- sapply(comparisons, function(comp) {
        group1 <- meta[meta[[input$alphaGroupingColumn]] == comp[1], metric]
        group2 <- meta[meta[[input$alphaGroupingColumn]] == comp[2], metric]

        # Check that both groups have data before running the test
        if (length(group1) > 0 && length(group2) > 0) {
          wilcox.test(group1, group2)$p.value
        } else {
          NA  # Return NA if any group lacks data
        }
      })

      # Apply BH correction to the p-values, handling any NA values
      adj_p_values <- p.adjust(p_values, method = "BH", n = sum(!is.na(p_values)))

      # Create a data frame for this metric
      stats_df <- data.frame(
        Metric = metric,
        Comparison = sapply(comparisons, function(comp) paste(comp, collapse = " vs ")),
        P_value = round(p_values, 4),  # Round p-values to 4 decimal places
        Adjusted_P_value = round(adj_p_values, 4),  # Round adjusted p-values to 4 decimal places
        stringsAsFactors = FALSE
      )

      print(stats_df)  # Debugging to check if stats_df has data for each metric

      stats_list[[metric]] <- stats_df
    }

    # Combine all metrics into one data frame
    full_stats <- do.call(rbind, stats_list)
    print("Full stats data frame:")  # Debugging statement
    print(full_stats)  # Check if full_stats has data before saving to alpha_stats

    # Update the reactive value
    alpha_stats(full_stats)
  })

  # Render the statistics table using DT
  output$statsTable <- DT::renderDT({
    req(alpha_stats())
    DT::datatable(alpha_stats(), options = list(pageLength = 10, autoWidth = TRUE),
                  editable = 'cell', rownames = FALSE,
                  caption = 'Table 1: Wilcoxon Rank-Sum Test Results: Pairwise comparisons between groups with p-values adjusted using the Benjamini-Hochberg method.')
  })


  # Download handler for the statistics in TSV format
  output$downloadStats <- downloadHandler(
    filename = function() {
      paste0("alpha_diversity_stats_", Sys.Date(), ".tsv")
    },
    content = function(file) {
      stats <- alpha_stats()
      req(stats)  # Ensure stats are available before downloading

      # Check if stats is empty
      if (nrow(stats) == 0) {
        showNotification("No data available for download.", type = "error")
        return()
      }

      # Debugging information
      print("Writing stats to file...")
      print(stats)  # Check if stats has data

      # Write the stats to a TSV file
      write.table(stats, file, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

      print("File written successfully")
    }
  )




  # Server-side code to show the explanation for Alpha Diversity
  observeEvent(input$show_alpha_explanation, {
    showModal(modalDialog(
      title = "Alpha Diversity Explanation",
      HTML("<p><strong>Alpha Diversity</strong> is a measure of the diversity within a single sample or environment.
           It provides insight into the richness (number of different taxa) and evenness (distribution of taxa abundances) of the microbial community.</p>
          <p><strong>Key Features:</strong></p>
          <ul>
            <li><strong>Diversity Metrics:</strong> Alpha function from microbiome package was used to calculate various measures of richness, evenness, diversity, dominance, and rarity with default parameters,
            including Shannon, Simpson, and Observed Species. Each metric emphasizes different aspects of diversity.</li>
            <li><strong>Grouping Column:</strong> Samples are grouped based on metadata, allowing for comparisons of alpha diversity across different conditions or treatments.</li>
            <li><strong>Statistical Test:</strong> Pairwise comparisons between groups are performed using the Wilcoxon Rank-Sum Test,
            and the p-values are adjusted using the Benjamini-Hochberg (BH) method to control the false discovery rate.</li>
            <li><strong>Plot Interpretation:</strong> The boxplots display the distribution of diversity metrics across groups.
            The p-values indicate whether the observed differences between groups are statistically significant.</li>
          </ul>
          <p>This analysis helps in understanding the complexity and structure of microbial communities within individual samples and across different conditions.</p>"),
      easyClose = TRUE,
      footer = modalButton("Close")
    ))
  })






  ######## Beta diversity


  # UI Updates for Taxonomic Ranks

  observe({
    req(current_physeq())
    pseq <- current_physeq()

    # Update taxRank choices
    updateSelectInput(session, "taxRank", choices = rank_names(pseq))

    # Update other dropdowns based on sample metadata
    updateSelectInput(session, "groupingColumnBeta_plotting", choices = colnames(phyloseq::sample_data(pseq)))
    updateSelectInput(session, "groupingColumnBeta_statistics", choices = colnames(phyloseq::sample_data(pseq)))
    updateSelectInput(session, "shapeColumnBeta", choices = colnames(phyloseq::sample_data(pseq)))
    updateSelectInput(session, "selectedVariables", choices = colnames(phyloseq::sample_data(pseq)))
    updateSelectInput(session, "strata", choices = colnames(phyloseq::sample_data(pseq)))
  })





  # Helper function to calculate distance and transform for each method (excluding gunifrac)
  get_transformed_pseq <- function(pseq, distanceMetric, normalizationMethod, taxRank) {
    if (distanceMetric == "robust.aitchison" && normalizationMethod != "identity") {
      stop("'Aitchison' distance requires count data. Please select 'identity' as the normalization method or choose another distance metric.")
    }

    if (distanceMetric == "jaccard") {
      pseq_transformed <- pseq %>%
        microViz::tax_transform(normalizationMethod, rank = taxRank) %>%
        microViz::dist_calc(dist = "jaccard", binary = TRUE)
    } else if (distanceMetric == "bray") {
      pseq_transformed <- pseq %>%
        microViz::tax_transform(normalizationMethod, rank = taxRank) %>%
        microViz::dist_calc(dist = "bray", binary = FALSE)
    } else if (distanceMetric == "robust.aitchison") {
      pseq_transformed <- pseq %>%
        microViz::tax_transform("identity", rank = taxRank) %>%
        microViz::dist_calc(dist = "robust.aitchison")
    } else {
      stop("Unsupported distance metric selected.")
    }

    return(pseq_transformed)
  }

  # Reactive values for plot objects
  ggplot_obj <- reactiveVal(NULL)
  ggplot_pca_obj <- reactiveVal(NULL)
  ggplot_nmds_obj <- reactiveVal(NULL)

  # Plotting: PCoA
  observeEvent(input$plotPCoA, {
    req(current_physeq(), input$distanceMetric, input$normalizationMethod, input$taxRank, input$groupingColumnBeta_plotting, input$shapeColumnBeta)
    pseq <- current_physeq()

    tryCatch({
      pseq_transformed <- get_transformed_pseq(pseq, input$distanceMetric, input$normalizationMethod, input$taxRank)
      ordination_res <- microViz::ord_calc(pseq_transformed, method = "PCoA")
      plot <- microViz::ord_plot(ordination_res, color = input$groupingColumnBeta_plotting, shape = input$shapeColumnBeta, plot_taxa = 1:5) +
        ggplot2::stat_ellipse(ggplot2::aes_string(fill = input$groupingColumnBeta_plotting, color = input$groupingColumnBeta_plotting), linetype = 1, segments = 10, lwd = 1.2, alpha = 0.25, level = 0.6) +
        ggplot2::scale_colour_brewer(palette = "Dark2", aesthetics = c("fill", "colour")) +
        ggplot2::theme_bw() +
        ggside::geom_xsideboxplot(ggplot2::aes_string(fill = input$groupingColumnBeta_plotting, y = input$groupingColumnBeta_plotting), orientation = "y") +
        ggside::geom_ysideboxplot(ggplot2::aes_string(fill = input$groupingColumnBeta_plotting, x = input$groupingColumnBeta_plotting), orientation = "x") +
        ggside::scale_xsidey_discrete(labels = NULL) +
        ggside::scale_ysidex_discrete(labels = NULL) +
        ggside::theme_ggside_void()

      ggplot_obj(plot)
      output$betaDiversityPlot <- renderPlot({ print(ggplot_obj()) })
    }, error = function(e) {
      output$betaDiversityNote <- renderText({ e$message })
      showModal(modalDialog(title = "Error", e$message, easyClose = TRUE, footer = modalButton("Close")))
    })
  })

  output$download_PCoAPlot <- downloadHandler(
    filename = function() { paste0("PCoA_beta_diversity_", Sys.Date(), ".", input$beta_filetype) },
    content = function(file) { req(ggplot_obj()); ggsave(file, plot = ggplot_obj(), width = input$plot_width, height = input$plot_height, units = "in", device = input$beta_filetype) }
  )

  # Plotting: PCA
  observeEvent(input$plotPCA, {
    req(current_physeq(), input$normalizationMethod, input$taxRank, input$groupingColumnBeta_plotting, input$shapeColumnBeta)
    pseq <- current_physeq()

    tryCatch({
      ordination_res <- pseq %>% microViz::tax_transform(input$normalizationMethod, rank = input$taxRank) %>% microViz::ord_calc(method = "PCA")
      plot <- microViz::ord_plot(ordination_res, color = input$groupingColumnBeta_plotting, shape = input$shapeColumnBeta, plot_taxa = 1:5) +
        ggplot2::stat_ellipse(ggplot2::aes_string(fill = input$groupingColumnBeta_plotting, color = input$groupingColumnBeta_plotting), linetype = 1, segments = 10, lwd = 1.2, alpha = 0.25, level = 0.6) +
        ggplot2::scale_colour_brewer(palette = "Dark2", aesthetics = c("fill", "colour")) +
        ggplot2::theme_bw() +
        ggside::geom_xsideboxplot(ggplot2::aes_string(fill = input$groupingColumnBeta_plotting, y = input$groupingColumnBeta_plotting), orientation = "y") +
        ggside::geom_ysideboxplot(ggplot2::aes_string(fill = input$groupingColumnBeta_plotting, x = input$groupingColumnBeta_plotting), orientation = "x") +
        ggside::scale_xsidey_discrete(labels = NULL) +
        ggside::scale_ysidex_discrete(labels = NULL) +
        ggside::theme_ggside_void()

      ggplot_pca_obj(plot)
      output$betaDiversityPlot <- renderPlot({ print(ggplot_pca_obj()) })
    }, error = function(e) {
      output$betaDiversityNote <- renderText({ e$message })
      showModal(modalDialog(title = "Error", e$message, easyClose = TRUE, footer = modalButton("Close")))
    })
  })

  output$download_PCAPlot <- downloadHandler(
    filename = function() { paste0("PCA_beta_diversity_", Sys.Date(), ".", input$beta_filetype) },
    content = function(file) { req(ggplot_pca_obj()); ggsave(file, plot = ggplot_pca_obj(), width = input$plot_width, height = input$plot_height, units = "in", device = input$beta_filetype) }
  )


  # Plotting: NMDS
  observeEvent(input$plotNMDS, {
    req(current_physeq(), input$distanceMetric, input$normalizationMethod, input$taxRank, input$groupingColumnBeta_plotting, input$shapeColumnBeta)
    pseq <- current_physeq()

    tryCatch({
      # Transform the data
      pseq_transformed <- get_transformed_pseq(pseq, input$distanceMetric, input$normalizationMethod, input$taxRank)

      # Perform NMDS ordination
      ordination_res <- microViz::ord_calc(pseq_transformed, method = "NMDS")

      # Print the stress value
      if ("stress" %in% names(ordination_res)) {  # Check if the result has a stress attribute
        print(ordination_res$stress)
      } else if (inherits(ordination_res, "metaMDS")) {  # Handle specific NMDS object types
        print(ordination_res$stress)
      } else {
        print("Unable to retrieve stress value. Object does not contain a recognizable stress attribute.")
      }

      # Generate NMDS plot
      plot <- microViz::ord_plot(
        ordination_res,
        color = input$groupingColumnBeta_plotting,
        shape = input$shapeColumnBeta,
        plot_taxa = 1:5
      ) +
        ggplot2::stat_ellipse(
          ggplot2::aes_string(fill = input$groupingColumnBeta_plotting, color = input$groupingColumnBeta_plotting),
          linetype = 1,
          segments = 10,
          lwd = 1.2,
          alpha = 0.25,
          level = 0.6
        ) +
        ggplot2::scale_colour_brewer(palette = "Dark2", aesthetics = c("fill", "colour")) +
        ggplot2::theme_bw() +
        ggside::geom_xsideboxplot(
          ggplot2::aes_string(fill = input$groupingColumnBeta_plotting, y = input$groupingColumnBeta_plotting),
          orientation = "y"
        ) +
        ggside::geom_ysideboxplot(
          ggplot2::aes_string(fill = input$groupingColumnBeta_plotting, x = input$groupingColumnBeta_plotting),
          orientation = "x"
        ) +
        ggside::scale_xsidey_discrete(labels = NULL) +
        ggside::scale_ysidex_discrete(labels = NULL) +
        ggside::theme_ggside_void()

      # Add custom shape palette if needed
      plot <- plot + ggplot2::scale_shape_manual(values = c(0:6, 15:20))  # Extend shapes if >6 groups


      # Render the plot
      ggplot_nmds_obj(plot)
      output$betaDiversityPlot <- renderPlot({ print(ggplot_nmds_obj()) })
    }, error = function(e) {
      output$betaDiversityNote <- renderText({ e$message })
      showModal(modalDialog(title = "Error", e$message, easyClose = TRUE, footer = modalButton("Close")))
    })
  })

  output$download_NMDSPlot <- downloadHandler(
    filename = function() { paste0("NMDS_beta_diversity_", Sys.Date(), ".", input$beta_filetype) },
    content = function(file) { req(ggplot_nmds_obj()); ggsave(file, plot = ggplot_nmds_obj(), width = input$plot_width, height = input$plot_height, units = "in", device = input$beta_filetype)
    }
  )




  # ########### statistics tab

  observeEvent(input$calculatePERMANOVA, {
    req(current_physeq(), input$groupingColumnBeta_statistics, input$distanceMetric, input$taxRank)

    pseq <- current_physeq()
    tryCatch({

      # Convert phyloseq object to relative abundance
      pseq_relabund <- microbiome::transform(pseq, "compositional")
      pseq_relabund <- aggregate_taxa(pseq_relabund, input$taxRank)
      otu_table <- as.matrix(phyloseq::otu_table(pseq_relabund))

      if (phyloseq::taxa_are_rows(pseq_relabund)) {
        otu_table <- t(otu_table)
      }

      meta <- data.frame(phyloseq::sample_data(pseq))  # Metadata as data frame

      # Validate grouping column exists in metadata
      validate(need(input$groupingColumnBeta_statistics %in% colnames(meta),
                    "Grouping column not found in metadata."))

      # Extract strata column if specified
      strata_var <- if (input$strata != "None") meta[[input$strata]] else NULL
      # Dynamically construct formula
      formula <- as.formula(paste("otu_table ~", input$groupingColumnBeta_statistics))


      # PERMANOVA Analysis
      set.seed(12346)
      permanova_res <- vegan::adonis2(
        formula = formula,
        data = meta,
        permutations = 999,
        method = input$distanceMetric,
        strata = strata_var
      )

      # Adjust p-values using Benjamini-Hochberg (FDR) method
      permanova_res$`Pr(>F)` <- p.adjust(permanova_res$`Pr(>F)`, method = "BH")
      # Render PERMANOVA results using `verbatimTextOutput`
      output$PERMANOVAResultsUI <- renderUI({
        tagList(
          h4("PERMANOVA Results"),
          verbatimTextOutput("permanovaResults")
        )
      })

      output$permanovaResults <- renderPrint({
        print(permanova_res)
      })

    }, error = function(e) {
      # Show error message
      showModal(modalDialog(
        title = "Error in PERMANOVA Analysis",
        e$message,
        easyClose = TRUE,
        footer = modalButton("Close")
      ))
    })
  })


  #### for dispersion test
  # Reactive value to track which button is pressed
  dispersion_output_type <- reactiveVal(NULL)

  # Dispersion Test Analysis
  observeEvent(input$calculateDispersionTest, {
    req(current_physeq(), input$distanceMetric, input$groupingColumnBeta_statistics, input$taxRank)

    pseq <- current_physeq()
    validate(need(!is.null(pseq), "Phyloseq object is missing."))
    validate(need(input$groupingColumnBeta_statistics %in% colnames(phyloseq::sample_data(pseq)), "Selected grouping column not found."))

    tryCatch({
      # Convert phyloseq object to relative abundance
      pseq_relabund <- microbiome::transform(pseq, "compositional")
      pseq_relabund <- aggregate_taxa(pseq_relabund, input$taxRank)
      otu_table <- as.matrix(phyloseq::otu_table(pseq_relabund))

      if (phyloseq::taxa_are_rows(pseq_relabund)) {
        otu_table <- t(otu_table)
      }

      # Distance calculation
      dist <- vegan::vegdist(otu_table, method = input$distanceMetric)
      grouping_var <- phyloseq::sample_data(pseq_relabund)[[input$groupingColumnBeta_statistics]]

      # Beta-dispersion
      dispersion <- vegan::betadisper(dist, grouping_var)

      # Update the reactive value
      dispersion_output_type("results")

      # Render ANOVA and permutation results
      output$dispersionTestResultsUI <- renderUI({
        tagList(
          h4("Dispersion Test Results"),
          verbatimTextOutput("dispersionAnova"),
          verbatimTextOutput("dispersionPermTest")
        )
      })

      output$dispersionAnova <- renderPrint({
        print(anova(dispersion))
      })
      output$dispersionPermTest <- renderPrint({
        print(vegan::permutest(dispersion))
      })
    }, error = function(e) {
      showModal(modalDialog(
        title = "Error in Dispersion Test Analysis",
        e$message,
        easyClose = TRUE, footer = modalButton("Close")
      ))
    })
  })



  # Centroid Plot Analysis
  observeEvent(input$plotCentroids, {
    req(current_physeq(), input$distanceMetric, input$groupingColumnBeta_statistics, input$taxRank)

    pseq <- current_physeq()
    validate(need(!is.null(pseq), "Phyloseq object is missing."))
    validate(need(input$groupingColumnBeta_statistics %in% colnames(phyloseq::sample_data(pseq)), "Selected grouping column not found."))

    tryCatch({
      # Convert phyloseq object to relative abundance
      pseq_relabund <- microbiome::transform(pseq, "compositional")
      #pseq_relabund <- aggregate_taxa(pseq_relabund, input$taxRank)
      otu_table <- as.matrix(phyloseq::otu_table(pseq_relabund))

      if (phyloseq::taxa_are_rows(pseq_relabund)) {
        otu_table <- t(otu_table)
      }

      # Distance calculation
      dist <- vegan::vegdist(otu_table, method = input$distanceMetric)
      grouping_var <- phyloseq::sample_data(pseq_relabund)[[input$groupingColumnBeta_statistics]]

      # Beta-dispersion
      dispersion <- vegan::betadisper(dist, grouping_var)

      # Update the reactive value
      dispersion_output_type("plot")

      # Render centroid plot
      output$centroidPlotUI <- renderUI({
        tagList(
          h4("Centroid Plot"),
          plotOutput("dispersionPlot", height = "400px")
        )
      })

      output$dispersionPlot <- renderPlot({
        plot(dispersion, ellipse = TRUE, hull = TRUE, main = "Centroids Plot")
      })
    }, error = function(e) {
      showModal(modalDialog(
        title = "Error in Dispersion Test Analysis",
        e$message,
        easyClose = TRUE, footer = modalButton("Close")
      ))
    })
  })


  # Server-side code to show the explanation for Beta Diversity
  observeEvent(input$show_beta_explanation, {
    showModal(modalDialog(
      title = "Beta Diversity Explanation",
      HTML("<p><strong>Beta Diversity</strong> is a measure of the differences between microbial communities from different environments or samples.</p>
         <p><strong>Key Concepts:</strong></p>
         <ul>
           <li><strong>Normalization Methods:</strong> Various normalization methods like CLR, log10, and compositional transformations are used to standardize the data before distance calculation.</li>
           <li><strong>Distance Metrics:</strong> Metrics like Bray-Curtis, Jaccard, and Aitchison measure the dissimilarity between samples, capturing different aspects of community variation.</li>
           <li><strong>Ordination Techniques:</strong> PCoA, PCA, and NMDS are ordination methods that reduce the dimensionality of the data, making it easier to visualize and interpret the differences between samples.</li>
           <li><strong>Group Comparisons:</strong> Samples are grouped by metadata columns, allowing comparisons of beta diversity across different conditions or treatments.</li>
         </ul>
         <p>Understanding beta diversity helps in exploring how microbial communities differ across samples and identifying potential factors influencing these differences.</p>"),
      easyClose = TRUE,
      footer = modalButton("Close")
    ))
  })




  ############# differential abundance


  observe({
    req(current_physeq())
    metadata_columns <- colnames(phyloseq::sample_data(current_physeq()))
    updateSelectizeInput(session, "fixed_variables", choices = metadata_columns, selected = NULL)
    updateSelectizeInput(session, "random_variables", choices = metadata_columns, selected = NULL)
  })


  # Define roots for folder navigation
  roots <- c(Home = fs::path_home())

  # Server logic for folder selection
  observe({
    shinyFiles::shinyDirChoose(input, "select_folder", roots = roots, session = session, allowDirCreate = TRUE) # Enable directory creation
  })

  observeEvent(input$select_folder, {
    # Parse the selected directory
    req(input$select_folder)
    selected_dir <- parseDirPath(roots, input$select_folder)

    # Update the text input with the selected folder path
    updateTextInput(session, "Maaslin2_output", value = selected_dir)

    # Display the final selected folder path
    output$selected_folder <- renderText({
      paste("Selected Folder: ", selected_dir)
    })
  })

  # Optional: Render content of the selected directory for debugging
  output$selected_folder_files <- renderPrint({
    req(input$select_folder)
    selected_dir <- parseDirPath(roots, input$select_folder)
    list.files(selected_dir)
  })



  # Reference Variable and Level Logic

  # Populate reference variable dropdown based on fixed variables
  observe({
    req(input$fixed_variables) # Ensure fixed variables are selected
    updateSelectInput(
      session,
      "reference_variable",
      choices = input$fixed_variables, # Use fixed variables for selection
      selected = NULL
    )
  })

  # Populate reference levels dropdown based on the selected reference variable
  observeEvent(input$reference_variable, {
    req(current_physeq(), input$reference_variable) # Ensure input exists
    ref_variable <- input$reference_variable

    # Extract unique levels for the selected variable
    ref_levels <- unique(as.character(phyloseq::sample_data(current_physeq())[[ref_variable]]))

    updateSelectInput(
      session,
      "reference_level",
      choices = ref_levels,
      selected = NULL
    )
  })

  # Combine selected reference variable and level into the reference vector for Maaslin2
  observeEvent(input$reference_level, {
    req(input$reference_variable, input$reference_level) # Ensure both inputs exist
    combined_reference <- paste(input$reference_variable, input$reference_level, sep = ",")

    # Update the text input for Maaslin2 reference
    updateTextInput(session, "reference", value = combined_reference)
  })

  ### Maaslin2 Execution Logic ###
  observeEvent(input$run_maaslin2, {
    req(current_physeq(), input$fixed_variables, input$Maaslin2_output)

    # Validate or create output folder
    selected_folder <- input$Maaslin2_output
    if (is.null(selected_folder) || selected_folder == "") {
      showModal(modalDialog(
        title = "Error",
        "Please specify a valid output directory for Maaslin2 results.",
        easyClose = TRUE,
        footer = modalButton("Close")
      ))
      return()
    }

    if (!dir.exists(selected_folder)) {
      tryCatch({
        dir.create(selected_folder, recursive = TRUE)
      }, error = function(e) {
        showModal(modalDialog(
          title = "Error",
          paste("Failed to create the output directory:", e$message),
          easyClose = TRUE,
          footer = modalButton("Close")
        ))
        return()
      })
    }

    # Aggregate physeq data at the selected taxonomy level
    physeq_agg <- microbiome::aggregate_taxa(current_physeq(), input$taxonomy_level)


    # Validate inputs
    Count_table <- data.frame(t(otu_table(physeq_agg))) # Transpose OTU table
    metadata <- data.frame(phyloseq::sample_data(physeq_agg)) # Extract metadata

    if (ncol(Count_table) == 0 || nrow(Count_table) == 0) {
      showModal(modalDialog(
        title = "Error",
        "The OTU table is empty or invalid. Please check your input data.",
        easyClose = TRUE,
        footer = modalButton("Close")
      ))
      return()
    }

    if (nrow(metadata) == 0) {
      showModal(modalDialog(
        title = "Error",
        "The sample metadata is empty or invalid. Please check your input data.",
        easyClose = TRUE,
        footer = modalButton("Close")
      ))
      return()
    }

    tryCatch({
      # Construct reference vector
      req(input$reference_variable, input$reference_level)
      reference_vector <- c(input$reference_variable, input$reference_level)

      if (length(reference_vector) != 2 || any(reference_vector == "")) {
        showModal(modalDialog(
          title = "Error",
          "Please specify a valid Reference Variable and Level.",
          easyClose = TRUE,
          footer = modalButton("Close")
        ))
        return()
      }

      # Debug: Print inputs
      print("Maaslin2 Inputs:")
      print(paste("Fixed Effects:", paste(input$fixed_variables, collapse = ", ")))
      print(paste("Random Effects:", paste(input$random_variables, collapse = ", ")))
      print(paste("Reference Vector:", paste(reference_vector, collapse = ", ")))
      print(paste("Output Folder:", selected_folder))

      # Run Maaslin2
      maaslin_da <- Maaslin2::Maaslin2(
        input_data = Count_table,
        input_metadata = metadata,
        normalization = input$normalization,
        standardize = input$standardize,
        transform = input$transform,
        analysis_method = input$analysis_method,
        max_significance = input$max_significance,
        output = selected_folder,
        fixed_effects = input$fixed_variables,
        random_effects = input$random_variables,
        reference = reference_vector,  # Pass formatted reference vector
        min_abundance = input$min_abundance,
        min_prevalence = input$min_prevalence,
        plot_heatmap = input$plot_heatmap,
        heatmap_first_n = input$heatmap_first_n,
        plot_scatter = input$plot_scatter,
        correction = input$correction,
        cores = input$cores
      )

      # Display significant results
      significant_results <- maaslin_da$results %>%
        dplyr::filter(qval < input$max_significance) %>%
        dplyr::arrange(qval)  # Sort by qval in ascending order

      if (nrow(significant_results) == 0) {
        showModal(modalDialog(
          title = "No Significant Results",
          "No taxa were found to be significantly differentially abundant with the given parameters.",
          easyClose = TRUE,
          footer = modalButton("Close")
        ))
      } else {
        output$maaslin2_results <- DT::renderDT({
          DT::datatable(significant_results, options = list(pageLength = 10, scrollX = TRUE))
        })
      }

    }, error = function(e) {
      # Handle errors in Maaslin2 execution
      showModal(modalDialog(
        title = "Error Running Maaslin2",
        paste("An error occurred while running Maaslin2:", e$message),
        easyClose = TRUE,
        footer = modalButton("Close")
      ))
    })
  })


  # Server-side code to show the explanation for Maaslin2
  observeEvent(input$show_maaslin2_explanation, {
    showModal(modalDialog(
      title = "Maaslin2 Explanation",
      HTML("<p><strong>Maaslin2</strong> is a tool used for determining associations between microbial features (e.g., taxa) and metadata (e.g., sample groupings, experimental conditions) using multivariate statistical models.</p>
         <p><strong>Key Concepts:</strong></p>
         <ul>
           <li><strong>Fixed Effects:</strong> Variables of interest (e.g., experimental group, treatment, or other metadata) that are tested for associations with microbial features.</li>
           <li><strong>Random Effects:</strong> Variables that introduce random variability (e.g., batch effects, repeated measures, or nested designs) and are included to account for this variability.</li>
           <li><strong>Normalization Methods:</strong> Transformations such as Total Sum Scaling (TSS), Centered Log Ratio (CLR), or Cumulative Sum Scaling (CSS) to standardize feature abundances.</li>
           <li><strong>Transformation Methods:</strong> Transformations such as log or logit are applied to the data to improve model performance.</li>
           <li><strong>Effect Size (coef):</strong> The estimated magnitude and direction of the association between the microbial feature and the metadata variable. Positive values indicate enrichment, while negative values indicate depletion.</li>
           <li><strong>Significance (pval and qval):</strong> The p-value indicates the statistical significance of the association, while the q-value accounts for multiple testing corrections. Significant results have q-values below the specified threshold (e.g., 0.05).</li>
           <li><strong>Prevalence and Abundance Thresholds:</strong> Filters to exclude features that do not meet minimum abundance or prevalence criteria, ensuring robust results.</li>
         </ul>
         <p><strong>How to Interpret the Results:</strong></p>
         <ul>
           <li><strong>Significant Features:</strong> Focus on features with q-values below the threshold. These are considered significantly associated with the metadata.</li>
           <li><strong>Effect Sizes:</strong> Positive coefficients indicate enrichment in the group of interest, while negative coefficients indicate depletion.</li>
           <li><strong>Prevalence (N.not.zero):</strong> Indicates the number of samples in which the feature is present. Features with very low prevalence may require further validation.</li>
         </ul>
         <p>Maaslin2 helps identify meaningful relationships between microbial communities and experimental conditions, providing insights into how microbiota composition is influenced by various factors.</p>"),
      easyClose = TRUE,
      footer = modalButton("Close")
    ))
  })





  ############# community typing(DMM)




  generateDriverPlot <- function(best, physeq, k) {
    d <- reshape2::melt(fitted(best))
    colnames(d) <- c("ASV", "cluster", "value")

    # Get the tax_table from the phyloseq object
    tax_df <- as.data.frame(phyloseq::tax_table(physeq))

    # Create a new column combining Genus, Species, and ASV
    d$Taxa <- with(tax_df[d$ASV, ], paste0(Genus, "_", Species, " (", d$ASV, ")"))

    # Replace NA values with "Unknown"
    d$Taxa[is.na(d$Taxa)] <- paste0("Unknown (", d$ASV[is.na(d$Taxa)], ")")

    # Filter and prepare the data for plotting
    d <- d %>% filter(cluster == k) %>%
      arrange(value) %>%
      mutate(Taxa = factor(Taxa, levels = unique(Taxa))) %>%
      filter(abs(value) > quantile(abs(value), 0.8))

    # Create the plot
    ggplot2::ggplot(d, aes(x = Taxa, y = value)) +
      ggplot2::theme_bw() +
      geom_bar(stat = "identity", fill = "darkblue") +
      coord_flip() +
      labs(title = paste("Top drivers: community type", k))
  }


  observe({
    req(current_physeq())
    pseq <- current_physeq()
    updateSelectInput(session, "taxRankDMM", choices = c(rank_names(pseq), "Genus", "Species"))
  })

  observeEvent(input$runDMM, {
    req(current_physeq(), input$taxRankDMM, input$detectionThreshold, input$prevalenceThreshold)

    # Data Transformation and Filtering
    pseq.comp <- microbiome::transform(current_physeq(), "compositional")
    taxa <- microbiome::core_members(pseq.comp, detection = input$detectionThreshold / 100, prevalence = input$prevalenceThreshold / 100)
    pseq <- phyloseq::prune_taxa(taxa, current_physeq())

    # Aggregate at selected taxonomic rank
    if (input$taxRankDMM %in% c("Genus", "Species")) {
      pseq <- phyloseq::tax_glom(pseq, taxrank = input$taxRankDMM)
    } else {
      pseq <- phyloseq::tax_glom(pseq, taxrank = input$taxRankDMM)
    }

    dat <- microbiome::abundances(pseq)
    count <- as.matrix(t(dat))
    count <- count[rowSums(count) > 0, ]

    if (nrow(count) > 0) {
      # Model Fitting
      fit <- lapply(1:input$numComponents, DirichletMultinomial::dmn, count = count, verbose = TRUE)
      lplc <- sapply(fit, DirichletMultinomial::laplace)
      aic <- sapply(fit, DirichletMultinomial::AIC)
      bic <- sapply(fit, DirichletMultinomial::BIC)

      # Render Model Fit Plot
      output$modelFitPlot <- renderPlot({
        plot(lplc, type = "b", xlab = "Number of Dirichlet Components", ylab = "Model Fit")
        lines(aic, type = "b", lty = 2)
        lines(bic, type = "b", lty = 3)
        legend("topright", legend = c("Laplace", "AIC", "BIC"), lty = 1:3, bty = "n")
      })

      best <- fit[[which.min(unlist(lplc))]]

      # Render Mixture Parameters
      output$mixtureParams <- renderPrint({
        DirichletMultinomial::mixturewt(best)
      })

      # Render Sample Assignments
      sample_assignments <- apply(DirichletMultinomial::mixture(best), 1, which.max)
      output$sampleAssignments <- renderPrint({
        sample_assignments
      })

      # Filter metadata to only include samples in the count matrix
      valid_samples <- rownames(count)
      metadata <- as(phyloseq::sample_data(current_physeq()), "data.frame")
      metadata_filtered <- metadata[valid_samples, ]

      # Update metadata with sample assignments
      metadata_filtered$DMM_Cluster <- factor(sample_assignments)

      # Re-create phyloseq object with filtered metadata
      updated_pseq <- phyloseq(otu_table(current_physeq()), phyloseq::tax_table(current_physeq()), phyloseq::sample_data(metadata_filtered), phy_tree(current_physeq()))

      # Render Driver Plots
      output$driverPlotsUI <- renderUI({
        plots <- lapply(1:ncol(fitted(best)), function(k) {
          plotname <- paste("driverPlot", k, sep = "")
          plotOutput(plotname)
        })
        do.call(tagList, plots)
      })

      # Using the function for rendering plots
      lapply(1:ncol(fitted(best)), function(k) {
        output[[paste("driverPlot", k, sep = "")]] <- renderPlot({
          generateDriverPlot(best, current_physeq(), k)
        })
      })

      # Download Handlers for Updated Phyloseq Object and Metadata
      output$downloadUpdatedPhyseq <- downloadHandler(
        filename = function() {
          paste0("updated_phyloseq_", Sys.Date(), ".rds")
        },
        content = function(file) {
          saveRDS(updated_pseq, file)
        }
      )

      output$downloadMetadata <- downloadHandler(
        filename = function() {
          paste0("metadata_", Sys.Date(), ".csv")
        },
        content = function(file) {
          write.csv(metadata_filtered, file, row.names = TRUE)
        }
      )

      # Download Handlers for Other Results
      output$downloadModelFitPlot <- downloadHandler(
        filename = function() {
          paste0("model_fit_", Sys.Date(), ".png")
        },
        content = function(file) {
          png(file)
          plot(lplc, type = "b", xlab = "Number of Dirichlet Components", ylab = "Model Fit")
          lines(aic, type = "b", lty = 2)
          lines(bic, type = "b", lty = 3)
          legend("topright", legend = c("Laplace", "AIC", "BIC"), lty = 1:3, bty = "n")
          dev.off()
        }
      )

      output$downloadMixtureParams <- downloadHandler(
        filename = function() {
          paste0("mixture_params_", Sys.Date(), ".txt")
        },
        content = function(file) {
          write(mixturewt(best), file)
        }
      )

      output$downloadSampleAssignments <- downloadHandler(
        filename = function() {
          paste0("sample_assignments_", Sys.Date(), ".csv")
        },
        content = function(file) {
          sample_assignments_df <- sample_assignments_reactive()  # Access the stored data
          write.csv(sample_assignments_df, file, row.names = FALSE)
        }
      )

      # Using the function for downloading driver plots
      output$downloadDriverPlots <- downloadHandler(
        filename = function() {
          paste0("driver_plots_", Sys.Date(), ".pdf")
        },
        content = function(file) {
          pdf(file)
          lapply(1:ncol(fitted(best)), function(k) {
            print(generateDriverPlot(best, current_physeq(), k))
          })
          dev.off()
        }
      )

    } else {
      # Handle cases where no valid samples are found
      output$modelFitPlot <- renderUI({
        tags$p("No valid samples with non-zero counts.", style = "color: red; text-align: center; font-size: 18px;")})
      output$mixtureParams <- renderPrint({"No valid samples with non-zero counts."})
      output$sampleAssignments <- renderDT({
        DT::datatable(data.frame(Sample = character(0), Cluster = character(0)), options = list(pageLength = 10, autoWidth = TRUE))
      })
      output$driverPlotsUI <- renderUI({ NULL })
    }
  })


  # Server-side code to show the explanation for Community Typing (DMM)
  observeEvent(input$show_community_typing_explanation, {
    showModal(modalDialog(
      title = "Community Typing (Dirichlet Multinomial Mixture) Explanation",
      HTML("<p><strong>Community Typing using Dirichlet Multinomial Mixture (DMM)</strong> is an approach to categorize or 'type' microbial communities into distinct groups based on their composition.</p>
         <p><strong>Purpose:</strong></p>
         <p>This analysis aims to identify and characterize different 'community types' in your dataset. By grouping samples based on the similarity of their microbial compositions, we can gain insights into the underlying ecological processes or environmental factors that influence these communities.</p>
         <p><strong>Algorithm and Method:</strong></p>
         <ul>
           <li><strong>Dirichlet Multinomial Mixture (DMM):</strong> The DMM model assumes that each sample is drawn from one of several possible 'community types', each represented by a distinct distribution of microbial taxa. The model estimates these distributions and assigns each sample to a community type based on the likelihood of its observed composition.</li>
           <li><strong>Steps in the Analysis:</strong>
             <ul>
               <li><em>Data Preprocessing:</em> The data is first normalized using compositional transformations to standardize the microbial abundances across samples.</li>
               <li><em>Model Fitting:</em> The DMM model is fitted to the data, with the number of community types (components) specified by the user. The model's fit is evaluated using criteria like Laplace, AIC, and BIC to determine the optimal number of community types.</li>
               <li><em>Community Assignment:</em> Each sample is assigned to a community type based on the highest posterior probability.</li>
               <li><em>Driver Analysis:</em> Key microbial taxa ('drivers') that differentiate between community types are identified and visualized.</li>
             </ul>
           </li>
           <li><strong>Applications:</strong> Understanding community types can help in identifying patterns related to health, disease, environmental conditions, or treatment responses. It provides a way to summarize complex microbial data into interpretable groups.</li>
         </ul>
         <p>This method is particularly useful in studies where microbial communities are expected to be heterogeneous, such as in gut microbiome research, environmental microbiology, or in comparing healthy vs. diseased states.</p>"),
      easyClose = TRUE,
      footer = modalButton("Close")
    ))
  })



  ############ RDA and dbRDA plot

  # Initialize reactive values to store plots and data
  plotly_plots_rda <- reactiveVal(NULL)
  ggplot_plots_rda <- reactiveVal(NULL)
  rounded_permanova <- reactiveVal(NULL)
  species_scores <- reactiveVal(NULL)

  # Function to impute missing values with the median
  impute_median <- function(pseq) {
    metadata <- phyloseq::sample_data(pseq)
    imputed_metadata <- metadata
    for (col in colnames(metadata)) {
      if (is.numeric(metadata[[col]])) {
        imputed_metadata[[col]][is.na(imputed_metadata[[col]])] <- median(imputed_metadata[[col]], na.rm = TRUE)
      }
    }
    phyloseq::sample_data(pseq) <- imputed_metadata
    return(pseq)
  }

  # Create a reactive value to store the pseq object
  reactive_pseq <- reactiveVal()

  # Impute missing values and update select inputs
  observe({
    req(current_physeq())
    pseq <- current_physeq()

    # Impute missing values before analysis
    pseq <- impute_median(pseq)

    # Store the imputed pseq in the reactive value
    reactive_pseq(pseq)

    # Update taxRankRDA choices based on the available taxonomic ranks
    if (!is.null(phyloseq::tax_table(reactive_pseq()))) {
      updateSelectInput(session, "taxRankRDA", choices = c("ASV", rank_names(reactive_pseq())))
    } else {
      updateSelectInput(session, "taxRankRDA", choices = "ASV")
      showNotification("Tax table is missing in the phyloseq object, using ASV level.", type = "warning")
    }

    # Update other inputs
    sample_data_cols <- colnames(phyloseq::sample_data(pseq))
    updateSelectInput(session, "constraints", choices = sample_data_cols)
    updateSelectInput(session, "conditions", choices = c("None", sample_data_cols))
    updateSelectInput(session, "colorBy", choices = sample_data_cols)
    updateSelectInput(session, "shapeBy", choices = sample_data_cols)
  })


  observeEvent(input$plotRDA, {
    req(reactive_pseq(), input$RDAdistanceMetric, input$taxRankRDA, input$normMethod, input$colorBy, input$shapeBy, input$constraints, input$conditions)

    tryCatch({
      # Retrieve the pseq from the reactive object
      pseq <- reactive_pseq()

      # Normalize and aggregate if needed
      if (input$taxRankRDA != "ASV") {
        pseq <- aggregate_taxa(pseq, input$taxRankRDA)
      }
      if (input$normMethod != "identity") {
        pseq <- microbiome::transform(pseq, input$normMethod)
      }

      # Construct the formula
      formula <- if (length(input$constraints) == 1) {
        as.formula(paste("otu_table ~", input$constraints))
      } else {
        as.formula(paste("otu_table ~", paste(input$constraints, collapse = " + ")))
      }
      if (!is.null(input$conditions) && input$conditions != "None") {
        formula <- update(formula, paste0(". ~ . + Condition(", paste(input$conditions, collapse = " + "), ")"))
      }

      # Extract OTU table and metadata
      otu_table <- as.matrix(phyloseq::otu_table(pseq))
      if (phyloseq::taxa_are_rows(pseq)) {
        otu_table <- t(otu_table)
      }
      metadata <- data.frame(phyloseq::sample_data(pseq))

      # Perform capscale analysis
      CAP_distance <- vegan::capscale(formula = formula, data = metadata, distance = input$RDAdistanceMetric, sqrt.dist = TRUE)



      # Extract species and site scores
      species_scores <- vegan::scores(CAP_distance, display = "species")
      site_scores <- vegan::scores(CAP_distance, display = "sites")
      env_vectors <- vegan::scores(CAP_distance, display = "bp")

      # Convert scores to data frames
      species_df <- as.data.frame(species_scores)
      species_df$Taxa <- rownames(species_df)

      site_df <- as.data.frame(site_scores)
      site_df$Group <- metadata[[input$colorBy]]

      env_df <- as.data.frame(env_vectors)
      env_df$Constraint <- rownames(env_df)

      # Define colors and shapes dynamically
      group_levels <- unique(site_df$Group)
      group_colors <- scales::hue_pal()(length(group_levels))
      group_shapes <- seq(1, length(group_levels))

      # Generate the plot
      rda_plot <- ggplot() +
        ggplot2::geom_point(data = species_df, aes(x = CAP1, y = CAP2), size = 2, shape = 17, alpha = 0.7, color = "red") +
        geom_text(data = species_df, aes(x = CAP1, y = CAP2, label = Taxa), size = 3, hjust = 1.1, vjust = 1.1, color = "darkred") +
        ggplot2::geom_point(data = site_df, aes(x = CAP1, y = CAP2, color = Group, shape = Group), size = 4, alpha = 0.7) +
        geom_segment(data = env_df, aes(x = 0, y = 0, xend = CAP1, yend = CAP2), arrow = arrow(length = unit(0.3, "cm")), color = "blue") +
        ggplot2::scale_color_manual(values = setNames(group_colors, group_levels)) +
        ggplot2::scale_shape_manual(values = setNames(group_shapes, group_levels)) +
        labs(title = "CAP (dbRDA) Plot", x = "CAP1", y = "CAP2", color = "Group", shape = "Group") +
        theme_minimal()

      # Render the ggplot
      ggplot_plots_rda(rda_plot)
      output$rdaPlot <- renderPlot({ ggplot_plots_rda() })

      # Enable plot download with dynamic file type
      output$download_rdaPlot <- downloadHandler(
        filename = function() {
          paste0("CAP_RDA_Plot_", Sys.Date(), ".", input$rda_filetype)
        },
        content = function(file) {
          ggsave(
            file,
            plot = ggplot_plots_rda(),
            width = 10,
            height = 8,
            units = "in",
            device = input$rda_filetype
          )
        }
      )


      ############## Extract statistics

      # Extract summary information
      cap_summary <- summary(CAP_distance)
      inertia_table <- as.data.frame(cap_summary$cont$importance)
      inertia_table <- cbind(
        Category = rownames(inertia_table),
        inertia_table
      )
      rownames(inertia_table) <- NULL  # Remove row names for cleaner output



      eigenvalues <- as.data.frame(vegan::eigenvals(CAP_distance))
      adjusted_r2 <- vegan::RsquareAdj(CAP_distance)
      regression_coefficients <- as.data.frame(coef(CAP_distance))
      site_scores <- as.data.frame(scores(CAP_distance, display = "sites"))
      species_scores <- as.data.frame(scores(CAP_distance, display = "species"))
      env_vectors <- scores(CAP_distance, display = "bp")
      anova_overall <- as.data.frame(anova(CAP_distance))
      anova_terms <- as.data.frame(anova(CAP_distance, by = "terms", permu = 200))
      anova_axes <- as.data.frame(anova(CAP_distance, by = "axis", perm.max = 500))

      # Render summary inertia table
      output$inertiaTable <- DT::renderDT({
        DT::datatable(inertia_table, options = list(pageLength = 5), rownames = FALSE) })

      # Render other tables
      output$eigenTable <- DT::renderDT({ DT::datatable(eigenvalues) })
      output$r2Table <- DT::renderDT({ DT::datatable(data.frame(Adjusted_R2 = adjusted_r2)) })
      output$coeffTable <- DT::renderDT({ DT::datatable(regression_coefficients) })
      output$siteScoresTable <- DT::renderDT({ DT::datatable(site_scores) })
      output$speciesScoresTable <- DT::renderDT({ DT::datatable(species_scores) })
      output$anovaOverallTable <- DT::renderDT({ DT::datatable(anova_overall) })
      output$anovaTermsTable <- DT::renderDT({ DT::datatable(anova_terms) })
      output$anovaAxesTable <- DT::renderDT({ DT::datatable(anova_axes) })

      # Download Handlers

      # Enable downloading of the inertia table
      output$downloadInertia <- downloadHandler(
        filename = function() { paste0("cap_inertia_summary_", Sys.Date(), ".csv") },
        content = function(file) { write.csv(inertia_table, file, row.names = FALSE) }
      )
      output$downloadEigen <- downloadHandler(
        filename = function() { "eigenvalues.csv" },
        content = function(file) { write.csv(eigenvalues, file) }
      )
      output$downloadR2 <- downloadHandler(
        filename = function() { "adjusted_r2.csv" },
        content = function(file) { write.csv(data.frame(Adjusted_R2 = adjusted_r2), file) }
      )
      output$downloadCoeff <- downloadHandler(
        filename = function() { "regression_coefficients.csv" },
        content = function(file) { write.csv(regression_coefficients, file) }
      )
      output$downloadSiteScores <- downloadHandler(
        filename = function() { "site_scores.csv" },
        content = function(file) { write.csv(site_scores, file) }
      )
      output$downloadSpeciesScores <- downloadHandler(
        filename = function() { "species_scores.csv" },
        content = function(file) { write.csv(species_scores, file) }
      )
      output$downloadAnovaOverall <- downloadHandler(
        filename = function() { "anova_overall.csv" },
        content = function(file) { write.csv(anova_overall, file) }
      )
      output$downloadAnovaTerms <- downloadHandler(
        filename = function() { "anova_terms.csv" },
        content = function(file) { write.csv(anova_terms, file) }
      )
      output$downloadAnovaAxes <- downloadHandler(
        filename = function() { "anova_axes.csv" },
        content = function(file) { write.csv(anova_axes, file) }
      )

    }, error = function(e) {
      showModal(modalDialog(title = "Error in CAP Analysis", e$message, easyClose = TRUE, footer = modalButton("Close")))
    })
  })

  observeEvent(input$show_rda_explanation, {
    showModal(modalDialog(
      title = "CAP and RDA Explanation",
      HTML("<p><strong>Constrained Analysis of Principal Coordinates (CAP)</strong> is an ordination method similar to Redundancy Analysis (RDA), but it allows the use of non-Euclidean distance metrics, such as Bray-Curtis or Jaccard distances. CAP is particularly useful for ecological and microbiome studies where non-Euclidean distances better represent community dissimilarities.</p>
         <p>When Euclidean distances are used, CAP produces results identical to RDA, but it is less computationally efficient for such cases. CAP is a constrained version of metric scaling (Principal Coordinates Analysis) that can also handle unconstrained ordination with extended dissimilarities.</p>

         <p><strong>Why Use CAP in This App?</strong></p>
         <ul>
           <li><strong>Flexibility in Dissimilarity Measures:</strong> CAP allows non-Euclidean distance metrics, which are common in microbiome analyses.</li>
           <li><strong>Constrained Analysis:</strong> CAP identifies how specific metadata variables (constraints) explain variations in the microbiome data.</li>
           <li><strong>Visual and Statistical Insights:</strong> CAP provides plots showing relationships among samples, taxa, and metadata, along with statistical summaries to evaluate the significance of constraints.</li>
         </ul>

         <p>In this app, CAP allows users to analyze microbiome data constrained by selected metadata variables (e.g., experimental conditions or environmental factors). Users can:</p>
         <ul>
           <li>Select taxonomic ranks (e.g., Genus or ASV) and normalize the data using methods like Hellinger or CLR transformations.</li>
           <li>Choose distance metrics, such as Bray-Curtis or Jaccard, for the analysis.</li>
           <li>Visualize relationships between samples, taxa, and metadata in an interactive plot.</li>
           <li>Access detailed statistical results, including inertia breakdowns, eigenvalues, regression coefficients, and ANOVA tests.</li>
         </ul>
         <p>This provides both visual and quantitative insights into how metadata explains variation in the microbial community.</p>"),
      easyClose = TRUE,
      footer = modalButton("Close")
    ))
  })
}
