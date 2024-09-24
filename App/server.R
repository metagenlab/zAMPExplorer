#server

shinyServer(function(input, output, session) {


  # Load and preprocess the phyloseq
  physeq <- reactive({
    req(input$physeqFile)
    pseq <- readRDS(input$physeqFile$datapath)

    # Preprocess the phyloseq object: remove samples or taxa with NA values and zero counts
    pseq <- prune_samples(sample_sums(pseq) > 0, pseq)
    pseq <- prune_taxa(taxa_sums(pseq) > 0, pseq)

    # Drop existing diversity metrics from sample data
    drop <- c("Observed", "Chao1", "se.chao1", "ACE", "se.ACE", "Shannon", "Simpson", "InvSimpson", "Fisher", "Observed_min_1")
    sample_data(pseq) <- sample_data(pseq)[, !(names(sample_data(pseq)) %in% drop)]

    return(pseq)
  })

  # make sure the UI is updated after the data is loaded
  observe({
    pseq <- physeq()
    req(pseq)


    # Dynamically get all possible diversity metrics from microbiome::alpha()
    diversity_metrics <- colnames(microbiome::alpha(physeq(), index = "all"))

    # Update bottomns based on user choice
    updateSelectInput(session, "phylum", choices = unique(tax_table(pseq)[, "Phylum"]))
    updateSelectInput(session, "color_by", choices = sample_variables(pseq))
    updateSelectInput(session, "taxaLevel", choices = rank_names(pseq))
    updateSelectInput(session, "sampleGroup", choices = sample_variables(pseq))
    updateSelectInput(session, "hueRank", choices = rank_names(pseq))
    updateSelectInput(session, "shadeRank", choices = rank_names(pseq))
    updateSelectInput(session, "facetBy", choices = c("ASV", sample_variables(pseq)))
    updateSelectInput(session, "heatmapRank", choices = rank_names(pseq))
    updateSelectInput(session, "annotationColumn1", choices = sample_variables(pseq))
    updateSelectInput(session, "annotationColumn2", choices = sample_variables(pseq))
    # updateSelectInput(session, "vennTaxonomicRank", choices = c("ASV", rank_names(physeq())))
    # updateSelectInput(session, "vennGroupingColumn", choices = sample_variables(physeq()))
    updateSelectInput(session, "alphaMetric", choices = diversity_metrics)
    updateSelectInput(session, "alphaGroupingColumn", choices = sample_variables(pseq))
    updateSelectInput(session, "groupingColumn", choices = sample_variables(pseq))
    updateSelectInput(session, "groupingColumnBeta", choices = sample_variables(pseq))
    updateSelectInput(session, "shapeColumnBeta", choices = sample_variables(pseq))
    updateSelectInput(session, "taxRank", choices = rank_names(pseq))
    updateSelectInput(session, "dmmTaxonomicRank", choices = c("ASV", rank_names(physeq())))
    updateSelectInput(session, "rdaTaxonomicRank", choices = c("ASV", rank_names(physeq())))
    updateSelectInput(session, "shapeColumn", choices = sample_variables(pseq))
    updateSelectInput(session, "explanatoryvariables", choices = sample_variables(pseq))
    updateSelectInput(session, "colorBy", choices = sample_variables(pseq))
  })

  # Populate the selectInput choices for phylogenetic tree
  observe({
    req(physeq())

    # Check if the Phylum rank exists in the taxonomy table
    if ("Phylum" %in% colnames(tax_table(physeq()))) {
      available_phyla <- unique(tax_table(physeq())[,"Phylum"])
      if (!is.null(available_phyla) && length(available_phyla) > 0) {
        updateSelectInput(session, "phylum", choices = available_phyla)
      } else {
        updateSelectInput(session, "phylum", choices = NULL)
        showNotification("Phylum rank is missing or empty in the physeq object", type = "warning")
      }
    } else {
      updateSelectInput(session, "phylum", choices = NULL)
      showNotification("Phylum rank is missing in the physeq object", type = "warning")
    }

    # Populate other dropdowns
    available_sample_columns <- colnames(sample_data(physeq()))
    updateSelectInput(session, "color_by", choices = available_sample_columns)
    available_tax_ranks <- colnames(tax_table(physeq()))
    updateSelectInput(session, "shape_by", choices = available_tax_ranks)
    updateSelectInput(session, "label_tips", choices = available_tax_ranks)
  })

  # Render Phylogenetic Tree when button is clicked
  observeEvent(input$plot_phylogenetic_tree, {
    req(input$phylum, input$color_by, input$shape_by, input$label_tips, physeq())

    # Check if the selected phylum is valid
    if (input$phylum %in% unique(tax_table(physeq())[,"Phylum"])) {
      pseq_filtered <- subset_taxa(physeq(), Phylum == input$phylum)

      output$phylogenetic_tree_plot <- renderPlot({
        plot_tree(
          pseq_filtered,
          color = input$color_by,
          shape = input$shape_by,
          label.tips = input$label_tips,
          size = "abundance",
          plot.margin = input$plot_size
        )
      })
    } else {
      showNotification("Selected phylum is not valid.", type = "error")
    }
  })




                                    ######### phyloseq components

  # Reactive value to track which button is pressed
  output_view <- reactiveVal("")

  # Observers for button clicks
  observeEvent(input$show_metadata, {
    output_view("metadata")
  })

  observeEvent(input$show_count_table, {
    output_view("count/Abundance_table")
  })

  observeEvent(input$show_combined_table, {
    output_view("combined_table")
  })

  observeEvent(input$show_phylogenetic_tree, {
    output_view("phylogenetic_tree")
  })

  observeEvent(input$show_summary_statistics, {
    output_view("summary_statistics")
  })





  # Show the selected data section
  output$dynamic_tables <- renderUI({
    if (output_view() == "metadata") {
      fluidRow(
        box(title = "Metadata Structure", width = 12, status = "primary",
            DT::dataTableOutput("metadata_structure")
        )
      )
    } else if (output_view() == "combined_table") {
      fluidRow(
        box(title = "Abundance/Taxonomy Table", width = 12, status = "primary",
            # Display combined table
            DT::dataTableOutput("combined_table"),
            # File format selection for download (under the table)
            selectInput("filetype", "Choose file type:", choices = c("csv", "xlsx", "tsv")),
            # Download button for the combined table
            downloadButton("download_combined_table", "Download Abundance/Taxonomy Table")
        )
      )
    } else if (output_view() == "phylogenetic_tree") {
      fluidRow(
        box(title = "Phylogenetic Tree", width = 12, status = "primary",
            # Add dynamic options for customization
            selectInput("phylum", "Choose Phylum:", choices = unique(tax_table(physeq())[, "Phylum"])),
            selectInput("color_by", "Choose Color By (metadata):", sample_variables(physeq())),
            selectInput("shape_by", "Choose Shape By (taxonomy):", choices = rank_names(physeq())),
            selectInput("label_tips", "Label Tips By (taxonomy):", choices = rank_names(physeq())),
            actionButton("plot_tree", "Plot Phylogenetic Tree"),
            div(style = 'overflow-x: auto; overflow-y: auto; width: 100%; height: 800px;',
                plotOutput("phylogenetic_tree_plot")
            ),
            # Add download button for the tree
            downloadButton("download_ggtree", "Download Phylogenetic Tree")
        )
      )
    } else if (output_view() == "summary_statistics") {
      fluidRow(
        box(title = "Summary Statistics", width = 12, status = "primary",
            DT::dataTableOutput("summary_statistics")
        )
      )
    }
  })



  #Render Metadata Structure when button is clicked
    output$metadata_structure <- DT::renderDataTable({
      req(output_view() == "metadata")
      as.data.frame(sample_data(physeq()))
    }, options = list(pageLength = 10, scrollX = TRUE))

    # Store the combined table in a reactive expression to avoid duplication
    combined_table <- reactive({
      # Extract the OTU/ASV count table
      otu_table <- as.data.frame(otu_table(physeq()))

      # Extract the taxonomy table
      taxonomy_table <- as.data.frame(tax_table(physeq()))

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
    output$combined_table <- DT::renderDataTable({
      req(output_view() == "combined_table")
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


    # Render Phylogenetic Tree using ggtree
    output$phylogenetic_tree_plot <- renderPlot({
      req(input$plot_tree)  # Ensure the user has clicked the plot button
      req(input$phylum, input$color_by, input$shape_by, input$label_tips)  # Ensure input exists


    # Subset the tree by selected Phylum
    GP <- prune_taxa(taxa_sums(physeq()) > 0, physeq())

    # Subset by phylum and check if the phylum has any taxa
    tree_subset <- subset_taxa(GP, Phylum == input$phylum)

    # Check if the phylum subset is empty
    if (is.null(tree_subset) || ntaxa(tree_subset) == 0) {
      showNotification("No taxa found for the selected phylum", type = "error")
      return(NULL)
    }

    # Check if color_by, shape_by, and label_tips are valid
    if (is.null(input$color_by) || !input$color_by %in% names(sample_data(physeq()))) {
      showNotification("Selected column for color is not found in the metadata", type = "error")
      return(NULL)
    }

    if (is.null(input$shape_by) || !input$shape_by %in% names(tax_table(physeq()))) {
      showNotification("Selected column for shape is not found in the taxonomy table", type = "error")
      return(NULL)
    }

    if (is.null(input$label_tips) || !input$label_tips %in% names(tax_table(physeq()))) {
      showNotification("Selected column for labeling tips is not found in the taxonomy table", type = "error")
      return(NULL)
    }

    # Generate ggtree plot
    tree_plot <- ggtree(tree_subset, ladderize = FALSE) +
      geom_text2(aes(subset = !isTip, label = label), hjust = -0.2, size = 4) +
      geom_tiplab(aes(label = get(input$label_tips)), hjust = -0.3) +
      geom_point(aes(x = x + hjust, color = get(input$color_by), shape = get(input$shape_by)), na.rm = TRUE) +
      scale_size_continuous(trans = log_trans(5)) +
      theme(legend.position = "right") +
      ggtitle(paste("Phylogenetic Tree (Phylum: ", input$phylum, ")", sep = ""))

    print(tree_plot)
  })

  # Download ggtree plot
  output$download_ggtree <- downloadHandler(
    filename = function() {
      paste0("phylogenetic_tree_", Sys.Date(), ".pdf")
    },
    content = function(file) {
      # Render the ggtree plot for download
      GP <- prune_taxa(taxa_sums(physeq()) > 0, physeq())
      tree_subset <- subset_taxa(GP, Phylum == input$phylum)

      tree_plot <- ggtree(tree_subset, ladderize = FALSE) +
        geom_text2(aes(subset = !isTip, label = label), hjust = -0.2, size = 4) +
        geom_tiplab(aes(label = get(input$label_tips)), hjust = -0.3) +
        geom_point(aes(x = x + hjust, color = get(input$color_by), shape = get(input$shape_by),
                       size = Abundance), na.rm = TRUE) +
        scale_size_continuous(trans = log_trans(5)) +
        theme(legend.position = "right") +
        ggtitle(paste("Phylogenetic Tree (Phylum: ", input$phylum, ")", sep = ""))

      # Save to file
      ggsave(file, tree_plot, width = 10, height = 8)
    }
  )


  # Render Summary Statistics
  output$summary_statistics <- DT::renderDataTable({
    req(output_view() == "summary_statistics")
    num_samples <- nsamples(physeq())
    num_asvs <- ntaxa(physeq())
    num_phylum <- length(unique(tax_table(physeq())[,"Phylum"]))
    num_family <- length(unique(tax_table(physeq())[,"Family"]))
    num_genus <- length(unique(tax_table(physeq())[,"Genus"]))
    num_species <- length(unique(tax_table(physeq())[,"Species"]))

    summary_df <- data.frame(
      "Metric" = c("Number of Samples", "Number of ASVs", "Number of Phyla", "Number of Families", "Number of Genera", "Number of Species"),
      "Count" = c(num_samples, num_asvs, num_phylum, num_family, num_genus, num_species)
    )

    DT::datatable(summary_df, options = list(pageLength = 10, scrollX = TRUE), rownames = FALSE)
  })


#
#                               ######### phyloseq components
#
#   # Reactive value to track which button is pressed
#   output_view <- reactiveVal("")
#
#   # Observers for button clicks
#   observeEvent(input$show_metadata, {
#     output_view("metadata")
#   })
#
#   observeEvent(input$show_count_table, {
#     output_view("count/Abundance_table")
#   })
#
#   observeEvent(input$show_summary_statistics, {
#     output_view("summary_statistics")
#   })
#
#   observeEvent(input$show_phylogenetic_tree, {
#     output_view("phylogenetic_tree")
#   })
#
#   observeEvent(input$show_combined_table, {
#     output_view("combined_table")
#   })
#
#   # Show the selected data section
#   output$dynamic_tables <- renderUI({
#     if (output_view() == "metadata") {
#       fluidRow(
#         box(title = "Metadata Structure", width = 12, status = "primary",
#             DT::dataTableOutput("metadata_structure")
#         )
#       )
#     } else if (output_view() == "combined_table") {
#       fluidRow(
#         box(title = "Abundance/Taxonomy Table", width = 12, status = "primary",
#             # Display combined table
#             DT::dataTableOutput("combined_table"),
#             # File format selection for download (under the table)
#             selectInput("filetype", "Choose file type:", choices = c("csv", "xlsx", "tsv")),
#             # Download button for the combined table
#             downloadButton("download_combined_table", "Download Abundance/Taxonomy Table")
#         )
#       )
#     } else if (output_view() == "summary_statistics") {
#       fluidRow(
#         box(title = "Summary Statistics", width = 12, status = "primary",
#             DT::dataTableOutput("summary_statistics")
#         )
#       )
#     } else if (output_view() == "phylogenetic_tree") {
#       fluidRow(
#         box(title = "Phylogenetic Tree", width = 12, status = "primary",
#             # Add dynamic options for customization
#             selectInput("phylum", "Choose Phylum:", choices = unique(tax_table(physeq())[, "Phylum"])),
#             selectInput("color_by", "Choose Color By (metadata):", sample_variables(physeq())),
#             selectInput("shape_by", "Choose Shape By (taxonomy):", choices = rank_names(physeq())),
#             selectInput("label_tips", "Label Tips By (taxonomy):", choices = rank_names(physeq())),
#             actionButton("plot_tree", "Plot Phylogenetic Tree"),
#             div(style = 'overflow-x: auto; overflow-y: auto; width: 100%; height: 800px;',
#                 plotOutput("phylogenetic_tree_plot")
#             ),
#             # Add download button for the tree
#             downloadButton("download_ggtree", "Download Phylogenetic Tree")
#         )
#       )
#     }
#   })
#
#
#
#   # Render Phylogenetic Tree using ggtree
#   output$phylogenetic_tree_plot <- renderPlot({
#     req(input$plot_tree)  # Wait for the button to be clicked
#
#     # Subset the tree by selected Phylum
#     GP <- prune_taxa(taxa_sums(physeq()) > 0, physeq())
#
#     # Subset by phylum and check if the phylum has any taxa
#     tree_subset <- subset_taxa(GP, Phylum == input$phylum)
#
#     # Check if the phylum subset is empty
#     if (is.null(tree_subset) || ntaxa(tree_subset) == 0) {
#       showNotification("No taxa found for the selected phylum", type = "error")
#       return(NULL)
#     }
#
#     # Check if color_by, shape_by, and label_tips are valid
#     if (is.null(input$color_by) || !input$color_by %in% names(sample_data(physeq()))) {
#       showNotification("Selected column for color is not found in the metadata", type = "error")
#       return(NULL)
#     }
#
#     if (is.null(input$shape_by) || !input$shape_by %in% names(tax_table(physeq()))) {
#       showNotification("Selected column for shape is not found in the taxonomy table", type = "error")
#       return(NULL)
#     }
#
#     if (is.null(input$label_tips) || !input$label_tips %in% names(tax_table(physeq()))) {
#       showNotification("Selected column for labeling tips is not found in the taxonomy table", type = "error")
#       return(NULL)
#     }
#
#     # Generate ggtree plot
#     tree_plot <- ggtree(tree_subset, ladderize = FALSE) +
#       geom_text2(aes(subset = !isTip, label = label), hjust = -0.2, size = 4) +
#       geom_tiplab(aes(label = get(input$label_tips)), hjust = -0.3) +
#       geom_point(aes(x = x + hjust, color = get(input$color_by), shape = get(input$shape_by)), na.rm = TRUE) +
#       scale_size_continuous(trans = log_trans(5)) +
#       theme(legend.position = "right") +
#       ggtitle(paste("Phylogenetic Tree (Phylum: ", input$phylum, ")", sep = ""))
#
#     print(tree_plot)
#   })
#
#
#
#   # Download ggtree plot
#   output$download_ggtree <- downloadHandler(
#     filename = function() {
#       paste0("phylogenetic_tree_", Sys.Date(), ".pdf")
#     },
#     content = function(file) {
#       # Render the ggtree plot for download
#       GP <- prune_taxa(taxa_sums(physeq()) > 0, physeq())
#       tree_subset <- subset_taxa(GP, Phylum == input$phylum)
#
#       tree_plot <- ggtree(tree_subset, ladderize = FALSE) +
#         geom_text2(aes(subset = !isTip, label = label), hjust = -0.2, size = 4) +
#         geom_tiplab(aes(label = get(input$label_tips)), hjust = -0.3) +
#         geom_point(aes(x = x + hjust, color = get(input$color_by), shape = get(input$shape_by),
#                        size=Abundance), na.rm = TRUE) +
#         scale_size_continuous(trans = log_trans(5)) +
#         theme(legend.position = "right") +
#         ggtitle(paste("Phylogenetic Tree (Phylum: ", input$phylum, ")", sep = ""))
#
#       # Save to file
#       ggsave(file, tree_plot, width = 10, height = 8)
#     }
#   )
#
#
#   # Render Metadata Structure when button is clicked
#   output$metadata_structure <- DT::renderDataTable({
#     req(output_view() == "metadata")
#     as.data.frame(sample_data(physeq()))
#   }, options = list(pageLength = 10, scrollX = TRUE))
#
#   # Store the combined table in a reactive expression to avoid duplication
#   combined_table <- reactive({
#     # Extract the OTU/ASV count table
#     otu_table <- as.data.frame(otu_table(physeq()))
#
#     # Extract the taxonomy table
#     taxonomy_table <- as.data.frame(tax_table(physeq()))
#
#     # Ensure that the OTU table and taxonomy table align by rows (taxa/ASVs)
#     if (nrow(taxonomy_table) != nrow(otu_table)) {
#       stop("Number of rows in the taxonomy table does not match the OTU table")
#     }
#
#     # Combine the OTU/ASV count table with the taxonomy table
#     combined <- merge(otu_table, taxonomy_table, by = 'row.names', all = TRUE)
#
#     # Rename the Row.names column to ASVs
#     colnames(combined)[1] <- "ASVs"
#
#     return(combined)
#   })
#
#   # Render Combined Abundance/Taxonomy Table
#   output$combined_table <- DT::renderDataTable({
#     req(output_view() == "combined_table")
#     DT::datatable(combined_table(), options = list(pageLength = 10, scrollX = TRUE), rownames = FALSE)
#   })
#
#   # Generate download handler for the Abundance/Taxonomy table
#   output$download_combined_table <- downloadHandler(
#     filename = function() {
#       paste0("abundance_taxonomy_table_", Sys.Date(), ".", input$filetype)
#     },
#     content = function(file) {
#       combined <- combined_table()  # Use the precomputed combined table
#
#       # Save the file based on the selected filetype
#       if (input$filetype == "csv") {
#         write.csv(combined, file, row.names = FALSE)
#       } else if (input$filetype == "xlsx") {
#         writexl::write_xlsx(combined, path = file)
#       } else if (input$filetype == "tsv") {
#         write.table(combined, file, sep = "\t", row.names = FALSE, quote = FALSE)
#       }
#     }
#   )
#
#   # Render Phylogenetic Tree with dynamically calculated plot size
#   output$phylogenetic_tree_plot <- renderPlot({
#     req(input$plot_tree)  # Wait for the button to be clicked
#
#     # Extract the phylogenetic tree from the phyloseq object
#     tree <- phy_tree(physeq())
#
#     # Dynamically calculate plot width and height based on the number of ASVs
#     num_asvs <- Ntip(tree)  # Get the number of tips (ASVs)
#
#     # Set label size proportional to the number of ASVs
#     label_size <- if (num_asvs > 500) 0.5 else 1  # Adjust label size if tree is too large
#
#     # Adjust the plot size dynamically
#     par(mar = c(5, 4, 4, 2) + 0.1)  # Adjust margins to avoid overlap
#
#     if (input$tree_type == "Unrooted") {
#       plot(tree, edge.width = 2, cex = label_size, main = "Unrooted Phylogenetic Tree")
#     } else if (input$tree_type == "Rooted") {
#       # Re-root the tree at a specific ASV (you can specify a default ASV, e.g., "ASV_1")
#       rerooted_tree <- root(tree, outgroup = "ASV_1", resolve.root = TRUE)
#       plot(rerooted_tree, edge.width = 2, cex = label_size, main = "Rooted Phylogenetic Tree")
#     }
#   },
#   # Dynamically calculate width and height within these arguments with a maximum limit
#   width = function() {
#     num_asvs <- Ntip(phy_tree(physeq()))  # Calculate the number of ASVs dynamically
#     min(max(400, num_asvs * 10), 3000)  # Cap the width at 3000 pixels
#   },
#   height = function() {
#     num_asvs <- Ntip(phy_tree(physeq()))  # Calculate the number of ASVs dynamically
#     min(max(400, num_asvs * 10), 3000)  # Cap the height at 3000 pixels
#   })
#
#
#   # Render Summary Statistics
#   output$summary_statistics <- DT::renderDataTable({
#     req(output_view() == "summary_statistics")
#     num_samples <- nsamples(physeq())
#     num_asvs <- ntaxa(physeq())
#     num_phylum <- length(unique(tax_table(physeq())[,"Phylum"]))
#     num_family <- length(unique(tax_table(physeq())[,"Family"]))
#     num_genus <- length(unique(tax_table(physeq())[,"Genus"]))
#     num_species <- length(unique(tax_table(physeq())[,"Species"]))
#
#     summary_df <- data.frame(
#       "Metric" = c("Number of Samples", "Number of ASVs", "Number of Phyla", "Number of Families", "Number of Genera", "Number of Species"),
#       "Count" = c(num_samples, num_asvs, num_phylum, num_family, num_genus, num_species)
#     )
#
#     DT::datatable(summary_df, options = list(pageLength = 10, scrollX = TRUE), rownames = FALSE)
#   })


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




                                                            ######### Reads QC

  # Explanation popup for Reads QC
  observeEvent(input$explanation_button_reads_qc, {
    showModal(modalDialog(
      title = "Reads QC Section Overview",
      p("This section provides quality control analysis for sequencing reads."),
      p("You can analyze the distribution of reads across samples or groups."),
      p("Choose a method to view the reads distribution, or filter samples based on the number of reads or their names. You can also download the plots and the filtered phyloseq object."),
      easyClose = TRUE,
      footer = modalButton("Close")
    ))
  })

  # Existing reactive value to track which button is pressed in Reads QC tab
  output_reads_qc_view <- reactiveVal("")

  # Observe button clicks and update reactive value accordingly
  observeEvent(input$show_reads_distribution_samples, {
    output_reads_qc_view("reads_distribution_samples")
  })

  observeEvent(input$show_reads_distribution_groups, {
    output_reads_qc_view("reads_distribution_groups")
  })

  observeEvent(input$show_number_of_reads, {
    output_reads_qc_view("number_of_reads")
  })

  # Dynamically render UI content based on the button clicked
  output$reads_qc_content <- renderUI({
    if (output_reads_qc_view() == "reads_distribution_samples") {
      fluidRow(
        box(title = "Read Distribution Across Samples", width = 12, status = "primary",
            plotlyOutput("reads_distribution_plot_samples"),
            selectInput("plot_filetype_samples", "Select File Type", choices = c("png", "pdf", "svg")),
            numericInput("plot_width_samples", "Plot Width (inches)", value = 8, min = 4),
            numericInput("plot_height_samples", "Plot Height (inches)", value = 6, min = 4),
            downloadButton("download_plot_samples", "Download Plot")
        )
      )
    } else if (output_reads_qc_view() == "reads_distribution_groups") {
      fluidRow(
        box(title = "Read Distribution Across Groups", width = 12, status = "primary",
            selectInput("group_column", "Select Grouping Column", choices = sample_variables(physeq())),
            plotlyOutput("reads_distribution_plot_groups"),
            selectInput("plot_filetype_groups", "Select File Type", choices = c("png", "pdf", "svg")),
            numericInput("plot_width_groups", "Plot Width (inches)", value = 8, min = 4),
            numericInput("plot_height_groups", "Plot Height (inches)", value = 6, min = 4),
            downloadButton("download_plot_groups", "Download Plot")
        )
      )
    } else if (output_reads_qc_view() == "number_of_reads") {
      fluidRow(
        box(title = "Number of Reads per Sample", width = 12, status = "primary",
            DT::dataTableOutput("number_of_reads_table"),
            numericInput("reads_threshold", "Reads Threshold for Sample Filtration", value = 1000, min = 0),
            selectInput("sample_name_filter", "Filter by Sample Name(s)", choices = sample_names(physeq()), multiple = TRUE),
            actionButton("filter_samples", "Filter Samples"),
            DT::dataTableOutput("filtered_samples_table"),
            downloadButton("download_filtered_physeq", "Download Filtered Phyloseq")
        )
      )
    }
  })

  # Plot for Reads Distribution Across Samples (Base R Plot)
  output$reads_distribution_plot_samples <- renderPlotly({
    req(physeq())
    dat <- data.frame(sample_data(physeq()), num_reads = sample_sums(physeq())) %>% rownames_to_column("sample_names")
    p <- ggplot(dat, aes(x = sample_names, y = num_reads)) +
      geom_bar(stat = 'identity') +
      labs(x='', y = 'Number of reads', colour = NULL) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    
    ggplotly(p)
  })

  # Plot for Reads Distribution Across Groups (Using plot_read_distribution and ggplotly)
  output$reads_distribution_plot_groups <- renderPlotly({
    req(physeq(), input$group_column)
    dat2 <- data.frame(sample_data(physeq()), num_reads = sample_sums(physeq()))
    p <- ggplot(dat2, aes(x = factor(!!sym(input$group_column)), y = num_reads, fill = factor(!!sym(input$group_column)))) +
      geom_boxplot(alpha=0.7) + ggtitle(paste0('Read distribution by ', input$group_column)) +
      labs(x='', y = 'Number of reads', fill = input$group_column) +
      theme_bw()
    ggplotly(p)
  })

  # Download handler for samples plot
  output$download_plot_samples <- downloadHandler(
    filename = function() {
      paste0("reads_distribution_samples_", Sys.Date(), ".", input$plot_filetype_samples)
    },
    content = function(file) {
      ggsave(file, plot = ggplot(data.frame(sample_data(physeq()), num_reads = sample_sums(physeq())) %>% rownames_to_column("sample_names"), 
                                 aes(x = sample_names, y = num_reads)) +
               geom_bar(stat = 'identity') +
               labs(x='', y = 'Number of reads', colour = NULL) +
               theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)),
             device = input$plot_filetype_samples, width = input$plot_width_samples, height = input$plot_height_samples)
    }
  )

  # Download handler for groups plot
  output$download_plot_groups <- downloadHandler(
    filename = function() {
      paste0("reads_distribution_groups_", Sys.Date(), ".", input$plot_filetype_groups)
    },
    content = function(file) {
      ggsave(file, plot =  ggplot(data.frame(sample_data(physeq()), num_reads = sample_sums(physeq())), 
                                  aes(x = factor(!!sym(input$group_column)), y = num_reads, fill = factor(!!sym(input$group_column)))) +
               geom_boxplot(alpha=0.7) + ggtitle(paste0('Read distribution by ', input$group_column)) +
               labs(x='', y = 'Number of reads', fill = input$group_column) +
               theme_bw(),
             device = input$plot_filetype_groups, width = input$plot_width_groups, height = input$plot_height_groups)
    }
  )

  # Table for Number of Reads per Sample
  output$number_of_reads_table <- DT::renderDataTable({
    req(physeq())
    reads_per_sample <- sample_sums(physeq())
    sample_names <- sample_names(physeq())
    reads_df <- data.frame(Sample = sample_names, Reads = reads_per_sample)

    DT::datatable(reads_df, options = list(pageLength = 10, scrollX = TRUE))
  })

  # Filtering samples based on reads threshold and sample names, triggered by button click
  observeEvent(input$filter_samples, {
    req(input$reads_threshold, physeq())
    pseq <- physeq()

    # Apply threshold filter first
    pseq <- prune_samples(sample_sums(pseq) >= input$reads_threshold, pseq)
    pseq <- prune_taxa(taxa_sums(pseq) > 0, pseq)

    # Apply sample name filter to remove selected samples
    if (!is.null(input$sample_name_filter) && length(input$sample_name_filter) > 0) {
      pseq <- prune_samples(!sample_names(pseq) %in% input$sample_name_filter, pseq)
      pseq <- prune_taxa(taxa_sums(pseq) > 0, pseq)
    }

    # Update the filtered samples table to show remaining samples after filtration
    output$filtered_samples_table <- DT::renderDataTable({
      filtered_samples <- sample_sums(pseq)
      filtered_samples_df <- data.frame(Sample = names(filtered_samples), Reads = filtered_samples)

      DT::datatable(filtered_samples_df, options = list(pageLength = 10, scrollX = TRUE))
    })

    # Download handler for filtered phyloseq object
    output$download_filtered_physeq <- downloadHandler(
      filename = function() {
        paste("filtered_phyloseq_", Sys.Date(), ".rds", sep = "")
      },
      content = function(file) {
        saveRDS(pseq, file)
      }
    )

    # Stay on the "Number of Reads per Sample" section after filtering
    output_reads_qc_view("number_of_reads")
  })




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

  observeEvent(input$show_dominant_taxa_groups, {
    output_taxa_overview_view("dominant_taxa_groups")
  })

  observeEvent(input$show_dominant_taxa_samples, {
    output_taxa_overview_view("dominant_taxa_samples")
  })

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
    } else if (output_taxa_overview_view() == "dominant_taxa_groups") {
      fluidRow(
        box(title = "Dominant Taxa Across Groups", width = 12, status = "primary",
            selectInput("dominant_taxa_level_groups", "Select Taxonomic Level", choices = rank_names(physeq())),
            selectInput("group_column", "Select Grouping Column", choices = sample_variables(physeq())),
            selectInput("dominance_method", "Select Dominance Calculation Method",
                        choices = c("Abundance-Based Dominance", "Prevalence-Based Dominance", "Combined Metrics"),
                        selected = "Abundance-Based Dominance"),
            DT::dataTableOutput("dominant_taxa_table_groups"),
            plotOutput("dominant_taxa_plot_groups"),
            downloadButton("download_dominant_taxa_groups", "Download Dominant Taxa Data")
        )
      )
    } else if (output_taxa_overview_view() == "dominant_taxa_samples") {
      fluidRow(
        box(title = "Dominant Taxa Per Sample", width = 12, status = "primary",
            selectInput("dominant_taxa_level_samples", "Select Taxonomic Level", choices = rank_names(physeq())),
            DT::dataTableOutput("dominant_taxa_table_samples"),
            downloadButton("download_dominant_taxa_samples", "Download Dominant Taxa Data")
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

  observeEvent(input$explain_dominant_taxa_groups, {
    showModal(modalDialog(
      title = "Dominant Taxa Across Groups - Explanation",
      p("This section allows you to identify the dominant taxa across different groups using selected criteria such as abundance, prevalence, or combined metrics."),
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
  output$taxa_prevalence_table_samples <- DT::renderDataTable({
    req(input$taxa_level_samples, physeq())
    taxa_level <- input$taxa_level_samples

    # Aggregate OTU table at the selected taxonomic level
    physeq_agg <- aggregate_taxa(physeq(), taxa_level)

    prevalence_df <- apply(as.data.frame(otu_table(physeq_agg)), 1, function(x) sum(x > 0) / nsamples(physeq_agg) * 100)
    prevalence_df <- round(prevalence_df, 3)  # Round to 3 decimal places

    absolute_abundance <- rowSums(otu_table(physeq_agg))

    prevalence_table <- data.frame(Species = tax_table(physeq_agg)[, taxa_level], Prevalence = prevalence_df, Abundance = absolute_abundance)
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
  output$taxa_prevalence_table_groups <- DT::renderDataTable({
    req(input$taxa_level_groups, input$group_column, physeq())
    taxa_level <- input$taxa_level_groups
    group_column <- input$group_column

    # Aggregate OTU table at the selected taxonomic level
    physeq_agg <- aggregate_taxa(physeq(), taxa_level)

    # Get the OTU table and the grouping variable
    otu_table_agg <- as.data.frame(otu_table(physeq_agg))
    groups <- as.factor(sample_data(physeq_agg)[[group_column]])

    # Calculate prevalence within each group
    prevalence_list <- lapply(levels(groups), function(group) {
      group_samples <- which(groups == group)
      apply(otu_table_agg[, group_samples], 1, function(x) sum(x > 0) / length(group_samples) * 100)
    })

    # Combine the prevalence results into a single data frame
    prevalence_df <- do.call(cbind, prevalence_list)
    colnames(prevalence_df) <- levels(groups)

    # Create the prevalence table
    prevalence_table <- data.frame(Species = tax_table(physeq_agg)[, taxa_level], prevalence_df)
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



  # Dominant Taxa Across Groups

  # Create a reactive value to store the dominant taxa table
  dominant_table_groups_reactive <- reactiveVal()

  # Calculate dominant taxa across groups based on selected method
  output$dominant_taxa_table_groups <- DT::renderDataTable({
    req(input$dominant_taxa_level_groups, input$group_column, physeq(), input$dominance_method)

    taxa_level <- input$dominant_taxa_level_groups
    group_column <- input$group_column
    method <- input$dominance_method

    # Aggregate OTU table at the selected taxonomic level
    physeq_agg <- microbiome::aggregate_taxa(physeq(), taxa_level)

    # Get OTU and taxonomic data frames
    abundance_df <- as.data.frame(otu_table(physeq_agg))
    taxa_df <- tax_table(physeq_agg)[, taxa_level, drop = FALSE]
    sample_data_df <- as.data.frame(sample_data(physeq_agg))

    # Combine abundance data with group information
    combined_df <- abundance_df %>%
      rownames_to_column(var = "Taxa") %>%
      pivot_longer(cols = -Taxa, names_to = "SampleID", values_to = "Abundance") %>%
      dplyr::left_join(sample_data_df %>% rownames_to_column(var = "SampleID"), by = "SampleID") %>%
      group_by(!!sym(group_column), Taxa) %>%
      summarise(
        MeanAbundance = mean(Abundance, na.rm = TRUE),
        Prevalence = sum(Abundance > 0) / n() * 100,
        .groups = "drop"
      )

    # Apply selected dominance method
    dominant_taxa <- combined_df %>%
      group_by(!!sym(group_column))

    if (method == "Abundance-Based Dominance") {
      dominant_taxa <- dominant_taxa %>% slice_max(MeanAbundance, n = 1)
    } else if (method == "Prevalence-Based Dominance") {
      dominant_taxa <- dominant_taxa %>% slice_max(Prevalence, n = 1)
    } else if (method == "Combined Metrics") {
      # Example of combined metrics (mean of normalized abundance and prevalence)
      dominant_taxa <- dominant_taxa %>%
        mutate(CombinedScore = (MeanAbundance / max(MeanAbundance)) + (Prevalence / max(Prevalence))) %>%
        slice_max(CombinedScore, n = 1)
    }

    dominant_taxa <- dominant_taxa %>% ungroup()

    # Store the table in the reactive value for later use in downloadHandler
    dominant_table_groups_reactive(dominant_taxa)

    # Render the table without row names
    DT::datatable(dominant_taxa, options = list(pageLength = 10, scrollX = TRUE), rownames = FALSE)
  })

  # Download handler for dominant taxa across groups
  output$download_dominant_taxa_groups <- downloadHandler(
    filename = function() {
      paste("dominant_taxa_groups_", Sys.Date(), ".tsv", sep = "")
    },
    content = function(file) {
      # Access the data stored in the reactive value
      dominant_table_groups <- dominant_table_groups_reactive()
      write.table(dominant_table_groups, file, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
    }
  )


  # Render the UI with the table and the note
  output$dominant_taxa_ui <- renderUI({
    tagList(
      DT::dataTableOutput("dominant_taxa_table_groups"),
      uiOutput("dominant_taxa_overview_note")
    )
  })


  #dominant taxa per samples

  # Create a reactive value to store the dominant taxa table
  dominant_table_reactive <- reactiveVal()

  # Calculate dominant taxa per sample
  output$dominant_taxa_table_samples <- DT::renderDataTable({
    req(input$dominant_taxa_level_samples, physeq())

    taxa_level <- input$dominant_taxa_level_samples

    # Aggregate OTU table at the selected taxonomic level
    physeq_agg <- aggregate_taxa(physeq(), taxa_level)

    # Get OTU and taxonomic data frames
    abundance_df <- as.data.frame(otu_table(physeq_agg))
    taxa_df <- tax_table(physeq_agg)[, taxa_level, drop = FALSE]

    # Initialize vectors to store results
    dominant_taxa <- vector("character", ncol(abundance_df))
    abundance_vals <- vector("numeric", ncol(abundance_df))

    # Loop through each sample to find the dominant taxa
    for (i in seq_along(dominant_taxa)) {
      sample_abundances <- abundance_df[, i]
      dominant_index <- which.max(sample_abundances)
      dominant_taxa[i] <- taxa_df[dominant_index, taxa_level]
      abundance_vals[i] <- sample_abundances[dominant_index]
    }

    # Create a data frame for the results
    dominant_table <- data.frame(
      Sample = colnames(abundance_df),
      Dominant_Taxa = dominant_taxa,
      Abundance = abundance_vals
    )

    # Store the table in the reactive value for later use in downloadHandler
    dominant_table_reactive(dominant_table)

    # Render the table without row names
    DT::datatable(dominant_table, options = list(pageLength = 10, scrollX = TRUE), rownames = FALSE)
  })


  # Download handler for dominant taxa per sample
  output$download_dominant_taxa_samples <- downloadHandler(
    filename = function() {
      paste("dominant_taxa_samples_", Sys.Date(), ".tsv", sep = "")
    },
    content = function(file) {
      # Access the data stored in the reactive value
      dominant_table <- dominant_table_reactive()
      write.table(dominant_table, file, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

    }
  )



                                                  ########### compositional barplot

  compositional_barplot_reactive <- reactive({
    req(physeq(), input$hueRank, input$shadeRank, input$facetBy)

    hueRank <- input$hueRank
    shadeRank <- input$shadeRank

    pseq2 <- physeq() %>%
      tax_sort(by = sum, at = shadeRank) %>%
      tax_sort(by = sum, at = hueRank) %>%
      tax_agg(rank = shadeRank)

    nHues <- input$nHues
    nShades <- input$nShades

    hierarchicalPalInfo <- data.frame(
      hue = as.vector(tt_get(pseq2)[, hueRank]),
      shade = as.vector(tt_get(pseq2)[, shadeRank]),
      counts = taxa_sums(otu_get(pseq2))
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
    palette_colors <- colorRampPalette(brewer.pal(12, "Set3"))(length(unique_taxa))
    names(palette_colors) <- unique_taxa

    facet_by <- if (input$facetBy == "None") NULL else input$facetBy

    p <- pseq2 %>%
      ps_get() %>%
      tax_mutate("Taxa" = hierarchicalPalInfo$Taxa, .keep = "none") %>%
      comp_barplot(
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
    req(physeq(), input$normalizationMethod, input$heatmapRank, input$annotationColumn1, input$annotationColumn2, input$topTaxa)

    # Conditionally aggregate the phyloseq object based on the user's selection
    aggregated_physeq <- reactive({
      if (input$heatmapRank == "ASV") {
        # Only normalize the physeq object without aggregation
        microbiome::transform(physeq(), input$normalizationMethod)
      } else {
        # Both normalize and aggregate at the selected taxonomic level
        tax_transform(physeq(), input$normalizationMethod, rank = input$heatmapRank)
      }
    })

    # Ensure aggregated_physeq has data
    pseq <- aggregated_physeq()
    req(nrow(otu_table(pseq)) > 0)

    # Select top taxa
    top_taxa <- names(sort(taxa_sums(pseq), decreasing = TRUE))[1:input$topTaxa]
    psq_normalized_pruned <- prune_taxa(top_taxa, pseq)

    # Define annotation colors dynamically based on the number of unique values
    unique_vals1 <- unique(sample_data(psq_normalized_pruned)[[input$annotationColumn1]])
    cols1 <- distinct_palette(n = length(unique_vals1), add = NA)
    names(cols1) <- unique_vals1

    unique_vals2 <- unique(sample_data(psq_normalized_pruned)[[input$annotationColumn2]])
    cols2 <- distinct_palette(n = length(unique_vals2), add = NA)
    names(cols2) <- unique_vals2

    # Prepare sample annotations
    sample_anno <- sampleAnnotation(
      State1 = anno_sample_cat(input$annotationColumn1, legend_title = input$annotationColumn1),
      col = list(State1 = cols1, State2 = cols2), border = FALSE,
      State2 = anno_sample_cat(input$annotationColumn2, col = cols2, legend_title = input$annotationColumn2),
      annotation_label = c(input$annotationColumn1, input$annotationColumn2)
    )

    # Generate the heatmap using comp_heatmap
    heatmap_obj <- comp_heatmap(
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
  # ht1 <- draw(heatmap_obj)
    # Render the interactive heatmap using InteractiveComplexHeatmap
    makeInteractiveComplexHeatmap(input, output, session, ht1,
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




                                                            ######## Venn diagram

  # observe({
  #   req(physeq())  # Ensure physeq is available
  #   updateSelectInput(session, "vennGroupingColumn", choices = colnames(sample_data(physeq())))
  # })
  #
  # venn_data <- reactive({
  #   # Ensure that physeq is available and that inputs are valid
  #   req(physeq(), input$vennGroupingColumn, input$vennTaxonomicRank)
  #
  #   # Attempt to assign the physeq object to pseq
  #   pseq <- physeq()  # Assign physeq object to pseq
  #
  #   # Debugging output: Ensure pseq is assigned and is valid
  #   if (is.null(pseq)) {
  #     stop("Error: pseq is NULL. The physeq() function did not return a valid phyloseq object.")
  #   }
  #   print("pseq has been assigned successfully.")
  #
  #   # Conditionally aggregate the phyloseq object based on the user's selection
  #   if (input$vennTaxonomicRank != "ASV") {
  #     pseq <- microbiome::aggregate_taxa(pseq, input$vennTaxonomicRank)
  #     print(paste("Aggregated pseq at rank:", input$vennTaxonomicRank))
  #   }
  #
  #   # Convert to relative abundances
  #   pseq <- microbiome::transform(pseq, "compositional")
  #   print("Transformed pseq to relative abundances.")
  #
  #   group_column <- input$vennGroupingColumn
  #
  #   # Debugging: Ensure that group_column exists in sample_data
  #   if (!group_column %in% colnames(sample_data(pseq))) {
  #     stop("Error: Selected grouping column not found in the sample data.")
  #   }
  #
  #   # Extract unique group states from the selected grouping column
  #   group_states <- unique(sample_data(pseq)[[group_column]])
  #
  #   # Handle the case where group_states is empty
  #   if (length(group_states) == 0) {
  #     stop("Error: No unique group states found in the selected grouping column.")
  #   }
  #   req(group_states)
  #
  #   list_core <- list()
  #   for (n in group_states) {
  #     ps_sub <- subset_samples(pseq, sample_data(pseq)[[group_column]] == n)
  #     core_m <- core_members(ps_sub, detection = 0.001, prevalence = 0.75)
  #     list_core[[n]] <- core_m
  #   }
  #
  #   return(list_core)
  # })
  #
  # output$vennDiagramPlot <- renderPlot({
  #   req(input$plotVenn)  # Wait until the button is clicked
  #   data <- venn_data()
  #   req(data)
  #
  #   plot(venn(data))
  # })
  #
  #
  #
  # observe({
  #   req(physeq())  # Ensure physeq is available
  #   updateSelectInput(session, "vennGroupingColumn", choices = colnames(sample_data(physeq())))
  # })
  #
  # venn_data <- reactive({
  #   # Ensure that physeq is available and that inputs are valid
  #   req(physeq(), input$vennGroupingColumn, input$vennTaxonomicRank)
  #
  #   # Attempt to assign the physeq object to pseq
  #   pseq <- physeq()  # Assign physeq object to pseq
  #
  #   # Debugging output: Ensure pseq is assigned and is valid
  #   if (is.null(pseq)) {
  #     stop("Error: pseq is NULL. The physeq() function did not return a valid phyloseq object.")
  #   }
  #   print("pseq has been assigned successfully.")
  #
  #   # Conditionally aggregate the phyloseq object based on the user's selection
  #   if (input$vennTaxonomicRank != "ASV") {
  #     pseq <- microbiome::aggregate_taxa(pseq, input$vennTaxonomicRank)
  #     print(paste("Aggregated pseq at rank:", input$vennTaxonomicRank))
  #   }
  #
  #   # Convert to relative abundances
  #   pseq <- microbiome::transform(pseq, "compositional")
  #   print("Transformed pseq to relative abundances.")
  #
  #   group_column <- input$vennGroupingColumn
  #
  #   # Debugging: Ensure that group_column exists in sample_data
  #   if (!group_column %in% colnames(sample_data(pseq))) {
  #     stop("Error: Selected grouping column not found in the sample data.")
  #   }
  #
  #   # Extract unique group states from the selected grouping column
  #   group_states <- unique(sample_data(pseq)[[group_column]])
  #
  #   # Handle the case where group_states is empty
  #   if (length(group_states) == 0) {
  #     stop("Error: No unique group states found in the selected grouping column.")
  #   }
  #   req(group_states)
  #
  #   list_core <- list()
  #   for (n in group_states) {
  #     ps_sub <- subset_samples(pseq, sample_data(pseq)[[group_column]] == n)
  #     core_m <- core_members(ps_sub, detection = 0.001, prevalence = 0.75)
  #     list_core[[n]] <- core_m
  #   }
  #
  #   return(list_core)
  # })
  #
  # output$vennDiagramPlot <- renderPlot({
  #   req(input$plotVenn)  # Wait until the button is clicked
  #   data <- venn_data()
  #   req(data)
  #
  #   plot(venn(data))
  # })
  #
  #



                                                              ###### Alpha diversity
  observe({
    req(physeq())
    pseq <- physeq()

    # Calculate diversity metrics
    tab <- microbiome::alpha(physeq(), index = "all")

    # Extract sample data and add diversity metrics
    meta <- as(sample_data(pseq), "data.frame")
    for (metric in colnames(tab)) {
      meta[[metric]] <- tab[[metric]]
    }

    # Update the physeq object with the new sample data
    updated_physeq <- physeq()
    sample_data(updated_physeq) <- sample_data(meta)

    # Update selectInput choices
    updateSelectInput(session, "alphaMetric", choices = colnames(tab))
    updateSelectInput(session, "alphaGroupingColumn", choices = sample_variables(updated_physeq))

    # Update the reactive physeq object
    isolate({
      physeq <<- reactive({ updated_physeq })
    })
  })

  # Reactive expression to store the ggplot2 and plotly plots
  plotly_plots <- reactiveVal(list())
  ggplot_plots <- reactiveVal(list())



  # Generate Alpha Diversity Plot
  observeEvent(input$plotAlpha, {
    req(physeq(), input$alphaMetric, input$alphaGroupingColumn)
    pseq <- physeq()

    meta <- as(sample_data(pseq), "data.frame")
    plotly_lst <- list()
    ggplot_lst <- list()

    clusters <- levels(factor(meta[[input$alphaGroupingColumn]]))
    comparisons <- combn(seq_along(clusters), 2, simplify = FALSE, FUN = function(i) clusters[i])

    custom_colors <- c("#0072B2", "#009E73", "#F0E442", "#56B4E9", "#D55E00", "#CC79A7",
                       "indianred", "darkblue", "seagreen", "skyblue", "lightsteelblue2",
                       "plum2")


    for (metric in input$alphaMetric) {
      p <- ggplot(meta, aes_string(x = input$alphaGroupingColumn, y = metric, color = input$alphaGroupingColumn)) +
        geom_boxplot(outlier.shape = NA, width = 0.4, alpha = 0.75) +
        geom_jitter(height = 0, width = 0.2, alpha = 0.5, size = 2) +
        stat_pwc(method = "wilcox_test",ref.group= "all",  label = "{p.adj.format}{p.adj.signif}",
                 p.adjust.method = "hochberg", p.adjust.by = "panel",remove.bracket=FALSE,
                 step.increase = 0.12, tip.length = 0.03, size = 0.3)+
        scale_color_manual(values = custom_colors) +  # Use only scale_color_manual
        theme_minimal(base_size = 15) +
        labs(y = metric, x = "Groups") +
        theme(
          legend.position = "right",  # Keep the legend position in ggplot
          plot.title = element_text(hjust = 0.5),
          axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1),
          plot.margin = margin(t = 40, r = 40, b = 40, l = 60, unit = "pt")
        )

      plotly_lst[[length(plotly_lst) + 1]] <- ggplotly(p) %>%
        layout(
          xaxis = list(tickangle = -45, title = "", automargin = TRUE),
          yaxis = list(title = metric, automargin = TRUE),
          showlegend = FALSE,  # Remove legends from individual plots
          margin = list(l = 120, r = 40, t = 40, b = 40)  # Increased left margin
        )

      ggplot_lst[[length(ggplot_lst) + 1]] <- p
    }

    plotly_plots(plotly_lst)  # Store the plotly plots in the reactiveVal
    ggplot_plots(ggplot_lst)  # Store the ggplot plots in the reactiveVal
  })

  # Render the Alpha Diversity Plot
  output$alphaDiversityPlot <- renderPlotly({
    plots <- plotly_plots()
    req(length(plots) > 0)  # Ensure plots are available before rendering


      if (length(plots) == 1) {
        plots[[1]] %>%
          layout(dragmode = "zoom", margin = list(l = 120, r = 40, t = 40, b = 40), showlegend = TRUE) %>%  # Show combined legend
          config(displayModeBar = TRUE, scrollZoom = TRUE)  # Enable interactive features
      } else {
        subplot(plots, nrows = ceiling(length(plots) / 2), shareX = TRUE, shareY = FALSE, titleX = TRUE, titleY = TRUE) %>%
          layout(dragmode = "zoom", margin = list(l = 120, r = 40, t = 40, b = 40), showlegend = TRUE) %>%  # Show combined legend
          config(displayModeBar = TRUE, scrollZoom = TRUE)
      }
    })


  # Reactive expression to store the Wilcoxon test results
  alpha_stats <- reactiveVal(data.frame())

  # Calculate Wilcoxon test statistics when the button is clicked
  observeEvent(input$showStats, {
    req(physeq(), input$alphaMetric, input$alphaGroupingColumn)
    pseq <- physeq()
    meta <- as(sample_data(pseq), "data.frame")

    clusters <- levels(factor(meta[[input$alphaGroupingColumn]]))
    comparisons <- combn(seq_along(clusters), 2, simplify = FALSE, FUN = function(i) clusters[i])

    stats_list <- list()

    for (metric in input$alphaMetric) {
      p_values <- sapply(comparisons, function(comp) {
        group1 <- meta[meta[[input$alphaGroupingColumn]] == comp[1], metric]
        group2 <- meta[meta[[input$alphaGroupingColumn]] == comp[2], metric]
        wilcox.test(group1, group2)$p.value
      })

      # Apply BH correction to the p-values
      adj_p_values <- p.adjust(p_values, method = "BH")

      # Create a data frame for this metric
      stats_df <- data.frame(
        Metric = metric,
        Comparison = sapply(comparisons, function(comp) paste(comp, collapse = " vs ")),
        P_value = round(p_values, 4),  # Round p-values to 4 decimal places
        Adjusted_P_value = round(adj_p_values, 4),  # Round adjusted p-values to 4 decimal places
        stringsAsFactors = FALSE
      )


      stats_list[[metric]] <- stats_df
    }

    # Combine all metrics into one data frame
    full_stats <- do.call(rbind, stats_list)

    # Update the reactive value
    alpha_stats(full_stats)
  })

  # Render the statistics table using DT
  output$statsTable <- renderDT({
    req(alpha_stats())
    datatable(alpha_stats(), options = list(pageLength = 10, autoWidth = TRUE),
              editable = 'cell', rownames = FALSE,
              caption = 'Table 1: Wilcoxon Rank-Sum Test Results: Pairwise comparisons between groups with p-values adjusted using the Benjamini-Hochberg method.')
  })

  output$downloadStats <- downloadHandler(
    filename = function() {
      paste0("alpha_diversity_stats_", Sys.Date(), ".tsv")
    },
    content = function(file) {
      stats <- alpha_stats()
      req(stats)  # Ensure stats are available before downloading

      print("Writing stats to file...")
      print(file)  # Print the file path for debugging

      # Write the stats to a TSV file
      write.table(stats, file, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

      print("File written successfully")
    }
  )

  # Download handler for the alpha diversity plot (html, pdf, svg, and png are possible)
  output$download_alphaDiversityPlot <- downloadHandler(
    filename = function() {
      paste0("alpha_diversity_", Sys.Date(), ".", input$alpha_filetype)
    },
    content = function(file) {
      if (input$alpha_filetype == "html") {
        # Re-generate the Plotly object to ensure it's updated
        plots <- plotly_plots()
        p <- do.call(subplot, c(plots, list(nrows = ceiling(length(plots) / 2), shareX = TRUE, shareY = FALSE)))
        htmlwidgets::saveWidget(as_widget(p), file)
      } else {
        # Re-generate the ggplot object to ensure it's updated
        plots <- ggplot_plots()
        p <- marrangeGrob(plots, ncol = 2, nrow = ceiling(length(plots) / 2))
        ggsave(file, plot = p, width = input$plot_width, height = input$plot_height, units = "in",
               device = input$alpha_filetype, limitsize = FALSE)  # Ensure width and height are from the latest inputs
      }
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

  # UI updates
  observe({
    req(physeq())
    pseq <- physeq()
    updateSelectInput(session, "taxRank", choices = rank_names(pseq))
    updateSelectInput(session, "groupingColumnBeta", choices = colnames(sample_data(physeq())))
    updateSelectInput(session, "shapeColumnBeta", choices = colnames(sample_data(physeq())))
  })

  # Initialize reactive values to store plots
  ggplot_obj <- reactiveVal(NULL)  # Use reactiveVal to store the ggplot object



  ### PCoA ###
  observeEvent(input$plotPCoA, {
    req(physeq(), input$distanceMetric, input$normalizationMethod, input$taxRank, input$groupingColumnBeta, input$shapeColumnBeta)

    print("PCoA button clicked")

    # Access the actual phyloseq object
    pseq <- physeq()

    tryCatch({
      if (input$distanceMetric == "aitchison" && input$normalizationMethod != "identity") {
        stop("'Aitchison' distance requires count data. Please select 'identity' as the normalization method or choose another distance metric.")
      }

      if (input$distanceMetric == "gunifrac") {
        pseq_transformed <- pseq %>%
          tax_transform(input$normalizationMethod, rank = "unique") %>%
          dist_calc("gunifrac", gunifrac_alpha = 0.5)
      } else {
        pseq_transformed <- pseq %>%
          tax_transform(input$normalizationMethod, rank = input$taxRank) %>%
          dist_calc(dist = input$distanceMetric)

        if (input$distanceMetric == "jaccard") {
          pseq_transformed <- dist_calc(pseq_transformed, dist = "jaccard", binary = TRUE)
        }
      }


      # Perform PCoA ordination
      ordination_res <- pseq_transformed %>%
        ord_calc("PCoA")

      # Generate the plot using ord_plot and additional ggplot layers
      plot <- ord_plot(ordination_res, color = input$groupingColumnBeta, shape = input$shapeColumnBeta,
                       plot_taxa = 1:5, size = 2) +
        stat_ellipse(linetype = 1, segments = 10, lwd = 1.2, alpha = 0.25, level = 0.6,
                     geom = "polygon", aes(fill = input$groupingColumnBeta)) +
        scale_colour_brewer(palette = "Dark2", aesthetics = c("fill", "colour")) +
        theme_bw() +
        ggside::geom_xsideboxplot(aes_string(fill = input$groupingColumnBeta, y = input$groupingColumnBeta), orientation = "y") +
        ggside::geom_ysideboxplot(aes_string(fill = input$groupingColumnBeta, x = input$groupingColumnBeta), orientation = "x") +
        ggside::scale_xsidey_discrete(labels = NULL) +
        ggside::scale_ysidex_discrete(labels = NULL) +
        ggside::theme_ggside_void()

      # Store the ggplot object in the reactive value
      ggplot_obj(plot)

      # Render the ggplot
      output$betaDiversityPlot <- renderPlot({
        print(ggplot_obj())
      })

    }, error = function(e) {
      output$betaDiversityNote <- renderText({ e$message })
      showModal(modalDialog(
        title = "Error",
        e$message,
        easyClose = TRUE,
        footer = modalButton("Close")
      ))
    })
  })

  # Download handler for PCoA plot
  output$download_PCoAPlot <- downloadHandler(
    filename = function() {
      paste0("PCoA_beta_diversity_", Sys.Date(), ".", input$beta_filetype)
    },
    content = function(file) {
      req(ggplot_obj())  # Ensure ggplot_obj contains a valid plot
      ggsave(file, plot = ggplot_obj(), width = input$plot_width, height = input$plot_height, units = "in", device = input$beta_filetype)
    }
  )



  # Initialize reactive value for ggplot object
  ggplot_pca_obj <- reactiveVal(NULL)  # Use reactiveVal to store the PCA ggplot object

  ### PCA ###
  observeEvent(input$plotPCA, {
    req(physeq(), input$normalizationMethod, input$taxRank, input$groupingColumnBeta, input$shapeColumnBeta)

    print("PCA button clicked")

    # Access the actual phyloseq object
    pseq <- physeq()

    tryCatch({
      # Perform PCA ordination directly on the transformed data
      ordination_res <- pseq %>%
        tax_transform(input$normalizationMethod, rank = input$taxRank) %>%
        ord_calc("PCA")

      # Generate the plot using ord_plot and additional ggplot layers
      plot <- ord_plot(ordination_res, color = input$groupingColumnBeta, shape = input$shapeColumnBeta,
                       plot_taxa = 1:5, size = 2) +
        scale_colour_brewer(palette = "Dark2", aesthetics = c("fill", "colour")) +
        theme_bw() +
        ggside::geom_xsideboxplot(aes_string(fill = input$groupingColumnBeta, y = input$groupingColumnBeta), orientation = "y") +
        ggside::geom_ysideboxplot(aes_string(fill = input$groupingColumnBeta, x = input$groupingColumnBeta), orientation = "x") +
        ggside::scale_xsidey_discrete(labels = NULL) +
        ggside::scale_ysidex_discrete(labels = NULL) +
        ggside::theme_ggside_void()

      # Store the ggplot object in the reactive value
      ggplot_pca_obj(plot)

      # Directly render the ggplot
      output$betaDiversityPlot <- renderPlot({
        print(ggplot_pca_obj())
      })

    }, error = function(e) {
      output$betaDiversityNote <- renderText({ e$message })
      showModal(modalDialog(
        title = "Error",
        e$message,
        easyClose = TRUE,
        footer = modalButton("Close")
      ))
    })
  })

  # Download handler for PCA plot
  output$download_PCAPlot <- downloadHandler(
    filename = function() {
      paste0("PCA_beta_diversity_", Sys.Date(), ".", input$beta_filetype)
    },
    content = function(file) {
      req(ggplot_pca_obj())  # Ensure ggplot_pca_obj contains a valid plot
      ggsave(file, plot = ggplot_pca_obj(), width = input$plot_width, height = input$plot_height, units = "in", device = input$beta_filetype)
    }
  )


  # Initialize a reactive value to store the NMDS plot object
  ggplot_nmds_obj <- reactiveVal(NULL)  # Use reactiveVal to store the NMDS ggplot object

  ### NMDS ###
  observeEvent(input$plotNMDS, {
    req(physeq(), input$distanceMetric, input$normalizationMethod, input$taxRank, input$groupingColumnBeta, input$shapeColumnBeta)

    print("NMDS button clicked")
    # Access the actual phyloseq object
    pseq <- physeq()

    # Filter out samples with NA or empty values in the grouping and shape columns
    valid_samples <- complete.cases(sample_data(pseq)[[input$groupingColumnBeta]],
                                    sample_data(pseq)[[input$shapeColumnBeta]])
    pseq <- prune_samples(valid_samples, pseq)

    # Transform the data and calculate distance
    pseq_transformed <- pseq %>%
      tax_transform(input$normalizationMethod, rank = input$taxRank)

    # Create a count table for NMDS
    count_table <- data.frame(t(otu_table(pseq_transformed)))

    # Perform NMDS ordination
    nmds <- metaMDS(count_table, distance = input$distanceMetric, binary = (input$distanceMetric == "jaccard"),
                    noshare = TRUE, autotransform = FALSE, k = 2, expand = TRUE,
                    maxit = 100, trymax = 100, wascores = TRUE)

    # Check if NMDS converged, stop if it didn't
    if (nmds$converged == FALSE) {
      showModal(modalDialog(
        title = "NMDS Warning",
        "NMDS did not converge. The results may not be reliable.",
        easyClose = TRUE,
        footer = modalButton("Close")
      ))
      return(NULL)
    }

    # Extract NMDS values
    nmds_values <- as.data.frame(scores(nmds, "sites"))
    nmds_values <- na.omit(nmds_values)
    colnames(nmds_values) <- c("NMDS1", "NMDS2")

    # Combine NMDS values with metadata
    meta <- as(sample_data(pseq), "data.frame")
    scores_df <- cbind(nmds_values, meta)

    # Ensure columns used in aesthetics are factors
    scores_df[[input$groupingColumnBeta]] <- as.factor(scores_df[[input$groupingColumnBeta]])
    scores_df[[input$shapeColumnBeta]] <- as.factor(scores_df[[input$shapeColumnBeta]])

    # Define axis labels for NMDS
    x_label <- "NMDS1"
    y_label <- "NMDS2"

    # Create ggplot with ellipses
    plot <- ggplot(data = scores_df, aes_string(x = "NMDS1", y = "NMDS2")) +
      theme_bw() +
      geom_point(aes(shape = input$shapeColumnBeta, color = input$groupingColumnBeta), alpha = 0.7, size = 2) +
      geom_vline(xintercept = c(0), color = "grey70", linetype = 1) +
      geom_hline(yintercept = c(0), color = "grey70", linetype = 1) +
      scale_color_manual(values = c("darkseagreen2","cornflowerblue", "indianred", "plum1", "skyblue", "yellow2")) +
      stat_ellipse(linetype = 1, segments = 10, lwd = 1.2, alpha = 0.25, level = 0.6,
                   geom = "polygon", aes(fill = input$groupingColumnBeta)) +
      theme(legend.title = element_text(size = 6),
            legend.text = element_text(size = 6),
            legend.key.size = unit(1, 'cm')) +
      coord_fixed(ratio = 1) +
      labs(title = paste("NMDS -", input$distanceMetric, "distance"), x = x_label, y = y_label)

    # Store the ggplot object in the reactive value
    ggplot_nmds_obj(plot)

    # Render the ggplot
    output$betaDiversityPlot <- renderPlot({
      print(ggplot_nmds_obj())
    })
  })

  # Download handler for NMDS plot
  output$download_NMDSPlot <- downloadHandler(
    filename = function() {
      paste0("NMDS_beta_diversity_", Sys.Date(), ".", input$beta_filetype)
    },
    content = function(file) {
      req(ggplot_nmds_obj())  # Ensure ggplot_nmds_obj contains a valid plot
      ggsave(file, plot = ggplot_nmds_obj(), width = input$plot_width, height = input$plot_height, units = "in", device = input$beta_filetype)
    }
  )



  ##for statistics
  observeEvent(input$Statistics, {
    req(physeq(), input$distanceMetric, input$normalizationMethod, input$taxRank, input$groupingColumnBeta)

    # Access the actual phyloseq object
    pseq <- physeq()
    variable_names <- input$groupingColumnBeta

    # Convert the phyloseq object to TreeSummarizedExperiment
    tse <- makeTreeSummarizedExperimentFromPhyloseq(pseq)
    tse <- agglomerateByRank(tse, input$taxRank)

    # Apply relative transform
    tse <- transformSamples(tse, method = "relabundance")

    # PERMANOVA Analysis
    set.seed(12346)
    assay_data <- t(assay(tse, "relabundance"))
    mod <- as.formula(paste("assay_data ~", paste(variable_names, collapse = "+")))

    permanova2 <- vegan::adonis2(mod,
                                 by = "margin",  # Each term analyzed individually
                                 data = colData(tse),
                                 method = input$distanceMetric,
                                 permutations = 999)  # More permutations for stability

    # Adjust p-values using Benjamini-Hochberg (FDR) method
    permanova2$`Pr(>F)` <- p.adjust(permanova2$`Pr(>F)`, method = "BH")

    output$permanovaTable <- DT::renderDataTable({
      DT::datatable(as.data.frame(permanova2), options = list(pageLength = 5))
    })

    output$download_permanovaTable <- downloadHandler(
      filename = function() {
        paste0("permanova_", input$distanceMetric, "_", Sys.Date(), ".csv")
      },
      content = function(file) {
        write.csv(as.data.frame(permanova2), file)
      }
    )

    # Dispersion Test Analysis
    dis <- vegan::vegdist(t(assays(tse)$counts), method = input$distanceMetric)
    b <- vegan::betadisper(dis, colData(tse)[[input$groupingColumnBeta]])
    dispersion_anova <- anova(b)

    output$betadisperTable <- DT::renderDataTable({
      DT::datatable(as.data.frame(dispersion_anova), options = list(pageLength = 5))
    })

    output$download_betadisperTable <- downloadHandler(
      filename = function() {
        paste0("betadisper_", input$distanceMetric, "_", Sys.Date(), ".csv")
      },
      content = function(file) {
        write.csv(as.data.frame(dispersion_anova), file)
      }
    )

    # Boxplot of distances to centroid
    p <- cbind(distance = as.numeric(b$distances),
               groupingColumnBeta = colData(tse)[[input$groupingColumnBeta]]) %>%
      as_tibble() %>%
      mutate(distance = as.numeric(distance)) %>%
      ggplot(aes(groupingColumnBeta, distance)) +
      geom_boxplot() +
      theme_light()

    output$dispersionBoxplot <- renderPlot({
      print(p)
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




                                            ############# community typing(DMM)


  generateDriverPlot <- function(best, physeq, k) {
    d <- melt(fitted(best))
    colnames(d) <- c("ASV", "cluster", "value")

    # Get the tax_table from the phyloseq object
    tax_df <- as.data.frame(tax_table(physeq))

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
    ggplot(d, aes(x = Taxa, y = value)) +
      theme_bw() +
      geom_bar(stat = "identity", fill = "darkblue") +
      coord_flip() +
      labs(title = paste("Top drivers: community type", k))
  }


  observe({
    req(physeq())
    pseq <- physeq()
    updateSelectInput(session, "taxRankDMM", choices = c(rank_names(pseq), "Genus", "Species"))
  })

  observeEvent(input$runDMM, {
    req(physeq(), input$taxRankDMM, input$detectionThreshold, input$prevalenceThreshold)

    # Data Transformation and Filtering
    pseq.comp <- microbiome::transform(physeq(), "compositional")
    taxa <- core_members(pseq.comp, detection = input$detectionThreshold / 100, prevalence = input$prevalenceThreshold / 100)
    pseq <- prune_taxa(taxa, physeq())

    # Aggregate at selected taxonomic rank
    if (input$taxRankDMM %in% c("Genus", "Species")) {
      pseq <- tax_glom(pseq, taxrank = input$taxRankDMM)
    } else {
      pseq <- tax_glom(pseq, taxrank = input$taxRankDMM)
    }

    dat <- abundances(pseq)
    count <- as.matrix(t(dat))
    count <- count[rowSums(count) > 0, ]

    if (nrow(count) > 0) {
      # Model Fitting
      fit <- lapply(1:input$numComponents, dmn, count = count, verbose = TRUE)
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
        mixturewt(best)
      })

      # Render Sample Assignments
      sample_assignments <- apply(mixture(best), 1, which.max)
      output$sampleAssignments <- renderPrint({
        sample_assignments
      })

      # Filter metadata to only include samples in the count matrix
      valid_samples <- rownames(count)
      metadata <- as(sample_data(physeq()), "data.frame")
      metadata_filtered <- metadata[valid_samples, ]

      # Update metadata with sample assignments
      metadata_filtered$DMM_Cluster <- factor(sample_assignments)

      # Re-create phyloseq object with filtered metadata
      updated_pseq <- phyloseq(otu_table(physeq()), tax_table(physeq()), sample_data(metadata_filtered), phy_tree(physeq()))

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
          generateDriverPlot(best, physeq(), k)
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
            print(generateDriverPlot(best, physeq(), k))
          })
          dev.off()
        }
      )

    } else {
      # Handle cases where no valid samples are found
      output$modelFitPlot <- renderPlot({ plot.new(); text(0.5, 0.5, "No valid samples with non-zero counts.") })
      output$mixtureParams <- renderPrint({"No valid samples with non-zero counts."})
      output$sampleAssignments <- renderDT({
        datatable(data.frame(Sample = character(0), Cluster = character(0)), options = list(pageLength = 10, autoWidth = TRUE))
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
    metadata <- sample_data(pseq)
    imputed_metadata <- metadata
    for (col in colnames(metadata)) {
      if (is.numeric(metadata[[col]])) {
        imputed_metadata[[col]][is.na(imputed_metadata[[col]])] <- median(imputed_metadata[[col]], na.rm = TRUE)
      }
    }
    sample_data(pseq) <- imputed_metadata
    return(pseq)
  }

  # Create a reactive value to store the pseq object
  reactive_pseq <- reactiveVal()

  # Impute missing values and update select inputs
  observe({
    req(physeq())
    pseq <- physeq()

    # Impute missing values before analysis
    pseq <- impute_median(pseq)

    # Store the imputed pseq in the reactive value
    reactive_pseq(pseq)

    # Update taxRankRDA choices based on the available taxonomic ranks
    if (!is.null(tax_table(reactive_pseq()))) {
      updateSelectInput(session, "taxRankRDA", choices = c("ASV", rank_names(reactive_pseq())))
    } else {
      updateSelectInput(session, "taxRankRDA", choices = "ASV")
      showNotification("Tax table is missing in the phyloseq object, using ASV level.", type = "warning")
    }

    # Update other inputs
    sample_data_cols <- colnames(sample_data(pseq))
    updateSelectInput(session, "constraints", choices = sample_data_cols)
    updateSelectInput(session, "conditions", choices = c("None", sample_data_cols))
    updateSelectInput(session, "colorBy", choices = sample_data_cols)
    updateSelectInput(session, "shapeBy", choices = sample_data_cols)
  })

  # Handle the RDA plot generation
  observeEvent(input$plotRDA, {
    req(reactive_pseq(), input$taxRankRDA, input$normMethod, input$colorBy, input$shapeBy, input$constraints, input$conditions)

    # Retrieve pseq from the reactive object
    pseq <- reactive_pseq()

    if (input$normMethod != "identity" & input$taxRankRDA != "ASV") {
      pseq <- aggregate_taxa(pseq, input$taxRankRDA)
    }

    # Handle constraints and conditions correctly for multiple selections
    constraints <- if (is.null(input$constraints) || length(input$constraints) == 0 || "None" %in% input$constraints) {
      NULL
    } else {
      input$constraints
    }

    conditions <- if (is.null(input$conditions) || length(input$conditions) == 0 || "None" %in% input$conditions) {
      NULL
    } else {
      input$conditions
    }

    # Convert the phyloseq object to TreeSummarizedExperiment
    tse <- makeTreeSummarizedExperimentFromPhyloseq(pseq)

    # Add a pseudocount to avoid zeros (this should be done before any log-ratio transformation)
    assay(tse, "counts") <- assay(tse, "counts") + 1  # Add 1 as a pseudocount to avoid zeros

    # Apply the normalization/transformation method
    tse <- transformAssay(
      tse,
      assay.type = "counts",  # Using the original counts assay
      method = input$normMethod,
      MARGIN = "samples"
    )

    # After transformation, you can proceed with RDA or other analyses
    transformed_assay_name <- switch(input$normMethod,
                                     "log1p" = "logcounts",
                                     "counts" = "counts",
                                     "relabundance"="relabundance",
                                     "hellinger"="hellinger",
                                     "standardize"="standardize",
                                     "normalize"="normalize",
                                     "total"="total",
                                     "z"="z")

    # Prepare the formula for RDA based on selected grouping columns
    rda_formula <- as.formula(paste("assay ~", paste(input$constraints, collapse = " + ")))

    # Run RDA on the transformed data
    tse <- runRDA(
      tse,
      assay.type = transformed_assay_name,  # Use the correct transformed assay name
      formula = rda_formula,
      distance = input$distanceMetric,
      na.action = na.exclude,
      normalize = FALSE
    )

    # Store PERMANOVA results for later use
    rda_info <- attr(reducedDim(tse, "RDA"), "significance")

    # Create the ggplot object
    ggplot_obj <- plotRDA(tse, "RDA", colour.by = input$colorBy) +
      stat_ellipse(aes(color = colData(tse)[[input$colorBy]]), level = 0.95, alpha = 0.3) +
      geom_point(aes(shape = colData(tse)[[input$shapeBy]]), size = 3) +
      theme_minimal() +
      labs(color = input$colorBy, shape = input$shapeBy, title = "RDA Plot")

    # Convert ggplot to plotly
    plotly_obj <- ggplotly(ggplot_obj, tooltip = c("x", "y", input$colorBy, input$shapeBy)) %>%
      layout(dragmode = "zoom") %>%
      config(displayModeBar = TRUE, scrollZoom = TRUE)

    # Store the ggplot and plotly objects in reactive values
    ggplot_plots_rda(ggplot_obj)
    plotly_plots_rda(plotly_obj)

    # Render the plotly plot
    output$rdaPlot <- renderPlotly({
      plotly_plots_rda()
    })

    # Render the rda_info table as an HTML table using DT
    output$rdaInfoTable <- DT::renderDataTable({
      rounded_permanova <- rda_info$permanova %>%
        mutate(across(everything(), round, digits = 4))

      DT::datatable(rounded_permanova, rownames = TRUE, options = list(pageLength = 5))
    })

    # Download handler for the PERMANOVA table
    output$download_rdaInfoTable <- downloadHandler(
      filename = function() {
        paste0("permanova_statistics_", Sys.Date(), ".csv")
      },
      content = function(file) {
        write.csv(rounded_permanova(), file)
      }
    )
  })





  # #### dbRDA Plot and Analysis
  #
  # observeEvent(input$plotDbRDA, {
  #   req(reactive_pseq(), input$normMethod, input$distanceMetric, input$constraints, input$conditions, input$taxRankRDA)
  #
  #   # Retrieve pseq from the reactive object
  #   pseq <- reactive_pseq()
  #   if (input$normMethod != "identity" & input$taxRankRDA != "ASV") {
  #     pseq <- aggregate_taxa(pseq, input$taxRankRDA)
  #   }
  #
  #   # Extract the OTU table (assume it's already normalized)
  #   otu_table <- as.data.frame(t(otu_table(pseq)))
  #   metadata <- as.data.frame(sample_data(pseq))
  #
  #   # Ensure that the variables used in the formula are in the metadata
  #   constraints <- if (!is.null(input$constraints) && length(input$constraints) > 0 && !"None" %in% input$constraints) {
  #     input$constraints
  #   } else {
  #     stop("Please select valid constraints.")
  #   }
  #
  #   conditions <- if (!is.null(input$conditions) && length(input$conditions) > 0 && !"None" %in% input$conditions) {
  #     input$conditions
  #   } else {
  #     NULL
  #   }
  #
  #   # Construct the formula for dbRDA
  #   formula_string <- paste("otu_table ~", paste(constraints, collapse = " + "))
  #   if (!is.null(conditions)) {
  #     formula_string <- paste(formula_string, "+ Condition(", paste(conditions, collapse = " + "), ")")
  #   }
  #   dbRDA_formula <- as.formula(formula_string)
  #
  #   # Perform the dbRDA
  #   dbrda_model <- dbrda(
  #     formula = dbRDA_formula,
  #     data = metadata,
  #     dist = input$distanceMetric,
  #     sqrt.dist = TRUE
  #   )
  #
  #   # Calculate species length (splen)
  #   splen <- scores(dbrda_model, display = "species", scaling = 2)$species %>%
  #     apply(1, function(x) sqrt(sum(x^2)))
  #
  #   # Plot the dbRDA with splen
  #   output$dbRdaPlot <- renderPlot({
  #     plot(dbrda_model, scaling = 2, main = "Triplot dbRDA", type = "none",
  #          xlab = "dbRDA1", ylab = "dbRDA2", xlim = c(-10, 20), ylim = c(-10, 10))
  #     points(dbrda_model, display = "sites", pch = 21, col = "black", cex = 1.2)
  #     text(dbrda_model, display = "species", select = splen > 0.5, arrow = TRUE, length = 0.05)
  #   })
  #
  #   # Store species scores in a reactive value for downloading
  #   species_scores <- data.frame(scores(dbrda_model, display = "species"))
  #   output$speciesScoresTable <- DT::renderDataTable({
  #     DT::datatable(species_scores, options = list(pageLength = 5))
  #   })
  #
  #   # Download handler for species scores
  #   output$download_speciesScoresTable <- downloadHandler(
  #     filename = function() {
  #       paste0("species_scores_", Sys.Date(), ".csv")
  #     },
  #     content = function(file) {
  #       write.csv(species_scores, file)
  #     }
  #   )
  #
  #   # Perform statistical tests and display in DT
  #   model_significance <- anova(dbrda_model)
  #   terms_significance <- anova(dbrda_model, by = "term", permutations = 199)
  #
  #   output$dbRdaInfoTable <- DT::renderDataTable({
  #     DT::datatable(
  #       as.data.frame(model_significance),
  #       options = list(pageLength = 5)
  #     )
  #   })
  #
  #   output$termsSignificanceTable <- DT::renderDataTable({
  #     DT::datatable(
  #       as.data.frame(terms_significance),
  #       options = list(pageLength = 5)
  #     )
  #   })
  #
  #   # Download handler for the statistical test results
  #   output$download_dbRdaInfoTable <- downloadHandler(
  #     filename = function() {
  #       paste0("dbRDA_statistics_", Sys.Date(), ".csv")
  #     },
  #     content = function(file) {
  #       write.csv(as.data.frame(model_significance), file)
  #     }
  #   )
  #
  #   output$download_termsSignificanceTable <- downloadHandler(
  #     filename = function() {
  #       paste0("dbRDA_terms_significance_", Sys.Date(), ".csv")
  #     },
  #     content = function(file) {
  #       write.csv(as.data.frame(terms_significance), file)
  #     }
  #   )
  # })

  # Download handler for the RDA plot
  output$download_rdaPlot <- downloadHandler(
    filename = function() {
      paste0("rda_plot_", Sys.Date(), ".", input$rda_filetype)
    },
    content = function(file) {
      if (input$rda_filetype == "html") {
        tempfile <- tempfile(fileext = ".html")
        htmlwidgets::saveWidget(as_widget(plotly_plots_rda()), tempfile)
        file.copy(tempfile, file)
      } else {
        ggsave(file, ggplot_plots_rda(), width = 16, height = 12, units = "in", device = input$rda_filetype)
      }
    }
  )

  # Download the updated phyloseq object
  output$downloadPhyseq <- downloadHandler(
    filename = function() {
      paste("imputed_phyloseq_object", Sys.Date(), ".rds", sep = "")
    },
    content = function(file) {
      saveRDS(reactive_pseq(), file)
    }
  )

  # Add the download button in the UI
  output$downloadPhyseqUI <- renderUI({
    downloadButton("downloadPhyseq", "Download Updated Phyloseq Object")
  })

  observeEvent(input$show_rda_explanation, {
    showModal(modalDialog(
      title = "Redundancy Analysis (RDA) Explanation",
      HTML("<p><strong>Redundancy Analysis (RDA)</strong> is a constrained ordination method used to understand how much variation in the data can be explained by specific environmental or experimental variables.</p>
         <p><strong>Why RDA?</strong></p>
         <ul>
           <li><strong>Identify Drivers:</strong> RDA helps identify which variables are the main drivers behind the differences in microbial communities across samples.</li>
           <li><strong>Visual Interpretation:</strong> The RDA plot provides a visual interpretation of the relationships between samples, taxa, and the chosen constraints/conditions.</li>
           <li><strong>Statistical Significance:</strong> It allows for testing the statistical significance of the relationships between the constraints and the community composition.</li>
         </ul>
         <p>In this app, RDA is used to explore how different metadata columns (constraints) can explain the variation in the microbiome data across samples.</p>"),
      easyClose = TRUE,
      footer = modalButton("Close")
    ))
  })



#   ############ RDA and dbRDA plot
#
#   # Initialize reactive values to store plots and data
#   plotly_plots_rda <- reactiveVal(NULL)
#   ggplot_plots_rda <- reactiveVal(NULL)
#   rounded_permanova <- reactiveVal(NULL)
#   species_scores <- reactiveVal(NULL)
#
#   # Function to impute missing values with the median
#   impute_median <- function(pseq) {
#     metadata <- sample_data(pseq)
#     imputed_metadata <- metadata
#     for (col in colnames(metadata)) {
#       if (is.numeric(metadata[[col]])) {
#         imputed_metadata[[col]][is.na(imputed_metadata[[col]])] <- median(imputed_metadata[[col]], na.rm = TRUE)
#       }
#     }
#     sample_data(pseq) <- imputed_metadata
#     return(pseq)
#   }
#
#   # Create a reactive value to store the pseq object
#   reactive_pseq <- reactiveVal()
#
#   # Impute missing values and update select inputs
#   observe({
#     req(physeq())
#     pseq <- physeq()
#
#     # Impute missing values before analysis
#     pseq <- impute_median(pseq)
#
#     # Store the imputed pseq in the reactive value
#     reactive_pseq(pseq)
#
#     # Update taxRankRDA choices based on the available taxonomic ranks
#     if (!is.null(tax_table(reactive_pseq()))) {
#       updateSelectInput(session, "taxRankRDA", choices = c("ASV", rank_names(reactive_pseq())))
#     } else {
#       updateSelectInput(session, "taxRankRDA", choices = "ASV")
#       showNotification("Tax table is missing in the phyloseq object, using ASV level.", type = "warning")
#     }
#
#     # Update other inputs
#     sample_data_cols <- colnames(sample_data(pseq))
#     updateSelectInput(session, "constraints", choices = sample_data_cols)
#     updateSelectInput(session, "conditions", choices = c("None", sample_data_cols))
#     updateSelectInput(session, "colorBy", choices = sample_data_cols)
#     updateSelectInput(session, "shapeBy", choices = sample_data_cols)
#   })
#
#   # Handle the RDA plot generation
#   observeEvent(input$plotRDA, {
#     req(reactive_pseq(), input$taxRankRDA, input$normMethod, input$colorBy, input$shapeBy, input$constraints, input$conditions)
#
#     # Retrieve pseq from the reactive object
#     pseq <- reactive_pseq()
#
#     if (input$normMethod != "identity" & input$taxRankRDA != "ASV") {
#       pseq <- aggregate_taxa(pseq, input$taxRankRDA)
#     }
#
#     # Handle constraints and conditions correctly for multiple selections
#     constraints <- if (is.null(input$constraints) || length(input$constraints) == 0 || "None" %in% input$constraints) {
#       NULL
#     } else {
#       input$constraints
#     }
#
#     conditions <- if (is.null(input$conditions) || length(input$conditions) == 0 || "None" %in% input$conditions) {
#       NULL
#     } else {
#       input$conditions
#     }
#
#     # Convert the phyloseq object to TreeSummarizedExperiment
#     tse <- makeTreeSummarizedExperimentFromPhyloseq(pseq)
#
#     # Add a pseudocount to avoid zeros (this should be done before any log-ratio transformation)
#     assay(tse, "counts") <- assay(tse, "counts") + 1  # Add 1 as a pseudocount to avoid zeros
#
#     # Apply the normalization/transformation method
#     tse <- transformAssay(
#       tse,
#       assay.type = "counts",  # Using the original counts assay
#       method = input$normMethod,
#       MARGIN = "samples"
#     )
#
#     # After transformation, you can proceed with RDA or other analyses
#     transformed_assay_name <- switch(input$normMethod,
#                                      "log1p" = "logcounts",
#                                      "counts" = "counts",
#                                      "relabundance"="relabundance",
#                                      "hellinger"="hellinger",
#                                      "standardize"="standardize",
#                                      "normalize"="normalize",
#                                      "total"="total",
#                                      "z"="z")
#
#     # Prepare the formula for RDA based on selected grouping columns
#     rda_formula <- as.formula(paste("assay ~", paste(input$constraints, collapse = " + ")))
#
#     # Run RDA on the transformed data
#     tse <- runRDA(
#       tse,
#       assay.type = transformed_assay_name,  # Use the correct transformed assay name
#       formula = rda_formula,
#       distance = input$distanceMetric,
#       na.action = na.exclude,
#       normalize = FALSE
#     )
#
#     # Store PERMANOVA results for later use
#     rda_info <- attr(reducedDim(tse, "RDA"), "significance")
#
#     # Create the ggplot object
#     ggplot_obj <- plotRDA(tse, "RDA", colour.by = input$colorBy) +
#       stat_ellipse(aes(color = colData(tse)[[input$colorBy]]), level = 0.95, alpha = 0.3) +
#       geom_point(aes(shape = colData(tse)[[input$shapeBy]]), size = 3) +
#       theme_minimal() +
#       labs(color = input$colorBy, shape = input$shapeBy, title = "RDA Plot")
#
#     # Convert ggplot to plotly
#     plotly_obj <- ggplotly(ggplot_obj, tooltip = c("x", "y", input$colorBy, input$shapeBy)) %>%
#       layout(dragmode = "zoom") %>%
#       config(displayModeBar = TRUE, scrollZoom = TRUE)
#
#     # Store the ggplot and plotly objects in reactive values
#     ggplot_plots_rda(ggplot_obj)
#     plotly_plots_rda(plotly_obj)
#
#     # Render the plotly plot
#     output$rdaPlot <- renderPlotly({
#       plotly_plots_rda()
#     })
#
#     # Render the rda_info table as an HTML table using DT
#     output$rdaInfoTable <- DT::renderDataTable({
#       rounded_permanova <- rda_info$permanova %>%
#         mutate(across(everything(), round, digits = 4))
#
#       DT::datatable(rounded_permanova, rownames = TRUE, options = list(pageLength = 5))
#     })
#
#     # Download handler for the PERMANOVA table
#     output$download_rdaInfoTable <- downloadHandler(
#       filename = function() {
#         paste0("permanova_statistics_", Sys.Date(), ".csv")
#       },
#       content = function(file) {
#         write.csv(rounded_permanova(), file)
#       }
#     )
#   })
#
#
# ####dbrda
#
#   observeEvent(input$plotDbRDA, {
#     req(reactive_pseq(), input$normMethod, input$distanceMetric, input$constraints, input$conditions, input$taxRankRDA)
#
#     # Retrieve pseq from the reactive object
#     pseq <- reactive_pseq()
#     if (input$normMethod != "identity" & input$taxRankRDA != "ASV") {
#       pseq <- aggregate_taxa(pseq, input$taxRankRDA)
#     }
#
#     # Extract the OTU table (assume it's already normalized)
#     otu_table <- as.matrix(otu_table(pseq))
#     metadata <- as.data.frame(sample_data(pseq))
#
#     # Ensure that the variables used in the formula are in the metadata
#     constraints <- if (!is.null(input$constraints) && length(input$constraints) > 0 && !"None" %in% input$constraints) {
#       input$constraints
#     } else {
#       stop("Please select valid constraints.")
#     }
#
#     conditions <- if (!is.null(input$conditions) && length(input$conditions) > 0 && !"None" %in% input$conditions) {
#       input$conditions
#     } else {
#       NULL
#     }
#
#     # Construct the formula for dbRDA
#     formula_string <- paste("otu_table ~", paste(constraints, collapse = " + "))
#     if (!is.null(conditions)) {
#       formula_string <- paste(formula_string, "+ Condition(", paste(conditions, collapse = " + "), ")")
#     }
#     dbRDA_formula <- as.formula(formula_string)
#
#     # Perform the dbRDA
#     dbrda_model <- dbrda(
#       formula = dbRDA_formula,
#       data = metadata,
#       dist = input$distanceMetric,
#       sqrt.dist = TRUE
#     )
#
#     # Calculate species length (splen)
#     splen <- scores(dbrda_model, display = "species", scaling = 2)$species %>%
#       apply(1, function(x) sqrt(sum(x^2)))
#
#     # Plot the dbRDA with splen
#     output$rdaPlot <- renderPlot({
#       plot(dbrda_model, scaling = 2, main = "Triplot dbRDA", type = "none",
#            xlab = "dbRDA1", ylab = "dbRDA2", xlim = c(-10, 20), ylim = c(-10, 10))
#       points(dbrda_model, display = "sites", pch = 21, col = "black", cex = 1.2)
#       text(dbrda_model, display = "species", select = splen > 0.5, arrow = TRUE, length = 0.05)
#     })
#
#     # Store species scores in a reactive value for downloading
#     species_scores <- data.frame(scores(dbrda_model, display = "species"))
#     output$speciesScoresTable <- DT::renderDataTable({
#       DT::datatable(species_scores, options = list(pageLength = 5))
#     })
#
#     # Download handler for species scores
#     output$download_speciesScoresTable <- downloadHandler(
#       filename = function() {
#         paste0("species_scores_", Sys.Date(), ".csv")
#       },
#       content = function(file) {
#         write.csv(species_scores, file)
#       }
#     )
#
#     # Perform statistical tests and display in DT
#     model_significance <- anova(dbrda_model)
#     terms_significance <- anova(dbrda_model, by = "term", permutations = 199)
#
#     output$rdaInfoTable <- DT::renderDataTable({
#       DT::datatable(
#         as.data.frame(model_significance),
#         options = list(pageLength = 5)
#       )
#     })
#
#     output$termsSignificanceTable <- DT::renderDataTable({
#       DT::datatable(
#         as.data.frame(terms_significance),
#         options = list(pageLength = 5)
#       )
#     })
#
#     # Download handler for the statistical test results
#     output$download_rdaInfoTable <- downloadHandler(
#       filename = function() {
#         paste0("dbRDA_statistics_", Sys.Date(), ".csv")
#       },
#       content = function(file) {
#         write.csv(as.data.frame(model_significance), file)
#       }
#     )
#
#     output$download_termsSignificanceTable <- downloadHandler(
#       filename = function() {
#         paste0("dbRDA_terms_significance_", Sys.Date(), ".csv")
#       },
#       content = function(file) {
#         write.csv(as.data.frame(terms_significance), file)
#       }
#     )
#   })
#
#   # Download handler for the RDA plot
#   output$download_rdaPlot <- downloadHandler(
#     filename = function() {
#       paste0("rda_plot_", Sys.Date(), ".", input$rda_filetype)
#     },
#     content = function(file) {
#       if (input$rda_filetype == "html") {
#         tempfile <- tempfile(fileext = ".html")
#         htmlwidgets::saveWidget(as_widget(plotly_plots_rda()), tempfile)
#         file.copy(tempfile, file)
#       } else {
#         ggsave(file, ggplot_plots_rda(), width = 16, height = 12, units = "in", device = input$rda_filetype)
#       }
#     }
#   )
#
#   # Download the updated phyloseq object
#   output$downloadPhyseq <- downloadHandler(
#     filename = function() {
#       paste("imputed_phyloseq_object", Sys.Date(), ".rds", sep = "")
#     },
#     content = function(file) {
#       saveRDS(reactive_pseq(), file)
#     }
#   )
#
#   # Add the download button in the UI
#   output$downloadPhyseqUI <- renderUI({
#     downloadButton("downloadPhyseq", "Download Updated Phyloseq Object")
#   })
#
#   observeEvent(input$show_rda_explanation, {
#     showModal(modalDialog(
#       title = "Redundancy Analysis (RDA) Explanation",
#       HTML("<p><strong>Redundancy Analysis (RDA)</strong> is a constrained ordination method used to understand how much variation in the data can be explained by specific environmental or experimental variables.</p>
#          <p><strong>Why RDA?</strong></p>
#          <ul>
#            <li><strong>Identify Drivers:</strong> RDA helps identify which variables are the main drivers behind the differences in microbial communities across samples.</li>
#            <li><strong>Visual Interpretation:</strong> The RDA plot provides a visual interpretation of the relationships between samples, taxa, and the chosen constraints/conditions.</li>
#            <li><strong>Statistical Significance:</strong> It allows for testing the statistical significance of the relationships between the constraints and the community composition.</li>
#          </ul>
#          <p>In this app, RDA is used to explore how different metadata columns (constraints) can explain the variation in the microbiome data across samples.</p>"),
#       easyClose = TRUE,
#       footer = modalButton("Close")
#     ))
#   })


  #works
  # ############ RDA plot
  #
  # plotly_plots_rda <- reactiveVal(NULL)
  # ggplot_plots_rda <- reactiveVal(NULL)
  # rounded_permanova <- reactiveVal(NULL)
  #
  #
  # # Function to impute missing values with the median
  # impute_median <- function(pseq) {
  #   metadata <- sample_data(pseq)
  #   imputed_metadata <- metadata
  #   for (col in colnames(metadata)) {
  #     if (is.numeric(metadata[[col]])) {
  #       imputed_metadata[[col]][is.na(imputed_metadata[[col]])] <- median(imputed_metadata[[col]], na.rm = TRUE)
  #     }
  #   }
  #   sample_data(pseq) <- imputed_metadata
  #   return(pseq)
  # }
  #
  # # Create a reactive value to store the pseq object
  # reactive_pseq <- reactiveVal()
  #
  # # Impute missing values and update select inputs
  # observe({
  #   req(physeq())
  #   pseq <- physeq()
  #
  #   # Impute missing values before analysis
  #   pseq <- impute_median(pseq)
  #
  #   # Store the imputed pseq in the reactive value
  #   reactive_pseq(pseq)
  #
  #   # Update taxRankRDA choices based on the available taxonomic ranks
  #   if (!is.null(tax_table(reactive_pseq()))) {
  #     updateSelectInput(session, "taxRankRDA", choices = c("ASV", rank_names(reactive_pseq())))
  #   } else {
  #     updateSelectInput(session, "taxRankRDA", choices = "ASV")
  #     showNotification("Tax table is missing in the phyloseq object, using ASV level.", type = "warning")
  #   }
  #
  #   # Update other inputs
  #   sample_data_cols <- colnames(sample_data(pseq))
  #   updateSelectInput(session, "constraints", choices = sample_data_cols)
  #   updateSelectInput(session, "conditions", choices = c("None", sample_data_cols))
  #   updateSelectInput(session, "colorBy", choices = sample_data_cols)
  #   updateSelectInput(session, "shapeBy", choices = sample_data_cols)
  # })
  #
  # # Handle the RDA plot generation
  # observeEvent(input$plotRDA, {
  #   req(reactive_pseq(), input$taxRankRDA, input$normMethod, input$colorBy, input$shapeBy, input$constraints, input$conditions)
  #
  #   # Retrieve pseq from the reactive object
  #   pseq <- reactive_pseq()
  #
  #
  #   if (input$normMethod != "identity" & input$taxRankRDA != "ASV") {
  #     pseq <- aggregate_taxa(pseq, input$taxRankRDA)
  #   }
  #
  #   # Handle constraints and conditions correctly for multiple selections
  #   constraints <- if (is.null(input$constraints) || length(input$constraints) == 0 || "None" %in% input$constraints) {
  #     NULL
  #   } else {
  #     input$constraints
  #   }
  #
  #   conditions <- if (is.null(input$conditions) || length(input$conditions) == 0 || "None" %in% input$conditions) {
  #     NULL
  #   } else {
  #     input$conditions
  #   }
  #
  #   # Convert the phyloseq object to TreeSummarizedExperiment
  #   tse <- makeTreeSummarizedExperimentFromPhyloseq(pseq)
  #
  #   # Add a pseudocount to avoid zeros (this should be done before any log-ratio transformation)
  #   assay(tse, "counts") <- assay(tse, "counts") + 1  # Add 1 as a pseudocount to avoid zeros
  #
  #   # Check for negative values (though unlikely in counts data, we check just in case)
  #   if (any(assay(tse, "counts") <= 0)) {
  #     stop("The count data contains negative values, which is not expected. Please check your data.")
  #   }
  #
  #   # Apply the normalization/transformation method
  #   tse <- transformAssay(
  #     tse,
  #     assay.type = "counts",  # Using the original counts assay
  #     method = input$normMethod,
  #     MARGIN = "samples"
  #   )
  #
  #   # After transformation, you can proceed with RDA or other analyses
  #   transformed_assay_name <- switch(input$normMethod,
  #                                    "log1p" = "logcounts",
  #                                    "counts" = "counts",
  #                                    "relabundance"="relabundance",
  #                                    "hellinger"="hellinger",
  #                                    "standardize"="standardize",
  #                                    "normalize"="normalize",
  #                                    "total"="total",
  #                                    "z"="z")
  #
  #
  #   # Prepare the formula for RDA based on selected grouping columns
  #   rda_formula <- as.formula(paste("assay ~", paste(input$constraints, collapse = " + ")))
  #
  #   # Run RDA on the transformed data
  #   tse <- runRDA(
  #     tse,
  #     assay.type = transformed_assay_name,  # Use the correct transformed assay name
  #     formula = rda_formula,
  #     distance = input$distanceMetric,
  #     na.action = na.exclude,
  #     normalize = FALSE
  #   )
  #
  #
  #   # Store PERMANOVA results for later use
  #   rda_info <- attr(reducedDim(tse, "RDA"), "significance")
  #
  #   # Create the ggplot object
  #   ggplot_obj <- plotRDA(tse, "RDA", colour.by = input$colorBy) +
  #     stat_ellipse(aes(color = colData(tse)[[input$colorBy]]), level = 0.95, alpha = 0.3) +
  #     geom_point(aes(shape = colData(tse)[[input$shapeBy]]), size = 3) +
  #     theme_minimal() +
  #     labs(color = input$colorBy, shape = input$shapeBy, title = "RDA Plot")
  #
  #   # Convert ggplot to plotly
  #   plotly_obj <- ggplotly(ggplot_obj, tooltip = c("x", "y", input$colorBy, input$shapeBy)) %>%
  #     layout(dragmode = "zoom") %>%
  #     config(displayModeBar = TRUE, scrollZoom = TRUE)
  #
  #   # Store the ggplot and plotly objects in reactive values
  #   ggplot_plots_rda(ggplot_obj)
  #   plotly_plots_rda(plotly_obj)
  #
  #   # Render the plotly plot
  #   output$rdaPlot <- renderPlotly({
  #     plotly_plots_rda()
  #   })
  #
  #   # Render the rda_info table as an HTML table using DT
  #   output$rdaInfoTable <- DT::renderDataTable({
  #     rounded_permanova <- rda_info$permanova %>%
  #       mutate(across(everything(), round, digits = 4))
  #
  #     DT::datatable(rounded_permanova, rownames = TRUE, options = list(pageLength = 5))
  #   })
  #
  #   # Download handler for the PERMANOVA table
  #   output$download_rdaInfoTable <- downloadHandler(
  #     filename = function() {
  #       paste0("permanova_statistics_", Sys.Date(), ".csv")
  #     },
  #     content = function(file) {
  #       write.csv(rounded_permanova(), file)
  #     }
  #   )
  # })
  #
  # # Download handler for the RDA plot
  # output$download_rdaPlot <- downloadHandler(
  #   filename = function() {
  #     paste0("rda_plot_", Sys.Date(), ".", input$rda_filetype)
  #   },
  #   content = function(file) {
  #     if (input$rda_filetype == "html") {
  #       tempfile <- tempfile(fileext = ".html")
  #       htmlwidgets::saveWidget(as_widget(plotly_plots_rda()), tempfile)
  #       file.copy(tempfile, file)
  #     } else {
  #       ggsave(file, ggplot_plots_rda(), width = 16, height = 12, units = "in", device = input$rda_filetype)
  #     }
  #   }
  # )
  #
  # # Download the updated phyloseq object
  # output$downloadPhyseq <- downloadHandler(
  #   filename = function() {
  #     paste("imputed_phyloseq_object", Sys.Date(), ".rds", sep = "")
  #   },
  #   content = function(file) {
  #     saveRDS(reactive_pseq(), file)
  #   }
  # )
  #
  # # Add the download button in the UI
  # output$downloadPhyseqUI <- renderUI({
  #   downloadButton("downloadPhyseq", "Download Updated Phyloseq Object")
  # })
  #
  # observeEvent(input$show_rda_explanation, {
  #   showModal(modalDialog(
  #     title = "Redundancy Analysis (RDA) Explanation",
  #     HTML("<p><strong>Redundancy Analysis (RDA)</strong> is a constrained ordination method used to understand how much variation in the data can be explained by specific environmental or experimental variables.</p>
  #        <p><strong>Why RDA?</strong></p>
  #        <ul>
  #          <li><strong>Identify Drivers:</strong> RDA helps identify which variables are the main drivers behind the differences in microbial communities across samples.</li>
  #          <li><strong>Visual Interpretation:</strong> The RDA plot provides a visual interpretation of the relationships between samples, taxa, and the chosen constraints/conditions.</li>
  #          <li><strong>Statistical Significance:</strong> It allows for testing the statistical significance of the relationships between the constraints and the community composition.</li>
  #        </ul>
  #        <p>In this app, RDA is used to explore how different metadata columns (constraints) can explain the variation in the microbiome data across samples.</p>"),
  #     easyClose = TRUE,
  #     footer = modalButton("Close")
  #   ))
  # })



  #   # Update choices based on the available taxonomic ranks and sample data columns
  #   updateSelectInput(session, "taxRankRDA", choices = rank_names(pseq))
  #   sample_data_cols <- colnames(sample_data(pseq))
  #   updateSelectInput(session, "constraints", choices = sample_data_cols)
  #   updateSelectInput(session, "conditions", choices = c("None", sample_data_cols))
  #   updateSelectInput(session, "colorBy", choices = sample_data_cols)
  #   updateSelectInput(session, "shapeBy", choices = sample_data_cols)
  # })
  #
  # # Handle the RDA plot generation
  # observeEvent(input$plotRDA, {
  #   req(reactive_pseq(), input$taxRankRDA, input$normMethod, input$colorBy, input$shapeBy, input$constraints, input$conditions)
  #
  #   # Retrieve pseq from the reactive object
  #   pseq <- reactive_pseq()
  #
  #   pseq <- phyloseq_validate(pseq, remove_undetected = TRUE)
  #   pseq <- tax_fix(pseq)
  #   pseq <- tax_agg(ps = pseq, input$taxRankRDA)
  #
  #   if (input$normMethod != "identity") {
  #     pseq <- tax_transform(pseq, trans = input$normMethod, rank = input$taxRankRDA)
  #   }
  #
  #   # Handle constraints and conditions correctly for multiple selections
  #   constraints <- if (is.null(input$constraints) || length(input$constraints) == 0 || "None" %in% input$constraints) {
  #     NULL
  #   } else {
  #     input$constraints
  #   }
  #
  #   conditions <- if (is.null(input$conditions) || length(input$conditions) == 0 || "None" %in% input$conditions) {
  #     NULL
  #   } else {
  #     input$conditions
  #   }
  #
  #   # RDA calculation
  #   ordination <- pseq %>%
  #     ord_calc(
  #       method = "RDA",
  #       constraints = constraints,
  #       conditions = conditions
  #     )
  #
  #
  #
  #   output$rdaPlot <- renderPlot({
  #     tryCatch({
  #       # Ensure correct data mapping and factor levels
  #       pseq_data <- sample_data(pseq)
  #       pseq_data[[input$colorBy]] <- factor(pseq_data[[input$colorBy]])
  #       pseq_data[[input$shapeBy]] <- factor(pseq_data[[input$shapeBy]])
  #
  #       # Create the ggplot object with ord_plot
  #       plot <- ord_plot(
  #         ordination,
  #         data = pseq_data,  # Explicitly passing the data
  #         colour = input$colorBy,
  #         size = 2,
  #         alpha = 0.5,
  #         shape = input$shapeBy,
  #         interactive = FALSE,  # Ensuring the plot is static
  #         plot_taxa = 1:8,
  #         tax_vec_length = 4.5,
  #         tax_lab_length = 4.6,
  #         tax_lab_style = tax_lab_style(
  #           type = "text", max_angle = 90, fontface = "bold.italic"
  #         ),
  #         constraint_vec_style = vec_constraint(1.5, alpha = 0.5),
  #         constraint_vec_length = 3, constraint_lab_length = 3.3,
  #         constraint_lab_style = constraint_lab_style(
  #           alpha = 0.8, size = 3, max_angle = 90, perpendicular = TRUE
  #         )
  #       ) +
  #         stat_ellipse(aes(linetype = .data[[input$colorBy]], colour = .data[[input$colorBy]])) +
  #         coord_fixed(ratio = 1, clip = "off", xlim = c(-6, 6)) +
  #         scale_colour_brewer(palette = "Set1") +
  #         scale_shape_manual(values = c(
  #           "ASPI" = 21, "BAL" = 22, "EXP" = 23, "MOCK" = 24, "NPS" = 25,
  #           "PCRneg" = 26, "TEneg" = 27, "TS" = 28  # Ensure these values match your data
  #         )) +
  #         labs(colour = input$colorBy, shape = input$shapeBy) +  # Add separate legend titles
  #         guides(
  #           colour = guide_legend(title = input$colorBy),
  #           shape = guide_legend(title = input$shapeBy)
  #         ) +  # Ensures separate legends
  #         theme_minimal()
  #
  #       # Print the static ggplot object
  #       print(plot)
  #
  #     }, error = function(e) {
  #       showNotification(paste("Error in rendering plot:", e$message), type = "error")
  #     })
  #   })
  #
  #   # Download handler for saving the plot in different formats
  #   output$download_rdaPlot <- downloadHandler(
  #     filename = function() {
  #       paste("rda_plot", Sys.Date(), ".", input$rda_filetype, sep = "")
  #     },
  #     content = function(file) {
  #       # Re-generate the ggplot object to ensure it's updated
  #       plot <- ord_plot(
  #         ordination,
  #         data = sample_data(pseq),  # Use the same data setup as in the renderPlot
  #         colour = input$colorBy,
  #         size = 2,
  #         alpha = 0.5,
  #         shape = input$shapeBy,
  #         interactive = FALSE,  # Static plot for saving
  #         plot_taxa = 1:8,
  #         tax_vec_length = 4.5,
  #         tax_lab_length = 4.6,
  #         tax_lab_style = tax_lab_style(
  #           type = "text", max_angle = 90, fontface = "bold.italic"
  #         ),
  #         constraint_vec_style = vec_constraint(1.5, alpha = 0.5),
  #         constraint_vec_length = 3, constraint_lab_length = 3.3,
  #         constraint_lab_style = constraint_lab_style(
  #           alpha = 0.8, size = 3, max_angle = 90, perpendicular = TRUE
  #         )
  #       ) +
  #         stat_ellipse(aes(linetype = .data[[input$colorBy]], colour = .data[[input$colorBy]])) +
  #         coord_fixed(ratio = 1, clip = "off", xlim = c(-6, 6)) +
  #         scale_colour_brewer(palette = "Set1") +
  #         scale_shape_manual(values = c(
  #           "ASPI" = 21, "BAL" = 22, "EXP" = 23, "MOCK" = 24, "NPS" = 25,
  #           "PCRneg" = 26, "TEneg" = 27, "TS" = 28  # Ensure these values match your data
  #         )) +
  #         labs(colour = input$colorBy, shape = input$shapeBy) +  # Add separate legend titles
  #         guides(
  #           colour = guide_legend(title = input$colorBy),
  #           shape = guide_legend(title = input$shapeBy)
  #         ) +  # Ensures separate legends
  #         theme_minimal()
  #
  #       # Save the plot to the desired format
  #       ggsave(file, plot = plot, device = input$rda_filetype, width = 10, height = 8, units = "in")
  #     }
  #   )
  #
  # })
  #
  # # Download the updated phyloseq object
  # output$downloadPhyseq <- downloadHandler(
  #   filename = function() {
  #     paste("imputed_phyloseq_object", Sys.Date(), ".rds", sep = "")
  #   },
  #   content = function(file) {
  #     saveRDS(reactive_pseq(), file)
  #   }
  # )
  #
  # # Add the download button in the UI
  # output$downloadPhyseqUI <- renderUI({
  #   downloadButton("downloadPhyseq", "Download Updated Phyloseq Object")
  # })
  #
  # observeEvent(input$show_rda_explanation, {
  #   showModal(modalDialog(
  #     title = "Redundancy Analysis (RDA) Explanation",
  #     HTML("<p><strong>Redundancy Analysis (RDA)</strong> is a constrained ordination method used to understand how much variation in the data can be explained by specific environmental or experimental variables.</p>
  #        <p><strong>Why RDA?</strong></p>
  #        <ul>
  #          <li><strong>Identify Drivers:</strong> RDA helps identify which variables are the main drivers behind the differences in microbial communities across samples.</li>
  #          <li><strong>Visual Interpretation:</strong> The RDA plot provides a visual interpretation of the relationships between samples, taxa, and the chosen constraints/conditions.</li>
  #          <li><strong>Statistical Significance:</strong> It allows for testing the statistical significance of the relationships between the constraints and the community composition.</li>
  #        </ul>
  #        <p>In this app, RDA is used to explore how different metadata columns (constraints) can explain the variation in the microbiome data across samples.</p>"),
  #     easyClose = TRUE,
  #     footer = modalButton("Close")
  #   ))
  # })

})

# # generates vector of ggplot2 shapes
# ggplot2_shapes <- function() {
#   c(
#     "circle", paste("circle", c("open", "filled", "cross", "plus", "small")),
#     "bullet",
#     "square", paste("square", c("open", "filled", "cross", "plus", "triangle")),
#     "diamond", paste("diamond", c("open", "filled", "plus")),
#     "triangle", paste("triangle", c("open", "filled", "square")),
#     paste("triangle down", c("open", "filled")),
#     "plus", "cross", "asterisk"
#   )
# }

  # # Function to impute missing values with the median
  # impute_median <- function(pseq) {
  #   metadata <- sample_data(pseq)
  #   imputed_metadata <- metadata
  #   for (col in colnames(metadata)) {
  #     if (is.numeric(metadata[[col]])) {
  #       imputed_metadata[[col]][is.na(imputed_metadata[[col]])] <- median(imputed_metadata[[col]], na.rm = TRUE)
  #     }
  #   }
  #   sample_data(pseq) <- imputed_metadata
  #   return(pseq)
  # }
  #
  # # Create a reactive value to store the pseq object
  # reactive_pseq <- reactiveVal()
  #
  # # Impute missing values and update select inputs
  # observe({
  #   req(physeq())
  #   pseq <- physeq()
  #
  #   # Impute missing values before analysis
  #   pseq <- impute_median(pseq)
  #
  #   # Store the imputed pseq in the reactive value
  #   reactive_pseq(pseq)
  #
  #   # Update choices based on the available taxonomic ranks and sample data columns
  #   updateSelectInput(session, "taxRankRDA", choices = rank_names(pseq))
  #   sample_data_cols <- colnames(sample_data(pseq))
  #   updateSelectInput(session, "constraints", choices = sample_data_cols)
  #   updateSelectInput(session, "conditions", choices = c("None", sample_data_cols))
  #   updateSelectInput(session, "colorBy", choices = sample_data_cols)
  #   updateSelectInput(session, "shapeBy", choices = sample_data_cols)
  # })
  #
  # # Handle the RDA plot generation
  # observeEvent(input$plotRDA, {
  #   req(reactive_pseq(), input$taxRankRDA, input$normMethod, input$colorBy, input$shapeBy, input$constraints, input$conditions)
  #
  #   # Retrieve pseq from the reactive object
  #   pseq <- reactive_pseq()
  #
  #   pseq <- phyloseq_validate(pseq, remove_undetected = TRUE)
  #   pseq <- tax_fix(pseq)
  #   pseq <- tax_agg(ps = pseq, input$taxRankRDA)
  #
  #   if (input$normMethod != "identity") {
  #     pseq <- tax_transform(pseq, trans = input$normMethod, rank = input$taxRankRDA)
  #   }
  #
  #   # Handle constraints and conditions correctly for multiple selections
  #   constraints <- if (is.null(input$constraints) || length(input$constraints) == 0 || "None" %in% input$constraints) {
  #     NULL
  #   } else {
  #     input$constraints
  #   }
  #
  #   conditions <- if (is.null(input$conditions) || length(input$conditions) == 0 || "None" %in% input$conditions) {
  #     NULL
  #   } else {
  #     input$conditions
  #   }
  #
  #   # RDA calculation
  #   ordination <- pseq %>%
  #     ord_calc(
  #       method = "RDA",
  #       constraints = constraints,
  #       conditions = conditions
  #     )
  #
  #   # Display RDA statistics
  #   rda_stats <- summary(ordination) # Or another method to extract statistics
  #   output$rdaInfoTable <- DT::renderDataTable({
  #     rda_stats$cont  # Adjust as necessary based on the structure of the statistics object
  #   })
  #
  #
  #
  #   output$rdaPlot <- renderPlotly({
  #     tryCatch({
  #       # Create the ggplot object with ord_plot
  #       interactive_plot <- ord_plot(
  #         ordination,
  #         colour = input$colorBy,
  #         size = 2,
  #         alpha = 0.5,
  #         shape = input$shapeBy,
  #         interactive = FALSE,  # Set to FALSE to generate ggplot for accurate legend display
  #         plot_taxa = 1:8,
  #         tax_vec_length = 4.5,
  #         tax_lab_length = 4.6,
  #         tax_lab_style = tax_lab_style(
  #           type = "text", max_angle = 90, fontface = "bold.italic"
  #         ),
  #         constraint_vec_style = vec_constraint(1.5, alpha = 0.5),
  #         constraint_vec_length = 3, constraint_lab_length = 3.3,
  #         constraint_lab_style = constraint_lab_style(
  #           alpha = 0.8, size = 3, max_angle = 90, perpendicular = TRUE
  #         )
  #       ) +
  #         stat_ellipse(aes(linetype = .data[[input$colorBy]], colour = .data[[input$colorBy]])) +
  #         coord_fixed(ratio = 1, clip = "off", xlim = c(-6, 6)) +
  #         scale_colour_brewer(palette = "Set1") +
  #         labs(colour = input$colorBy, shape = input$shapeBy) +  # Add separate legend titles
  #         theme_minimal()
  #
  #       # Convert the ggplot object to plotly for interactive rendering
  #       plotly_plot <- ggplotly(interactive_plot) %>%
  #         layout(legend = list(orientation = "v", x = 1.02, y = 1))
  #
  #       plotly_plot
  #
  #     }, error = function(e) {
  #       showNotification(paste("Error in rendering plot:", e$message), type = "error")
  #     })
  #   })
  #
  #
  #   # Download plot as selected file type using orca for plotly
  #   output$download_rdaPlot <- downloadHandler(
  #     filename = function() {
  #       paste("rda_plot", Sys.Date(), ".", input$rda_filetype, sep = "")
  #     },
  #     content = function(file) {
  #       if (input$rda_filetype %in% c("png", "jpeg", "pdf", "svg")) {
  #         orca(ggplotly(interactive_plot), file)  # Save using orca for plotly
  #       } else {
  #         saveWidget(as_widget(ggplotly(interactive_plot)), file, selfcontained = TRUE)  # Save as HTML
  #       }
  #     }
  #   )
  # })
  #
  #
  #
  #
  # # Download the updated phyloseq object
  # output$downloadPhyseq <- downloadHandler(
  #   filename = function() {
  #     paste("imputed_phyloseq_object", Sys.Date(), ".rds", sep = "")
  #   },
  #   content = function(file) {
  #     saveRDS(reactive_pseq(), file)
  #   }
  # )
  #
  #
  # # Add the download button in the UI
  # output$downloadPhyseqUI <- renderUI({
  #   downloadButton("downloadPhyseq", "Download Updated Phyloseq Object")
  # })
  #
  # observeEvent(input$show_rda_explanation, {
  #   showModal(modalDialog(
  #     title = "Redundancy Analysis (RDA) Explanation",
  #     HTML("<p><strong>Redundancy Analysis (RDA)</strong> is a constrained ordination method used to understand how much variation in the data can be explained by specific environmental or experimental variables.</p>
  #          <p><strong>Why RDA?</strong></p>
  #          <ul>
  #            <li><strong>Identify Drivers:</strong> RDA helps identify which variables are the main drivers behind the differences in microbial communities across samples.</li>
  #            <li><strong>Visual Interpretation:</strong> The RDA plot provides a visual interpretation of the relationships between samples, taxa, and the chosen constraints/conditions.</li>
  #            <li><strong>Statistical Significance:</strong> It allows for testing the statistical significance of the relationships between the constraints and the community composition.</li>
  #          </ul>
  #          <p>In this app, RDA is used to explore how different metadata columns (constraints) can explain the variation in the microbiome data across samples.</p>"),
  #     easyClose = TRUE,
  #     footer = modalButton("Close")
  #   ))
  # })

