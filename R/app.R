#' @title vcf2ploidy Shiny Application
#'
#' @description Open a shiny gadget for running vcf2ploidy functions
#'
#' @details Drag and drop VCF file to run vcf2had(), vcf2ploidy(), or vcf2colony() on and select parameters in graphical interface
#'
#' @examples
#' \dontrun{vcf2ploidy_app()}
#'
#' @importFrom shiny icon
#' @importFrom shiny fileInput
#' @importFrom shiny textInput
#' @importFrom shiny selectInput
#' @importFrom shiny numericInput
#' @importFrom shiny radioButtons
#' @importFrom shiny checkboxInput
#' @importFrom shiny checkboxGroupInput
#' @importFrom shiny conditionalPanel
#' @importFrom shiny actionButton
#' @importFrom shiny observeEvent
#' @importFrom shiny runGadget
#' @importFrom shiny stopApp
#' @importFrom miniUI miniPage
#' @importFrom miniUI gadgetTitleBar
#' @importFrom miniUI miniContentPanel
#' @importFrom miniUI miniTabstripPanel
#' @importFrom miniUI miniTabPanel
#' @export
vcf2ploidy_app <- function(){

    ui <- miniPage(
        gadgetTitleBar("vcf2ploidy", right=NULL),
        miniTabstripPanel(
            # vcf2had panel
            miniTabPanel("Convert to HAD", icon = icon("list-ul"),
                         miniContentPanel(
                             fileInput(inputId = "vcfFile_vcf2had",
                                      label = "Drag and drop VCF file here",
                                      multiple = FALSE,
                                      buttonLabel = "Browse...",
                                      placeholder = "No file selected"),
                            numericInput(inputId="skip_lines_vcf2had",
                                         label="Number of metadata lines to skip (leave blank to automatically count metadata lines):",
                                         value=NULL),
                            radioButtons(inputId="remove_double_hets_vcf2had",
                                         label="Remove Double Heterozygotes?",
                                         choiceNames=c("Yes", "No"),
                                         choiceValues=c(TRUE, FALSE)),

                            # Save options for output
                            radioButtons(inputId="save_output_vcf2had",
                                         label="Save output to:",
                                         choiceNames = c("Data frame in global env.",
                                                         "Text file",
                                                         "RDS file"),
                                         choiceValues = c("save_dataframe",
                                                          "save_textfile",
                                                          "save_rsdfile")),
                            # Option for naming a data frame
                            conditionalPanel(condition = "input.save_output_vcf2had == 'save_dataframe'",
                                             textInput(inputId="out_df_name_vcf2had",
                                                       label = "Data frame name:",
                                                       value = "had_df")),
                            # Option for naming a text file
                            conditionalPanel(condition = "input.save_output_vcf2had == 'save_textfile'",
                                             textInput(inputId="out_textfilename_vcf2had",
                                                       label = "File name:",
                                                       value = "converted.txt")),
                            # Option for naming an RDS file
                            conditionalPanel(condition = "input.save_output_vcf2had == 'save_rsdfile'",
                                             textInput(inputId="out_RDSfilename_vcf2had",
                                                       label = "File name:",
                                                       value = "converted.rds")),
                            actionButton(inputId = "convert_to_HAD",
                                         label = "Convert to HAD",
                                         style="color: #fff;
                                     background-color: #337ab7;
                                     border-color: #2e6da4"))
            ),


            # vcf2ploidy panel
            miniTabPanel("Estimate Ploidy", icon=icon("dna"),
                         miniContentPanel(
                             fileInput(inputId = "vcfFile_vcf2ploidy",
                                       label = "Drag and drop VCF file here",
                                       multiple = FALSE,
                                       buttonLabel = "Browse...",
                                       placeholder = "No file selected"),
                             numericInput(inputId="skip_lines_vcf2ploidy",
                                          label="Number of metadata lines to skip (leave blank to automatically count metadata lines):",
                                          value=NULL),
                             radioButtons(inputId="remove_double_hets_vcf2ploidy",
                                          label="Remove Double Heterozygotes?",
                                          choiceNames=c("Yes", "No"),
                                          choiceValues = c(TRUE, FALSE)),
                             checkboxInput(inputId = "allow_custom_props",
                                           label = "Custom proportions?",
                                           value = FALSE),
                             conditionalPanel(condition = "input.allow_custom_props == false",
                                              checkboxGroupInput(inputId = "checkbox_props",
                                                                 label = "Allowed proportions:",
                                                                 choices=c(0.25, 0.33, 0.5, 0.66, 0.75),
                                                                 selected=c(0.25, 0.33, 0.5, 0.66, 0.75))),
                             conditionalPanel(condition = "input.allow_custom_props == true",
                                              textInput(inputId = "custom_props",
                                                        label = "Allowed proportions (separate values by commas, no spaces):",
                                                        value = "0.25,0.33,0.5,0.66,0.75")),
                             numericInput(inputId = "mcmc_nchain",
                                          label = "Number of chains for Markov chain Monte Carlo:",
                                          value = 2),
                             numericInput(inputId = "mcmc_steps",
                                          label = "Number of post burnin iterations for each chain:",
                                          value = 10000),
                             numericInput(inputId = "mcmc_burnin",
                                          label = "Number of iterations to discard from each chain as a burnin:",
                                          value = 1000),
                             numericInput(inputId = "mcmc_thin",
                                          label = "Thinning interval for MCMC:",
                                          value = 2),
                             radioButtons(inputId = "train",
                                          label = "Use training set?",
                                          choices = c("Yes"=TRUE, "No"=FALSE),
                                          selected = "FALSE"),
                             conditionalPanel(condition = "input.train == 'TRUE'",
                                              textInput(inputId = "known_individuals_set",
                                                        label = "Which individuals have known ploidies? (list column indices, separated by commas, no spaces)",
                                                        value = "2,3,4"),
                                              textInput(inputId = "known_ploidies",
                                                        label = "What are their respective ploidies? List them, separated by commas, no spaces.",
                                                        value = "Diploid,Triploid,Tetraploid"),
                                              numericInput(inputId = "nclasses",
                                                           label = "Number of cytotypes:",
                                                           value = 2)),
                             textInput(inputId = "pcs",
                                       label = "pcs: \"A vector giving the PC to use for DA\"",
                                       value = "1,2"),

                             # Save options for output
                             radioButtons(inputId="save_output_vcf2ploidy",
                                          label="Save output to:",
                                          choiceNames = c("List object in global env.",
                                                          "RDS file"),
                                          choiceValues = c("save_list",
                                                           "save_rsdfile")),
                             # Option for naming a list
                             conditionalPanel(condition = "input.save_output_vcf2ploidy == 'save_list'",
                                              textInput(inputId="out_list_name_vcf2ploidy",
                                                        label = "Output list object name:",
                                                        value = "vcf2ploidy_out")),
                             # Option for naming an RDS file
                             conditionalPanel(condition = "input.save_output_vcf2ploidy == 'save_rsdfile'",
                                              textInput(inputId="out_RDSfilename_vcf2ploidy",
                                                        label = "File name:",
                                                        value = "vcf2ploidy_out.rds")),
                             actionButton(inputId = "estimate_ploidy",
                                          label = "Estimate Ploidy",
                                          style="color: #fff;
                                     background-color: #337ab7;
                                     border-color: #2e6da4"))
            ),


            # vcf2colony panel
            miniTabPanel("Convert to COLONY", icon=icon("sitemap"),
                         miniContentPanel(
                             fileInput(inputId = "vcfFile_vcf2colony",
                                       label = "Drag and drop VCF file here",
                                       multiple = FALSE,
                                       buttonLabel = "Browse...",
                                       placeholder = "No file selected"),
                             numericInput(inputId="skip_lines_vcf2colony",
                                          label="Number of metadata lines to skip (leave blank to automatically count metadata lines):",
                                          value=NULL),
                             textInput(inputId = "out_filename_vcf2colony",
                                       label = "Name of converted file (saves in current working directory:",
                                       value = "converted.txt"),
                             actionButton(inputId = "convert_to_colony",
                                          label = "Convert to Colony",
                                          style="color: #fff;
                                     background-color: #337ab7;
                                     border-color: #2e6da4"))
            )
        )
    )

    server <- function(input, output, session){
        observeEvent(input$convert_to_HAD,
                     {shiny::req(input$vcfFile_vcf2had)
                      had_df <- vcf2had(input$vcfFile_vcf2had$datapath,
                                    skip_lines = input$skip_lines_vcf2had,
                                    remove_double_hets = as.logical(input$remove_double_hets_vcf2had))

                      ## Options for saving output
                      # Data frame
                      if(input$save_output_vcf2had == "save_dataframe"){
                        assign(input$out_df_name_vcf2had, had_df, envir=.GlobalEnv)
                        message(paste("Output assigned to data frame named:", input$out_df_name_vcf2had))
                      }
                      # Text file
                      else if(input$save_output_vcf2had == "save_textfile"){
                          # Write the converted data to a text file
                          write.table(had_df, file=input$out_textfilename_vcf2had,
                                      quote=FALSE, sep="\t",
                                      row.names=FALSE, col.names=TRUE)
                          message(paste("Saved to text file:", input$out_textfilename_vcf2had))
                      }
                      # RDS file
                      else{
                          saveRDS(had_df, input$out_RDSfilename_vcf2had)
                          message(paste("Saved to RDS file:", input$out_RDSfilename_vcf2had))
                      }
                      stopApp();
                     })
        observeEvent(input$estimate_ploidy,
                     {shiny::req(input$vcfFile_vcf2ploidy)

                      # Clean inputs for vcf2ploidy()
                      clean_args <- handle_vcf2ploidy_args(input)


                      vcf2ploidy_out <-
                          vcf2ploidy(filename = input$vcfFile_vcf2ploidy$datapath,
                                     skip_lines = input$skip_lines_vcf2ploidy,
                                     remove_double_hets = input$remove_double_hets_vcf2ploidy,
                                     props = clean_args$props,
                                     mcmc.nchain = input$mcmc_nchain,
                                     mcmc.steps = input$mcmc_steps,
                                     mcmc.burnin = input$mcmc_burnin,
                                     mcmc.thin = input$mcmc_thin,
                                     train = clean_args$train,
                                     pl = clean_args$pl,
                                     set = clean_args$set,
                                     nclasses = input$nclasses,
                                     pcs = clean_args$pcs)

                      ## Options for saving output
                      # List
                      if(input$save_output_vcf2ploidy == "save_list"){
                          assign(input$out_list_name_vcf2ploidy, vcf2ploidy_out, envir=.GlobalEnv)
                          message(paste("Output assigned to list named:", input$out_df_name_vcf2ploidy))
                      }
                      # RDS file
                      else{
                          saveRDS(vcf2ploidy_out, input$out_RDSfilename_vcf2ploidy)
                          message(paste("Saved to RDS file:", input$out_RDSfilename_vcf2ploidy))
                      }



                     stopApp()
                     })
        observeEvent(input$convert_to_colony,
                     {shiny::req(input$vcfFile_vcf2colony)
                      vcf2colony(input$vcfFile_vcf2colony$datapath,
                                 skip_lines = input$skip_lines_vcf2colony,
                                 out_filename = input$out_filename_vcf2colony)
                      stopApp()
                      })
    }

    runGadget(ui, server)
}



#' @title Handle Arguments for vcf2ploidy from Shiny Input
#'
#' @description Convert text input values from vcf2ploidy_app() into appropriate vectors of values for vcf2ploidy arguments
#'
#' @return A list containing cleaned arguments for vcf2ploidy()
handle_vcf2ploidy_args <- function(input){

    # Handle proportions inputs
    if(input$allow_custom_props == "TRUE"){
        # Split text by a comma and cast to numeric
        props <- strsplit(input$custom_props, ",")[[1]]
        props <- as.numeric(props)
    }
    else{
        # Use checkboxed proportions if no custom proportions
        props <- input$checkbox_props
    }

    # Handle training data inputs
    if(input$train == "TRUE"){
        # Fix train value
        train <- TRUE
        # Split text by comma and cast to numeric
        set <- strsplit(input$known_individuals_set, ",")[[1]]
        set <- as.numeric(set)

        # Split ploidies by comma
        pl <- strsplit(input$known_ploidies, ",")[[1]]
    }
    else{
        # Fix train values
        train <- FALSE
        pl <- NA
        set <- NA
    }

    # Handle pcs
    pcs <- strsplit(input$pcs, ",")[[1]]

    # Combine arguments into list to be returned
    args_out <- list("props"=props, "train"=train, "set"=set, "pl"=pl, "pcs"=pcs)

    return(args_out)
}
