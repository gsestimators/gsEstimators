# Load packages ----------------------------------------------------------------
library(shiny)
library(bslib)
library(ggplot2)
library(tools)
library(thematic)
library(dplyr)
library(DT)
library(mvtnorm)
library(gsDesign)
library(rootSolve)
library(shinythemes)

# Load data --------------------------------------------------------------------

#load("")
#thematic::thematic_shiny()

#load r files containing estimators --------------------------------------------

source("R/estimator_inputs/binary_estimators_single_arm.R")
source("R/estimator_inputs/binary_estimators_two_arm.R")
source("R/estimator_inputs/normal_estimators_single_arm.R")
source("R/estimator_inputs/normal_estimators_parallel.R")
source("R/estimator_inputs/normal_estimators_paired.R")
source("R/estimator_inputs/normal_estimators_crossover.R")
source("R/estimator_inputs/time_to_event_estimators.R")
source("R/operating_characteristics/binary_data_operating_characteristics.R")
source("R/operating_characteristics/normal_data_operating_characteristics.R")
source("R/operating_characteristics/time_to_event_data_operating_characteristics.R")
source("R/operating_characteristics/correlation_matrix.R")
source("R/dists/dsnorm.R")
source("R/dists/psnorm.R")
source("R/dists/z1_conditional_dsnorm.R")
source("R/dists/zk_conditional_dsnorm.R")
source("R/error_checking/vector_input_check.R")
source("R/error_checking/stage_input_check.R")
source("R/error_checking/limit_input_check.R")
source("R/error_checking/binary_single_arm_input_check.R")
source("R/error_checking/binary_two_arm_input_check.R")
source("R/error_checking/normal_single_arm_input_check.R")
source("R/error_checking/normal_parallel_input_check.R")
source("R/error_checking/normal_paired_input_check.R")
source("R/error_checking/normal_crossover_input_check.R")
source("R/error_checking/time_to_event_check.R")
source("R/estimators/mue.R")
source("R/estimators/cmue.R")
source("R/estimators/umvue.R")
source("R/estimators/umvcue.R")
source("R/estimators/ubc_mle.R")
source("R/estimators/cbc_mle.R")

# Define UI --------------------------------------------------------------------
ui <- page_sidebar(
  theme = bs_theme(bootswatch = "cerulean"),
  title = "Group Sequential Estimators",
  sidebar= sidebar(
#SIDEBAR INPUTS ----------------------------------------------------------------

#SelectInput for the type of endpoint used in the trial
    tags$h5("Estimators"),
    selectInput(inputId = "endpoint_type",
                label = HTML("Select type of endpoint:"),
                choices = c( "Binary",
                             "Normal",
                             "Time to Event"),
                selected = "Binary"),

#selectInput for the type of Binary trial
    conditionalPanel(condition = "input.endpoint_type == 'Binary'",
                 selectInput(inputId = "binary_trial_type",
                             label = HTML("Choose type of trial:"),
                             choices = c("Single-arm",
                                         "Two-arm"),
                             selected = "Two-arm"
                            )
    ),

#selectInput for the type of Normal trial
    conditionalPanel(condition = "input.endpoint_type == 'Normal'",
                 selectInput(inputId = "normal_trial_type",
                             label = HTML("Choose type of trial:"),
                             choices = c("Single-arm",
                                         "Parallel",
                                         "Paired",
                                         "Two Period Crossover"
                                         ),
                             selected = "Parallel"
                 )
    ),

#selectInput for the type of Time to Event trial
    conditionalPanel(condition = "input.endpoint_type == 'Time to Event'",
                     selectInput(inputId = "time_to_event_trial_type",
                                 label = HTML("Choose type of trial:"),
                                 choices = c("Single-arm",
                                             "Two-arm"),
                                 selected = "Two-arm"
                                )
    ),




#CHOOSE WHICH ESTIMATOR TO CALCULATE -------------------------------------------

#SelectInput for Binary estimators
    conditionalPanel(condition = "input.endpoint_type == 'Binary'",
                     selectInput(inputId = "binary_estimator",
                                 label = HTML("Choose estimator:"),
                                 choices = c("Overall MLE" ,
                                             "Stage 1 MLE" ,
                                             "MUE",
                                             "UMVUE",
                                             "ubc-MLE",
                                             "cMLE",
                                             "cMUE",
                                             "UMVCUE",
                                             "cbc-MLE"),
                                 selected = "Overall MLE"
                     )
    ),

#SelectInput for Time to Event estimators
    conditionalPanel(condition = "input.endpoint_type == 'Time to Event'",
                     selectInput(inputId = "time_to_event_estimator",
                                 label = HTML("Choose estimator:"),
                                 choices = c("Overall MLE",
                                             "Stage 1 MLE",
                                             "MUE",
                                             "UMVUE",
                                             "ubc-MLE",
                                             "cMLE",
                                             "cMUE",
                                             "UMVCUE",
                                             "cbc-MLE"),
                                 selected = "Overall MLE"
                     )
    ),

#SelectInput for Normal estimators
    conditionalPanel(condition = "input.endpoint_type == 'Normal'",
                     selectInput(inputId = "normal_estimator",
                                 label = HTML("Choose estimator:"),
                                 choices = c("Overall MLE",
                                             "Stage 1 MLE",
                                             "MUE",
                                             "UMVUE",
                                             "ubc-MLE",
                                             "cMLE",
                                             "cMUE",
                                             "UMVCUE",
                                             "cbc-MLE"),
                                 selected = "Overall MLE"
                     )
    ),




#INPUT INFORMATION ABOUT TRIAL STAGES ------------------------------------------
    tags$hr(),
    tags$h5("Stages"),

#numericInput for maximum number of stages
      numericInput(inputId = "maximum_stages",
             label = HTML("Maximum number of trial stages:"),
             min = 1,
             max = NA,
             value = 2,
             step = 1),

#numericInput for stage realized
    numericInput(inputId = "stage_reached",
                 label = HTML("Trial stage reached:"),
                 min = 1,
                 max = NA,
                 value = 2,
                 step = 1),




#DEFINE UPPER AND LOWER STOPPING BOUNDS ----------------------------------------
    tags$hr(),
    tags$h5("Stoppping bounds"),

#textInput the futility bounds for the trial
    textInput(inputId = "lower_bound",
              label = HTML("Futility bounds (csv):"),
              value = "-Inf"),

#textInput the superiority bounds for the trial
    textInput(inputId = "upper_bound",
              label = HTML("Superiority bound (csv):"),
              value = "2.797")
  ),






#CONTROL PANEL INPUTS ----------------------------------------------------------
  layout_columns(
    card(tags$h5("Control Panel"),

#numericInput for epsilon
fluidRow(
    column( width = 6,
           tags$h6("Arguments for Density Approximation"),
           numericInput(inputId = "epsilon",
                        label = HTML("Epsilon (argument for probability density):"),
                        max = 1,
                        min = 1e-10,
                        value = 1e-5)
    )
),

#ARGUMENTS FOR OPTIMIZE SEARCH REGION-------------------------------------------

#Binary inputs
conditionalPanel(condition = "input.endpoint_type == 'Binary' &
                 input.binary_estimator == 'ubc-MLE' ||
                 input.endpoint_type == 'Binary' &
                 input.binary_estimator == 'cbc-MLE'",
                  fluidRow( column( width = 6,
                       tags$h6("Arguments for optimize"),
                       textInput(inputId = "binary_search_region_optimize",
                                 label = HTML("Upper and lower search limits
                                            (csv)"),
                                 value = "-1,1")
                                  )
                          )
                ),

#Normal inputs
conditionalPanel(condition = "input.endpoint_type == 'Normal' &
                 input.normal_estimator == 'ubc-MLE' ||
                 input.endpoint_type == 'Normal' &
                 input.normal_estimator == 'cbc-MLE'",
                  fluidRow( column(width = 6,
                            tags$h6("Arguments for optimize"),
                            textInput(inputId = "normal_search_region_optimize",
                                      label = HTML("Upper and lower search limits
                                            (csv)"),
                                      value = "-10,10")
                                  )
                          )
                ),


#NUMERICINPUT FOR STAGE TO CONDITION ON-----------------------------------------
      fluidRow(
          column( width = 6,
            conditionalPanel(condition = "input.endpoint_type == 'Binary' &
                             input.binary_estimator == 'cMLE' ||
                             input.endpoint_type == 'Normal' &
                             input.normal_estimator == 'cMLE'",
                          tags$h6("Conditional Stopping Stage"),
                          numericInput(inputId = "conditional_stop_stage",
                                       label = HTML("Choose stage to condition on
                                       (must be greater than total stages):"),
                                       min = 2,
                                       max = 1000000,
                                       value = 2)
                            )
                )
              ),

#ARGUMENTS FOR CALCULATION OF THE UMVUE AND UMVCUE -----------------------------
         conditionalPanel(condition = "input.endpoint_type == 'Binary' &
                          input.binary_estimator == 'UMVUE' ||
                          input.endpoint_type == 'Binary' &
                          input.binary_estimator == 'UMVCUE' ||
                          input.endpoint_type == 'Binary' &
                          input.binary_estimator == 'ubc-MLE' ||
                          input.endpoint_type == 'Binary' &
                          input.binary_estimator == 'cbc-MLE' ||
                          input.endpoint_type == 'Normal' &
                          input.normal_estimator == 'UMVUE' ||
                          input.endpoint_type == 'Normal' &
                          input.normal_estimator == 'UMVCUE' ||
                          input.endpoint_type == 'Normal' &
                          input.normal_estimator == 'ubc-MLE' ||
                          input.endpoint_type == 'Normal' &
                          input.normal_estimator == 'cbc-MLE'",
                          tags$h6("Arguments for Integrate Function"),
                          fluidRow(

#numericInput for relative tolerance
                          column( width = 4,
                            numericInput(inputId = "relative_tolerance",
                                         label = HTML("Relative Tolerance
                                         (argument for integrate):"),
                                         min = 0.0001220703,
                                         max = 1,
                                         value = 0.0001220703)
                                ),

#numericInput for absolute tolerance
                          column( width = 4,
                            numericInput(inputId = "absolute_tolerance",
                                         label = HTML("Absolute Tolerance
                                         (argument for integrate):"),
                                         max = 1,
                                         min = 0.0001220703,
                                         value = 0.0001220703)
                                ),

#numericInput for number of subdivisions
                          column( width = 4,
                            numericInput(inputId = "subdivisions",
                                         label = HTML("Subdivisions
                                         (argument for integrate):"),
                                         max = 10000,
                                         min = 100,
                                         value = 100)
                                )
                          )
         ),

#ARGUMENTS FOR CALCULATION OF THE MUE ------------------------------------------

#textInput for upper and lower values of uniroot search region
conditionalPanel( condition = "input.endpoint_type == 'Binary' &
                  input.binary_estimator == 'MUE'||
                  input.endpoint_type == 'Binary' &
                  input.binary_estimator == 'cMUE'",
                  tags$h6("Uniroot Search Region"),
                  fluidRow( column( width = 6,
                                    textInput(inputId = "binary_search_region_uniroot",
                                              label = HTML("Upper and lower search limits
                                               (csv):"),
                                                value = "-1,1")
                                  )
                          )
                ),

conditionalPanel( condition = "input.endpoint_type == 'Normal' &
                  input.normal_estimator == 'MUE'||
                  input.endpoint_type == 'Normal'&
                  input.normal_estimator == 'cMUE'",
                  tags$h6("Uniroot Search Region"),
                  fluidRow( column( width = 6,
                                    textInput(inputId = "normal_search_region_uniroot",
                                              label = HTML("Upper and lower search limits
                                               (csv):"),
                                              value = "-10,10")
                                   )
                           )
                ),


         conditionalPanel(condition = "input.endpoint_type == 'Binary' &
                          input.binary_estimator == 'MUE' ||
                          input.endpoint_type == 'Binary' &
                          input.binary_estimator == 'cMUE' ||
                          input.endpoint_type == 'Normal' &
                          input.normal_estimator == 'MUE' ||
                          input.endpoint_type == 'Normal' &
                          input.normal_estimator == 'cMUE'",
                          tags$h6("Arguments for Uniroot"),
                          fluidRow(

#numericInput for tolerance
                          column( width = 6,
                            numericInput(inputId = "tolerance",
                                         label = HTML("Tolerance:"),
                                         max = 1,
                                         min = 0.0001220703,
                                         value = 0.0001220703)
                                ),

#numericInput for maximum number of iterations
                          column( width = 6,
                            numericInput(inputId = "maximum_iterations",
                                         label = HTML("Maximum iterations"),
                                         min = 1000,
                                         max = 10000,
                                         value = 1000)
                          ),

#numericInput for maximum number of sub-intervals
                          column( width = 6,
                            numericInput(inputId = "maximum_sub_intervals",
                                         label = HTML("Maximum sub-intervals:"),
                                         min = 100,
                                         max = 10000,
                                         value = 100)
                                ),


#radioButtons for which tail to search for estimator in
                          column( width = 6,
                            radioButtons(inputId = "tail_type",
                                         label = HTML("Tail type for median unbiased
                                         estimation:"),
                                         choices = list(
                                           "Two-tailed" = "two",
                                           "Upper tail only" = "upper",
                                           "Lower tail only" = "lower"
                                         ),
                                         selected = "upper")
                                )
                          )
         )
    ),






#TRIAL DATA CARD ---------------------------------------------------------------
    card( tags$h5("Trial Data"),

#BINARY DATA (SINGLE_ARM) ------------------------------------------------------
conditionalPanel(condition = "input.endpoint_type == 'Binary' &
                               input.binary_trial_type == 'Single-arm'",
                 tags$h6("", style = "margin-top: -8px;"),

  fluidRow(

#textInput panel for cumulative sample size (experimental arm)
  column( width = 6,
                           textInput(inputId = "single_arm_sample_size_binary",
                                     label = HTML("Cumulative stagewise sample size
                                         in the Experimental arm (comma separated
                                         input):"),
                                     value = "101, 143")
          ),

  #textInput panel for cumulative number of events (experimental arm)
  column( width = 6,
                           textInput(inputId = "single_arm_events",
                                     label = HTML("Cumulative number of events in
                                         the Experimental arm (comma separated
                                         input):"),
                                     value = "20, 36")
          ),

  #textInput panel for the assumed control event rate
  column( width = 6,
                           textInput(inputId = "event_rate_null",
                                     label = HTML("Assumed event rate under the null
                                     hypothesis (comma separated values):"),
                                     value = "0.125")
          )
  )
),


#BINARY DATA (TWO-ARM) ---------------------------------------------------------
conditionalPanel(condition = "input.endpoint_type == 'Binary' &
                               input.binary_trial_type == 'Two-arm'",
                 tags$h6("", style = "margin-top: -8px;"),
  fluidRow(

#textInput panel for cumulative sample size (control arm)
  column( width = 6,
                               textInput(inputId = "control_sample_size_binary",
                                         label = HTML("Cumulative stagewise sample size
                                         in the Control arm (comma separated
                                         input):"),
                                         value = "97, 134")
    ),

#textInput panel for cumulative number of events (control arm)
  column( width = 6,
          textInput(inputId = "control_events",
                    label = HTML("Cumulative number of events in
                                         the Control arm (comma separated input):"),
                    value = "12, 21")
  ),

#textInput panel for cumulative sample size (experimental arm)
  column( width = 6,
                               textInput(inputId = "experimental_sample_size_binary",
                                         label = HTML("Cumulative stagewise sample size
                                         in the Experimental arm (comma separated
                                         input):"),
                                         value = "101, 143")
    ),

#textInput panel for cumulative number of events (experimental arm)
  column( width = 6,
                               textInput(inputId = "experimental_events",
                                         label = HTML("Cumulative number of events in
                                         the Experimental arm (comma separated
                                         input):"),
                                         value = "27, 42")
    )
  )
),

#NORMAL DATA (SINGLE-ARM) ------------------------------------------------------
conditionalPanel(condition = "input.endpoint_type == 'Normal' &
                       input.normal_trial_type == 'Single-arm'",
                 tags$h6("", style = "margin-top: -8px;"),
  fluidRow(

#textInput panel for mean at each stage (experimental arm)
  column( width = 6,
                           textInput(inputId = "means_normal_single_arm",
                                     label = HTML("Experimental arm mean per
                                 stage (comma separated input):"),
                                     value = "5.2, 6")
          ),

#textInput panel for sample size at each stage (experimental arm)
  column( width = 6,
                           textInput(inputId = "sample_size_normal_single_arm",
                                     label = HTML("Cumulative number of events in
                                     the Experimental arm (comma separated input):"),
                                     value = "34, 78")
          ),

#textInput panel for known variance in the experimental arm
  column(width = 6,
                          textInput(inputId = "variance_normal_single_arm",
                                    label = HTML("Known variance in the experimental arm:)"),
                                    value = "25")
         ),

#textInput panel for mean under null hypothesis
  column(width = 6,
                          textInput(inputId = "mean_null_normal_single_arm",
                                    label = HTML("Mean under the null hypothesis"),
                                    value = "3")
         )
  )
),

#NORMAL DATA (PARALLEL) --------------------------------------------------------
conditionalPanel(condition = "input.endpoint_type == 'Normal' &
                       input.normal_trial_type == 'Parallel'",
                 tags$h6("", style = "margin-top: -8px;"),
  fluidRow(

#textInput panel for mean at each stage (control arm)
  column( width = 6,
                       textInput(inputId = "parallel_control_arm_means",
                                 label = HTML("Control arm mean per
                                 stage (comma separated input):"),
                                 value = "3, 2.8")
   ),

#textInput panel for mean at each stage(experimental arm)
  column( width = 6,
                       textInput(inputId = "parallel_experimental_arm_means",
                                 label = HTML("Experimental arm mean per
                                 stage (comma separated input):"),
                                 value = "5.2, 6")
    ),

#textInput panel for sample size at each stage (control arm)
  column( width = 6,
                       textInput(inputId = "parallel_control_sample_size_normal",
                                 label = HTML("Cumulative stagewise sample size
                                         in the Control arm (comma separated
                                         input):"),
                                 value = "34, 78")
    ),

#textInput panel for sample size at each stage (experimental arm)
  column( width = 6,
                     textInput(inputId = "parallel_experimental_sample_size_normal",
                               label = HTML("Cumulative stagewise sample size in the
                               Experimental arm (comma separated
                                       input):"),
                               value = "34, 78")
    ),

#textInput panel for known variance (control arm)
  column( width = 6,
                     textInput(inputId = "parallel_variance_control",
                               label = HTML("Known variance in the control arm:"),
                               value = "25")
    ),

#textInput panel for known variance (control arm)
  column( width = 6,
                   textInput(inputId = "parallel_variance_experimental",
                             label = HTML("Known variance in the experimental arm:"),
                             value = "25")
    )
  )
),

#NORMAL DATA (PAIRED) ----------------------------------------------------------
conditionalPanel(condition = "input.endpoint_type == 'Normal' &
                       input.normal_trial_type == 'Paired'",
                 tags$h6("", style = "margin-top: -8px;"),
 fluidRow(

#textInput panel for mean at each stage (control arm)
  column( width = 6,
                           textInput(inputId = "paired_control_arm_means",
                                     label = HTML("Control arm mean per
                                 stage (comma separated input):"),
                                     value = "3, 2.8")
  ),

#textInput panel for mean at each stage(experimental arm)
  column( width = 6,
                           textInput(inputId = "paired_experimental_arm_means",
                                     label = HTML("Experimental arm mean per
                                 stage (comma separated input):"),
                                     value = "5.2, 6")
  ),

#textInput panel for the number of pairs of observations
  column(width = 6,
                          textInput(inputId = "paired_sample_size",
                                    label = HTML("Number of pairs of observations
                                    at each stage:"),
                                    value = "34, 78")
         ),

#textInput panel for the variance of pairs
  column(width = 6,
                          textInput(inputId = "paired_variance_normal",
                                    label = HTML("Known variance:"),
                                    value = "25")
         )
  )
),

#NORMAL DATA (CROSSOVER) -------------------------------------------------------
conditionalPanel(condition = "input.endpoint_type == 'Normal' &
                       input.normal_trial_type == 'Two Period Crossover'",
                 tags$h6("", style = "margin-top: -8px;"),
  fluidRow(

#textInput panel for means in the AB period group
  column(width = 6,
                          textInput(inputId = "crossover_means_ab",
                                    label = HTML("Mean in the AB period group at
                                    each stage:"),
                                    value = "5.2, 6")
         ),

#textInput panel for means in the BA period group
  column(width = 6,
                          textInput(inputId = "crossover_means_ba",
                                    label = HTML("Mean in the BA period group at
                                      each stage:"),
                                    value = "3, 2.8")
          ),

#textInput panel for sample size at each stage in the AB period group
  column(width = 6,
                          textInput(inputId = "sample_size_ab_crossover",
                                    label = HTML("Sample size in the AB period group
                                    at each stage:"),
                                    value = "34, 78")
         ),

#textInput panel for sample size at each stage in the AB period group
  column(width = 6,
                          textInput(inputId = "sample_size_ba_crossover",
                                    label = HTML("Sample size in the BA period group
                                      at each stage:"),
                                    value = "34, 78")
          ),

#textInput for the variance in the crossover trial
  column(width = 6,
                          textInput(inputId = "variance_crossover_ab",
                                           label = HTML("Variance in AB period arm:"),
                                           value = "25")
         ),

#textInput for the variance in the crossover trial
  column(width = 6,
                          textInput(inputId = "variance_crossover_ba",
                                    label = HTML("Variance in BA period arm:"),
                                    value = "25")
         )
  )
),

#TIME TO EVENT DATA ------------------------------------------------------------

#fileinput for Time to Event data
  conditionalPanel(condition = "input.endpoint_type == 'Time to Event' &
                    input.time_to_event_trial_type == 'Two-arm'",
                   tags$h6("", style = "margin-top: -8px;"),
                    fluidRow(column( width = 12,
                             fileInput(inputId = "time_to_event_data_csv",
                                       label = HTML("Please input a csv file containing the
                                       survival analysis data:"),
                                       accept = c(".csv"))
                                  )
                          )
    )
),





#OUTPUT PANEL ------------------------------------------------------------------
    card( tags$h5("Output"),

#actionButton that calculates estimate when pressed
          actionButton(inputId = "calc_butt",
                       label = "Calculate estimator"),

#uiOutput displaying either the estimate or any errors that occur
          uiOutput(outputId = "text")
    ),
    col_widths = c(6,6,12),
    row_heights = c(4,2,2)
  )
)




#DEFINE SERVER -----------------------------------------------------------------
server <- function(input, output, session) {


#Set reactive estimator value
  estimator <- reactiveVal(0)

#OBSERVER FUNCTION FOR THE BINARY (SINGLE-ARM) ESTIMATORS ----------------------
  observeEvent(c(input$calc_butt, input$endpoint_type, input$binary_estimator,
                 input$binary_trial_type),
               {

                 #Overall MLE
                 if(input$binary_estimator == "Overall MLE" &
                    input$endpoint_type == "Binary" &
                    input$binary_trial_type == "Single-arm"){
                   overall_mle_binary = b_ovr_mle_single(input$single_arm_events,
                                                  input$single_arm_sample_size_binary,
                                                  input$event_rate_null,
                                                  input$stage_reached,
                                                  input$maximum_stages)
                   estimator(overall_mle_binary)
                 }

                 #Stage 1 MLE
                 if(input$binary_estimator == "Stage 1 MLE" &
                    input$endpoint_type == "Binary" &
                    input$binary_trial_type == "Single-arm"){
                   stage_1_mle_binary = b_stg1_mle_single(input$single_arm_events,
                                                   input$single_arm_sample_size_binary,
                                                   input$event_rate_null,
                                                   input$stage_reached,
                                                   input$maximum_stages)
                   estimator(stage_1_mle_binary)
                 }

                 #MUE
                 if(input$binary_estimator == "MUE" &
                    input$endpoint_type == "Binary" &
                    input$binary_trial_type == "Single-arm"){
                   mue_binary = b_mue_single(input$single_arm_events,
                                      input$single_arm_sample_size_binary,
                                      input$event_rate_null,
                                      input$lower_bound,
                                      input$upper_bound,
                                      input$binary_search_region_uniroot,
                                      input$stage_reached,
                                      input$tolerance,
                                      input$maximum_iterations,
                                      input$tail_type,
                                      input$maximum_sub_intervals,
                                      input$maximum_stages)
                   estimator(mue_binary)
                 }

                 #UMVUE
                 if(input$binary_estimator == "UMVUE" &
                    input$endpoint_type == "Binary" &
                    input$binary_trial_type == "Single-arm"){
                   umvue_binary = b_UMVUE_single(input$single_arm_events,
                                          input$single_arm_sample_size_binary,
                                          input$event_rate_null,
                                          input$lower_bound,
                                          input$upper_bound,
                                          input$epsilon,
                                          input$absolute_tolerance,
                                          input$relative_tolerance,
                                          input$stage_reached,
                                          input$subdivisions,
                                          input$maximum_stages)
                   estimator(umvue_binary)
                 }

                 #ubc-MLE
                 if(input$binary_estimator == "ubc-MLE" &
                    input$endpoint_type == "Binary" &
                    input$binary_trial_type == "Single-arm"){
                   ubc_mle_binary = b_ubc_MLE_single(input$single_arm_events,
                                              input$single_arm_sample_size_binary,
                                              input$event_rate_null,
                                              input$lower_bound,
                                              input$upper_bound,
                                              input$epsilon,
                                              input$absolute_tolerance,
                                              input$relative_tolerance,
                                              input$binary_search_region_optimize,
                                              input$stage_reached,
                                              input$maximum_stages)
                   estimator(ubc_mle_binary)
                 }

                 #cMLE
                 if(input$binary_estimator == "cMLE" &
                    input$endpoint_type == "Binary" &
                    input$binary_trial_type == "Single-arm"){
                   conditional_mle_binary = b_cMLE_single(input$single_arm_events,
                                                   input$single_arm_sample_size_binary,
                                                   input$event_rate_null,
                                                   input$stage_reached,
                                                   input$conditional_stop_stage,
                                                   input$maximum_stages)
                   estimator(conditional_mle_binary)
                 }

                 #cMUE
                 if(input$binary_estimator == "cMUE" &
                    input$endpoint_type == "Binary" &
                    input$binary_trial_type == "Single-arm"){
                   cMUE_binary = b_cMUE_single(input$single_arm_events,
                                        input$single_arm_sample_size_binary,
                                        input$event_rate_null,
                                        input$lower_bound,
                                        input$upper_bound,
                                        input$binary_search_region_uniroot,
                                        input$stage_reached,
                                        input$tolerance,
                                        input$maximum_iterations,
                                        input$tail_type,
                                        input$maximum_sub_intervals,
                                        input$maximum_stages)
                   estimator(cMUE_binary)
                 }

                 #UMVCUE
                 if(input$binary_estimator == "UMVCUE" &
                    input$endpoint_type == "Binary" &
                    input$binary_trial_type == "Single-arm"){
                   umvcue_binary = b_umvcue_single(input$single_arm_events,
                                            input$single_arm_sample_size_binary,
                                            input$event_rate_null,
                                            input$lower_bound,
                                            input$upper_bound,
                                            input$epsilon,
                                            input$absolute_tolerance,
                                            input$relative_tolerance,
                                            input$stage_reached,
                                            input$subdivisions,
                                            input$maximum_stages)
                   estimator(umvcue_binary)
                 }

                 #cbc-MLE
                 if(input$binary_estimator == "cbc-MLE"&
                    input$endpoint_type == "Binary" &
                    input$binary_trial_type == "Single-arm"){
                   cbc_mle_binary = b_cbc_MLE_single(input$single_arm_events,
                                              input$single_arm_sample_size_binary,
                                              input$event_rate_null,
                                              input$lower_bound,
                                              input$upper_bound,
                                              input$epsilon,
                                              input$absolute_tolerance,
                                              input$relative_tolerance,
                                              input$binary_search_region_optimize,
                                              input$stage_reached,
                                              input$conditional_stop_stage,
                                              input$maximum_stages)
                   estimator(cbc_mle_binary)
                 }
               }
  )





#OBSERVER FUNCTION FOR THE BINARY (TWO-ARM) ESTIMATORS -------------------------
  observeEvent(c(input$calc_butt, input$endpoint_type, input$binary_estimator,
                 input$binary_trial_type),
               {

#Overall MLE
                 if(input$binary_estimator == "Overall MLE" &
                    input$endpoint_type == "Binary" &
                    input$binary_trial_type == "Two-arm"){
                   overall_mle_binary = b_ovr_mle(input$control_events,
                                       input$experimental_events,
                                       input$control_sample_size_binary,
                                       input$experimental_sample_size_binary,
                                       input$stage_reached,
                                       input$maximum_stages)
                   estimator(overall_mle_binary)
                 }

#Stage 1 MLE
                  if(input$binary_estimator == "Stage 1 MLE" &
                         input$endpoint_type == "Binary" &
                     input$binary_trial_type == "Two-arm"){
                   stage_1_mle_binary = b_stg1_mle(input$control_events,
                                         input$experimental_events,
                                         input$control_sample_size_binary,
                                         input$experimental_sample_size_binary,
                                         input$stage_reached,
                                         input$maximum_stages)
                   estimator(stage_1_mle_binary)
                 }

#MUE
                  if(input$binary_estimator == "MUE" &
                         input$endpoint_type == "Binary" &
                     input$binary_trial_type == "Two-arm"){
                   mue_binary = b_mue(input$control_events,
                               input$experimental_events,
                               input$control_sample_size_binary,
                               input$experimental_sample_size_binary,
                               input$lower_bound,
                               input$upper_bound,
                               input$binary_search_region_uniroot,
                               input$stage_reached,
                               input$tolerance,
                               input$maximum_iterations,
                               input$tail_type,
                               input$maximum_sub_intervals,
                               input$maximum_stages)
                   estimator(mue_binary)
                 }

#UMVUE
                  if(input$binary_estimator == "UMVUE" &
                         input$endpoint_type == "Binary" &
                     input$binary_trial_type == "Two-arm"){
                   umvue_binary = b_UMVUE(input$control_events,
                                   input$experimental_events,
                                   input$control_sample_size_binary,
                                   input$experimental_sample_size_binary,
                                   input$lower_bound,
                                   input$upper_bound,
                                   input$epsilon,
                                   input$absolute_tolerance,
                                   input$relative_tolerance,
                                   input$stage_reached,
                                   input$subdivisions,
                                   input$maximum_stages)
                   estimator(umvue_binary)
                 }

#ubc-MLE
                  if(input$binary_estimator == "ubc-MLE" &
                         input$endpoint_type == "Binary" &
                     input$binary_trial_type == "Two-arm"){
                   ubc_mle_binary = b_ubc_MLE(input$control_events,
                                       input$experimental_events,
                                       input$control_sample_size_binary,
                                       input$experimental_sample_size_binary,
                                       input$lower_bound,
                                       input$upper_bound,
                                       input$epsilon,
                                       input$absolute_tolerance,
                                       input$relative_tolerance,
                                       input$binary_search_region_optimize,
                                       input$stage_reached,
                                       input$maximum_stages)
                   estimator(ubc_mle_binary)
                 }

#cMLE
                  if(input$binary_estimator == "cMLE" &
                         input$endpoint_type == "Binary" &
                     input$binary_trial_type == "Two-arm"){
                   conditional_mle_binary = b_cMLE(input$control_events,
                                  input$experimental_events,
                                  input$control_sample_size_binary,
                                  input$experimental_sample_size_binary,
                                  input$stage_reached,
                                  input$conditional_stop_stage,
                                  input$maximum_stages)
                   estimator(conditional_mle_binary)
                 }

#cMUE
                  if(input$binary_estimator == "cMUE" &
                         input$endpoint_type == "Binary" &
                     input$binary_trial_type == "Two-arm"){
                  cMUE_binary = b_cMUE(input$control_events,
                                input$experimental_events,
                                input$control_sample_size_binary,
                                input$experimental_sample_size_binary,
                                input$lower_bound,
                                input$upper_bound,
                                input$binary_search_region_uniroot,
                                input$stage_reached,
                                input$tolerance,
                                input$maximum_iterations,
                                input$tail_type,
                                input$maximum_sub_intervals,
                                input$maximum_stages)
                  estimator(cMUE_binary)
                 }

 #UMVCUE
                  if(input$binary_estimator == "UMVCUE" &
                         input$endpoint_type == "Binary" &
                     input$binary_trial_type == "Two-arm"){
                   umvcue_binary = b_umvcue(input$control_events,
                                     input$experimental_events,
                                     input$control_sample_size_binary,
                                     input$experimental_sample_size_binary,
                                     input$lower_bound,
                                     input$upper_bound,
                                     input$epsilon,
                                     input$absolute_tolerance,
                                     input$relative_tolerance,
                                     input$stage_reached,
                                     input$subdivisions,
                                     input$maximum_stages)
                   estimator(umvcue_binary)
                 }

#cbc-MLE
                  if(input$binary_estimator == "cbc-MLE"&
                          input$endpoint_type == "Binary" &
                     input$binary_trial_type == "Two-arm"){
                   cbc_mle_binary = b_cbc_MLE(input$control_events,
                                       input$experimental_events,
                                       input$control_sample_size_binary,
                                       input$experimental_sample_size_binary,
                                       input$lower_bound,
                                       input$upper_bound,
                                       input$epsilon,
                                       input$absolute_tolerance,
                                       input$relative_tolerance,
                                       input$binary_search_region_optimize,
                                       input$stage_reached,
                                       input$conditional_stop_stage,
                                       input$maximum_stages)
                   estimator(cbc_mle_binary)
                 }
               }
 )




#OBSERVER FUNCTION FOR THE NORMAL (SINGLE-ARM) ESTIMATORS ----------------------
  observeEvent(c(input$calc_butt, input$endpoint_type, input$normal_estimator,
                 input$normal_trial_type),
               {

                 #Overall MLE
                 if(input$normal_estimator == "Overall MLE" &
                    input$endpoint_type == "Normal" &
                    input$normal_trial_type == "Single-arm"){
                   overall_mle = n_ovr_mle_single_arm(input$means_normal_single_arm,
                                                    input$mean_null_normal_single_arm,
                                                    input$sample_size_normal_single_arm,
                                                    input$variance_normal_single_arm,
                                                    input$stage_reached,
                                                    input$maximum_stages)
                   estimator(overall_mle)
                 }

                 #Stage 1 MLE
                 if(input$normal_estimator == "Stage 1 MLE" &
                    input$endpoint_type == "Normal" &
                    input$normal_trial_type == "Single-arm"){
                   stage_1_mle = n_stage1_mle_single_arm(input$means_normal_single_arm,
                                                        input$mean_null_normal_single_arm,
                                                        input$sample_size_normal_single_arm,
                                                        input$variance_normal_single_arm,
                                                        input$stage_reached,
                                                        input$maximum_stages)
                   estimator(stage_1_mle)
                 }

                 #MUE
                 if(input$normal_estimator == "MUE" &
                    input$endpoint_type == "Normal" &
                    input$normal_trial_type == "Single-arm"){
                   parallel_mue = n_mue_single_arm(input$means_normal_single_arm,
                                                 input$mean_null_normal_single_arm,
                                                 input$sample_size_normal_single_arm,
                                                 input$variance_normal_single_arm,
                                                 input$lower_bound,
                                                 input$upper_bound,
                                                 input$normal_search_region_uniroot,
                                                 input$stage_reached,
                                                 input$tolerance,
                                                 input$maximum_iterations,
                                                 input$tail_type,
                                                 input$maximum_sub_intervals,
                                                 input$maximum_stages)
                   estimator(parallel_mue)
                 }

                 #UMVUE
                 if(input$normal_estimator == "UMVUE" &
                    input$endpoint_type == "Normal" &
                    input$normal_trial_type == "Single-arm"){
                   umvue = n_umvue_single_arm(input$means_normal_single_arm,
                                            input$mean_null_normal_single_arm,
                                            input$sample_size_normal_single_arm,
                                            input$variance_normal_single_arm,
                                            input$lower_bound,
                                            input$upper_bound,
                                            input$epsilon,
                                            input$absolute_tolerance,
                                            input$relative_tolerance,
                                            input$stage_reached,
                                            input$subdivisions,
                                            input$maximum_stages)
                   estimator(umvue)
                 }

                 #ubc-MLE
                 if(input$normal_estimator == "ubc-MLE" &
                    input$endpoint_type == "Normal" &
                    input$normal_trial_type == "Single-arm"){
                   ubc_mle_normal = n_ubc_mle_single_arm(input$means_normal_single_arm,
                                                       input$mean_null_normal_single_arm,
                                                       input$sample_size_normal_single_arm,
                                                       input$variance_normal_single_arm,
                                                       input$lower_bound,
                                                       input$upper_bound,
                                                       input$epsilon,
                                                       input$absolute_tolerance,
                                                       input$relative_tolerance,
                                                       input$normal_search_region_optimize,
                                                       input$stage_reached,
                                                       input$maximum_stages)
                   estimator(ubc_mle_normal)
                 }

                 #cMLE
                 if(input$normal_estimator == "cMLE" &
                    input$endpoint_type == "Normal" &
                    input$normal_trial_type == "Single-arm"){
                   cMLE = n_cmle_single_arm(input$means_normal_single_arm,
                                          input$mean_null_normal_single_arm,
                                          input$sample_size_normal_single_arm,
                                          input$variance_normal_single_arm,
                                          input$stage_reached,
                                          input$conditional_stop_stage,
                                          input$maximum_stages)
                   estimator(cMLE)
                 }

                 #cMUE
                 if(input$normal_estimator == "cMUE" &
                    input$endpoint_type == "Normal" &
                    input$normal_trial_type == "Single-arm"){
                   cMUE = n_cmue_single_arm(input$means_normal_single_arm,
                                          input$mean_null_normal_single_arm,
                                          input$sample_size_normal_single_arm,
                                          input$variance_normal_single_arm,
                                          input$lower_bound,
                                          input$upper_bound,
                                          input$normal_search_region_uniroot,
                                          input$stage_reached,
                                          input$tolerance,
                                          input$maximum_iterations,
                                          input$tail_type,
                                          input$maximum_sub_intervals,
                                          input$maximum_stages)
                   estimator(cMUE)
                 }

                 #UMVCUE
                 if(input$normal_estimator == "UMVCUE" &
                    input$endpoint_type == "Normal" &
                    input$normal_trial_type == "Single-arm"){
                   umvcue = n_umvcue_single_arm(input$means_normal_single_arm,
                                              input$mean_null_normal_single_arm,
                                              input$sample_size_normal_single_arm,
                                              input$variance_normal_single_arm,
                                              input$lower_bound,
                                              input$upper_bound,
                                              input$epsilon,
                                              input$absolute_tolerance,
                                              input$relative_tolerance,
                                              input$stage_reached,
                                              input$subdivisions,
                                              input$maximum_stages)
                   estimator(umvcue)
                 }

                 #cbc-MLE
                 if(input$normal_estimator == "cbc-MLE" &
                    input$endpoint_type == "Normal" &
                    input$normal_trial_type == "Single-arm"){
                   cbc_mle = n_cbc_mle_single_arm(input$means_normal_single_arm,
                                                input$mean_null_normal_single_arm,
                                                input$sample_size_normal_single_arm,
                                                input$variance_normal_single_arm,
                                                input$lower_bound,
                                                input$upper_bound,
                                                input$epsilon,
                                                input$absolute_tolerance,
                                                input$relative_tolerance,
                                                input$normal_search_region_optimize,
                                                input$stage_reached,
                                                input$maximum_stages)
                   estimator(cbc_mle)
                 }
               }
  )




#OBSERVER FUNCTION FOR THE NORMAL (PARALLEL) ESTIMATORS ------------------------
  observeEvent(c(input$calc_butt, input$endpoint_type, input$normal_estimator,
                 input$normal_trial_type),
               {

#Overall MLE
                 if(input$normal_estimator == "Overall MLE" &
                    input$endpoint_type == "Normal" &
                    input$normal_trial_type == "Parallel"){
                    overall_mle = n_ovr_mle_parallel(input$parallel_control_arm_means,
                                                     input$parallel_experimental_arm_means,
                                                     input$parallel_control_sample_size_normal,
                                                     input$parallel_experimental_sample_size_normal,
                                                     input$parallel_variance_control,
                                                     input$parallel_variance_experimental,
                                                     input$stage_reached,
                                                     input$maximum_stages)
                   estimator(overall_mle)
                 }

#Stage 1 MLE
                  if(input$normal_estimator == "Stage 1 MLE" &
                     input$endpoint_type == "Normal" &
                     input$normal_trial_type == "Parallel"){
                     stage_1_mle = n_stage1_mle_parallel(input$parallel_control_arm_means,
                                                         input$parallel_experimental_arm_means,
                                                         input$parallel_control_sample_size_normal,
                                                         input$parallel_experimental_sample_size_normal,
                                                         input$parallel_variance_control,
                                                         input$parallel_variance_experimental,
                                                         input$stage_reached,
                                                         input$maximum_stages)
                     estimator(stage_1_mle)
                 }

#MUE
                  if(input$normal_estimator == "MUE" &
                     input$endpoint_type == "Normal" &
                     input$normal_trial_type == "Parallel"){
                     parallel_mue = n_mue_parallel(input$parallel_control_arm_means,
                                          input$parallel_experimental_arm_means,
                                          input$parallel_control_sample_size_normal,
                                          input$parallel_experimental_sample_size_normal,
                                          input$parallel_variance_control,
                                          input$parallel_variance_experimental,
                                          input$lower_bound,
                                          input$upper_bound,
                                          input$normal_search_region_uniroot,
                                          input$stage_reached,
                                          input$tolerance,
                                          input$maximum_iterations,
                                          input$tail_type,
                                          input$maximum_sub_intervals,
                                          input$maximum_stages)
                   estimator(parallel_mue)
                 }

#UMVUE
                  if(input$normal_estimator == "UMVUE" &
                     input$endpoint_type == "Normal" &
                     input$normal_trial_type == "Parallel"){
                   umvue = n_umvue_parallel(input$parallel_control_arm_means,
                                            input$parallel_experimental_arm_means,
                                            input$parallel_control_sample_size_normal,
                                            input$parallel_experimental_sample_size_normal,
                                            input$parallel_variance_control,
                                            input$parallel_variance_experimental,
                                            input$lower_bound,
                                            input$upper_bound,
                                            input$epsilon,
                                            input$absolute_tolerance,
                                            input$relative_tolerance,
                                            input$stage_reached,
                                            input$subdivisions,
                                            input$maximum_stages)
                   estimator(umvue)
                 }

#ubc-MLE
                  if(input$normal_estimator == "ubc-MLE" &
                     input$endpoint_type == "Normal" &
                     input$normal_trial_type == "Parallel"){
                   ubc_mle_normal = n_ubc_mle_parallel(input$parallel_control_arm_means,
                                                input$parallel_experimental_arm_means,
                                                input$parallel_control_sample_size_normal,
                                                input$parallel_experimental_sample_size_normal,
                                                input$parallel_variance_control,
                                                input$parallel_variance_experimental,
                                                input$lower_bound,
                                                input$upper_bound,
                                                input$epsilon,
                                                input$absolute_tolerance,
                                                input$relative_tolerance,
                                                input$normal_search_region_optimize,
                                                input$stage_reached,
                                                input$maximum_stages)
                   estimator(ubc_mle_normal)
                 }

#cMLE
                  if(input$normal_estimator == "cMLE" &
                     input$endpoint_type == "Normal" &
                     input$normal_trial_type == "Parallel"){
                   cMLE = n_cmle_parallel(input$parallel_control_arm_means,
                                          input$parallel_experimental_arm_means,
                                          input$parallel_control_sample_size_normal,
                                          input$parallel_experimental_sample_size_normal,
                                          input$parallel_variance_control,
                                          input$parallel_variance_experimental,
                                          input$stage_reached,
                                          input$conditional_stop_stage,
                                          input$maximum_stages)
                   estimator(cMLE)
                 }

#cMUE
                  if(input$normal_estimator == "cMUE" &
                     input$endpoint_type == "Normal" &
                     input$normal_trial_type == "Parallel"){
                   cMUE = n_cmue_parallel(input$parallel_control_arm_means,
                                          input$parallel_experimental_arm_means,
                                          input$parallel_control_sample_size_normal,
                                          input$parallel_experimental_sample_size_normal,
                                          input$parallel_variance_control,
                                          input$parallel_variance_experimental,
                                          input$lower_bound,
                                          input$upper_bound,
                                          input$normal_search_region_uniroot,
                                          input$stage_reached,
                                          input$tolerance,
                                          input$maximum_iterations,
                                          input$tail_type,
                                          input$maximum_sub_intervals,
                                          input$maximum_stages)
                   estimator(cMUE)
                 }

#UMVCUE
                  if(input$normal_estimator == "UMVCUE" &
                     input$endpoint_type == "Normal" &
                     input$normal_trial_type == "Parallel"){
                   umvcue = n_umvcue_parallel(input$parallel_control_arm_means,
                                              input$parallel_experimental_arm_means,
                                              input$parallel_control_sample_size_normal,
                                              input$parallel_experimental_sample_size_normal,
                                              input$parallel_variance_control,
                                              input$parallel_variance_experimental,
                                              input$lower_bound,
                                              input$upper_bound,
                                              input$epsilon,
                                              input$absolute_tolerance,
                                              input$relative_tolerance,
                                              input$stage_reached,
                                              input$subdivisions,
                                              input$maximum_stages)
                   estimator(umvcue)
                 }

#cbc-MLE
                 if(input$normal_estimator == "cbc-MLE" &
                    input$endpoint_type == "Normal" &
                    input$normal_trial_type == "Parallel"){
                   cbc_mle = n_cbc_mle_parallel(input$parallel_control_arm_means,
                                                input$parallel_experimental_arm_means,
                                                input$parallel_control_sample_size_normal,
                                                input$parallel_experimental_sample_size_normal,
                                                input$parallel_variance_control,
                                                input$parallel_variance_experimental,
                                                input$lower_bound,
                                                input$upper_bound,
                                                input$epsilon,
                                                input$absolute_tolerance,
                                                input$relative_tolerance,
                                                input$normal_search_region_optimize,
                                                input$stage_reached,
                                                input$maximum_stages)
                   estimator(cbc_mle)
                 }
               }
 )





#OBSERVER FUNCTION FOR THE NORMAL (PAIRED) ESTIMATORS --------------------------
  observeEvent(c(input$calc_butt, input$endpoint_type, input$normal_estimator,
                 input$normal_trial_type),
               {

                 #Overall MLE
                 if(input$normal_estimator == "Overall MLE" &
                    input$endpoint_type == "Normal" &
                    input$normal_trial_type == "Paired"){
                   overall_mle = n_ovr_mle_paired(input$paired_experimental_arm_means,
                                                  input$parallel_control_arm_means,
                                                  input$paired_variance_normal,
                                                  input$paired_sample_size,
                                                  input$stage_reached,
                                                  input$maximum_stages)
                   estimator(overall_mle)
                 }

                 #Stage 1 MLE
                 if(input$normal_estimator == "Stage 1 MLE" &
                    input$endpoint_type == "Normal" &
                    input$normal_trial_type == "Paired"){
                   stage_1_mle = n_stage1_mle_paired(input$paired_experimental_arm_means,
                                                     input$parallel_control_arm_means,
                                                     input$paired_variance_normal,
                                                     input$paired_sample_size,
                                                     input$stage_reached,
                                                     input$maximum_stages)
                   estimator(stage_1_mle)
                 }

                 #MUE
                 if(input$normal_estimator == "MUE" &
                    input$endpoint_type == "Normal" &
                    input$normal_trial_type == "Paired"){
                   parallel_mue = n_mue_paired(input$paired_experimental_arm_means,
                                               input$parallel_control_arm_means,
                                               input$paired_variance_normal,
                                               input$paired_sample_size,
                                               input$lower_bound,
                                               input$upper_bound,
                                               input$normal_search_region_uniroot,
                                               input$stage_reached,
                                               input$tolerance,
                                               input$maximum_iterations,
                                               input$tail_type,
                                               input$maximum_sub_intervals,
                                               input$maximum_stages)
                   estimator(parallel_mue)
                 }

                 #UMVUE
                 if(input$normal_estimator == "UMVUE" &
                    input$endpoint_type == "Normal" &
                    input$normal_trial_type == "Paired"){
                   umvue = n_umvue_paired(input$paired_experimental_arm_means,
                                          input$parallel_control_arm_means,
                                          input$paired_variance_normal,
                                          input$paired_sample_size,
                                          input$lower_bound,
                                          input$upper_bound,
                                          input$epsilon,
                                          input$absolute_tolerance,
                                          input$relative_tolerance,
                                          input$stage_reached,
                                          input$subdivisions,
                                          input$maximum_stages)
                   estimator(umvue)
                 }

                 #ubc-MLE
                 if(input$normal_estimator == "ubc-MLE" &
                    input$endpoint_type == "Normal" &
                    input$normal_trial_type == "Paired"){
                   ubc_mle_normal = n_ubc_mle_paired(input$paired_experimental_arm_means,
                                                     input$parallel_control_arm_means,
                                                     input$paired_variance_normal,
                                                     input$paired_sample_size,
                                                     input$lower_bound,
                                                     input$upper_bound,
                                                     input$epsilon,
                                                     input$absolute_tolerance,
                                                     input$relative_tolerance,
                                                     input$normal_search_region_optimize,
                                                     input$stage_reached,
                                                     input$maximum_stages)
                   estimator(ubc_mle_normal)
                 }

                 #cMLE
                 if(input$normal_estimator == "cMLE" &
                    input$endpoint_type == "Normal" &
                    input$normal_trial_type == "Paired"){
                   cMLE = n_cmle_paired(input$paired_experimental_arm_means,
                                        input$parallel_control_arm_means,
                                        input$paired_variance_normal,
                                        input$paired_sample_size,
                                        input$stage_reached,
                                        input$conditional_stop_stage,
                                        input$maximum_stages)
                   estimator(cMLE)
                 }

                 #cMUE
                 if(input$normal_estimator == "cMUE" &
                    input$endpoint_type == "Normal" &
                    input$normal_trial_type == "Paired"){
                   cMUE = n_cmue_paired(input$paired_experimental_arm_means,
                                        input$parallel_control_arm_means,
                                        input$paired_variance_normal,
                                        input$paired_sample_size,
                                        input$lower_bound,
                                        input$upper_bound,
                                        input$normal_search_region_uniroot,
                                        input$stage_reached,
                                        input$tolerance,
                                        input$maximum_iterations,
                                        input$tail_type,
                                        input$maximum_sub_intervals,
                                        input$maximum_stages)
                   estimator(cMUE)
                 }

                 #UMVCUE
                 if(input$normal_estimator == "UMVCUE" &
                    input$endpoint_type == "Normal" &
                    input$normal_trial_type == "Paired"){
                   umvcue = n_umvcue_paired(input$paired_experimental_arm_means,
                                            input$parallel_control_arm_means,
                                            input$paired_variance_normal,
                                            input$paired_sample_size,
                                            input$lower_bound,
                                            input$upper_bound,
                                            input$epsilon,
                                            input$absolute_tolerance,
                                            input$relative_tolerance,
                                            input$stage_reached,
                                            input$subdivisions,
                                            input$maximum_stages)
                   estimator(umvcue)
                 }

                 #cbc-MLE
                 if(input$normal_estimator == "cbc-MLE" &
                    input$endpoint_type == "Normal" &
                    input$normal_trial_type == "Paired"){
                   cbc_mle = n_cbc_mle_paired(input$paired_experimental_arm_means,
                                              input$parallel_control_arm_means,
                                              input$paired_variance_normal,
                                              input$paired_sample_size,
                                              input$lower_bound,
                                              input$upper_bound,
                                              input$epsilon,
                                              input$absolute_tolerance,
                                              input$relative_tolerance,
                                              input$normal_search_region_optimize,
                                              input$stage_reached,
                                              input$maximum_stages)
                   estimator(cbc_mle)
                 }
               }
  )



#OBSERVER FUNCTION FOR THE NORMAL (CROSSOVER) ESTIMATORS -----------------------
  observeEvent(c(input$calc_butt, input$endpoint_type, input$normal_estimator,
                 input$normal_trial_type),
               {

                 #Overall MLE
                 if(input$normal_estimator == "Overall MLE" &
                    input$endpoint_type == "Normal" &
                    input$normal_trial_type == "Two Period Crossover"){
                   overall_mle = n_ovr_mle_crossover(input$crossover_means_ab,
                                                     input$crossover_means_ba,
                                                     input$variance_crossover_ab,
                                                     input$variance_crossover_ba,
                                                     input$sample_size_ab_crossover,
                                                     input$sample_size_ba_crossover,
                                                     input$stage_reached,
                                                     input$maximum_stages)
                   estimator(overall_mle)
                 }

                 #Stage 1 MLE
                 if(input$normal_estimator == "Stage 1 MLE" &
                    input$endpoint_type == "Normal" &
                    input$normal_trial_type == "Two Period Crossover"){
                   stage_1_mle = n_stage1_mle_crossover(input$crossover_means_ab,
                                                        input$crossover_means_ba,
                                                        input$variance_crossover_ab,
                                                        input$variance_crossover_ba,
                                                        input$sample_size_ab_crossover,
                                                        input$sample_size_ba_crossover,
                                                        input$stage_reached,
                                                        input$maximum_stages)
                   estimator(stage_1_mle)
                 }

                 #MUE
                 if(input$normal_estimator == "MUE" &
                    input$endpoint_type == "Normal" &
                    input$normal_trial_type == "Two Period Crossover"){
                    mue = n_mue_crossover(input$crossover_means_ab,
                                          input$crossover_means_ba,
                                          input$variance_crossover_ab,
                                          input$variance_crossover_ba,
                                          input$sample_size_ab_crossover,
                                          input$sample_size_ba_crossover,
                                          input$lower_bound,
                                          input$upper_bound,
                                          input$normal_search_region_uniroot,
                                          input$stage_reached,
                                          input$tolerance,
                                          input$maximum_iterations,
                                          input$tail_type,
                                          input$maximum_sub_intervals,
                                          input$maximum_stages)
                   estimator(mue)
                 }

                 #UMVUE
                 if(input$normal_estimator == "UMVUE" &
                    input$endpoint_type == "Normal" &
                    input$normal_trial_type == "Two Period Crossover"){
                   umvue = n_umvue_crossover(input$crossover_means_ab,
                                             input$crossover_means_ba,
                                             input$variance_crossover_ab,
                                             input$variance_crossover_ba,
                                             input$sample_size_ab_crossover,
                                             input$sample_size_ba_crossover,
                                             input$lower_bound,
                                             input$upper_bound,
                                             input$epsilon,
                                             input$absolute_tolerance,
                                             input$relative_tolerance,
                                             input$stage_reached,
                                             input$subdivisions,
                                             input$maximum_stages)
                   estimator(umvue)
                 }

                 #ubc-MLE
                 if(input$normal_estimator == "ubc-MLE" &
                    input$endpoint_type == "Normal" &
                    input$normal_trial_type == "Two Period Crossover"){
                   ubc_mle = n_ubc_mle_crossover(input$crossover_means_ab,
                                                        input$crossover_means_ba,
                                                 input$variance_crossover_ab,
                                                 input$variance_crossover_ba,
                                                        input$sample_size_ab_crossover,
                                                        input$sample_size_ba_crossover,
                                                        input$lower_bound,
                                                        input$upper_bound,
                                                        input$epsilon,
                                                        input$absolute_tolerance,
                                                        input$relative_tolerance,
                                                        input$normal_search_region_optimize,
                                                        input$stage_reached,
                                                        input$maximum_stages)
                   estimator(ubc_mle)
                 }

                 #cMLE
                 if(input$normal_estimator == "cMLE" &
                    input$endpoint_type == "Normal" &
                    input$normal_trial_type == "Two Period Crossover"){
                   cMLE = n_cmle_crossover(input$crossover_means_ab,
                                           input$crossover_means_ba,
                                           input$variance_crossover_ab,
                                           input$variance_crossover_ba,
                                           input$sample_size_ab_crossover,
                                           input$sample_size_ba_crossover,
                                           input$stage_reached,
                                           input$conditional_stop_stage,
                                           input$maximum_stages)
                   estimator(cMLE)
                 }

                 #cMUE
                 if(input$normal_estimator == "cMUE" &
                    input$endpoint_type == "Normal" &
                    input$normal_trial_type == "Two Period Crossover"){
                   cMUE = n_cmue_crossover(input$crossover_means_ab,
                                           input$crossover_means_ba,
                                           input$variance_crossover_ab,
                                           input$variance_crossover_ba,
                                           input$sample_size_ab_crossover,
                                           input$sample_size_ba_crossover,
                                           input$lower_bound,
                                           input$upper_bound,
                                           input$normal_search_region_uniroot,
                                           input$stage_reached,
                                           input$tolerance,
                                           input$maximum_iterations,
                                           input$tail_type,
                                           input$maximum_sub_intervals,
                                           input$maximum_stages)
                   estimator(cMUE)
                 }

                 #UMVCUE
                 if(input$normal_estimator == "UMVCUE" &
                    input$endpoint_type == "Normal" &
                    input$normal_trial_type == "Two Period Crossover"){
                   umvcue = n_umvcue_crossover(input$crossover_means_ab,
                                               input$crossover_means_ba,
                                               input$variance_crossover_ab,
                                               input$variance_crossover_ba,
                                               input$sample_size_ab_crossover,
                                               input$sample_size_ba_crossover,
                                               input$lower_bound,
                                               input$upper_bound,
                                               input$epsilon,
                                               input$absolute_tolerance,
                                               input$relative_tolerance,
                                               input$stage_reached,
                                               input$subdivisions,
                                               input$maximum_stages)
                   estimator(umvcue)
                 }

                 #cbc-MLE
                 if(input$normal_estimator == "cbc-MLE" &
                    input$endpoint_type == "Normal" &
                    input$normal_trial_type == "Two Period Crossover"){
                   cbc_mle = n_cbc_mle_crossover(input$crossover_means_ab,
                                                 input$crossover_means_ba,
                                                 input$variance_crossover_ab,
                                                 input$variance_crossover_ba,
                                                 input$sample_size_ab_crossover,
                                                 input$sample_size_ba_crossover,
                                                 input$lower_bound,
                                                 input$upper_bound,
                                                 input$epsilon,
                                                 input$absolute_tolerance,
                                                 input$relative_tolerance,
                                                 input$normal_search_region_optimize,
                                                 input$stage_reached,
                                                 input$maximum_stages)
                   estimator(cbc_mle)
                 }
               }
  )




#OBSERVER FUNCTION FOR THE TIME TO EVENT ESTIMATORS ----------------------------
  observeEvent(c(input$calc_butt, input$endpoint_type,
                 input$time_to_event_estimator),
               {

#Overall MLE
                 if(input$time_to_event_estimator == "Overall MLE" &
                    input$endpoint_type == "Time to event" &
                    input$time_to_event_trial_type == "Two-arm"){

                 }

#Stage 1 MLE
                  if(input$time_to_event_estimator == "Stage 1 MLE" &
                     input$endpoint_type == "Time to event" &
                     input$time_to_event_trial_type == "Two-arm"){

                 }

#MUE
                  if(input$time_to_event_estimator == "MUE" &
                     input$endpoint_type == "Time to event" &
                     input$time_to_event_trial_type == "Two-arm"){

                 }

#UMVUE
                  if(input$time_to_event_estimator == "UMVUE" &
                     input$endpoint_type == "Time to event" &
                     input$time_to_event_trial_type == "Two-arm"){

                 }

#ubc-MLE
                  if(input$time_to_event_estimator == "ubc-MLE" &
                     input$endpoint_type == "Time to event" &
                     input$time_to_event_trial_type == "Two-arm"){

                 }

#cMLE
                  if(input$time_to_event_estimator == "cMLE" &
                     input$endpoint_type == "Time to event" &
                     input$time_to_event_trial_type == "Two-arm"){

                 }

#cMUE
                  if(input$time_to_event_estimator == "cMUE" &
                     input$endpoint_type == "Time to event" &
                     input$time_to_event_trial_type == "Two-arm"){

                 }

#UMVCUE
                  if(input$time_to_event_estimator == "UMVCUE" &
                     input$endpoint_type == "Time to event" &
                     input$time_to_event_trial_type == "Two-arm"){

                 }

#cbc-MLE
                 if(input$time_to_event_estimator == "cbc-MLE" &
                    input$endpoint_type == "Time to event" &
                    input$time_to_event_trial_type == "Two-arm"){

                 }
               }
  )




#Render output text (either estimate or any errors to be handled) --------------
  output$text <- renderUI({

#req statement to ensure reactive value estimator() is set ---------------------
    req(estimator())

#Locally set result to be estimator() ------------------------------------------
    result <- estimator()

#Begin loop to detect whether estimator should be outputted or error list ------
    if (is.numeric(result)) {

#Check to see if there is a single numeric result ------------------------------
      if (length(result) == 1) {

#print single numeric result ---------------------------------------------------
        if(input$endpoint_type == "Binary"){
        return(p(paste0("The ", input$binary_estimator,
                        " for the ", input$endpoint_type,
                        " endpoint is ", round(result, 5))))
        }

        if(input$endpoint_type == "Normal"){
          return(p(paste0("The ", input$normal_estimator,
                          " for the ", input$endpoint_type,
                          " endpoint is ", round(result, 5))))
        }
      }

#More than one numeric result so print multiple values -------------------------
      else {
        result_lines <- paste0("The ", input$endpoint_type," result ",
                               seq_along(result), " is: ", round(result, 5))
        return(tagList(lapply(result_lines, function(line) p(line))))
      }
    }

#If result isn't numeric then print errors -------------------------------------
    else {

#Calculate number of errors detected -------------------------------------------
      error_num <- length(result)

#Print list and total number of errors -----------------------------------------
      result <- c(paste0("You have ", error_num, " Error(s):"), result)
      return(tagList(lapply(result, function(line) p(line))))
    }
  })
}




#CREATE SHINY APP---------------------------------------------------------------
shinyApp(ui = ui, server = server)
