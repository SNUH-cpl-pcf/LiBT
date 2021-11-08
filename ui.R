library(shiny)
library(shinyjs)
library(shinydashboard)
library(shinydashboardPlus)
library(shinyWidgets)
library(shinyalert)
library(shinycssloaders)
library(shinyjqui)
library(shinybusy)
library(shinyBS)
library(httr)
library(sortable)
library(DT)

ui <- function(request) {shinyUI(
  dashboardPage(
    header = dashboardHeader(
      title = "LiBT-Analyst",
      titleWidth = 300
    ),
    sidebar = dashboardSidebar(
      width=300,
      sidebarMenu(
        menuItem("User Guide", selected = F, icon=icon("info"),tabName = "user_gide",
                 menuSubItem("Tutorial", tabName = "tutorial")),
        menuItem("Analysis", selected = T, icon=icon("chart-pie"),tabName = "analysis"),
        menuItem("TimeLine", selected = T, icon=icon("history"),
                 tabName = "timeline", uiOutput("timeline"))
      )
    ), # End of dashboardSidebar
    body = dashboardBody(
      tags$head(tags$link(rel="stylesheet",type="text/css",href="body.css"),
                tags$link(rel="stylesheet",type="text/css",href="rightsidebar.css"),
                tags$link(rel="stylesheet",type="text/css",href="timeline.css")),
      tags$head(tags$script(src="body.js")),
      tags$style(HTML("
                      ul{
                        list-style:none;
                      }
                      ")),
      useShinyjs(),
      tabItems(
        tabItem(
          tabName = "analysis",
          fluidRow(
            hidden(textInput("text", "Text")),
            hidden(downloadButton("download_log2_exp_btn","Save log2-transformed Expression Data", "download_log2_exp_btn")),
            hidden(downloadButton("download_preprocessed_exp_btn","Save preprocessed Expression Data", "download_preprocessed_exp_btn")),
            hidden(downloadButton("download_dep_info_btn","Save DEP info Data","download_dep_info_btn")),
            box(
              id="data_table",
              solidHeader = T,
              width = 12,
              use_busy_spinner(spin = "flower",
                               color = "#112446",
                               position = c("top-right"),
                               margins = c("50%", "50%"),
                               spin_id = NULL,
                               height = "50px",
                               width = "50px"),
              DT::dataTableOutput("uploaded_file_header")
            ), # End of uploaded file data table
            tabBox(
              id="plot_tabBox",
              width = 12,
              tabPanel("Results of DEA",
                       fluidRow(
                         use_busy_spinner(spin = "flower",
                                          color = "#112446",
                                          position = c("top-right"),
                                          margins = c("50%", "50%"),
                                          spin_id = NULL,
                                          height = "50px",
                                          width = "50px"),
                         box(
                           id = "pca_plot_box",
                           solidHeader = T,
                           width = 6,
                           hidden(prettySwitch(inputId = "show_sampleID", label = "Show SampleID", value = F, width=2, status = "success", inline=T)),
                           hidden(plotOutput("pca_plot")),
                           hidden(downloadButton("download_pca", "Save PCA Result"))
                         ),
                         box(
                           id = "correlation_matrix_box",
                           solidHeader = T,
                           width = 6,
                           br(),br(),
                           hidden(plotOutput("correlation_matrix")),
                           hidden(downloadButton("download_correlation", "Save Correlation Result"))
                         ),
                         br(),
                         box(
                           id = "volcano_box",
                           solidHeader = T,
                           width = 6,
                           hidden(plotOutput("volcano_plot", brush = "volcano_brush")),
                           DT::dataTableOutput("volcano_info"),
                           hidden(downloadButton("download_volcano", "Save Volcano Result","download_volcano"))
                         ),
                         box(
                           id = "heatmap_box",
                           solidHeader = T,
                           width = 6,
                           hidden(plotOutput("heatmap")),
                           hidden(downloadButton("download_heatmap", "Save Heatmap Result")),
                           hidden(downloadButton("download_gene_cluster","Save Clustering Result"))
                         )
                         ,useShinyalert()
                       )
              ),#tab1 end
              tabPanel("Results of GSA",
                       fluidRow(
                         use_busy_spinner(spin = "flower",
                                          color = "#112446",
                                          position = c("top-right"),
                                          margins = c("50%", "50%"),
                                          spin_id = NULL,
                                          height = "50px",
                                          width = "50px"),
                         box(
                           id = "gobp_gsa_box",
                           solidHeader = T,
                           width = 12,
                           #withSpinner(plotOutput("gobp_gsa_plot")),
                           hidden(plotOutput("gobp_gsa_plot")),
                           hidden(downloadButton("download_gsa_gobp_csv", "Save Table Result")),
                           hidden(downloadButton("download_gsa_gobp_png", "Save Plot Result"))
                         ),
                         box(
                           id = "gocc_gsa_box",
                           solidHeader = T,
                           width = 12,
                           #withSpinner(plotOutput("gocc_gsa_plot")),
                           hidden(plotOutput("gocc_gsa_plot")),
                           hidden(downloadButton("download_gsa_gocc_csv", "Save Table Result")),
                           hidden(downloadButton("download_gsa_gocc_png", "Save Plot Result"))
                         ),
                         box(
                           id = "gomf_gsa_box",
                           solidHeader = T,
                           width = 12,
                           #withSpinner(plotOutput("gomf_gsa_plot")),
                           hidden(plotOutput("gomf_gsa_plot")),
                           hidden(downloadButton("download_gsa_gomf_csv", "Save Table Result")),
                           hidden(downloadButton("download_gsa_gomf_png", "Save Plot Result"))
                         ),
                         box(
                           id = "kegg_gsa_box",
                           solidHeader = T,
                           width = 12,
                           #withSpinner(plotOutput("kegg_gsa_plot")),
                           hidden(plotOutput("kegg_gsa_plot")),
                           hidden(downloadButton("download_gsa_kegg_csv", "Save Table Result")),
                           hidden(downloadButton("download_gsa_kegg_png", "Save Plot Result"))
                         )
                         ,useShinyalert()
                       )
              ),#tab2 end
              tabPanel("Results of GSA Pathview",
                       fluidRow(
                         use_busy_spinner(spin = "flower",
                                          color = "#112446",
                                          position = c("top-right"),
                                          margins = c("50%", "50%"),
                                          spin_id = NULL,
                                          height = "50px",
                                          width = "50px"),
                         box(
                           title="Top 10 of KEGG pathway",
                           id = "top10_pathview_box",
                           solidHeader = T,
                           width = 12,
                           hidden(DT::dataTableOutput("topOfKeggDT"))
                         ),
                         box(
                           id = "pathview_select_box",
                           solidHeader = T,
                           width=6,
                           selectInput("pathID_selector", label="Choose a pathway ID",
                                       choices=""),
                           actionButton("render_pathway_btn", "Render pathway"),
                           actionButton("zoom_pathway_btn", "Zoom"),
                           hidden(downloadButton("download_pathview", "Save GSA Pathview zip"))
                         ),
                         box(
                           id = "pathview_result_box",
                           solidHeader = T,
                           width = 6,
                           hidden(plotOutput("pathview_result"))
                         )
                         ,useShinyalert()
                       )
              ),#tab3 end
              tabPanel("Results of GSEA",
                       fluidRow(
                         use_busy_spinner(spin = "flower",
                                          color = "#112446",
                                          position = c("top-right"),
                                          margins = c("50%", "50%"),
                                          spin_id = NULL,
                                          height = "50px",
                                          width = "50px"),
                         box(
                           id = "gobp_gsea_box",
                           solidHeader = T,
                           width = 12,
                           #withSpinner(plotOutput("gobp_gsea_plot")),
                           bsPopover("gobp_gsea_plot","Tip!","When click the red or blue point, you can show the enrichment plot of term",placement = "top"),
                           hidden(plotOutput("gobp_gsea_plot", click="gsea_gobp_click")),
                           hidden(downloadButton("download_gsea_gobp_csv", "Save Table Result")),
                           hidden(downloadButton("download_gsea_gobp_png", "Save Plot Result"))
                         ),
                         box(
                           id = "gocc_gsea_box",
                           solidHeader = T,
                           width = 12,
                           #withSpinner(plotOutput("gocc_gsea_plot")),
                           bsPopover("gocc_gsea_plot","Tip!","When click the red or blue point, you can show the enrichment plot of term",placement = "top"),
                           hidden(plotOutput("gocc_gsea_plot", click="gsea_gocc_click")),
                           hidden(downloadButton("download_gsea_gocc_csv", "Save Table Result")),
                           hidden(downloadButton("download_gsea_gocc_png", "Save Plot Result"))
                         ),
                         box(
                           id = "gomf_gsea_box",
                           solidHeader = T,
                           width = 12,
                           #withSpinner(plotOutput("gomf_gsea_plot")),
                           bsPopover("gomf_gsea_plot","Tip!","When click the red or blue point, you can show the enrichment plot of term",placement = "top"),
                           hidden(plotOutput("gomf_gsea_plot", click="gsea_gomf_click")),
                           hidden(downloadButton("download_gsea_gomf_csv", "Save Table Result")),
                           hidden(downloadButton("download_gsea_gomf_png", "Save Plot Result"))
                         ),
                         box(
                           id = "kegg_gsea_box",
                           solidHeader = T,
                           width = 12,
                           #withSpinner(plotOutput("kegg_gsea_plot")),
                           bsPopover("kegg_gsea_plot","Tip!","When click the red or blue point, you can show the enrichment plot of term",placement = "top"),
                           hidden(plotOutput("kegg_gsea_plot", click="gsea_kegg_click")),
                           hidden(downloadButton("download_gsea_kegg_csv", "Save Table Result")),
                           hidden(downloadButton("download_gsea_kegg_png", "Save Plot Result"))
                         )
                         ,useShinyalert()
                       )
              ),#End of GSEA
              tabPanel("Results of GSEA Pathview",
                       fluidRow(
                         use_busy_spinner(spin = "flower",
                                          color = "#112446",
                                          position = c("top-right"),
                                          margins = c("50%", "50%"),
                                          spin_id = NULL,
                                          height = "50px",
                                          width = "50px"),
                         box(
                           id = "gsea_pathview_box",
                           solidHeader = T,
                           width = 12,
                           hidden(imageOutput("gsea_pathview_image",height="auto")),
                           h3(id="gsea_up_title","Up-regulated:"),
                           hidden(DT::dataTableOutput("result_of_gsea_up_regulated")),
                           h3(id="gsea_down_title","Down-regulated:"),
                           hidden(DT::dataTableOutput("result_of_gsea_down_regulated")),
                           hidden(downloadButton("download_gsea_pathview","Save GSEA Pathview Zip", "download_gsea_pathview"))
                         )
                         ,useShinyalert()
                       )
              ),
              tabPanel("Result of PPI Network Analysis",
                       fluidRow(
                         use_busy_spinner(spin = "flower",
                                          color = "#112446",
                                          position = c("top-right"),
                                          margins = c("50%", "50%"),
                                          spin_id = NULL,
                                          height = "50px",
                                          width = "50px"),
                         box( 
                           id = "ppi_box",
                           solidHeader = T,
                           width = 12,
                           hidden(downloadButton("ppi_image_download_btn", "Save PPI Network Image")),
                           hidden(downloadButton("ppi_tsv_download_btn", "Save PPI Network Info")), 
                           uiOutput("ppi_image")
                         )
                         ,useShinyalert()
                       )
              )
            )#End of tabBox
          )#End of fluidRow
        ),#End of tabItem_analysis
        tabItem(
          tabName = "tutorial",
          fluidRow(
            h1("Tutorial"),
            downloadButton("lfq_example_file_download_btn", "Download LFQ Example File","lfq_example_file_download_btn"),
            downloadButton("ibaq_example_file_download_btn", "Download iBaq Example File","ibaq_example_file_download_btn"),
            downloadButton("tmt_example_file_download_btn", "Download TMT Example File","tmt_example_file_download_btn"),
            tags$br(),
            tags$br(),
            h3("1. Input Data preparation"),tags$br(),
            tags$div(id="input_data_box",
                     tags$ul(
                       tags$li(
                         "1) LiBT enables isobaric label-based Tandem Mass Tags (TMT) and label-free iBAQ and LFQ data analysis, and support txt format files. 
                         To perform the analysis properly, the data must have the following columns.The workflow varies slightly depending on the data type. 
                         A brief description for each type is as follows.",tags$br(),tags$br(),
                         tags$ul(
                           tags$li("a. Common",tags$br(),tags$br(),
                                   tags$table(border=1,
                                     tags$tr(
                                       tags$th(style="padding:5px; text-align:center","Columns"),
                                       tags$th(style="padding:5px; text-align:center","Description")
                                     ),
                                     tags$tr(
                                       tags$td(style="padding:5px; text-align:center","GeneName"),
                                       tags$td(style="padding:5px;","Name(s) of the gene(s) associated with the protein(s) contained within the group.")
                                     ),
                                     tags$tr(
                                       tags$td(style="padding:5px; text-align:center","ProteinID"),
                                       tags$td(style="padding:5px;","Identifier(s) of protein(s) contained in the protein group.")
                                     )
                                   )
                          ),tags$br(),tags$br(),
                          tags$li("b. Label-Free",
                                  tags$div("- iBAQ (the theoretical number of tryptic peptides to normalize each protein within the same sample): 
                                             The columns of the expression value of each sample must be the sample name including “iBAQ”."),
                                  tags$div("- LFQ (Normalized intensity resulted from maxLFQ algorithm): 
                                             LFQ input data file does not need the normalization step."),
                                  tags$div("- As long as the required columns are in the input file, 
                                     it is acceptable to use a custom-made file. However, it is recommended to use the original file."),
                                  tags$br(),
                                  tags$table(border=1,
                                             tags$tr(
                                               tags$th(style="padding:5px; text-align:center","Columns"),
                                               tags$th(style="padding:5px; text-align:center","Description")
                                             ),
                                             tags$tr(
                                               tags$td(style="padding:5px; text-align:center","Intensity"),
                                               tags$td(style="padding:5px;","column name must contain iBAQ or LFQ. ex) “iBAQ” + sample name, “LFQ” + sample name")
                                             ),
                                             tags$tr(
                                               tags$td(style="padding:5px; text-align:center","Unique peptides"),
                                               tags$td(style="padding:5px;","Number of peptides associated with each protein in protein group, occuring in the order as the protein IDs occur in the 'Protein IDs' column.")
                                             ),
                                             tags$tr(
                                               tags$td(style="padding:5px; text-align:center","Sequence coverage [%]"),
                                               tags$td(style="padding:5px;","Percentage of the sequence that is covered by the identified peptides of the best protein sequence contained in the group.")
                                             ),
                                             tags$tr(
                                               tags$td(style="padding:5px; text-align:center","Only identified by site"),
                                               tags$td(style="padding:5px;","When marked with '+', this particular protein group was identified only by a modification site.")
                                             ),
                                             tags$tr(
                                               tags$td(style="padding:5px; text-align:center","Reverse"),
                                               tags$td(style="padding:5px;","When marked with '+', this particular protein group contains no protein, made up of at least 50% of the peptides of the leading protein, with a peptide derived from the reversed part of the decoy database.")
                                             ),
                                             tags$tr(
                                               tags$td(style="padding:5px; text-align:center","Potential contaminant"),
                                               tags$td(style="padding:5px;","When marked with '+', this particular protein group was found to be a commonly occurring contaminant. These should be removed for further data analysis.")
                                             )
                                   )
                             
                          ),tags$br(),tags$br(),
                          tags$li("c. Label-Based",tags$br(),tags$br(),
                                  tags$div("- TMT : - Primarily, TMT input format is produced by Proteome Discoverer (PD). 
                                  However, recently Maxquant also supports TMT quantification data. 
                                  Therefore, LiBT accepts TMT quantification result files processed from both tools. 
                                  If the result is from a PD research tool, users are required to check if a normalization process was done. 
                                  TMT data processed by Maxquant needs to be original as it is from Maxquant to be used in LiBT. 
                                  Following columns can be used for identified filtering. Use the result from research tools as it is (Do not use custom-made files).")
                                  )
                        )
                       ),tags$br(),
                       tags$li("Input data can be filtered by selecting “Only identified by site”, “Reverse” and/or “Potential contaminant” options. Also, a log2-transformed input data is used for further analysis."),
                       tags$br(),
                       tags$li("2) Experimental design file is not required, 
                               and the sample list used in the experiment is provided to the user by referring to the “Intensity” column of the input data.
                               If the user selects a case group from the list, 
                               he can select the samples to be used as a control group among the remaining samples.
                               However, since only pairwise analysis is possible, care must be taken in selecting case and control samples.
                               After selecting submit, the following design matrix is automatically created."),
                               tags$br(),
                               tags$table(border=1, style="position:relative; left:50%; margin-left:-126px",
                                 tags$tr(
                                   tags$th(style="padding:5px; text-align:center",""),
                                   tags$th(style="padding:5px; text-align:center","condition"),
                                   tags$th(style="padding:5px; text-align:center","replicate")
                                 ),
                                 tags$tr(
                                   tags$td(style="padding:5px; text-align:center","A_1"),
                                   tags$td(style="padding:5px; text-align:center","case"),
                                   tags$td(style="padding:5px; text-align:center","1")
                                 ),
                                 tags$tr(
                                   tags$td(style="padding:5px; text-align:center","A_2"),
                                   tags$td(style="padding:5px; text-align:center","case"),
                                   tags$td(style="padding:5px; text-align:center","2")
                                 ),
                                 tags$tr(
                                   tags$td(style="padding:5px; text-align:center","A_3"),
                                   tags$td(style="padding:5px; text-align:center","case"),
                                   tags$td(style="padding:5px; text-align:center","3")
                                 ),
                                 tags$tr(
                                   tags$td(style="padding:5px; text-align:center","A_4"),
                                   tags$td(style="padding:5px; text-align:center","case"),
                                   tags$td(style="padding:5px; text-align:center","4")
                                 ),
                                 tags$tr(
                                   tags$td(style="padding:5px; text-align:center","B_1"),
                                   tags$td(style="padding:5px; text-align:center","control"),
                                   tags$td(style="padding:5px; text-align:center","1")
                                 ),
                                 tags$tr(
                                   tags$td(style="padding:5px; text-align:center","B_2"),
                                   tags$td(style="padding:5px; text-align:center","control"),
                                   tags$td(style="padding:5px; text-align:center","2")
                                 ),
                                 tags$tr(
                                   tags$td(style="padding:5px; text-align:center","B_3"),
                                   tags$td(style="padding:5px; text-align:center","control"),
                                   tags$td(style="padding:5px; text-align:center","3")
                                 )
                               ),
                              tags$br(),
                              tags$div("To perform downstream analysis, replicates are randomly numbered from 1 to the number of samples in each group, 
                                       but analysis using replicate is not supported. ")
                     )),
            tags$br(),
            h3("2. Data preprocessing"),
            tags$br(),
            tags$ul(
              tags$li("- Preprocessing can be done in the second tab of the right sidebar. Users can designate the preprocessing steps by drag and drop. 
                      The order of steps is important to take note of because the process proceeds in that set order. 
                      Users can exclude it by moving it to “Not to Use” at the bottom. Otherwise, at least one preprocessing must be performed."),
              tags$li("The following is a list of preprocessing options provided by LiBT.",
                      tags$br(),tags$br(),
                      tags$ul(
                        tags$li("1) Valid value This option is used to determine whether to use the protein according to the ratio of missing values. 
                        The number of valid values is calculated based on the case group, and when less than the number of samples is expressed 
                        in at least one group, the corresponding protein is removed. The default value is 70%."),
                        tags$br(),
                        tags$li("2) Imputation",
                                tags$br(),
                                tags$ul(
                                  tags$li("a. QRILC : A missing data imputation method that performs the imputation of left-censored missing data
                                          using random draws from a truncated distribution with parameters estimated using quantile regression."),
                                  tags$br(),
                                  tags$li("b. MinProb : Performs the imputation of left-censored missing data by random draws from a Gaussian distribution
                                  centered to a minimal value. Considering an expression data matrix with n samples and p features, 
                                  for each sample, the mean value of the Gaussian distribution is set to a minimal observed value in that sample. 
                                  The minimal value observed is estimated as being the q-th quantile (default q = 0.01) of the observed values 
                                  in that sample. The standard deviation is estimated as the median of the feature standard deviations. 
                                  Note that when estimating the standard deviation of the Gaussian distribution, 
                                  only the peptides/proteins which present more than 50% recorded values are considered."),
                                  tags$br(),
                                  tags$li("c. Man (Perseus-type, default) : This method is based on the popular missing value imputation procedure implemented in the Perseus software. 
                                          The missing values are replaced by random numbers drawn from a normal distribution of 1.8 standard deviation downshift 
                                          and with a width of 0.3 of each sample."),
                                  tags$br(),
                                  tags$li("d. Constant : Replaces the missing values by 0")
                                )),
                        tags$br(),
                        tags$li("3) Normalization",
                                tags$div("- This is a step that should be used differently depending on the data type. Users can select “No” or “Not to Use” to exclude it."),
                                tags$div("- Normalization is the variant stabilizing normalization (vsn) method on the protein intensity distribution in each sample, which is the default method of the DEP package.")
                                ),
                        tags$br(),
                        tags$li(
                                tags$div("4) Once preprocessing is completed, the user can check the quality control (QC) results in 4 types of plots using the DEP package. 
                                The bar plot shows the number of identified proteins for each sample, and can always be checked regardless of the user-selected step. 
                                However, only after the normalization step or when there are missing values, 
                                users can check the density plot that can show the result of distribution before or after imputation. 
                                Users can also check a box plot per sample that shows distributions before and after normalization as well as density plot to see missing value and bias as a result of density and cumulative fraction."),
                                tags$br(),
                                imageOutput("preprocessing4th", width="100%", height="520px")
                                )
                        
                      ))
            ),
            tags$br(),
            h3("3. Differential Expression Analysis"),
            tags$ul(
              tags$li(
                "1) Statistical analysis",
                tags$br(),
                tags$ul(
                  tags$li("- The following 4 options are provided. Statistical test results are calculated without considering replicates."),
                  tags$li("a. T-test (ver 4.0.5, default)"),
                  tags$li("b. Wilcoxon Rank Sum (ver 4.0.5)"),
                  tags$li("c. edgeR (ver 3.32.1)"),
                  tags$li("d. limma(ver 3.46.0)"),
                  tags$br(),
                  tags$li("- p-value correction option (Multiple testing)"),
                  tags$li("a. Benjamini-Hochberg (BH)"),
                  tags$li("b. Bonferroni"),
                  tags$br(),
                  tags$li("- Significant protein filtering criteria default to p.value 0.05 and FoldChange 1.5.")
                )
              ),
              tags$li(
                "2) Result plots",
                tags$br(),
                tags$ul(
                  tags$li(
                    tags$div("a. Principal Component Analysis (PCA) plot"),
                    tags$div("By default, the sample name is not shown. If the user needs a plot with the sample name, this can be adjusted with the toggle button."),
                    tags$br(),
                    imageOutput("pcaplotImg", width="100%", height="520px")
                  ),
                  tags$li(
                    tags$div("b. Volcano plot"),
                    tags$div("Each point indicated in the plot is a single protein, and if the protein (DEP) satisfies the significant protein filtering criteria, it is highlighted in red. If the user drags a desired area on the plot, the results of the statistical analysis of proteins in the area are displayed in the table below."),
                    tags$br(),
                    imageOutput("vcplotImg", width="100%", height="700px")
                  ),
                  tags$li(
                    tags$div("c. Correlation heatmap"),
                    tags$div("Visualize the result of correlation analysis (Pearson coefficient correlation) between samples as a heatmap."),
                    tags$br(),
                    imageOutput("corrplotImg", width="100%", height="520px")
                  ),
                  tags$li(
                    tags$div("d. Heatmap with clusters"),
                    tags$div("Heatmap with clusters visualizes the Z-score normalized log 2-transformed intensity value of proteins along with DEP Clustering results. The default value of k for gene clustering is set to 1. If the default value is 1, the optimal number of clusters k is automatically calculated, 
                             and clustering is possible even based on a specific value of k. The optimal k is calculated by the fviz_nbclust method of the factoextra library. Clustering results can be downloaded."),
                    tags$br(),
                    imageOutput("heatmapImg", width="100%", height="550px")
                  )
                )
              )
            ),
            tags$br(),
            h3("4. Further analysis "),
            tags$br(),
            tags$div("After DEP is completed, further analysis can be performed from the last tab in the right sidebar. Analysis provided by LiBT includes GSA, GSEA, KEGG pathway analysis, and PPI analysis. This section describes the options used for each analysis."),
            tags$br(),
            tags$ul(
              tags$li(
                tags$div("1) Gene Set Analysis (GSA)"),
                tags$br(),
                tags$ul(
                   tags$li("a. Basic analysis",
                           tags$div("GSA uses the “enrichR” package and is based on a total of 4 databases:"),
                           tags$ul(
                             tags$li("- GO_Molecular_Function_2021"),
                             tags$li("- GO_Cellular_Component_2021"),
                             tags$li("- GO_Biological_Process_2021"),
                             tags$li("- KEGG_2021_Human")
                           ),
                           tags$div("The GSA section has two input data options and one output data option. In order to proceed with the GSA, it is first necessary to select which group to conduct the analysis on. Select Case-up protein group or Ctrl-up (Case-down) protein group and input the corresponding gene symbols. 
                                    Next, users need to decide which proteins to base the GSA on. Since GSA analysis is usually performed based on DEP, the default value is selected as the DEP level. These GSA results are displayed as barplot, and the user can directly set how many barplot results to output. 
                                    All plots can be directly downloaded with the download button. The file that is downloaded as a csv has the entire GSA result that is not filtered."),
                           tags$br(),
                           imageOutput("barplotGSA", width="100%", height="280px")
                          ),
                   tags$li("b. KEGG pathway mapping",
                          tags$div("Based on the KEGG pathway analysis result of GSA, a KEGG pathway map can be generated. It provides integrated results by downloading the pathway graph using the pathview library and mapping the log2Foldchange value of the data input by the user to the graph. 
                                   The visualized result can be obtained by clicking the desired term in the table or clicking the desired term in the drop-down menu and pressing the render button. 
                                   The result obtained in this way is printed at the bottom of the table, and if the zoom button is selected, the result of a large image in the modal window can be checked. 
                                   If clicked at least once to check the graph, users can download it as a zip file by clicking the download button."),
                          tags$br(),
                          imageOutput("pathviewImgGSA", width="100%", height="600px")
                   )
                )
              ),
              tags$li(
                tags$div("2) Gene Set Enrichment Analysis (GSEA)"),
                tags$br(),
                tags$ul(
                  tags$li("a. Basic analysis",
                          tags$div("GSEA is done by the “fgsea” package using msigDB including:"),
                          tags$ul(
                            tags$li("- c2.cp.kegg.v7.4.symbols.gmt"),
                            tags$li("- c5.bp.v7.4.symbols.gmt"),
                            tags$li("- c5.cc.v7.4.symbols.gmt"),
                            tags$li("- c5.mf.v7.4.symbols.gmt")
                          ),
                          tags$div("Before performing GSEA, four gene ranking value options should be selected: log2FoldChange, P.value, P.adj, log2(FC)*-log10 (P.adj). 
                                   GSEA results are divided into UP-regulated and Down-regulated depending on whether the Normalized ES (NES) value is greater than or less than 0. 
                                   A plot is drawn after being cut off based on p.adj 0.05. Also, if you click the red and blue points of the segment plot, 
                                   you can check the enrichment plot of each result. All results and plots (except for enrichment plot) can be directly downloaded with the download button."),
                          tags$br(),
                          imageOutput("erplotGSEA", width="100%", height="450px")
                  ),
                  tags$li("b. KEGG pathway mapping",
                          tags$div("GSEA also performs KEGG pathway mapping based on the KEGG pathway analysis result. Using the pathview package, upregulated and downregulated pathways are separately provided by a table. 
                                   If each pathway name is selected, a modal window appears and the user can check the KEGG pathway graph to which the log2FoldChange value is mapped. 
                                   If the graph is reviewed by selecting it once, the user can download it as a zip file by clicking the download button."),
                          tags$br(),
                          imageOutput("pathviewImgGSEA", width="100%", height="440px")
                  )
                )
              ),
              tags$li(
                tags$div("3) Protein-Protein Interaction (PPI) network analysis"),
                tags$br(),
                tags$div("StringDB application programming interface (API) was used for PPI network analysis. 
                         The number of default input proteins is the highest value as the number of DEPs. A value higher than the number cannot be entered. 
                         If a number lower than the number of DEPs is entered, the sort is based on one selected among P.adj (default), 
                         P.value, and log2FoldChange, and the number of inputs is selected and inputted. If the user selects Protein Name, 
                         an area for text input appears, where the user enters the text separated by “,” as described in the example. 
                         If the Protein Name option is selected, PPI network analysis is possible even if the analysis of the previous workflow has not been performed. 
                         The resulting network can be downloaded as an image file, and network information can also be downloaded."),
                tags$br(),
                imageOutput("ppiImg", width="100%", height="520px")
              )
            ),
            tags$br(),
            h3("5. Download options"),
            tags$br(),
            tags$div("Once LiBT finishes preprocessing the input data, the plot results and all the table results can be directly downloaded through the download button located on each result panel."),
            tags$br(),
            tags$ul(
              tags$li(
                tags$div("1) Download the results table"),
                tags$br(),
                tags$ul(
                  tags$li("a) log2 transformed data matrix: In the input stage, the original data must undergo log2 transformation for the next analysis stage. This transformed data can be downloaded in csv format."),
                  tags$br(),
                  tags$li("b) Preprocessed data matrix: A user-defined preprocessing step is performed based on the data on which log2 transformation has been completed, and the preprocessed data can be downloaded in csv format."),
                  tags$br(),
                  tags$li("c) DEP info table: When DEP analysis is completed, p-value, fold-change, and corrected q-value can be obtained for each gene. This information can be downloaded in csv format."),
                  tags$br(),
                  tags$li("d) clustering info table: The k-means clustering result of DEP can be downloaded. k is determined according to the option setting, and protein name, accessionID, and cluster class information can be obtained as a result."),
                  tags$br(),
                  tags$li("e) GSA/GSEA result: The entire GSA/GSEA result can be downloaded. Only terms that meet the cutoff criteria are shown in the results shown in LiBT, but when downloaded to a table, the entire GSA/GSEA analysis result is downloaded. In the downloaded result, calculated values such as p-value, q-value, log2 fold-change, etc. can be checked, and users can manually change the cutoff criteria to check other results in addition to the cutoff results specified by LiBT."),
                  tags$br(),
                  tags$li("f) PPI network result: Users can download connection information between proteins. It includes which proteins are linked and what the linkage score is. By adding log2FoldChange of the corresponding protein, style customization is possible in Cytoscape, etc.")
                )
              ),
              tags$br(),
              tags$li(
                tags$div("2) Download the results plot"),
                tags$br(),
                tags$div("All plots shown as a result in LiBT can be downloaded using the respective download button. PPI analysis results are downloaded in JPEG format, and all result plots except PPI are downloaded in PNG format.")
              )
            )
          )
        )#End fo tabItems_tutorial
      )#tepItems
    ),# End of dashboardBody
    controlbar = dashboardControlbar(
      skin = "dark",
      collapsed = F,
      width = 350,
      controlbarMenu(
        id = "entireRightsidebar",
        controlbarItem(
          id="File_upload",
          active = T,
          title="",
          icon=icon("file-upload"),
          fluidRow(
            box(
              title = "Upload Data",
              width = 12,
              status = "danger",
              solidHeader=TRUE, 
              closable = F,
              radioButtons("file_type", label="",
                           choices = list("LFQ" = "LFQ", "iBAQ" = "iBAQ", "TMT" = "TMT"), 
                           selected = "LFQ"),
              tags$hr(),
              fileInput("fileBrowser", "Choose txt File",
                        multiple = FALSE,
                        accept = c("text/csv",
                                   "text/comma-separated-values,text/plain",
                                   ".csv")),
              hidden(radioButtons("TMT_tool_option", label="Choose TMT data research tool", 
                                  choices = list("Proteome Discover" = "PD", "MaxQuant" = "MQ"), 
                                  selected = "PD")),
              hidden(radioButtons("TMT_input_option", label="Get Normalized TMT data", 
                                  choices = list("YES" = "T", "NO" = "F"), 
                                  selected = "F")),
              checkboxGroupInput("nonTMT_input_option", label="Filtering Options (Remove all)", 
                                 choices = list("Potential contaminant" = "potential", 
                                                "Reverse" = "reverse", 
                                                "Only identified by site" = "identified")),
              actionButton("select_all_filtering_btn", "Select All"),
              actionButton("deselect_all_filtering_btn", "Deselect All"),
              tags$hr(),
              actionButton("file_upload_btn", "Upload file")
            ),
            box(
              title = "Experiment Design",
              width = 12,
              status = "danger",
              solidHeader=TRUE, 
              closable = F,
              footer = fluidRow(
                pickerInput(inputId = "case_group_selection", label="Case samples", multiple=TRUE,
                            choices=c(), options=list(`actions-box` = TRUE, `selected-text-format` = "count > 1")),
                pickerInput(inputId = "control_group_selection", label="Control samples", multiple=TRUE,
                            choices=c(), options=list(`actions-box` = TRUE, `selected-text-format` = "count > 1")),
                tags$hr(),
                actionButton("exp_design_submit_btn", "Submit"),
                actionButton("exp_design_reset_btn", "Reset")
              )
            ) 
          )
        ),# End of Upload Data box
        controlbarItem(
          id="preprocesseing",
          title="",
          icon=icon("dna"),
          active = F,
          box(
            title = "Preprocessing",
            width = 12,
            status = "primary",
            solidHeader=TRUE, 
            closable = F,
            hidden(prettyToggle(inputId="exp_design_check", 
                                label_on="submitted", label_off="not_submitted", value=FALSE)),
            hidden(prettyToggle(inputId="preprocess_check", 
                                label_on="preprocessed", label_off="not_preprocessed", value=FALSE)),
            bucket_list(
              header = "Drag to select options and set the order",
              add_rank_list(
                text = "Use",
                labels = list(
                  "Use_valid_value" = radioButtons("valid_value", label="Choose % of valid value",
                                                   choices = list("30%" = 0.3, "50%" = 0.5, "70%"= 0.7, "100%"=0), selected = 0.7),
                  "Use_imputation" = radioButtons("imputation", label="Choose Imputation",
                                                  choiceNames = list(HTML("QRILC<br/>(quantile regression-based<br/>left-censored function)"),
                                                                     HTML("MinProb<br/>(left-shifted Gaussian distribution)"),
                                                                     HTML("Man<br/>(manually defined<br/>left-shifted Gaussian distribution)"),
                                                                     HTML("Constant<br/>(replace with zero value)")),
                                                  choiceValues = list("QRILC", "MinProb", "man", "zero"), selected = "man"),
                  # choices = list("Normal distribution" = "normal_distribution",
                  #                "Constant" = "constant",
                  #                "NaN"="nan")),
                  "Use_normalization" = radioButtons("normalization", label="Choose method of normalization",
                                                     choices = list("Yes" = "Yes",
                                                                    "No" = "No"))
                  # choices = list("Width distribution" = "quantile",
                  #                "Z-score"="zscore", "none" = "none"))
                ),
                input_id = "use_options"
              ),
              add_rank_list(
                text = "Not to Use",
                labels = NULL,
                input_id = "not_to_use_options"
              )
            ),
            tags$hr(),
            actionButton("preprocess_btn", "Start preprocessing")
          ) 
        ),# End of preprocessing box
        controlbarItem(
          id="dea",
          title="",
          icon=icon("chart-bar"),
          active = F,
          box(
            title = HTML("Differential <br/>experimental analysis <br/><br/>_Test"),
            width = 12,
            status = "success",
            solidHeader=TRUE, 
            closable = F,
            uiOutput("dea_case"),
            uiOutput("dea_control"),
            radioButtons("test_method", label="Choose test method", 
                         choices = list("T-test" = "T-test", "Wilcoxon Rank Sum" = "Wilcoxon-Rank-Sum", "edgeR" = "edgeR", "Limma" = "Limma")),
            radioButtons("padj_method", label="Choose p.adj method", 
                         choices = list("Benjamini-Hocherg" = "BH", "Bonferroni" = "bonferroni")),
            actionButton("test_btn", "Start DEA")
          ),
          box(
            title = HTML("Differential <br/>experimental analysis <br/><br/>_Visulization"),
            width = 12,
            status = "success",
            solidHeader=TRUE, 
            closable = F,
            radioButtons("thres_type", label="Choose threshold type", 
                         choices = list("P.adj" = "P.adj", "P.value" = "P.value", "None" = "none"),selected = "P.value"),
            numericInput("dea_pvalue", label = HTML("Set the threshold <br/>for the adjusted P-value or P-value"),
                         value=0.05),
            numericInput("dea_log2fc", label = HTML("Set the threshold <br/>for the log<sub>2</sub> Fold Change"),
                         value=1.5),
            numericInput("dea_clusterNum", label = HTML("Set a number of cluster <br/>for heatmap<br/>(between 1 and number of samples)"), value=1),
            tags$hr(),
            actionButton("dea_btn", "Start Draw")
          ) 
        ),# End of DEP box
        controlbarItem(
          id="etc",
          title="",
          icon=icon("chart-pie"),
          active = F,
          box(
            title = HTML("Gene Set Analysis <br/><br/>"),
            width = 12,
            status = "warning",
            solidHeader=TRUE, 
            closable = F,
            uiOutput("dep_up"),
            uiOutput("dep_down"),
            tags$div(
              id="gsa_set", class="form-group shiny-input-radiogroup shiny-input-container",
              tags$label(class="control-label", `for`="gsa_set", "Choose Data set"),
              tags$div(class="shiny-options-group",
                       tags$div(class="radio",
                                tags$label(
                                  tags$input(type="radio", name="gsa_set", value="case_up", checked="checked",
                                             tags$span(HTML("Case-UP")))
                                )
                       ),
                       tags$div(class="radio",
                                tags$label(
                                  tags$input(type="radio", name="gsa_set", value="case_down",
                                             tags$span(HTML("Ctrl-UP")))
                                )
                       )
              )
            ),
            tags$div(
              id="gsa_input_set", class="form-group shiny-input-radiogroup shiny-input-container",
              tags$label(class="control-label", `for`="gsa_input_set", "Choose stats of gene level"),
              tags$div(class="shiny-options-group",
                       tags$div(class="radio",
                                tags$label(
                                  tags$input(type="radio", name="gsa_input_set", value="total", 
                                             tags$span(HTML("Total")))
                                )
                       ),
                       tags$div(class="radio",
                                tags$label(
                                  tags$input(type="radio", name="gsa_input_set", value="dep",checked="checked",
                                             tags$span(HTML("DEP")))
                                )
                       )
              )
            ),
            numericInput("set_nterm", label = HTML("Set a number of showed Term <br/>for barplot"),
                         value=10),
            actionButton("gsa_btn", "Start GSA")
          ),
          box(
            title = HTML("Gene Set <br/>Enrichment Analysis <br/><br/>"),
            width = 12,
            status = "warning",
            solidHeader=TRUE, 
            closable = F,
            tags$div(
              id="select_genelevel_stats", class="form-group shiny-input-radiogroup shiny-input-container",
              tags$label(class="control-label", `for`="select_genelevel_stats", "Choose stats of gene rank"),
              tags$div(class="shiny-options-group",
                       tags$div(class="radio",
                                tags$label(
                                  tags$input(type="radio", name="select_genelevel_stats", value="log2fc", checked="checked",
                                             tags$span(HTML("log<sub>2</sub>FoldChange")))
                                )
                       ),
                       tags$div(class="radio",
                                tags$label(
                                  tags$input(type="radio", name="select_genelevel_stats", value="P.value",
                                             tags$span(HTML("P.value")))
                                )
                       ),
                       tags$div(class="radio",
                                tags$label(
                                  tags$input(type="radio", name="select_genelevel_stats", value="P.adj",
                                             tags$span(HTML("P.adj")))
                                )
                       ),
                       tags$div(class="radio",
                                tags$label(
                                  tags$input(type="radio", name="select_genelevel_stats", value="log2fcxmlog10padj",
                                             tags$span(HTML("log<sub>2</sub>(FC)*-log<sub>10</sub>(P.adj)")))
                                )
                       )
              )
            ),
            actionButton("gsea_btn", "Start GSEA")
          ), # End of GSEA box
          box(
            title = HTML("Protein-Protein Interaction <br/>Network Analysis <br/><br/>_StringDB"),
            width = 12,
            status = "warning",
            solidHeader=TRUE, 
            closable = F,
            pickerInput("input_organism_toppi", label=HTML("Choose organism"), choices = c("Homo Sapiens")),
            numericInput("input_num_toppi", label = HTML("Set a number of input proteins <br/> (Defualt[=max] : # of DEP)"),value=100),
            tags$div(
              id="select_ppi_condition", class="form-group shiny-input-radiogroup shiny-input-container",
              tags$label(class="control-label", `for`="select_ppi_condition", "Choose condition of gene select"),
              tags$div(class="shiny-options-group",
                       tags$div(class="radio",
                                tags$label(
                                  tags$input(type="radio", name="select_ppi_condition", value="padj", checked="checked",
                                             tags$span(HTML("P.adj")))
                                )
                       ),
                       tags$div(class="radio",
                                tags$label(
                                  tags$input(type="radio", name="select_ppi_condition", value="pval",
                                             tags$span(HTML("P.value")))
                                )
                       ),
                       tags$div(class="radio",
                                tags$label(
                                  tags$input(type="radio", name="select_ppi_condition", value="log2fc",
                                             tags$span(HTML("log<sub>2</sub>FoldChange")))
                                )
                       ),
                       tags$div(class="radio",
                                tags$label(
                                  tags$input(type="radio", name="select_ppi_condition", value="Gene_Name",
                                             tags$span(HTML("Protein Name")))
                                )
                       )
              )
            ),
            hidden(textInput("input_gene_toppi", "Input gene name", placeholder = "Ex)TP53 or TP53, PLK1, RHOB")),
            actionButton("ppi_btn", "Start PPIA")
          ) # End of STRING box
        )
      )#End of controlbarMenu
    )
  ))}
