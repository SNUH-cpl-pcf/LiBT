library(dplyr)
library(readr)
library(stringr)
library(DEP)
library(SummarizedExperiment)
library(ggplot2)
library(ggrepel)
library(gridExtra)
library(edgeR)
library(factoextra)
library(enrichR)
library(tibble)
library(fgsea)
library(pathview)
library(zip)
library(data.table)
library(readxl)
library(BiocParallel)
library(promises)
library(future)
plan(multisession)


source("function.R")
# 2019.12.30

HEARTBEAT_INTERVAL_MILLIS = 1000  # 1 second

shinyServer(function(input,output, session){
  
  options(shiny.maxRequestSize=100*1024^2)
  
  #### Heartbeat ####
  # Define reactive variable
  cnt <- reactiveVal(0)
  # Define time dependent trigger
  autoInvalidate <- reactiveTimer(HEARTBEAT_INTERVAL_MILLIS)
  # Time dependent change of variable
  observeEvent(autoInvalidate(), {  cnt(cnt() + 1)  })
  #### Spark job ####
  result <- reactiveVal() # the result of the spark job
  busy <- reactiveVal(0)  # whether the spark job is running
  
  step <- c()
  info <- c()
  sample_num <- c()
  time <- c()
  user_id <- make_id()
  base_dir <- getwd()
  
  ############################
  ##      RENDERING         ##
  ############################
  
  observeEvent(input$file_type, {
    if(input$file_type=="TMT"){
      shinyjs::show("TMT_tool_option")
      shinyjs::show("TMT_input_option")
      shinyjs::hide("nonTMT_input_option")
      shinyjs::hide("select_all_filtering_btn")
      shinyjs::hide("deselect_all_filtering_btn")
    }else{
      shinyjs::hide("TMT_tool_option")
      shinyjs::hide("TMT_input_option")
      shinyjs::show("nonTMT_input_option")
      shinyjs::show("select_all_filtering_btn")
      shinyjs::show("deselect_all_filtering_btn")
    }
  })
  
  observeEvent(input$TMT_tool_option, {
    if(input$file_type=="TMT" & input$TMT_tool_option == "PD"){
      shinyjs::show("TMT_input_option")
      shinyjs::hide("nonTMT_input_option")
      shinyjs::hide("select_all_filtering_btn")
      shinyjs::hide("deselect_all_filtering_btn")
    }else{
      shinyjs::show("nonTMT_input_option")
      shinyjs::show("select_all_filtering_btn")
      shinyjs::show("deselect_all_filtering_btn")
      shinyjs::hide("TMT_input_option")
    }
  })
  
  observeEvent(input$fileBrowser, {
      temp <- file_input()
      if(is.null(temp)){
        shinyalert("Check your file type!", type="error", timer = 10000,
                   closeOnClickOutside = T, closeOnEsc = T)
        reset("fileBrowser")
      } else {
        if(input$file_type=="TMT"){
          if(input$TMT_tool_option == "PD"){
            shinyjs::show("TMT_input_option")
            state <- colnames(temp)[grep("normalized", colnames(temp), ignore.case = T)]
            updateRadioButtons(session, "TMT_input_option", label="Get Normalized TMT data",
                               choices = list("YES" = "T", "NO" = "F"), selected="F")
            if(length(state)==0){
              shinyjs::disable("TMT_input_option")
            }
          }
        }
        
        info = paste0("* File Type : ", input$file_type,"\n")
        timeLine <<- data.frame(step="Start!:)",info=info,
                                sample_num=as.numeric(nrow(temp)),
                                time=as.character(Sys.time()),color="red",icon="file-upload")
        addTimeLine(timeLine)
      }
  })
  
  observeEvent(input$select_all_filtering_btn, {
    if(input$select_all_filtering_btn == 0) {
      return(NULL)
    } else if(input$select_all_filtering_btn >0){
      updateCheckboxGroupInput(session, "nonTMT_input_option",
                               choices = list("Potential contaminant" = "potential",
                                              "Reverse" = "reverse",
                                              "Only identified by site" = "identified"),
                               selected = c("potential","reverse","identified"))
    }
  })
  
  observeEvent(input$deselect_all_filtering_btn, {
    if(input$deselect_all_filtering_btn == 0) {
      return(NULL)
    } else if(input$deselect_all_filtering_btn >0){
      updateCheckboxGroupInput(session, "nonTMT_input_option",
                               choices = list("Potential contaminant" = "potential",
                                              "Reverse" = "reverse",
                                              "Only identified by site" = "identified"),
                               selected = c())
    }
  })
  
  observeEvent(input$file_upload_btn, {
    setwd(base_dir)
    updatePrettyToggle(session, "preprocess_check", label=NULL, value=F)
    updatePrettyToggle(session, "exp_design_check", label=NULL, value=F)
    
    output$pca_plot <- NULL
    output$correlation_matrix <- NULL
    output$volcano_plot <- NULL
    output$heatmap <- NULL
    output$gobp_gsa_plot <- NULL
    output$gocc_gsa_plot <- NULL
    output$gomf_gsa_plot <- NULL
    output$kegg_gsa_plot <- NULL
    output$topOfKeggDT <- NULL
    output$pathID_selector <- NULL
    output$pathview_result <- NULL
    output$gobp_gsea_plot <- NULL
    output$gocc_gsea_plot <- NULL
    output$gomf_gsea_plot <- NULL
    output$kegg_gsea_plot <- NULL
    output$gsea_pathview_image <- NULL
    output$result_of_gsea_up_regulated <- NULL
    output$result_of_gsea_down_regulated <- NULL
    output$ppi_image <- NULL
    
    # || is.null(input$nonTMT_input_option)
    if((is.null(input$fileBrowser))){
      shinyalert("Choose option!", type="error", timer = 10000,
                 closeOnClickOutside = T, closeOnEsc = T)
    } else {
      shinyalert("Data will be converted to log2 scale in the following order!", type = "info")
      show_spinner()
      render_df <- main_data()
      hide_spinner()
      
      output$uploaded_file_header <- DT::renderDataTable({render_df},
        options = list(scrollX = TRUE, pageLength = 5,lengthMenu = c(5, 10, 15)))
      shinyjs::show("download_log2_exp_btn")
      
      samples <- total_samples()
      updatePickerInput(session, "case_group_selection", choices = samples)
      if(length(render_df) != 0){
        if(input$file_type != "TMT"){
          option <- c()
          for(i in 1:length(input$nonTMT_input_option)){
            if(i!=length(input$nonTMT_input_option)){
              option <- paste0(option,input$nonTMT_input_option[i]," / ")  
            } else{
              option <- paste0(option,input$nonTMT_input_option[i])  
            }
          }
          info <- paste0("* Numerical Filter\n  : Peptides = 0\n    Intensity = 0\n",
                         "* Categorical Filter\n", "  : ", option)
          newTL <- data.frame(step="Data Input",info=info,
                              sample_num=as.numeric(nrow(render_df)),
                              time=as.character(Sys.time()),color="red",icon="file-upload")
          
          if(nrow(timeLine)!=2){
            timeLine <<- rbind(timeLine,newTL)
          } else{
            timeLine[2,] <<- newTL
          }
          
          addTimeLine(timeLine)
        } else{
          info <- paste0("*Research tool : " , input$TMT_tool_option,"\n")
          if(input$TMT_tool_option=="PD"){
            info <- paste0(info,"* Is normalized ? : ", input$TMT_input_option)
          }else{
            option <- c()
            for(i in 1:length(input$nonTMT_input_option)){
              if(i!=length(input$nonTMT_input_option)){
                option <- paste0(option,input$nonTMT_input_option[i]," / ")  
              } else{
                option <- paste0(option,input$nonTMT_input_option[i])  
              }
            }
            info <- paste0(info,"* Numerical Filter\n  : Peptides = 0\n    Intensity = 0\n",
                           "* Categorical Filter\n", "  : ", option)
          }
          
          newTL <- data.frame(step="Data Input",info=info,
                              sample_num=as.numeric(nrow(render_df)),
                              time=as.character(Sys.time()),color="red",icon="file-upload")
          
          if(nrow(timeLine)!=2){
            timeLine <<- rbind(timeLine,newTL)
          } else{
            timeLine[2,] <<- newTL
          }
          
          addTimeLine(timeLine)
        }
      }
    }
  })
  
  
  observeEvent(input$case_group_selection, {
    choices <- control_toUpdate()
    updatePickerInput(session, "control_group_selection",
                      choices = choices, selected = choices)
  })
  
  observeEvent(input$exp_design_submit_btn, {
    updatePickerInput(session, "control_group_selection",
                      label=paste0("Control samples (n=", length(control_samples()), ")"),
                      choices = control_samples(), selected=control_samples())
    updatePickerInput(session, "case_group_selection",
                      label=paste0("Case samples (n=", length(case_samples()), ")"),
                      choices = case_samples(), selected=case_samples())
    
    updatePrettyToggle(session, "exp_design_check", label=NULL, value=TRUE)

    if(length(case_samples()) == 0 | length(control_samples()) == 0){
      shinyalert("There is no case & control sample!", "Please choice case & control samples", type="error", timer = 10000,
                 closeOnClickOutside = T, closeOnEsc = T)
    }else{
      if(input$file_type == "TMT"){
        condition <- make_condition(case_samples(),control_samples(), input$TMT_tool_option)
      } else{
        condition <- make_condition(case_samples(),control_samples(), "MQ")
      }
      
      shinyalert("Complete Submit", type="success", timer = 10000,
                 closeOnClickOutside = T, closeOnEsc = T)
      
      group_name <- unique(condition)
      
      tmp <- paste0("* Case Group : ", group_name[1], "\n",
                    "* Control Group : ", group_name[2], "\n",
                    "* # of Case : ", length(case_samples()),"\n",
                    "* # of Control : ", length(control_samples()))
      info <- paste0(info,tmp,"\n")
      newTL <- data.frame(step="Exp Design Submitted",
                          info=info,
                          sample_num=as.numeric(nrow(main_data())),
                          time=as.character(Sys.time()),color="red",icon="file-upload")
      timeLine <<- rbind(timeLine,newTL)
      addTimeLine(timeLine)
    }
  })
  
  observeEvent(input$exp_design_reset_btn, {
    updatePickerInput(session, "case_group_selection", label="Case samples",
                      choices = total_samples())
    updatePickerInput(session, "control_group_selection", label="Control samples",
                      choices = total_samples())
    updatePrettyToggle(session, "exp_design_check", label=NULL, value=FALSE)
  })
  
  observeEvent(input$preprocess_btn, {
    req(ready_for_dea())
    show_spinner()
    if(input$exp_design_check) {
      updatePrettyToggle(session, "preprocess_check", label=NULL, value=TRUE)
      data_se <- ready_for_dea()
      res <- assay(data_se)
      
      output$uploaded_file_header <- DT::renderDataTable({
        DT::datatable(res,options = list(scrollX = TRUE, pageLength = 5,lengthMenu = c(5, 10, 15)), selection="none") %>% 
          DT::formatRound(colnames(res), digits=2)
      })
      

      tmp <- c()
      preprocessing_options <- input$use_options
      for(i in 1:length(preprocessing_options)) {
        if(i != 1){
          tmp <- paste0(tmp,"\n")
        }
        switch(preprocessing_options[i],
               "Use_valid_value" = {
                 vv <- as.character(as.numeric(input$valid_value)*100)
                 tmp <- paste0(tmp,i,". Valid value : ", paste0(vv,"%"))},
               "Use_imputation" = {
                  if(!is.null(preprocessed_data()[[4]])){
                  tmp <- paste0(tmp,i,". Imputation : " , str_to_title(input$imputation))
                  }else{
                    tmp <- paste0(tmp,i,". Imputation : there is no\n    missing value")
                  }
                },
               "Use_normalization" = {
                 tmp <- paste0(tmp,i,". Normalization : ", str_to_title(input$normalization))},
               NULL = {
                 tmp <- tmp
               }
        )
      }
      
      info <- paste0(info,tmp,"\n")
      newTL <- data.frame(step="Preprocessing",
                          info=info,
                          sample_num=as.numeric(nrow(data_se)),
                          time=as.character(Sys.time()),color="blue",icon="dna")
      timeLine <<- rbind(timeLine, newTL)
      addTimeLine(timeLine)
      
      output$dea_case <- renderUI({
        box(
          title = "Case samples",
          closable = F,
          collapsed = T,
          enable_label = T,
          label_text = length(case_samples()),
          label_status = "danger",
          width = 12,
          solidHeader = F,
          collapsible = T,
          footer = HTML(paste(case_samples(), collapse=",<br/>"))
        )
      })
      
      output$dea_control <- renderUI({
        box(
          title = "Control samples",
          closable = F,
          collapsed = T,
          enable_label = T,
          label_text = length(control_samples()),
          label_status = "danger",
          width = 12,
          solidHeader = F,
          collapsible = T,
          footer = HTML(paste(control_samples(), collapse=",<br/>"))
        )
      })
      
      shinyjs::show("download_preprocessed_exp_btn") 
      
      if(length(preprocessing_options) == 3){
        showModal(modalDialog(
          renderUI({plotOutput("plotQC", width = "800px",height="700px")}),
          footer=actionButton("close_modal", label="Close")
        ))
      } else{
        showModal(modalDialog(
          renderUI({plotOutput("plotQC", width = "800px",height="350px")}),
          footer=actionButton("close_modal", label="Close")
        ))
      }
     
      
      Sys.sleep(1)
      hide_spinner()
      
    } else {
      shinyalert("Please submit experiment design", type="error", timer = 10000,
                 closeOnClickOutside = T, closeOnEsc = T)
    }
  })
  
  
  output$plotQC <- renderPlot({
    req(preprocessed_data())
    data_filt <- preprocessed_data()[[2]]
    data_norm <- preprocessed_data()[[3]]
    data_imp <- preprocessed_data()[[4]]

    if(is.null(data_imp)){
      grid.arrange(plot_numbers(data_filt), ncol=1)
    }else if(is.null(data_norm)){
      grid.arrange(plot_numbers(data_filt), plot_detect(data_filt), ncol=2)
    }else{
      grid.arrange(plot_numbers(data_filt), plot_detect(data_filt), plot_normalization(data_filt, data_norm), plot_imputation(data_norm,data_imp), ncol=2, nrow=2)
    }

  })
  
  observeEvent(input$test_btn, {
    if(input$preprocess_check){
      show_spinner()
      if(!is.null(res_test())){
        hide_spinner()
        shinyalert("Complete Tests!", type="success", timer = 10000,
                   closeOnClickOutside = T, closeOnEsc = T)
        info <- paste0("* Test Method : ", input$test_method,"\n",
                       "* P.adj Method : ", input$padj_method)
        newTL <- data.frame(step="DEA_Test",info=info,
                            sample_num=as.numeric(nrow(assay(res_test()))),
                            time=as.character(Sys.time()),color="green",icon="chart-bar")
        timeLine <<- rbind(timeLine,newTL)
        addTimeLine(timeLine)
      }
    } else {
      shinyalert("Missing preprocessing step", type="error", timer = 10000,
                 closeOnClickOutside = T, closeOnEsc = T)
    }
  })
  
  
  observeEvent(input$dea_btn, {
    req(dep())
    updatePrettySwitch(session, "show_sampleID", label = "Show SampleID", value = F)
    
    shinyalert("Start Differential\nexperimental analysis!","It will be start after calculate\nthe best k of kmeans cluster for heatmap", type="info", timer = 10000,
               closeOnClickOutside = T, closeOnEsc = T)
    sig <- which(rowData(dep())$significant==T)
    output$volcano_plot <- NULL
    output$pca_plot <- NULL
    output$correlation_matrix <- NULL
    output$heatmap <- NULL
    
    if(!is.null(ready_for_dea())){
      shinyjs::show("show_sampleID")
      shinyjs::show("pca_plot")
      shinyjs::show("download_pca")
      shinyjs::show("volcano_plot")
      shinyjs::show("download_volcano")
      
      
      output$volcano_plot <- renderPlot({
        volcano_input()
      })
      output$pca_plot <- renderPlot({
        pca_input_noSample()
      })
    }
    
    if(length(sig) == 0){
      shinyalert("There is no DEP!", "Change threshold value or threshold type", type="error", timer = 10000,
                 closeOnClickOutside = T, closeOnEsc = T)
 
      
      if(input$thres_type == "none"){
        info <- paste0("* Threshold type : ", input$thres_type,"\n",
                       "* There is no DEP")
      } else{
        info <- paste0("* Threshold type : ", input$thres_type,"\n",
                       "* Threshold value for\n  ", 
                       input$thres_type, " : ", input$dea_pvalue,"\n",
                       "  log2FC : ", input$dea_log2fc,"\n",
                       "* There is no DEP")
      }
    }else{
      best_k <- length(sig)
      try(best_k <- optimize_k_input())

      if((input$dea_clusterNum==1) & is.null(best_k)){
        if(length(sig)>2){
          best_k <- 2
        }else{
          best_k <- 1
        }
        updateNumericInput(session, "dea_clusterNum", value=as.numeric(best_k))
      } else if((input$dea_clusterNum==1) & !is.null(best_k)){
        updateNumericInput(session, "dea_clusterNum", value=as.numeric(best_k))
      } 
      
      updateNumericInput(session, "input_num_toppi", value=length(sig))
      
      output$correlation_matrix <- renderPlot({
        correlation_input()
      })
      output$heatmap <- renderPlot({
        heatmap_input()
      })
      
      withProgress(message = 'Plots calculations are in progress',
                   detail = 'Please wait for a while', value = 0, {
                     for (i in 1:20) {
                       incProgress(1/20)
                       Sys.sleep(0.25)
                     }
                   })
      
      shinyalert("Complete Visualization!", type="success", timer = 10000,
                 closeOnClickOutside = T, closeOnEsc = T)
      
      
      shinyjs::show("correlation_matrix")
      shinyjs::show("download_correlation")
      shinyjs::show("heatmap")
      shinyjs::show("download_heatmap")
      shinyjs::show("download_gene_cluster")
      shinyjs::show("download_dep_info_btn")
      
      
      if(input$thres_type == "none"){
        info <- paste0("* Threshold type : ", input$thres_type,"\n",
                       "* The number of cluster : ", best_k)
      } else{
        info <- paste0("* Threshold type : ", input$thres_type,"\n",
                       "* Threshold value for\n  ", 
                       input$thres_type, " : ", input$dea_pvalue,"\n",
                       "  log2FC : ", input$dea_log2fc,"\n",
                       "* The number of cluster : ", best_k)
      }
    }
    
    
    newTL <- data.frame(step="DEA_Visualization",info=info,
                        sample_num=paste0(length(sig)," / ",as.numeric(nrow(assay(res_test())))),
                        time=as.character(Sys.time()),color="green",icon="chart-bar")
    timeLine <<- rbind(timeLine,newTL)
    addTimeLine(timeLine)
  })
  
  output$download_log2_exp_btn <- downloadHandler(
    filename = function() {paste0(input$file_type,"_log2transform_Expression_data_", Sys.Date(), ".csv")},
    content = function(file) {
      if(!is.null(main_data())){
        data <- main_data()
        write.csv(data,file,row.names=F,quote=F)
      }
    }
  )
  
  output$download_preprocessed_exp_btn <- downloadHandler(
    filename = function() {paste0(input$file_type,"_proprocessed_Expression_data_", Sys.Date(), ".csv")},
    content = function(file) {
      if(!is.null(ready_for_dea())){
        data_se <- ready_for_dea()
        res <- assay(data_se)
        write.csv(res,file,quote=F)
      }
    }
  )
  
  output$download_dep_info_btn <- downloadHandler(
    filename = function() {paste0(input$file_type,"_DEP_info_data_", Sys.Date(), ".csv")},
    content = function(file) {
      if(!is.null(dep())){
        write.csv(rowData(dep()),file,row.names=F,quote=F)
      }
    }
  )
  
  observeEvent(input$show_sampleID ,{
    if(input$show_sampleID){
      output$pca_plot <- renderPlot({
        pca_input_Sample()
      })
    } else {
      output$pca_plot <- renderPlot({
        pca_input_noSample()
      })
    }
  })
  
  output$download_pca <- downloadHandler(
    filename = function() {paste0("PCA_plot_", Sys.Date(), ".png")},
    content = function(file) {
      if(input$show_sampleID){
        if(!is.null(pca_input_Sample())){
          png(file,width=2000,height=2000,res=300)
          print(pca_input_Sample())
          dev.off()
        } 
      }else{
        if(!is.null(pca_input_noSample())){
          png(file,width=2000,height=2000,res=300)
          print(pca_input_noSample())
          dev.off()
        } 
      }
    }
  )
  
  output$download_heatmap <- downloadHandler(
    filename = function() {paste0("Heatmap_", Sys.Date(), ".png")},
    content = function(file) {
      if(!is.null(heatmap_input())){
        png(file,width=2000,height=2000,res=300)
        print(heatmap_input())
        dev.off()
      }
    }
  )
  
  output$download_gene_cluster <- downloadHandler(
    filename = function() {paste0("Heatmap_gene_cluster_info_",Sys.Date(),".csv")},
    content = function(file){
      filtered = subset(rowData(dep()), cut = significant == T)
      if(!is.null(filtered)){
        res <- get_gene_cluster(assay(dep()), rowData(dep()), input$dea_clusterNum)
        write.csv(res,file,row.names = F, quote = F)
      } 
    }
  )
  
  output$download_correlation <- downloadHandler(
    filename = function() {paste0("Correlation_plot_", Sys.Date(), ".png")},
    content = function(file) {
      if(!is.null(correlation_input())){
        png(file,width=2000,height=2000,res=300)
        print(correlation_input())
        dev.off()
      }
    }
  )
  
  output$download_volcano <- downloadHandler(
    filename = function() {paste0("Volcano_plot_", Sys.Date(), ".png")},
    content = function(file) {
      if(!is.null(volcano_input())){
        png(file,width=2000,height=2000,res=300)
        print(volcano_input())
        dev.off()
      }
    }
  )
  
  ################gsa_rendering###################
  
  observeEvent(input$gsa_btn,{
    setwd(base_dir)
    output$gobp_gsa_plot <- NULL
    output$gocc_gsa_plot <- NULL
    output$gomf_gsa_plot <- NULL
    output$kegg_gsa_plot <- NULL
    
    if(!is.null(dep())){
      shinyalert("Start GSA!","Please wait for a while", type="info", timer = 10000,
                 closeOnClickOutside = T, closeOnEsc = T)
      show_spinner()
      result_gsa()
      
      output$gobp_gsa_plot <- renderPlot({
        Sys.sleep(0.1)
        plot_gsa_gobp()
      })
      
      output$gocc_gsa_plot <- renderPlot({
        Sys.sleep(0.1)
        plot_gsa_gocc()
      })
      
      output$gomf_gsa_plot <- renderPlot({
        Sys.sleep(0.1)
        plot_gsa_gomf()
      })
      
      output$kegg_gsa_plot <- renderPlot({
        Sys.sleep(0.1)
        plot_gsa_kegg()
      })
      
      kegg_info <- reverted_gsa_kegg()
      pathway_choices <- kegg_info$Term
      updateSelectInput(session, "pathID_selector",
                        choices = pathway_choices, selected = "")
      
      output$topOfKeggDT <- DT::renderDataTable({
        DT::datatable(top_of_kegg(), options = list(scrollX = TRUE, pageLength = 5,lengthMenu = c(5, 10, 15)),
                      selection ="single") 
      }, server=T)
      
      hide_spinner()
      shinyalert("Complete GSA!","", type="success", timer = 10000,
                 closeOnClickOutside = T, closeOnEsc = T)
      
      shinyjs::show("gobp_gsa_plot")
      shinyjs::show("download_gsa_gobp_csv")
      shinyjs::show("download_gsa_gobp_png")
      shinyjs::show("gocc_gsa_plot")
      shinyjs::show("download_gsa_gocc_csv")
      shinyjs::show("download_gsa_gocc_png")
      shinyjs::show("gomf_gsa_plot")
      shinyjs::show("download_gsa_gomf_csv")
      shinyjs::show("download_gsa_gomf_png")
      shinyjs::show("kegg_gsa_plot")
      shinyjs::show("download_gsa_kegg_csv")
      shinyjs::show("download_gsa_kegg_png")
      
      shinyjs::show("topOfKeggDT")
      shinyjs::show("download_pathview")
      shinyjs::show("pathview_result")
      
      info <- paste0("* Data set : ", input$gsa_set,"\n",
                     "* stats of gene level : ", input$gsa_input_set)
      newTL <- data.frame(step="GSA",info=info,
                          sample_num=paste0(length(sig <- which(rowData(dep())$significant==T))," / ",as.numeric(nrow(assay(dep())))),
                          time=as.character(Sys.time()),color="yellow",icon="chart-bar")
      timeLine <<- rbind(timeLine,newTL)
      addTimeLine(timeLine)
      
    }else{
      shinyalert("There is no data!", "Check your data or change stats of gene level option", type="error", timer = 10000,
                 closeOnClickOutside = T, closeOnEsc = T)
    }
    
    
  })
  
  output$download_gsa_gobp_csv <- downloadHandler(
    filename = function() {paste0("GSA_EnrichR_GOBP_result_", Sys.Date(), ".csv")},
    content = function(file) {
      if(!is.null(result_gsa_gobp())){
        write.csv(result_gsa_gobp(),file,row.names=F,quote=F)
      }
    }
  )
  
  output$download_gsa_gobp_png <- downloadHandler(
    filename = function() {paste0("GSA_EnrichR_GOBP_result_", Sys.Date(), ".png")},
    content = function(file) {
      if(!is.null(result_gsa_gobp())){
        res_gobp <- result_gsa_gobp()
        res_gobp$P.value2 <- -log10(res_gobp$P.value)
        res_gobp <- res_gobp[order(-res_gobp$P.value2),]
        res_gobp <- res_gobp[c(1:input$set_nterm),]
        x_num <- max(nchar(res_gobp$Term))
        width <- x_num*40
        y_num <- input$set_nterm
        height <- y_num*80
        png(file,width=width,height=height,res=150)
        print(plot_gsa_gobp())
        dev.off()
      }
    }
  )
  
  output$download_gsa_gocc_csv <- downloadHandler(
    filename = function() {paste0("GSA_EnrichR_GOCC_result_", Sys.Date(), ".csv")},
    content = function(file) {
      if(!is.null(result_gsa_gocc())){
        write.csv(result_gsa_gocc(),file,row.names=F,quote=F)
      }
    }
  )
  
  output$download_gsa_gocc_png <- downloadHandler(
    filename = function() {paste0("GSA_EnrichR_GOCC_result_", Sys.Date(), ".png")},
    content = function(file) {
      if(!is.null(result_gsa_gocc())){
        res_gocc <- result_gsa_gocc()
        res_gocc$P.value2 <- -log10(res_gocc$P.value)
        res_gocc <- res_gocc[order(-res_gocc$P.value2),]
        res_gocc <- res_gocc[c(1:input$set_nterm),]
        x_num <- max(nchar(res_gocc$Term))
        width <- x_num*40
        y_num <- input$set_nterm
        height <- y_num*80
        png(file,width=width,height=height,res=150)
        print(plot_gsa_gocc())
        dev.off()
      }
    }
  )
  
  output$download_gsa_gomf_csv <- downloadHandler(
    filename = function() {paste0("GSA_EnrichR_GOMF_result_", Sys.Date(), ".csv")},
    content = function(file) {
      if(!is.null(result_gsa_gomf())){
        write.csv(result_gsa_gomf(),file,row.names=F,quote=F)
      }
    }
  )
  
  output$download_gsa_gomf_png <- downloadHandler(
    filename = function() {paste0("GSA_EnrichR_GOMF_result_", Sys.Date(), ".png")},
    content = function(file) {
      if(!is.null(result_gsa_gomf())){
        res_gomf <- result_gsa_gomf()
        res_gomf$P.value2 <- -log10(res_gomf$P.value)
        res_gomf <- res_gomf[order(-res_gomf$P.value2),]
        res_gomf <- res_gomf[c(1:input$set_nterm),]
        x_num <- max(nchar(res_gomf$Term))
        width <- x_num*40
        y_num <- input$set_nterm
        height <- y_num*80
        png(file,width=width,height=height,res=150)
        print(plot_gsa_gomf())
        dev.off()
      }
    }
  )
  
  output$download_gsa_kegg_csv <- downloadHandler(
    filename = function() {paste0("GSA_EnrichR_KEGG_result_", Sys.Date(), ".csv")},
    content = function(file) {
      if(!is.null(result_gsa_kegg())){
        write.csv(result_gsa_kegg(),file,row.names=F,quote=F)
      }
    }
  )
  
  output$download_gsa_kegg_png <- downloadHandler(
    filename = function() {paste0("GSA_EnrichR_KEGG_result_", Sys.Date(), ".png")},
    content = function(file) {
      if(!is.null(result_gsa_kegg())){
        res_kegg <- result_gsa_kegg()
        res_kegg$P.value2 <- -log10(res_kegg$P.value)
        res_kegg <- res_kegg[order(-res_kegg$P.value2),]
        res_kegg <- res_kegg[c(1:input$set_nterm),]
        x_num <- max(nchar(res_kegg$Term))
        width <- x_num*40
        y_num <- input$set_nterm
        height <- y_num*80
        png(file,width=width,height=height,res=150)
        print(plot_gsa_kegg())
        dev.off()
      }
    }
  )

  observeEvent(input$pathID_selector, {
    if(!is.null(input$topOfKeggDT_rows_selected)){
      dtProxy = DT::dataTableProxy("topOfKeggDT", session=session)
      DT::reloadData(dtProxy, clearSelection=c("all"))
    }
  })
  
  observeEvent(input$render_pathway_btn,{
    setwd(base_dir)
    if(input$pathID_selector != "") {
      show_spinner()
      
      dir <- getwd()
      dir <- strsplit(dir,"/",fixed = T)
      dir <- as.character(unlist(dir))
      dir <- dir[length(dir)]
      if(dir != "GSA"){
        setwd("GSA")
      }
      
      outfile <- pathway_graph()
      hide_spinner()
      output$pathview_result <- renderImage({
        list(src=paste0("./",outfile), contentType="image/png",
             width="100%", height="100%",
             alt="Pathview_graph")}, deleteFile=F)
    } 
  })
  
  observeEvent(input$topOfKeggDT_rows_selected, {
    if(input$pathID_selector!=""){
      req(input$pathID_selector)
      updateSelectInput(session, "pathID_selector",
                        selected = "")
      dtProxy = DT::dataTableProxy("topOfKeggDT", session=session)
      DT::selectRows(dtProxy, selected=input$topOfKeggDT_rows_selected)
    }
    show_spinner()
    outfile <- pathway_graph()
    hide_spinner()
    output$pathview_result <- renderImage({
      list(src=outfile, contentType="image/png",
           width="100%", height="100%",
           alt="Pathview_graph")}, deleteFile=F)
  })
  
  observeEvent(input$zoom_pathway_btn,{
    setwd(base_dir)
    
    dir <- getwd()
    dir <- strsplit(dir,"/",fixed = T)
    dir <- as.character(unlist(dir))
    dir <- dir[length(dir)]
    if(dir != "GSA"){
      setwd("GSA")
    }
    
    showModal(modalDialog(
      renderImage({
        outfile <- pathway_graph()
        list(src=outfile, contentType="image/png",
             alt="Pathview_graph")}, deleteFile=F),
      footer=actionButton("close_modal", label="Close")
    ))
  })
  
  observeEvent(input$close_modal, {
    removeModal()
  })
 
  output$download_pathview <- downloadHandler(
    filename = function() {
      paste0("GSA_Pathview_",user_id,".zip")},
    content = function(file) {
      setwd(base_dir)
      dir <- getwd()
      dir <- strsplit(dir,"/",fixed = T)
      dir <- as.character(unlist(dir))
      dir <- dir[length(dir)]
      if(dir != "GSA"){
        setwd("GSA")
      }
      
      png_list <- list.files(path="./", pattern=".png")
      png_list <- png_list[grep("gsa",png_list,fixed = T)]
      file_set <- c()
      for(i in 1:length(png_list)){
        file_path <- paste0(png_list[i])
        file_set <- c(file_set, file_path)
      }
      zip(file,file_set)
    }
  )
  
  
  ########################GSEA_rendering###########################
  observeEvent(input$gsea_btn, {
    setwd(base_dir)
    output$gobp_gsea_plot <- NULL
    output$gocc_gsea_plot <- NULL
    output$gomf_gsea_plot <- NULL
    output$kegg_gsea_plot <- NULL
    
    if(!is.null(res_test())){
      shinyalert("Start GSEA!","Please wait for a while", type="info", timer = 10000,
                 closeOnClickOutside = T, closeOnEsc = T)
      show_spinner()
      result_gsea()
      
      output$gobp_gsea_plot <- renderPlot({
        Sys.sleep(0.5)
        plot_gsea_gobp()
      })
      
      output$gocc_gsea_plot <- renderPlot({
        Sys.sleep(0.5)
        plot_gsea_gocc()
      })
      
      output$gomf_gsea_plot <- renderPlot({
        Sys.sleep(0.5)
        plot_gsea_gomf()
      })
      
      output$kegg_gsea_plot <- renderPlot({
        Sys.sleep(0.5)
        plot_gsea_kegg()
      })
      
      shinyalert("Complete GSEA!","", type="success", timer = 10000,
                 closeOnClickOutside = T, closeOnEsc = T)
      
      gsea_res <- reverted_gsea_kegg()
      
      up_res <- gsea_res[gsea_res$Enrichment=="Up-regulated",]
      up_res <- up_res[,-ncol(up_res)]
      print(up_res)
      down_res <- gsea_res[gsea_res$Enrichment=="Down-regulated",]
      down_res <- down_res[,-ncol(down_res)]
      print(down_res)
      
      up_res$kegg_id <- as.character(up_res$kegg_id)
      down_res$kegg_id <- as.character(down_res$kegg_id)
      
      output$result_of_gsea_up_regulated <- DT::renderDataTable({
        DT::datatable(up_res,options = list(scrollX = TRUE, pageLength = 5,lengthMenu = c(5, 10, 15)), 
                      escape = T, selection ="single") %>% DT::formatRound(colnames(up_res)[3:ncol(up_res)], digits=3)
      },server = T) 
      
      output$result_of_gsea_down_regulated <- DT::renderDataTable({
        DT::datatable(down_res,options = list(scrollX = TRUE, pageLength = 5,lengthMenu = c(5, 10, 15)), 
                      escape = T, selection ="single") %>% DT::formatRound(colnames(down_res)[3:ncol(down_res)], digits=3)
      },server = T) 
      
      hide_spinner()
      
      shinyjs::show("gobp_gsea_plot")
      shinyjs::show("download_gsea_gobp_csv")
      shinyjs::show("download_gsea_gobp_png")
      shinyjs::show("gocc_gsea_plot")
      shinyjs::show("download_gsea_gocc_csv")
      shinyjs::show("download_gsea_gocc_png")
      shinyjs::show("gomf_gsea_plot")
      shinyjs::show("download_gsea_gomf_csv")
      shinyjs::show("download_gsea_gomf_png")
      shinyjs::show("kegg_gsea_plot")
      shinyjs::show("download_gsea_kegg_csv")
      shinyjs::show("download_gsea_kegg_png")
      
      shinyjs::show("result_of_gsea_up_regulated")
      shinyjs::show("result_of_gsea_down_regulated")
      shinyjs::show("download_gsea_pathview")
      
      gsea_select <- input$select_genelevel_stats
      if(gsea_select == "log2fcxmlog10padj"){
        gsea_select <- "log<sub>2</sub>(FC)*-log<sub>10</sub>(P.adj)"
      }
      info <- paste0("* stats of gene level : ", gsea_select)
      newTL <- data.frame(step="GSEA",info=info,
                          sample_num=as.numeric(nrow(assay(dep()))),
                          time=as.character(Sys.time()),color="yellow",icon="chart-bar")
      timeLine <<- rbind(timeLine,newTL)
      addTimeLine(timeLine)
    } else{
      shinyalert("There is no Data!", "Chack your data again", type="error", timer = 10000,
                 closeOnClickOutside = T, closeOnEsc = T)
    }
  })
  
  output$download_gsea_gobp_csv <- downloadHandler(
    filename = function() {paste0("GSEA_fgsea_GOBP_result_", Sys.Date(), ".csv")},
    content = function(file) {
      if(!is.null(result_gsea_gobp())){
        write.csv(result_gsea_gobp(),file,row.names=F,quote=F)
      }
    }
  )
  
  output$download_gsea_gobp_png <- downloadHandler(
    filename = function() {paste0("GSEA_fgsea_GOBP_result_", Sys.Date(), ".png")},
    content = function(file) {
      if(!is.null(result_gsea_gobp())){
        res_gobp <- result_gsea_gobp()
        x_num <- max(nchar(res_gobp$pathway))
        width <- x_num*40
        if(nrow(res_gobp) >= 20){
          png(file,width=width,height=1800,res=150)
        }else{
          height <- nrow(res_gobp)*90
          png(file,width=2000,height=height,res=150)
        }
        print(plot_gsea_gobp())
        dev.off()
      }
    }
  )
  
  output$download_gsea_gocc_csv <- downloadHandler(
    filename = function() {paste0("GSEA_fgsea_GOCC_result_", Sys.Date(), ".csv")},
    content = function(file) {
      if(!is.null(result_gsea_gocc())){
        write.csv(result_gsea_gocc(),file,row.names=F,quote=F)
      }
    }
  )
  
  output$download_gsea_gocc_png <- downloadHandler(
    filename = function() {paste0("GSEA_fgsea_GOCC_result_", Sys.Date(), ".png")},
    content = function(file) {
      if(!is.null(result_gsea_gocc())){
        res_gocc <- result_gsea_gocc()
        x_num <- max(nchar(res_gocc$pathway))
        width <- x_num*40
        if(nrow(res_gocc) >= 20){
          png(file,width=width,height=1800,res=150)
        }else{
          height <- nrow(res_gocc)*90
          png(file,width=2000,height=height,res=150)
        }
        print(plot_gsea_gocc())
        dev.off()
      }
    }
  )
  
  output$download_gsea_gomf_csv <- downloadHandler(
    filename = function() {paste0("GSEA_fgsea_GOMF_result_", Sys.Date(), ".csv")},
    content = function(file) {
      if(!is.null(result_gsea_gomf())){
        write.csv(result_gsea_gomf(),file,row.names=F,quote=F)
      }
    }
  )
  
  output$download_gsea_gomf_png <- downloadHandler(
    filename = function() {paste0("GSEA_fgsea_GOMF_result_", Sys.Date(), ".png")},
    content = function(file) {
      if(!is.null(result_gsea_gomf())){
        res_gomf <- result_gsea_gomf()
        x_num <- max(nchar(res_gomf$pathway))
        width <- x_num*40
        if(nrow(res_gomf) >= 20){
          png(file,width=width,height=1800,res=150)
        }else{
          height <- nrow(res_gomf)*90
          png(file,width=2000,height=height,res=150)
        }
        print(plot_gsea_gomf())
        dev.off()
      }
    }
  )
  
  output$download_gsea_kegg_csv <- downloadHandler(
    filename = function() {paste0("GSEA_fgsea_KEGG_result_", Sys.Date(), ".csv")},
    content = function(file) {
      if(!is.null(result_gsea_kegg())){
        write.csv(result_gsea_kegg(),file,row.names=F,quote=F)
      }
    }
  )
  
  output$download_gsea_kegg_png <- downloadHandler(
    filename = function() {paste0("GSEA_fgsea_KEGG_result_", Sys.Date(), ".png")},
    content = function(file) {
      if(!is.null(result_gsea_kegg())){
        res_kegg <- result_gsea_kegg()
        x_num <- max(nchar(res_kegg$pathway))
        width <- x_num*40
        if(nrow(res_kegg) >= 20){
          png(file,width=width,height=1800,res=150)
        }else{
          height <- nrow(res_kegg)*90
          png(file,width=2000,height=height,res=150)
        }
        print(plot_gsea_kegg())
        dev.off()
      }
    }
  )
  
  observeEvent(input$gsea_gobp_click,{
      res_gobp <- result_gsea_gobp()
      if(nrow(res_gobp) >= 20){
        filt_gobp <-  rbind(head(res_gobp, n = 10),tail(res_gobp, n = 10 ))  
      } else{
        filt_gobp <-  res_gobp
      }
      
      filt_gobp$"0" <- filt_gobp$NES
      filt_gobp$"reorder(pathway, NES)" <- filt_gobp$pathway
      df <- nearPoints(filt_gobp, input$gsea_gobp_click)
      
      pathway <- df$pathway
      myGO <- make_myGO("gobp")
      gene_list <- make_geneRank(rowData(res_test()),input$select_genelevel_stats)

      output$plotEnrichment <- renderPlot({
        plotEnrichment(myGO[[pathway]], gene_list) + labs(title=pathway)
      })
      
      showModal(modalDialog(
        renderUI({plotOutput("plotEnrichment", width = "800px",height="700px")}),
        footer=fluidRow(actionButton("close_modal", label="Close"))
      ))
  })
  
  observeEvent(input$gsea_gocc_click,{
    res_gocc <- result_gsea_gocc()
    if(nrow(res_gocc) >= 20){
      filt_gocc <-  rbind(head(res_gocc, n = 10),tail(res_gocc, n = 10 ))  
    } else{
      filt_gocc <-  res_gocc
    }
    
    filt_gocc$"0" <- filt_gocc$NES
    filt_gocc$"reorder(pathway, NES)" <- filt_gocc$pathway
    df <- nearPoints(filt_gocc, input$gsea_gocc_click)
    
    pathway <- df$pathway
    myGO <- make_myGO("gocc")
    gene_list <- make_geneRank(rowData(res_test()),input$select_genelevel_stats)
    
    output$plotEnrichment <- renderPlot({
      plotEnrichment(myGO[[pathway]], gene_list) + labs(title=pathway) 
    })
    
    showModal(modalDialog(
      renderUI({plotOutput("plotEnrichment", width = "800px",height="700px")}),
      footer=fluidRow(actionButton("close_modal", label="Close"))
    ))
  })
  
  observeEvent(input$gsea_gomf_click,{
    res_gomf <- result_gsea_gomf()
    if(nrow(res_gomf) >= 20){
      filt_gomf <-  rbind(head(res_gomf, n = 10),tail(res_gomf, n = 10 ))  
    } else{
      filt_gomf <-  res_gomf
    }
    
    filt_gomf$"0" <- filt_gomf$NES
    filt_gomf$"reorder(pathway, NES)" <- filt_gomf$pathway
    df <- nearPoints(filt_gomf, input$gsea_gomf_click)
    
    pathway <- df$pathway
    myGO <- make_myGO("gomf")
    gene_list <- make_geneRank(rowData(res_test()),input$select_genelevel_stats)
    
    output$plotEnrichment <- renderPlot({
      plotEnrichment(myGO[[pathway]], gene_list) + labs(title=pathway)
    })
    
    showModal(modalDialog(
      renderUI({plotOutput("plotEnrichment", width = "800px",height="700px")}),
      footer=fluidRow(actionButton("close_modal", label="Close"))
    ))
  })
  
  observeEvent(input$gsea_kegg_click,{
    res_kegg <- result_gsea_kegg()
    if(nrow(res_kegg) >= 20){
      filt_kegg <-  rbind(head(res_kegg, n = 10),tail(res_kegg, n = 10 ))  
    } else{
      filt_kegg <-  res_kegg
    }
    
    filt_kegg$"0" <- filt_kegg$NES
    filt_kegg$"reorder(pathway, NES)" <- filt_kegg$pathway
    df <- nearPoints(filt_kegg, input$gsea_kegg_click)
    
    pathway <- df$pathway
    myGO <- make_myGO("kegg")
    gene_list <- make_geneRank(rowData(res_test()),input$select_genelevel_stats)
    
    output$plotEnrichment <- renderPlot({
      plotEnrichment(myGO[[pathway]], gene_list) + labs(title=pathway)
    })
    
    showModal(modalDialog(
      renderUI({plotOutput("plotEnrichment", width = "800px",height="700px")}),
      footer=fluidRow(actionButton("close_modal", label="Close"))
    ))
  })
  
  observeEvent(input$result_of_gsea_up_regulated_rows_selected, {
    dtProxy = DT::dataTableProxy("result_of_gsea_up_regulated", session=session)
    DT::selectRows(dtProxy, selected=input$result_of_gsea_up_regulated_rows_selected)
    
    show_spinner()
    gsea_outfile <- gsea_up_pathway_graph()
    hide_spinner()
    
    setwd(base_dir)
    dir <- getwd()
    dir <- strsplit(dir,"/",fixed = T)
    dir <- as.character(unlist(dir))
    dir <- dir[length(dir)]
    if(dir != "GSEA"){
      setwd("GSEA")
    }
    
    showModal(modalDialog(
      renderImage({
        list(src=gsea_outfile, contentType="image/png",
             alt="Pathview_gesa_graph")}, deleteFile=F),
      footer=actionButton("close_modal", label="Close")
    ))
    
  })
  
  observeEvent(input$result_of_gsea_down_regulated_rows_selected, {
    dtProxy = DT::dataTableProxy("result_of_gsea_down_regulated", session=session)
    DT::selectRows(dtProxy, selected=input$result_of_gsea_down_regulated_rows_selected)
    
    show_spinner()
    gsea_outfile <- gsea_down_pathway_graph()
    hide_spinner()
    
    setwd(base_dir)
    dir <- getwd()
    dir <- strsplit(dir,"/",fixed = T)
    dir <- as.character(unlist(dir))
    dir <- dir[length(dir)]
    if(dir != "GSEA"){
      setwd("GSEA")
    }
    
    showModal(modalDialog(
      renderImage({
        list(src=gsea_outfile, contentType="image/png",
             alt="Pathview_gesa_graph")}, deleteFile=F),
      footer=actionButton("close_modal", label="Close")
    ))
  })
    
  observeEvent(input$close_modal, {
    removeModal()
  })
  
  output$download_gsea_pathview <- downloadHandler(
    filename = function() {
      paste0("GSEA_Pathview_",user_id,".zip")},
    content = function(file) {
      setwd(base_dir)
      dir <- getwd()
      dir <- strsplit(dir,"/",fixed = T)
      dir <- as.character(unlist(dir))
      dir <- dir[length(dir)]
      if(dir != "GSEA"){
        setwd("GSEA")
      }
      
      png_list <- list.files(path="./", pattern=".png")
      png_list <- png_list[grep("gsea",png_list,fixed = T)]
      file_set <- c()
      for(i in 1:length(png_list)){
        file_path <- paste0("./",png_list[i])
        file_set <- c(file_set, file_path)
      }
      zip(file,file_set)
    }
  )
  

  ####################################ppi_Rendering###############################################

  observeEvent(input$input_num_toppi,{
    req(dep())
    dep <- rowData(dep())
    dep <- dep[dep$significant==T,]
    num <- input$input_num_toppi
    if(num > nrow(dep)| is.na(input$input_num_toppi)){
      shinyalert("Input wrong number!","You can input a number between 1 and the number of DEP.", type="error", timer = 10000,
                 closeOnClickOutside = T, closeOnEsc = T)
    }
    output$ppi_image <- NULL
  })
  
  observeEvent(input$select_ppi_condition, {
    ppi_condition <- input$select_ppi_condition
    if(ppi_condition == "Gene_Name"){
      shinyjs::show("input_gene_toppi")
      shinyjs::hide("input_num_toppi")
    }else{
      shinyjs::hide("input_gene_toppi")
      shinyjs::show("input_num_toppi")
    }
    output$ppi_image <- NULL
  })
  
  observeEvent(input$ppi_btn, {
    if(input$input_num_toppi == 0 | is.na(input$input_num_toppi)){
      shinyjs::hide("ppi_box")
      shinyalert("There is no DEP!", "You should input DEP. Change threshold value or threshold type at DEA step.", type="error", timer = 10000,
                 closeOnClickOutside = T, closeOnEsc = T)
      hide_spinner()
    } else{
      shinyjs::show("ppi_box")
      shinyalert("Start PPI Network anlaysis!","please wait for for a while", type="info", timer = 10000,
                 closeOnClickOutside = T, closeOnEsc = T)
      show_spinner()
      output$ppi_image <- renderUI(div(
        tags$img(src=string_url(),
                 id="ppiImage",
                 align="left",
                 style = "position: center;"
        ))
      )
      Sys.sleep(5)
      shinyjs::show("ppi_image_download_btn")
      shinyjs::show("ppi_tsv_download_btn")
      
      hide_spinner()
      shinyalert("Compelte PPI Network anlaysis!","please wait for redering image", type="success", timer = 10000,
                 closeOnClickOutside = T, closeOnEsc = T)
      
      
      
      info <- paste0("* number of input gene : ", input$input_num_toppi,"\n",
                     "* condition of gene select\n  : ", input$select_ppi_condition)
      if(input$select_ppi_condition == "Gene_Name"){
        info <- paste0(info,"\n",
                       "  => ", input$input_gene_toppi)
      }
      
      newTL <- data.frame(step="PPI network",info=info,
                          sample_num=paste0(length(sig <- which(rowData(dep())$significant==T))," / ",as.numeric(nrow(assay(dep())))),
                          time=as.character(Sys.time()),color="yellow",icon="chart-bar")
      timeLine <<- rbind(timeLine,newTL)
      addTimeLine(timeLine)
    }
  })
  
  output$ppi_image_download_btn <- downloadHandler(
    filename = paste0("PPI_network_STRINGDB_", Sys.Date(), ".jpg"),
    content = function(file) {
      GET(string_image_download_url(), write_disk(file))
    }
  )
  
  output$ppi_tsv_download_btn <- downloadHandler(
    filename = paste0("PPI_network_STRINGDB_", Sys.Date(), ".tsv"),
    content = function(file) {
      
      GET(string_tsv_download_url(), write_disk(tmp <- tempfile(fileext = ".tsv")))
      tmp_tsv <- read.delim(tmp)
      colnames(tmp_tsv)[3] <- "name"

      sig_data <- rowData(dep())[rowData(dep())$significant==T,]

      tmp_tsv <- merge(tmp_tsv, sig_data, by="name")
      colnames(tmp_tsv)[1] <- "preferredName_A"
      tmp_tsv <- tmp_tsv[,c(2,3,1,4:13,19)]
      colnames(tmp_tsv)[ncol(tmp_tsv)] <- "log2FC_A"
      
      colnames(tmp_tsv)[4] <- "name"
      tmp_tsv <- merge(tmp_tsv, sig_data, by="name")
      colnames(tmp_tsv)[1] <- "preferredName_B"
      tmp_tsv <- tmp_tsv[,c(2,3,4,1,5:14,20)]
      colnames(tmp_tsv)[ncol(tmp_tsv)] <- "log2FC_B"
      
      write.table(tmp_tsv,file, row.names = F, quote = F, sep="\t")
    }
  )
  
  
  ####example file downlaod
  output$lfq_example_file_download_btn <- downloadHandler(
    filename = function() {paste0("LiBT_LFQ_example_file_", Sys.Date(), ".txt")},
    content = function(file) {
      data <- read.delim("base/LiBT_EXAMPLE_LFQ&iBaq_Proteins.txt")
      if(!is.null(data)){
        write.table(data,file,row.names=F,quote=F,sep="\t")
      }
    }
  )
  
  output$ibaq_example_file_download_btn <- downloadHandler(
    filename = function() {paste0("LiBT_iBAQ_example_file_", Sys.Date(), ".txt")},
    content = function(file) {
      data <- read.delim("base/LiBT_EXAMPLE_LFQ&iBaq_Proteins.txt")
      if(!is.null(data)){
        write.table(data,file,row.names=F,quote=F,sep="\t")
      }
    }
  )
  
  output$tmt_example_file_download_btn <- downloadHandler(
    filename = function() {paste0("LiBT_TMT_example_file_", Sys.Date(), ".txt")},
    content = function(file) {
      data <- read.delim("base/LiBT_EXAMPLE_TMT10_F_Proteins.txt")
      if(!is.null(data)){
        write.table(data,file,row.names=F,quote=F,sep="\t")
      }
    }
  )
  
  output$preprocessing4th <- renderImage({
    list(src = "base/figure1.png",
         alt = "preprocessing4th"
    )
  })
  
  output$pcaplotImg <- renderImage({
    list(src = "base/figure2.png",
         alt = "pcaplotImg"
    )
  })
  output$vcplotImg <- renderImage({
    list(src = "base/figure3.png",
         alt = "vcplotImg"
    )
  })
  output$corrplotImg <- renderImage({
    list(src = "base/figure4.png",
         alt = "corrplotImg"
    )
  })
  output$heatmapImg <- renderImage({
    list(src = "base/figure5.png",
         alt = "heatmapImg"
    )
  })
  output$barplotGSA <- renderImage({
    list(src = "base/figure6.png",
         alt = "barplotGSA"
    )
  })
  output$pathviewImgGSA <- renderImage({
    list(src = "base/figure7.png",
         alt = "pathviewImgGSA"
    )
  })
  output$erplotGSEA <- renderImage({
    list(src = "base/figure8.png",
         alt = "erplotImgGSA"
    )
  })
  output$pathviewImgGSEA <- renderImage({
    list(src = "base/figure9.png",
         alt = "pathviewImgGSEA"
    )
  })
  output$ppiImg <- renderImage({
    list(src = "base/figure10.png",
         alt = "ppiImg"
    )
  })
  ####################################################################################################################################
  ##--------------------------------------------------- reactive/EventReactive Section
  ####################################################################################################################################
  
  file_input <- reactive({NULL})
  file_input <- eventReactive(input$fileBrowser, {
    req(input$fileBrowser)
    if(is.null(input$fileBrowser)){
      return(NULL)
    }
    temp_df <- readLines(input$fileBrowser$datapath, n=1)
    
    if(grepl("\t", temp_df)){
      sep <- c("\t")
    } else if(grepl(";", temp_df)) {
      sep <- c(";")
    } else if(grepl(",", temp_df)) {
      sep <- c(",")
    } else {
      sep <- c(" ")
    }
    
    temp_df <- read.csv(input$fileBrowser$datapath,
                          header = T, fill = T,
                          sep = sep)

    state <- file_input_test(temp_df, input$file_type)
    if(length(state)==0){
      return(NULL)
    } else {
      return(temp_df)
    }
  })
  
  main_data <- reactive({NULL})
  main_data <- eventReactive(input$file_upload_btn, {
    temp_df <- file_input()
    file_type <- input$file_type
    tool_type <- input$TMT_tool_option
    if(file_type=="TMT" & tool_type=="PD"){
      temp_df <- get_main_data_PD(temp_df,input$TMT_input_option)
    }
    else {
      checked_option <- input$nonTMT_input_option
      temp_df <- filter_with_option(checked_option, temp_df)
      temp_df <- get_main_data_MQ(temp_df, file_type)
    }
  })
  
  total_samples <- reactive({
    df <- main_data()
    file_type <- input$file_type
    if(file_type=="TMT"){
      if(input$TMT_tool_option == "MQ"){
        samples <- make_case_samples_MQ(df, file_type)
      }else{
        samples <- make_case_samples_PD(df)
      }
    } else{
      samples <- make_case_samples_MQ(df, file_type)
    }
    return(samples)
  })
  
  case_samples <- reactive({
    case_samples <- input$case_group_selection
    return(case_samples)
  })
  
  control_toUpdate <- eventReactive(input$case_group_selection, {
    control_toUpdate <- setdiff(total_samples(), case_samples())
    return(control_toUpdate)
  })
  
  control_samples <- reactive({
    control_samples <- input$control_group_selection
    return(control_samples)
  })
  
  summarized_Data <- reactive({
    case <- case_samples()
    ctrl <- control_samples()
    
    if(input$file_type=="TMT"){
      design <- make_expDesignData(case, ctrl, input$TMT_tool_option)
    } else{
      design <- make_expDesignData(case, ctrl, "MQ")
    }
    
    main <- main_data()
    file_type <- input$file_type
    summary <- make_summarizedData(main,file_type,design)
    return(summary)
  })
  
  preprocessed_data <- reactive({
    file_type <- input$file_type
    transformation_item <- input$transformation
    filter_item <- input$valid_value
    imputation_item <- input$imputation
    normalization_item <- input$normalization
    
    case <- case_samples()
    control <- control_samples()
    
    data_se <- summarized_Data()
    
    data_filt <- NULL
    data_norm <- NULL
    data_imp <- NULL
    
    preprocessing_options <- input$use_options
    for(i in 1:length(preprocessing_options)) {
      print(preprocessing_options[i])
      switch(preprocessing_options[i],
             "Use_valid_value" = {
               data_se <- use_valid_option(data_se, case, control, input$valid_value)
               data_filt <- data_se},
             "Use_imputation" = {
               tryCatch({data_se <- use_imputation_option(data_se, case, control, input$imputation)
                          data_imp <- data_se},
                        warning = function(e) data_imp <- NULL)
             },
             "Use_normalization" = {
               data_se <- use_normalization_option(data_se,input$normalization)
               data_norm <- data_se},
             NULL = {
               data_se <- data_se
             }
      )
    }
    preprocessed_data <- list(data_se, data_filt, data_norm, data_imp)
    return (preprocessed_data) # format : SummarizedExperiment
  })
  
  ready_for_dea <- reactive({
    req(preprocessed_data())
    preprocessed_data()[[1]]
  })
  
  data_filt <- reactive({
    req(preprocessed_data())
    preprocessed_data()[[2]]
  })
  
  data_norm <- reactive({
    req(preprocessed_data())
    preprocessed_data()[[3]]
  })
  
  data_imp <- reactive({
    req(preprocessed_data())
    preprocessed_data()[[4]]
  })

  res_test <- eventReactive(input$test_btn,{
    req(preprocessed_data())
    data_diff <- test_diff(ready_for_dea(), type="all")
    if(input$test_method != "Limma"){
      diff_rowData <- rowData(data_diff)
      diff_colData <- colData(data_diff)
      
      data <- assay(data_diff)
      df <- test(data,diff_colData,input$test_method,input$padj_method)
      pval_pos <- grep("p.val",colnames(diff_rowData),fixed = T)
      padj_pos <- grep("_p.adj",colnames(diff_rowData),fixed = T)
      
      diff_rowData[,pval_pos] <- df$pval
      diff_rowData[,padj_pos] <- df$padj
      
      data_diff_edit <- SummarizedExperiment(assays = list(assay(data_diff)), rowData = diff_rowData, colData = diff_colData)
    } else{
      data_diff <- data_diff
    }
  })
  
  dep <- eventReactive(input$dea_btn,{
    req(res_test())
    type <- input$thres_type
    pvalue <- input$dea_pvalue
    log2fc <- input$dea_log2fc
    data_rejection <- c()
    
    if(type == "P.adj"){
      data_rejection <- add_rejections(res_test(), alpha=pvalue, lfc=log2fc)
      dep_rowData <- rowData(data_rejection)
      dep_rowData$name <- as.character(dep_rowData$name)
      data_rejection <- SummarizedExperiment(assays = list(assay(data_rejection)), rowData = dep_rowData, colData = colData(data_rejection))
    }
    else if(type == "P.value"){
      data_rejection <- add_rejections(res_test(), alpha=pvalue, lfc=log2fc)
      dep_rowData <- change_Sig(data_rejection,pvalue,log2fc)
      data_rejection <- SummarizedExperiment(assays = list(assay(data_rejection)), rowData = dep_rowData, colData = colData(data_rejection))
    }
    else {
      data_rejection <- add_rejections(res_test(), alpha=1, lfc=0)
      dep_rowData <- rowData(data_rejection)
      dep_rowData$name <- as.factor(dep_rowData$name)
      data_rejection <- SummarizedExperiment(assays = list(assay(data_rejection)), rowData = dep_rowData, colData = colData(data_rejection))
    } 
    return(data_rejection)
  })
  
  optimize_k_input <- reactive({
    req(dep())
    assay <- assay(dep())
    rowData <- rowData(dep())
    sig_gene <- as.character(rowData$name[rowData$significant==T])
    best_k <- get_optimize_k(assay, sig_gene)
    return(best_k)
  })
  
  volcano_input <- reactive({
    req(dep())
    condition <- dep()$condition
    case_name <- condition[1]
    ctrl_name <- condition[length(condition)]
    contrast <- paste0(case_name,"_vs_", ctrl_name)
    
    dep_rowData <- rowData(dep())
    pv_pos <- grep("p.val",colnames(dep_rowData),fixed = T)
    lfc_pos <- grep("diff",colnames(dep_rowData),fixed = T)
    
    input_vc <- data_frame(name=dep_rowData$name, lfc=dep_rowData[,lfc_pos], 
                           p=-log10(as.numeric(dep_rowData[,pv_pos])), sig=dep_rowData$significant)
    if(input$thres_type != "none"){
      plot_volcano(dep(), contrast=contrast, label_size=2, add_names=F) + 
        geom_point(data = filter(input_vc,sig), aes(lfc, p), color = "red", size= 2)+
        labs(x=expression(paste(log[2],FoldChange)),y=expression(paste(-log[10],P.value)))
    } else{
      plot_volcano(dep(), contrast=contrast, label_size=2, add_names=F) + 
        geom_point(data = filter(input_vc,sig), aes(lfc, p), color = "black", size= 2)+
        labs(x=expression(paste(log[2],FoldChange)),y=expression(paste(-log[10],P.value)))
    }
  })
  
  output$volcano_info <- DT::renderDataTable(DT::datatable({
    req(dep())
    dep_rowData <- rowData(dep())
    pv_pos <- grep("p.val",colnames(dep_rowData),fixed = T)
    lfc_pos <- grep("diff",colnames(dep_rowData),fixed = T)
    
    input_vc <- data.frame(name=dep_rowData$name, log2FC=dep_rowData[,lfc_pos], 
                           `-log10P.val`=-log10(as.numeric(dep_rowData[,pv_pos])), sig=dep_rowData$significant)
    input_vc <- data.frame(input_vc)
    colnames(input_vc) <- c("name", "log<sub>2</sub>FC", "-log<sub>10</sub>P.value", "sig")
    df <- brushedPoints(input_vc, input$volcano_brush, xvar="log<sub>2</sub>FC", yvar = "-log<sub>10</sub>P.value")
  }, options = list(scrollX = TRUE, pageLength = 5,lengthMenu = c(5, 10, 15)), escape = F, selection="none") %>%  DT::formatRound(c(2:3),digits=2))
  
  pca_input_noSample <- reactive({
    req(dep())
    set.seed(1234)
    n = nrow(assay(dep()))
    plot_pca(dep(), x=1, y=2, n=n, point_size=3, indicate = c("condition"))+
      ggtitle(paste("PCA plot -", n, "variable proteins", sep=" "))
  })
  
  pca_input_Sample <- reactive({
    req(dep())
    set.seed(1234)
    n = nrow(assay(dep()))
    pca_df <- get_pca_df(assay(dep()), colData(dep()), n)
    plot_pca(dep(), x=1, y=2, n=n, point_size=3, indicate = c("condition"))+
      ggtitle(paste("PCA plot -", n, "variable proteins", sep=" "))+
      geom_text_repel(data=pca_df,aes(label=rowname),nudge_y = 1)
  })
  
  correlation_input <- reactive({
    req(dep())
    plot_cor(dep(), significant = T, pal="Blues", lower=-1, upper=1, indicate=c("condition"))
  })
  
  heatmap_input <- reactive({
    req(dep())
    set.seed(1234)
    plot_heatmap(dep(), type="centered", kmeans=T, k=as.numeric(input$dea_clusterNum),
                 col_limit=2, show_row_names=F, indicate=c("condition"))
  })
  
  data_results <- reactive({
    req(dep())
    data_results <- get_results(dep())
  })
  
  ########################### GSA ######################################
  result_gsa <- reactive({
    data <- rowData(dep())
    res_gsa <- gsa(data,input$gsa_input_set,input$gsa_set)
    return(res_gsa)
  })
  
  result_gsa_gobp <- reactive({
    res_gsa <- result_gsa()
    if(!is.null(res_gsa)){
      gobp<-res_gsa[["GO_Biological_Process_2021"]]
    }
    return(gobp)
  })
  
  plot_gsa_gobp <- reactive({
    gobp <- result_gsa_gobp()
    gobp$P.value2 <- -log10(gobp$P.value)
    gobp <- gobp[order(-gobp$P.value2),]
    gobp <- gobp[c(1:input$set_nterm),]
    ggplot(data=gobp, aes(x=`P.value2`,y=reorder(`Term`,`P.value2`)))+
      geom_bar(stat = "identity",fill="#3c8dbc")+
      labs(title="GO_BP",x=expression(paste(-log[10],P.value)),y="")+
      theme_bw()+
      theme(axis.text=element_text(size=13),title = element_text(size=15,face="bold"), aspect.ratio = 0.5)
  })
  
  result_gsa_gocc <- reactive({
    res_gsa <- result_gsa()
    if(!is.null(res_gsa)){
      gocc<-res_gsa[["GO_Cellular_Component_2021"]]
    }
    return(gocc)
  })
  
  plot_gsa_gocc <- reactive({
    gocc <- result_gsa_gocc()
    gocc$P.value2 <- -log10(gocc$P.value)
    gocc <- gocc[order(-gocc$P.value2),]
    gocc <- gocc[c(1:input$set_nterm),]
    ggplot(data=gocc, aes(x=`P.value2`,y=reorder(`Term`,`P.value2`)))+
      geom_bar(stat = "identity",fill="#3c8dbc")+
      labs(title="GO_CC",x=expression(paste(-log[10],P.value)),y="")+
      theme_bw()+
      theme(axis.text=element_text(size=13),title = element_text(size=15,face="bold"), aspect.ratio = 0.5)
  })
  
  result_gsa_gomf <- reactive({
    res_gsa <- result_gsa()
    if(!is.null(res_gsa)){
      gomf<-res_gsa[["GO_Molecular_Function_2021"]]
    }
    return(gomf)
  })
  
  plot_gsa_gomf <- reactive({
    gomf <- result_gsa_gomf()
    gomf$P.value2 <- -log10(gomf$P.value)
    gomf <- gomf[order(-gomf$P.value2),]
    gomf <- gomf[c(1:input$set_nterm),]
    ggplot(data=gomf, aes(x=`P.value2`,y=reorder(`Term`,`P.value2`)))+
      geom_bar(stat = "identity",fill="#3c8dbc")+
      labs(title="GO_MF",x=expression(paste(-log[10],P.value)),y="")+
      theme_bw()+
      theme(axis.text=element_text(size=13),title = element_text(size=15,face="bold"), aspect.ratio = 0.5)
  })
  
  result_gsa_kegg <- reactive({
    res_gsa <- result_gsa()
    if(!is.null(res_gsa)){
      kegg<-res_gsa[["KEGG_2021_Human"]]
    }
    return(kegg)
  })
  
  reverted_gsa_kegg <- reactive({
    kegg <- result_gsa_kegg()
    kegg <- gsa_changePathwayID(kegg)
    return(kegg)
  })
  
  
  plot_gsa_kegg <- reactive({
    kegg <- result_gsa_kegg()
    kegg$P.value2 <- -log10(kegg$P.value)
    kegg <- kegg[order(-kegg$P.value2),]
    kegg <- kegg[c(1:input$set_nterm),]
    ggplot(data=kegg, aes(x=`P.value2`,y=reorder(`Term`,`P.value2`)))+
      geom_bar(stat = "identity",fill="#3c8dbc")+
      labs(title="KEGG",x=expression(paste(-log[10],P.value)),y="")+
      theme_bw()+
      theme(axis.text=element_text(size=13),title = element_text(size=15,face="bold"), aspect.ratio = 0.5)
  })
  
  selected_pathway <- reactive({
    if(!is.null(input$topOfKeggDT_rows_selected)){
      kegg_info <- top_of_kegg()
      selected_pathway <- kegg_info[input$topOfKeggDT_rows_selected, 1]
    } else if(input$pathID_selector !="") {
      selected_pathway <- input$pathID_selector  
    } 
    return(selected_pathway)
  })
  
  top_of_kegg <- reactive({
    kegg_info <- reverted_gsa_kegg()
    top_kegg <- kegg_info[order(kegg_info$Adjusted.P.value, decreasing=F),]
    top_kegg <- head(top_kegg, 10)
    top_kegg <- dplyr::select(top_kegg, -matches("Old|Odds|Score"))
    top_kegg <- top_kegg[,-ncol(top_kegg)]
    return(top_kegg)
  })

  pathway_graph <- eventReactive(list(input$render_pathway_btn, input$topOfKeggDT_rows_selected),{
    setwd(base_dir)
    
    pathway_name <- selected_pathway()
    kegg_info <- reverted_gsa_kegg()

    rowdt <- rowData(dep())
    rowdt <- rowdt[rowdt$significant == T,]

    fc_cols <- colnames(rowdt)[grepl("diff",colnames(rowdt))]
    fc <- rowdt[,fc_cols]
    names(fc) <- rowdt$name
    
    pathid <- as.character(kegg_info[kegg_info$Term==pathway_name, "kegg_id"])

    dir <- getwd()
    dir <- strsplit(dir,"/",fixed = T)
    dir <- as.character(unlist(dir))
    dir <- dir[length(dir)]
    if(dir != "GSA"){
      setwd("GSA")
    }
  
    pathview_dir <- "./xml"
    outfile <- paste0("hsa", pathid,".gsa_pathview_",user_id,".png")
    pathview(fc, pathway.id=pathid, gene.idtype="SYMBOL", species = "hsa",
             kegg.dir=pathview_dir,out.suffix=paste0("gsa_pathview_",user_id))
    
    return(outfile)
  })

  ########################### GSEA ######################################
  
  output$download_gsea_gobp <- downloadHandler(
    filename = function() {paste0("GSEA_GOBP_result_", Sys.Date(), ".csv")},
    content = function(file) {
      if(!is.null(result_gsea_gobp())){
        write.csv(result_gsea_gobp(),file,row.names=F,quote=F)
      }
    }
  )
  
  output$download_gsea_gocc <- downloadHandler(
    filename = function() {paste0("GSEA_GOCC_result_", Sys.Date(), ".csv")},
    content = function(file) {
      if(!is.null(result_gsea_gocc())){ 
        write.csv(result_gsea_gocc(),file,row.names=F,quote=F)
      }
    }
  )
  
  output$download_gsea_gomf <- downloadHandler(
    filename = function() {paste0("GSEA_GOMF_result_", Sys.Date(), ".csv")},
    content = function(file) {
      if(!is.null(result_gsea_gomf())){
        write.csv(result_gsea_gomf(),file,row.names=F,quote=F)
      }
    }
  )
  
  output$download_gsea_kegg <- downloadHandler(
    filename = function() {paste0("GSEA_Kegg_result_", Sys.Date(), ".csv")},
    content = function(file) {
      if(!is.null(result_gsea_kegg())){
        write.csv(result_gsea_kegg(),file,row.names=F,quote=F)
      }
    }
  )

  result_gsea <- reactive({
    req(res_test())
    rowdata <- rowData(res_test())
    res_gsea <- gsea(rowdata, input$select_genelevel_stats)
    return(res_gsea)
  })
  
  result_gsea_gobp <- reactive({
    res_gsea <- result_gsea()
    res_gobp <- res_gsea$result_GO_BP
    temp <- mapply(function(x){unlist(x)},res_gobp$leadingEdge)
    temp <- mapply(function(x){paste(x,collapse = ";")},temp)
    res_gobp$leadingEdge <- temp
    remove(temp)
    return(res_gobp)
  })
  
  plot_gsea_gobp <- reactive({
    res_gobp <- result_gsea_gobp()
    if(nrow(res_gobp) >= 20){
      filt_gobp <-  rbind(head(res_gobp, n = 10),tail(res_gobp, n = 10 ))  
    } else{
      filt_gobp <-  res_gobp
    }
    
    ggplot(filt_gobp, aes(reorder(pathway, NES), NES)) +
      geom_segment(aes(reorder(pathway, NES), xend=pathway, y=0, yend=NES)) +
      geom_point(size=5, aes( fill = Enrichment), shape=21, stroke=2) +
      scale_fill_manual(values = c("Down-regulated" = "dodgerblue", "Up-regulated" = "firebrick") ) +
      coord_flip() +
      labs(x="Pathway", y="Normalized Enrichment Score", title=paste0("GSEA - GOBP")) +
      theme(axis.text=element_text(size=13),title = element_text(size=15,face="bold"))
  })
  
  result_gsea_gocc <- reactive({
    res_gsea <- result_gsea()
    res_gocc <- res_gsea$result_GO_CC
    temp <- mapply(function(x){unlist(x)},res_gocc$leadingEdge)
    temp <- mapply(function(x){paste(x,collapse = ";")},temp)
    res_gocc$leadingEdge <- temp
    remove(temp)
    return(res_gocc)
  })
  
  plot_gsea_gocc <- reactive({
    res_gocc <- result_gsea_gocc()
    if(nrow(res_gocc) >= 20){
      filt_gocc <-  rbind(head(res_gocc, n = 10),tail(res_gocc, n = 10 ))
    }else{
      filt_gocc <-  res_gocc
    }
    
    ggplot(filt_gocc, aes(reorder(pathway, NES), NES)) +
      geom_segment(aes(reorder(pathway, NES), xend=pathway, y=0, yend=NES)) +
      geom_point(size=5, aes( fill = Enrichment), shape=21, stroke=2) +
      scale_fill_manual(values = c("Down-regulated" = "dodgerblue", "Up-regulated" = "firebrick") ) +
      coord_flip() + 
      labs(x="Pathway", y="Normalized Enrichment Score", title=paste0("GSEA - GOCC")) +
      theme(axis.text=element_text(size=13),title = element_text(size=15,face="bold"))
  })
  
  result_gsea_gomf <- reactive({
    res_gsea <- result_gsea()
    res_gomf <- res_gsea$result_GO_MF
    temp <- mapply(function(x){unlist(x)},res_gomf$leadingEdge)
    temp <- mapply(function(x){paste(x,collapse = ";")},temp)
    res_gomf$leadingEdge <- temp
    remove(temp)
    return(res_gomf)
  })
  
  plot_gsea_gomf <- reactive({
    res_gomf <- result_gsea_gomf()
    if(nrow(res_gomf) >= 20){
      filt_gomf <-  rbind(head(res_gomf, n = 10),tail(res_gomf, n = 10 ))
    }else{
      filt_gomf <- res_gomf
    }
   
    ggplot(filt_gomf, aes(reorder(pathway, NES), NES)) +
      geom_segment(aes(reorder(pathway, NES), xend=pathway, y=0, yend=NES)) +
      geom_point(size=5, aes( fill = Enrichment), shape=21, stroke=2) +
      scale_fill_manual(values = c("Down-regulated" = "dodgerblue", "Up-regulated" = "firebrick") ) +
      coord_flip() +
      labs(x="Pathway", y="Normalized Enrichment Score", title=paste0("GSEA - GOMF")) +
      theme(axis.text=element_text(size=13),title = element_text(size=15,face="bold"))
  })
  
  result_gsea_kegg <- reactive({
    res_gsea <- result_gsea()
    res_kegg<- res_gsea$result_Kegg
    temp <- mapply(function(x){unlist(x)},res_kegg$leadingEdge)
    temp <- mapply(function(x){paste(x,collapse = ";")},temp)
    res_kegg$leadingEdge <- temp
    remove(temp)
    return(res_kegg)
  })
  
  plot_gsea_kegg <- reactive({
    res_kegg <- result_gsea_kegg()
    if(nrow(res_kegg) >= 20){
      filt_kegg <-  rbind(head(res_kegg, n = 10),tail(res_kegg, n = 10 ))
    }else{
      filt_kegg <-  res_kegg
    }
    
    ggplot(filt_kegg, aes(reorder(pathway, NES), NES)) +
      geom_segment(aes(reorder(pathway, NES), xend=pathway, y=0, yend=NES)) +
      geom_point(size=5, aes( fill = Enrichment), shape=21, stroke=2) +
      scale_fill_manual(values = c("Down-regulated" = "dodgerblue", "Up-regulated" = "firebrick") ) +
      coord_flip() + 
      labs(x="Pathway", y="Normalized Enrichment Score", title=paste0("GSEA - KEGG")) +
      theme(axis.text=element_text(size=13),title = element_text(size=15,face="bold"))
  })
  
  reverted_gsea_kegg <- reactive({
    kegg <- result_gsea_kegg()
    kegg <- gsea_changePathwayID(kegg)
    return(kegg)
  })
  
  selected_gsea_up_pathway <- reactive({
    if(!is.null(input$result_of_gsea_up_regulated_rows_selected)){
      res <- reverted_gsea_kegg()
      up_res <- res[res$Enrichment=="Up-regulated",]
      up_res <- up_res[,-ncol(up_res)]
      up_res$kegg_id <- as.character(up_res$kegg_id)
      selected_up_pathway <- up_res[input$result_of_gsea_up_regulated_rows_selected, 1]
    } 
    return(selected_up_pathway)
  })
  
  selected_gsea_down_pathway <- reactive({
    if(!is.null(input$result_of_gsea_down_regulated_rows_selected)){
      res <- reverted_gsea_kegg()
      down_res <- res[res$Enrichment=="Down-regulated",]
      down_res <- down_res[,-ncol(down_res)]
      down_res$kegg_id <- as.character(down_res$kegg_id)
      selected_down_pathway <- down_res[input$result_of_gsea_down_regulated_rows_selected, 1]
    }
    
    return(selected_down_pathway)
  })
  
  gsea_up_pathway_graph <- eventReactive(input$result_of_gsea_up_regulated_rows_selected,{
    setwd(base_dir)
    
    pathid <- selected_gsea_up_pathway()
    gsea_rowdt <- rowData(dep())

    gsea_fc_cols <- colnames(gsea_rowdt)[grepl("diff",colnames(gsea_rowdt))]
    gsea_fc <- gsea_rowdt[,gsea_fc_cols]
    names(gsea_fc) <- gsea_rowdt$name
    
    dir <- getwd()
    dir <- strsplit(dir,"/",fixed = T)
    dir <- as.character(unlist(dir))
    dir <- dir[length(dir)]
    if(dir != "GSEA"){
      setwd("GSEA")
    }
    
    pathview_dir <- "./xml"
    gsea_outfile <- paste0("hsa", pathid,".gsea_pathview_",user_id,".png")
    pathview(gsea_fc, pathway.id=pathid, gene.idtype="SYMBOL", species = "hsa",
             kegg.dir=pathview_dir,out.suffix=paste0("gsea_pathview_",user_id))
    
    return(gsea_outfile)
  })
  
  gsea_down_pathway_graph <- eventReactive(input$result_of_gsea_down_regulated_rows_selected ,{
    setwd(base_dir)
    
    pathid <- selected_gsea_down_pathway()
    gsea_rowdt <- rowData(dep())

    gsea_fc_cols <- colnames(gsea_rowdt)[grepl("diff",colnames(gsea_rowdt))]
    gsea_fc <- gsea_rowdt[, gsea_fc_cols]
    names(gsea_fc) <- gsea_rowdt$name
    
    dir <- getwd()
    dir <- strsplit(dir,"/",fixed = T)
    dir <- as.character(unlist(dir))
    dir <- dir[length(dir)]
    if(dir != "GSEA"){
      setwd("GSEA")
    }
    
    pathview_dir <- "./xml"
    gsea_outfile <- paste0("hsa", pathid,".gsea_pathview_",user_id,".png")
    pathview(gsea_fc, pathway.id=pathid, gene.idtype="SYMBOL", species = "hsa",
             kegg.dir=pathview_dir,out.suffix=paste0("gsea_pathview_",user_id))
    
    return(gsea_outfile)
  })
  
  ########################### PPI network ######################################
  get_input_gene_toppi <- reactive({NULL})
  get_input_gene_toppi <- eventReactive(input$ppi_btn,{
    data <- rowData(dep())
    data <- data[data$significant==T,]
    gene <- as.character(data$name)
    fc <- abs(data[,grep("diff",colnames(data))])
    pval <- data[,grep("p.val",colnames(data))]
    qval <- data[,grep("p.adj",colnames(data))]
    dep <- data.frame(gene,fc,pval,qval)
    
    ppi_condition <- input$select_ppi_condition
    if(ppi_condition == "padj"){
      dep <- dep[order(dep$qval),]  
      gene <- as.character(dep$gene[c(1:input$input_num_toppi)])
    } else if (ppi_condition == "pval"){
      dep <- dep[order(dep$pval),]
      gene <- as.character(dep$gene[c(1:input$input_num_toppi)])
    } else if (ppi_condition == "log2fc") {
      dep <- dep[order(-dep$fc),]
      gene <- as.character(dep$gene[c(1:input$input_num_toppi)])
    } else{
      gene <- input$input_gene_toppi
      print(gene)
      gene <- gsub(" ","",gene, fixed=T)
      gene <- as.character(do.call('rbind', strsplit(as.character(gene), split = ',', fixed = TRUE)))
      print(gene)
    }
    
    return(gene)
  })
  
  string_url <- reactive({
    input_gene <- get_input_gene_toppi()
    URL <- string_url_builder(input$input_organism_toppi,input_gene)
    url <- gsub("http://","https://",URL)
    url_check <- readLines(url, encoding="UTF-8", skipNul = T)
    length_url <- length(url_check)
    slash_check <- grep("/",url_check,fixed = T)
    if(length(slash_check) != 0){
      length_url <- length_url - length(slash_check)
    }
    if(length_url < 7){
      shinyalert("There is no conection!", "Change PPI analysis option", type="error", timer = 10000,
                 closeOnClickOutside = T, closeOnEsc = T)
    }
   
    return(URL)
  })
  
  string_image_download_url <- reactive({
    input_gene <- get_input_gene_toppi()
    URL <- string_image_download_url_builder(input$input_organism_toppi,input_gene)
   
    return(URL)
  })
  
  string_tsv_download_url <- reactive({
    input_gene <- get_input_gene_toppi()
    URL <- string_tsv_download_url_builder(input$input_organism_toppi,input_gene)

    return(URL)
  })
  
  
  #############################################################################################################
  
  
  # Define the long Spark job here
  run_spark <- function(x) {
    # Environment setting
    library("SparkR", lib.loc = "/databricks/spark/R/lib")
    sparkR.session()
    
    irisDF <- createDataFrame(iris)
    collect(irisDF)
    Sys.sleep(3)
    x + 1
  }
  
  run_spark_sparklyr <- function(x) {
    # Environment setting
    library(sparklyr)
    library(dplyr)
    library("SparkR", lib.loc = "/databricks/spark/R/lib")
    sparkR.session()
    sc <- spark_connect(method = "databricks")
    
    iris_tbl <- copy_to(sc, iris, overwrite = TRUE)
    collect(iris_tbl)
    x + 1
  }
  
  
  addTimeLine <-  function(timeLine){
    output$timeline <- renderUI({
      timelineBlock(
        reversed = F,
        timelineEnd(icon=icon("hourglass-start"), color = "gray"),
        lapply(1:nrow(timeLine), FUN = function(i){
          tagList(
            timelineItem(
              icon = icon(timeLine$icon[i]),
              color = timeLine$color[i],
              time = timeLine$time[i],
              h4(timeLine$step[i]),
              footer= timeLine$info[i]
            )
            ,timelineLabel(paste0("# of proteins : ",timeLine$sample_num[i]), color = "purple")
          )
        }),
        timelineStart(icon=icon("hourglass-end"), color = "gray")
      )
    })
  } # End of addTimeLine

  
  session$onSessionEnded(function() {
    setwd(base_dir)
    file_list <- list.files(path="./GSA/", pattern=".png")
    file_list <- paste0("./GSA/",file_list)
    file.remove(file_list)
    file_list <- list.files(path="./GSEA/", pattern=".png")
    file_list <- paste0("./GSEA/",file_list)
    file.remove(file_list)
  })
})