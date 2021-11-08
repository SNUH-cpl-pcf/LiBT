make_id <- function(){
  user_sys_time <- as.character(Sys.time())
  user_sys_time <- gsub(" ","",user_sys_time,fixed=T)
  user_sys_time <- gsub("-","",user_sys_time,fixed=T)
  user_sys_time <- gsub(":","",user_sys_time,fixed=T)
  
  ran_char <- sample(letters,4)
  ran_num <- sample(0:9,4)
  user_id <- c()
  for(i in 1:4){
    user_id <- paste0(user_id,ran_char[i])
  }
  for(i in 1:4){
    user_id <- paste0(user_id,ran_num[i])
  }
  
  user_id <- paste0(user_id,"_",user_sys_time)
  print(paste0("user_id:",user_id))
  return(user_id)
}


file_input_test <- function(file, file_type) {
  switch(file_type,
         "TMT" = 
           {cols <- colnames(file)[grep("Abundance", colnames(file))]
           if(length(cols)==0){
             cols <- colnames(file)[grep("Reporter", colnames(file))]
           }},
         "iBAQ" = 
           {cols <- colnames(file)[grep("iBAQ", colnames(file))]},
         "LFQ" = 
           {cols <- colnames(file)[grep("LFQ", colnames(file))]})
  
  # if(file_type=="TMT"){
  #   cols <- colnames(file)[grep("Abundance", colnames(file))]
  # } else {
  #   cols <- colnames(file)[grep(file_type, colnames(file))]
  # }
  
  return(cols)
}


filter_with_option <- function(option, data) {
  data <- dplyr::filter(data, Peptides != 0)#4903
  data <- dplyr::filter(data, Intensity != 0)#4898
  
  option <- paste0(option, collapse=",")
  switch(option,
         "potential" = 
           {data <- dplyr::filter(data, Potential.contaminant != "+")},
         "reverse" = 
           {data <- dplyr::filter(data, Reverse != "+")},
         "identified" = 
           {data <- dplyr::filter(data, Only.identified.by.site != "+")},
         "potential,reverse" = 
           {data <- dplyr::filter(data, Potential.contaminant != "+" &
                                    Reverse != "+")},
         "potential,identified" = 
           {data <- dplyr::filter(data, Potential.contaminant != "+" &
                                    Only.identified.by.site != "+")},
         "reverse,identified" = 
           {data <- dplyr::filter(data, Reverse != "+" &
                                    Only.identified.by.site != "+")},
         "potential,reverse,identified" = 
           {data <- dplyr::filter(data, Potential.contaminant != "+" &
                                    Only.identified.by.site != "+" & Reverse != "+")},
         {data <- data}
  )
  return(data)
}



get_main_data_MQ <- function(data, file_type) {
  rv_db <- read.delim("base/uniprot-filtered-organism_Human(9606)_rv_20200114ver.txt")
  nrv_db <- read.delim("base/uniprot-filtered-organism_Human(9606)_nrv_20200114ver.txt")
  # data <- read.csv("./TMT_MQ.csv")
  
  ProteinID <- as.character(data$Majority.protein.IDs)
  Protein_count <- as.character(data$Peptide.counts..all.)
  GeneName <- as.character(data$Gene.names)

  multi_pos <- grep(";",Protein_count,fixed = T)
  ProteinID_multi <- ProteinID[multi_pos] 
  Protein_count_multi <- Protein_count[multi_pos]
  
  pt_multi <- c()
  i=1
  for(i in 1:length(ProteinID_multi)){
    pt_id <- ProteinID_multi[i]
    pt_id <- data.frame(do.call('rbind', strsplit(as.character(pt_id), split = ';', fixed = TRUE)))
    pt_id <- as.character(t(pt_id))
    
    pt_count <- Protein_count_multi[i]
    pt_count <- data.frame(do.call('rbind', strsplit(as.character(pt_count), split = ';', fixed = TRUE)))
    pt_count <- as.numeric(t(pt_count))
    pt_count <- pt_count[1:length(pt_id)]
    
    pt <- data.frame(Entry=pt_id,count=pt_count)
    id_pos <- as.numeric(grep("-",pt$Entry,fixed = T))
    if(length(id_pos)!=0){
      edit_id <- as.character(pt$Entry[id_pos])
      edit_id <- data.frame(do.call('rbind', strsplit(edit_id, split = '-', fixed = TRUE)))
      pt_id[id_pos] <- as.character(edit_id[,1])
      pt$Entry <- pt_id
    }
    
    pt_db <- merge(pt,rv_db,by="Entry")
    if(nrow(pt_db)!=0){
      pt_db <- pt_db[order(pt_db$count,decreasing = T),]
      pt_db <- pt_db[pt_db$count==pt_db$count[1],]
      pt_db <- pt_db[!duplicated(pt_db$Entry),]
    }else{
      pt_nrv <- merge(pt,nrv_db,by="Entry")
      if(nrow(pt_nrv) != 0){
        pt_nrv <- pt_nrv[pt_nrv$count==pt_nrv$count[1],]
        pt_db <- pt_nrv
        pt_db <- pt_db[!duplicated(pt_db$Entry),]
      }else{
        data_pos <- which(data$Majority.protein.IDs==ProteinID_multi[i])
        pt_db <- data.frame(Entry=data$Majority.protein.IDs[data_pos],Gene.symbol=data$Gene.names[data_pos])
      }
    }
    
    id <- c()
    gn <- c()
    if(nrow(pt_db)>1){
      for(j in 1:nrow(pt_db)){
        id <- paste0(id,";",pt_db$Entry[j])
        gn <- paste0(gn,";",pt_db$Gene.symbol[j])
      }
      pt_id <- gsub("^;","",id)
      pt_gn <- gsub("^;","",gn)
      gn_multi <- data.frame(ProteinID=pt_id,GeneName=pt_gn)
    } else{
      pt_gn <- pt_db$Gene.symbol
      gn_multi <- data.frame(ProteinID=pt_db$Entry,GeneName=pt_gn)
    }
    
    pt_multi <- rbind(pt_multi,gn_multi)
  }
  
  ProteinID[multi_pos] <- as.character(pt_multi$ProteinID)
  GeneName[multi_pos] <- as.character(pt_multi$GeneName)
  
  ProteinID_solo <- ProteinID[-multi_pos]
  for(i in 1:length(ProteinID_solo)){
    solo_pos <- which(data$Majority.protein.IDs==ProteinID_solo[i])
    pt <- data.frame(Entry=ProteinID_solo[i])
    pt_db <- merge(pt,rv_db,by="Entry")
    if(nrow(pt_db)!=0){
      pt_db <- pt_db[!duplicated(pt_db$Entry),]
    }else{
      pt_nrv <- merge(pt,nrv_db,by="Entry")
      if(nrow(pt_nrv) != 0){
        pt_db <- pt_nrv
        pt_db <- pt_db[!duplicated(pt_db$Entry),]
      }else{
        pt_db <- data.frame(Entry=data$Majority.protein.IDs[solo_pos],Gene.symbol=data$Gene.names[solo_pos])
      }
    }
    
    id <- c()
    gn <- c()
    if(nrow(pt_db)>1){
      for(j in 1:nrow(pt_db)){
        id <- paste0(id,";",pt_db$Entry[j])
        gn <- paste0(gn,";",pt_db$Gene.symbol[j])
      }
      pt_id <- gsub("^;","",id)
      pt_gn <- gsub("^;","",gn)
      pt_solo <- data.frame(ProteinID=pt_id,GeneName=pt_gn)
    } else{
      pt_gn <- pt_db$Gene.symbol
      pt_solo <- data.frame(ProteinID=pt_db$Entry,GeneName=pt_gn)
    }
    
    ProteinID[solo_pos] <- as.character(pt_solo$ProteinID)
    GeneName[solo_pos] <- as.character(pt_solo$GeneName)
  }
  
  cn_data <- colnames(data)
  
  if(file_type =="LFQ"){
    main_pos <- grep("LFQ",cn_data)
  } else if(file_type =="iBAQ") {
    main_pos <- grep("iBAQ.",cn_data)
  } else{
    total_pos <- grep("Reporter", cn_data)
    real_pos <- grep("corrected",cn_data)
    real_pos <- intersect(total_pos, real_pos)
    except_pos <- grep("FL_SET",cn_data)
    main_pos <- setdiff(real_pos,except_pos)
  }
  
  # total <- data[,total_pos]
  # real <- data[,real_pos]
  # except <- data[,except_pos]
  
  main_data <- data[,main_pos]
  
  unique_pos <- grep("Unique.peptides",cn_data)
  Unique.peptides <- data[,unique_pos]
  Sequence.coverage <- data$Sequence.coverage....
  
  main_df <- cbind(ProteinID,GeneName,main_data,Unique.peptides,Sequence.coverage)
  
  isDup_Gene <- main_df$GeneName %>% duplicated() %>% any()
  if(isDup_Gene==T){
    data_unique <- make_unique(main_df, "GeneName", "ProteinID")
  } else{
    data_col <- colnames(data)
    data_pos <- which("ProteinID"==data_col | "GeneName"==data_col)
    colnames(data)[data_pos] <- c("ID","name")
    data_unique <- data
  }
  
  isDup_name <- main_df$name %>% duplicated() %>% any()
  if(isDup_name == F){
    for(i in 1:nrow(data_unique)){
      if(is.na(data_unique$GeneName[i])){
        data_unique$name[i] <- paste0(data_unique$ID[i],"|?")
      }
    }
    data_unique <- data.frame(ID=data_unique$ID,name=data_unique$name,data_unique[,-c(1:2,(ncol(data_unique)-1):ncol(data_unique))])
    return(data_unique)
  }
}



get_main_data_PD <- function(data,normalization) {
  rv_db <- read.delim("base/uniprot-filtered-organism_Human(9606)_rv_20200114ver.txt")
  nrv_db <- read.delim("base/uniprot-filtered-organism_Human(9606)_nrv_20200114ver.txt")
  
  cn <- colnames(data)
  
  data$Accession <- as.character(data$Accession)
  data_class <- sapply(data,class)
  
  rm_pos <- grep("factor",data_class) 
  rm_pos <- c(rm_pos,grep("Ratio",cn))
  rm_pos <- c(rm_pos,grep("Grouped",cn))
  rm_pos <- c(rm_pos,grep("Scaled",cn))
  rm_pos <- c(rm_pos,grep("Count",cn))

  data <- data[,-rm_pos]
  cn <- cn[-rm_pos]
  main_pos <- grep("Accession",cn)
  
  cn_tmp <- grep("127C",cn)
  if(length(cn_tmp) != 0){
    TMT_marker <- c("126","127C","127N","128C","128N","129C","129N","130C","130N","131")
  } else{
    TMT_marker <- c("126","127","128","129","130","131")
  }
  
  for(i in 1:length(TMT_marker)){
    main_pos <- c(main_pos,grep(TMT_marker[i],cn))  
  }
  main_pos <- unique(main_pos)
  
  main_data <- data[,main_pos]
  cn <- colnames(main_data)
  
  if(normalization == "T"){
    nor_pos <- grep("Normalized",cn)
    main_pos <- c(1,nor_pos)
  } else{
    origin_pos <- grep("Abundance", cn)
    nor_pos <- grep("Normalized",cn)
    nor_pos <- setdiff(origin_pos, nor_pos)
    main_pos <- c(1,nor_pos)
  }
  
  main_data <- main_data[,main_pos]
  #main_data <- na.omit(main_data)
  
  acc <- main_data$Accession
  id_pos <- grep("-",acc)
  acc_id <- acc[id_pos]
  
  acc_new <- c()
  for(j in 1:length(acc_id)){
    tmp <- acc_id[j]
    tmp <- unlist(strsplit(as.character(tmp), split = '-', fixed = TRUE))
    acc_new <- c(acc_new,tmp[1])
  }
  acc[id_pos] <- acc_new
  main_data$Entry <- acc
  
  pt_data <- merge(main_data,rv_db,by="Entry",all.x = T)
  pt_rv <- pt_data[!is.na(pt_data$Entry.name),]
  pt_rv <- pt_rv[,-c(((ncol(pt_rv)-6):(ncol(pt_rv)-2)),ncol(pt_rv))]
  pt_rv <- pt_rv[,c(2,ncol(pt_rv),3:(ncol(pt_rv)-1))]
  
  pt_nrv <- pt_data[is.na(pt_data$Entry.name),]
  pt_nrv <- pt_nrv[,c(1:(ncol(pt_nrv)-7))]
  pt_nrv <- merge(pt_nrv,nrv_db,by="Entry",all.x = T)
  
  pt_none <- pt_nrv[is.na(pt_nrv$Entry.name),]
  pt_none <- pt_none[,-c(((ncol(pt_none)-6):(ncol(pt_none)-2)),ncol(pt_none))]
  pt_none <- pt_none[,c(2,ncol(pt_none),3:(ncol(pt_none)-1))]
  
  pt_nrv <- pt_nrv[!is.na(pt_nrv$Entry.name),]
  pt_nrv <- pt_nrv[,-c(((ncol(pt_nrv)-6):(ncol(pt_nrv)-2)),ncol(pt_nrv))]
  pt_nrv <- pt_nrv[,c(2,ncol(pt_nrv),3:(ncol(pt_nrv)-1))]

  main_df <- rbind(pt_rv,pt_nrv,pt_none)
  
  isDup_Gene <- main_df$Gene.symbol %>% duplicated() %>% any()
  if(isDup_Gene==T){
    data_unique <- make_unique(main_df, "Gene.symbol", "Accession")
  } else{
    data_col <- colnames(data)
    data_pos <- which("Accession"==data_col | "Gene.symbol"==data_col)
    colnames(data)[data_pos] <- c("ID","name")
    data_unique <- data
  }
  
  #isDup_name <- data$name %>% duplicated() %>% any()
  isDup_name <- main_df$name %>% duplicated() %>% any()
  if(isDup_name==F){
    for(i in 1:nrow(data_unique)){
      if(is.na(data_unique$Gene.symbol[i])){
        data_unique$name[i] <- paste0(data_unique$ID[i],"|?")
      }
    }
    data_unique <- data.frame(ID=data_unique$ID,name=data_unique$name,data_unique[,-c(1:2,(ncol(data_unique)-1):ncol(data_unique))])
    return(data_unique)
  }
}

make_case_samples_MQ <- function(data, file_type) {
  temp_col <- colnames(data)
  if(file_type == "TMT"){
    temp_col <- temp_col[grep("Reporter", temp_col)]
  }else{
    temp_col <- temp_col[grep(file_type, temp_col)]
  }
  
  return(temp_col)
}

make_case_samples_PD <- function(data) {
  # write.csv(data,"./main_data.csv")
  # data <- read.csv("./main_data.csv")
  data <- data[,c(3:ncol(data))]
  temp_col <- colnames(data)
  
  cn_tmp <- grep("127C",temp_col)
  if(length(cn_tmp) != 0){
    TMT_marker <- c("126","127C","127N","128C","128N","129C","129N","130C","130N","131")
  } else{
    TMT_marker <- c("126","127","128","129","130","131")
  }
  
  temp <- c()
  for(i in 1:length(TMT_marker)){
    col_pos <- str_locate(temp_col[grep(TMT_marker[i],temp_col)],TMT_marker[i]) 
    temp <- c(temp,substr(temp_col[i],col_pos[1,1],(nchar(temp_col[i])+1)))
  }

  col <- data.frame(id=temp)
  
  col1 <- col[1:(nrow(col)/2),]
  col2 <- col[((nrow(col)/2)+1):nrow(col),]
  
  col <- union(col1, col2)
  
  return(col)
}

make_case_samples_diffPD <- function(data,normalization) {
  cn <- colnames(data)
  data$ID <- as.character(data$ID)
  data_class <- sapply(data,class)
  rm_pos <- grep("factor",data_class) 
  rm_pos <- c(rm_pos,grep("Ratio",cn))
  
  data <- data[,-rm_pos]
  cn <- cn[-rm_pos]
  
  main_pos <- grep("ID",cn)
  
  TMT_marker <- c("126","127C","127N","128C","128N","129C","129N","130C","130N","131")
  for(i in 1:length(TMT_marker)){
    main_pos <- c(main_pos,grep(TMT_marker[i],cn))  
  }
  main_pos <- unique(main_pos)
  
  main_data <- data[,main_pos]
  cn <- colnames(main_data)
  
  nor_pos <- grep("Normalized",cn)
  main_pos <- c(1,nor_pos)
  main_data <- main_data[,main_pos]
  
  main_data <- na.omit(main_data)
  main_data <- main_data[,c(2:ncol(main_data))]
  temp_col <- colnames(main_data)
  
  temp <- c()
  for(j in 1:length(TMT_marker)){
    col_pos <- str_locate(temp_col[grep(TMT_marker[j],temp_col)],TMT_marker[j]) 
    temp <- c(temp,substr(temp_col[j],col_pos[1,1],(nchar(temp_col)+1)))
  }
  
  col <- data.frame(ID=temp)
  
  col1 <- col[1:(nrow(col)/2),]
  col2 <- col[((nrow(col)/2)+1):nrow(col),]
  
  col <- union(col1, col2)
  
  return(col)
}

make_expDesignData <- function(case,ctrl,discover_tool) {
  # discover_tool <- "PD"
  label <- c(case,ctrl)
  condition <- make_condition(case,ctrl,discover_tool)

  # replicate <- c()
  # replicate[c(1:length(label))] <- 1
  
  replicate_case <- c()
  for(i in 1:length(case)){
    replicate_case[i] <- i
  }
  
  replicate_ctrl <- c()
  for(j in 1:length(ctrl)){
    replicate_ctrl[j] <- j
  }
  
  replicate <- c(replicate_case,replicate_ctrl)
  
  
  exp_design <- data.frame(label=label, condition=condition, replicate=replicate)

  return(exp_design)
}

make_condition <- function(case,ctrl,discover_tool){
  # case <- c( "126..Sample..A2..n.a",  "127C..Sample..A4..n.a", "127N..Sample..A3..n.a", "128N..Sample..A5..n.a", "129N..Sample..A6..n.a")
  # ctrl <- c("128C..Sample..n.a..P1", "129C..Sample..n.a..P2", "130C..Sample..n.a..P5", "130N..Sample..n.a..P4", "131..Sample..n.a..P6" )
  # case <- c("LFQ.intensity.Total_309B", "LFQ.intensity.Total_445B", "LFQ.intensity.Total_555B", "LFQ.intensity.Total_588B", "LFQ.intensity.Total_636B",
  #            "LFQ.intensity.Total_667B", "LFQ.intensity.Total_741B", "LFQ.intensity.Total_764B", "LFQ.intensity.Total_876B", "LFQ.intensity.Total_883B")
  # ctrl <- c("LFQ.intensity.Total_309M", "LFQ.intensity.Total_445M", "LFQ.intensity.Total_555M", "LFQ.intensity.Total_588M", "LFQ.intensity.Total_636M",
  #            "LFQ.intensity.Total_667M", "LFQ.intensity.Total_741M", "LFQ.intensity.Total_764M", "LFQ.intensity.Total_876M", "LFQ.intensity.Total_883M")
  # case <- c("Reporter.intensity.corrected.1", "Reporter.intensity.corrected.4", "Reporter.intensity.corrected.5", "Reporter.intensity.corrected.6", "Reporter.intensity.corrected.7", "Reporter.intensity.corrected.9")
  # ctrl <- c("Reporter.intensity.corrected.0", "Reporter.intensity.corrected.2", "Reporter.intensity.corrected.3", "Reporter.intensity.corrected.8")
  if(discover_tool != "PD"){
    label_case <- strsplit(case,".",fixed = T)
    label_ctrl <- strsplit(ctrl,".",fixed = T)
    case_name <- c()
    i=1
    for(i in 1:length(label_case)){
      list_case <- label_case[[i]]
      case_name <- c(case_name,list_case[length(list_case)])
    }
    
    ctrl_name <- c()
    for(i in 1:length(label_ctrl)){
      list_ctrl <- label_ctrl[[i]]
      ctrl_name <- c(ctrl_name,list_ctrl[length(list_ctrl)])
    }
    
    case_name_end <- substr(case_name,nchar(case_name),nchar(case_name))
    ctrl_name_end <- substr(ctrl_name,nchar(ctrl_name),nchar(ctrl_name))
    
    case_name_end <- gsub("\\d+","",case_name_end)
    ctrl_name_end <- gsub("\\d+","",ctrl_name_end)
    
    for(i in 1:length(case_name_end)){
      if(case_name_end[i] == ""){
        case_name[i] <- substr(case_name[i],1,(nchar(case_name[i])-1))
      }else{
        case_name[i] <- gsub("\\d+","",case_name[i])
      }
    }
    
    
    for(i in 1:length(ctrl_name_end)){
      if(ctrl_name_end[i] == ""){
        ctrl_name[i] <- substr(ctrl_name[i],1,(nchar(ctrl_name[i])-1))
      }else{
        ctrl_name[i] <- gsub("\\d+","",ctrl_name[i])
      }
    }
    
    case_name <- gsub(".","",case_name, fixed = T)
    ctrl_name <- gsub(".","",ctrl_name, fixed = T)
    
    if(nchar(case_name)[1]==0){
      case_name[1:length(case_name)] <- "Case"
    }
    
    if(nchar(ctrl_name)[1]==0){
      ctrl_name[1:length(ctrl_name)] <- "Control"
    }
    
    condition <- c(case_name,ctrl_name) 
  } else{
    ####PD####
    label_case <- strsplit(case,"..",fixed = T)
    label_ctrl <- strsplit(ctrl,"..",fixed = T)
    
    case_name <- c()
    for(i in 1:length(label_case)){
      list_case <- label_case[[i]]
      list_case <- list_case[-c(1:2)]
      list_case <- list_case[-grep("n.a",list_case)]
      case_name <- c(case_name,list_case)
    }
    
    ctrl_name <- c()
    for(i in 1:length(label_ctrl)){
      list_ctrl <- label_ctrl[[i]]
      list_ctrl <- list_ctrl[-c(1:2)]
      list_ctrl <- list_ctrl[-grep("n.a",list_ctrl)]
      ctrl_name <- c(ctrl_name,list_ctrl)
    }
    
    case_name_end <- substr(case_name,nchar(case_name),nchar(case_name))
    ctrl_name_end <- substr(ctrl_name,nchar(ctrl_name),nchar(ctrl_name))
  
    case_name_end <- gsub("\\d+","",case_name_end)
    ctrl_name_end <- gsub("\\d+","",ctrl_name_end)
    
    i=1
    for(i in 1:length(case_name_end)){
      if(case_name_end[i] == ""){
        case_name[i] <- substr(case_name[i],1,(nchar(case_name[i])-1))
      }else{
        case_name[i] <- gsub("\\d+","",case_name[i])
      }
    }
    
    
    for(i in 1:length(ctrl_name_end)){
      if(ctrl_name_end[i] == ""){
        ctrl_name[i] <- substr(ctrl_name[i],1,(nchar(ctrl_name[i])-1))
      }else{
        ctrl_name[i] <- gsub("\\d+","",ctrl_name[i])
      }
    }
    
    case_name <- gsub(".","",case_name, fixed = T)
    ctrl_name <- gsub(".","",ctrl_name, fixed = T)
    
    if(nchar(case_name)[1]==0){
      case_name[1:length(case_name)] <- "Case"
    }
    
    if(nchar(ctrl_name)[1]==0){
      ctrl_name[1:length(ctrl_name)] <- "Control"
    }
    
    condition <- c(case_name,ctrl_name)
  }
  
  return(condition)
}


make_summarizedData <- function(main,file_type,design){
  data_col <- c()
  
  for(i in 1:nrow(design)){
    tmp <- design$label[i]
    if(file_type == "TMT"){
      data_col <- c(data_col,grep(tmp,colnames(main)))
      if(length(data_col)==0){
        data_col <- c(data_col,grep(tmp,colnames(main)))
      }
    } else if(file_type == "LFQ"){
      data_col <- grep("LFQ",colnames(main))
    } else{
      data_col <- grep("iBAQ",colnames(main))
    }
  }
  
  design$label <- as.character(design$label)
  design$condition <- as.character(design$condition)
  design$replicate <- as.numeric(design$replicate)
  
  # print(data_col)
  # main <- main[,c(1,2,data_col)]
  # print(colnames(main))
  # data_col <- c(3:ncol(main))
  
  summary_data <- make_se(main,data_col,design)
  return(summary_data)
}

# use_transformation_option <- function(data, case, control, options) {
#   samples <- c(case, control)
# 
#   if(options=="none") {
#     for(i in 1:nrow(data)) {
#       data[i,c(samples)] <- apply(data[i,c(samples)],1,function(x){2^(x)})
#     }
#   } else {
#     data <- data
#   }
# 
#   return(data)
# }

use_valid_option <- function(data, case, control, option) {
  option <- as.numeric(as.character(option))
  valid_num <- length(case) * option
  valid_num <- round(valid_num, digits = 0)
  data_filt <- filter_missval(data, thr=valid_num)
  
  # res <- assay(data_filt)
  # res <- round(res,digits = 2)
  # data_filt <- SummarizedExperiment(assays = list(as.data.frame(res)), rowData = rowData(data_filt), colData = colData(data_filt))
  
  return(data_filt)
}

use_normalization_option <- function(data, option) {
  if(option=="Yes") {
    print(option)
    data_norm <- normalize_vsn(data)
  }
  # } 
  # else if(option == "Z-score"){
  #   val_z <- assay(data)
  #   res <- c()
  #   for(i in 1:nrow(val_z)){
  #     tmp <- as.numeric(val_z[i,])
  #     res <- rbind(res,(tmp-mean(tmp))/sd(tmp))
  #   }
  #   data_norm <- SummarizedExperiment(assays = list(as.data.frame(res)), rowData = rowData(data), colData = colData(data))
  # }
  else {
    data_norm <- data
  }
  return(data_norm)
}


use_imputation_option <- function(data, case, control, option) {
  switch(option,
         "QRILC" =
           {data_imp <- impute(data, fun="QRILC")},
         "MinProb" =
           {data_imp <- impute(data, fun="MinProb", q=0.01)},
         "man" =
           {data_imp <- impute(data, fun = "man", shift = 1.8, scale = 0.3)},
         "zero" =
           {data_imp <- impute(data, fun="zero")}
  )
  
  return(data_imp)
}

test <- function(data,design,p_option,q_option){
  total_sample_num <- nrow(design)
  condition <- unique(design$condition)
  case_num <- length(grep(condition[1],design$condition))
  ctrl_num <- length(grep(condition[2],design$condition))
  case_dt <- data[,c(1:case_num)]
  ctrl_dt <- data[,c((case_num+1):(case_num+ctrl_num))]
  pval <- c()
  if(p_option == "T-test"){
    case_df <- as.data.frame(t(case_dt))
    ctrl_df <-  as.data.frame(t(ctrl_dt))
    pval <- mapply(function(x,y){t.test(x,y,alternative = c("two.sided"))},case_df,ctrl_df)
    pval <- as.data.frame(t(pval))
    pval <- unlist(pval$p.value)
  }
  else if(p_option == "Wilcoxon-Rank-Sum"){
    case_df <- t(case_dt)
    case_df <- as.data.frame(case_df)
    ctrl_df <- t(ctrl_dt)
    ctrl_df <- as.data.frame(ctrl_df)
    pval <- mapply(function(x,y){wilcox.test(x,y,alternative = c("two.sided"))},case_df,ctrl_df)
    pval <- t(pval)
    pval <- as.data.frame(pval)
    pval <- unlist(pval$p.value)
  }
  else{
    #edgeR
    condition <- design[,c(2,3)]
    condition$condition <- as.factor(condition$condition)
    design<-model.matrix(~condition,data=condition)
    
    edr<-DGEList(counts=data, genes=rownames(data))
    noredr<-calcNormFactors(edr, method="TMM")
    
    disedr<-estimateDisp(noredr, design)
    fitq<-glmQLFit(disedr,design)
    glf<-glmQLFTest(fitq)
    res_edgeR<-glf$table
    pval <- res_edgeR$PValue
  }
  
  qval <- p.adjust(pval,method = q_option)
  df <- data.frame(pval=pval,padj=qval)
  return(df)
}

change_Sig <- function(data_rejection,pvalue,log2fc){
  dep_rowData <- rowData(data_rejection)
  pval_pos <- grep("p.val",colnames(dep_rowData),fixed = T)
  log2fc_pos <- grep("diff",colnames(dep_rowData),fixed = T)
  sig_pos <- grep("significant",colnames(dep_rowData),fixed = T)
  
  pval_df <- as.numeric(dep_rowData[,pval_pos])
  log2fc_df <- abs(as.numeric(dep_rowData[,log2fc_pos]))
  res_pval <- as.logical(pval_df < pvalue)
  res_lfc <-  as.logical(log2fc_df > log2fc)
  TT_pos <- which(res_pval=="TRUE" & res_lfc=="TRUE")
  
  res <- c()
  res[1:nrow(dep_rowData)] <- "FALSE"
  res[TT_pos] <- "TRUE"
  res <- as.logical(res)
  
  dep_rowData[,sig_pos] <- res
  dep_rowData$name <- as.character(dep_rowData$name)
  return(dep_rowData)
}

get_optimize_k <- function(assay, sig_gene){
  data <- as.data.frame(assay)
  gene <- rownames(data)
  data <- cbind(gene,data)
  
  sig_gene <- data.frame(gene=sig_gene)
  
  data <- merge(sig_gene,data,by="gene")
  rownames(data) <- data$gene
  data <- data[,-1]
  
  # write.csv(data,"./input_heatmap.csv", row.names = F)
  # data <- read.csv("./input_heatmap.csv")
  data <- t(data)
  if(ncol(data) > 2){
    max_col = nrow(data)-1
    fiz <- fviz_nbclust(data, kmeans, method="silhouette", diss=dist(data,method="euclidean"), k.max=max_col)
    y <- fiz[["data"]][["y"]]
    best_k <- as.numeric(which(y==max(y)))

  }else{
    best_k <- 1
  }
  return(best_k)
}

get_gene_cluster <- function(assay, rowData, k){
  assay <- as.data.frame(assay)
  # write.csv(assay,"./assay.csv",row.names = F, quote = F)
  rowData <- as.data.frame(rowData)
  # write.csv(rowData,"./rowData.csv",row.names = F, quote = F)
  # assay <- read.csv("./assay.csv")
  # rowData <- read.csv("./rowData.csv")
  # k <- 2
  
  data <- cbind(assay,rowData)
  data <- data[data$significant==T,]
  
  geneName <- data$name
  ID <- data$ID
  
  data <- data[,c(1:ncol(assay))]
  mean <- rowMeans(data, na.rm=T)
  df <- data-mean
  
  set.seed(1234)
  df_kmeans <- kmeans(df, k)
  order <- data.frame(df) %>%
    cbind(., cluster = df_kmeans$cluster) %>%
    mutate(row = apply(.[, seq_len(ncol(.) - 1)], 1, function(x) max(x))) %>%
    group_by(cluster) %>%
    summarize(index = sum(row)/n()) %>%
    arrange(desc(index)) %>%
    pull(cluster) %>%
    match(seq_len(k), .)
  
  result <- data.frame(geneName, ID, cluster = df_kmeans$cluster)
  result <-result[order(result$cluster),]
  return(result)
}

get_pca_df <- function(assay, colData, n){
  assasy <- as.data.frame(assay)
  colData <- as.data.frame(colData)
  var <- apply(assay, 1, sd)
  df <- assay[order(var, decreasing = TRUE)[seq_len(n)],]
  
  # Calculate PCA
  pca <- prcomp(t(df), scale = FALSE)
  pca_df <- pca$x %>%
    data.frame() %>%
    rownames_to_column() %>%
    left_join(., data.frame(colData), by = c("rowname" = "ID"))
  
  return(pca_df)
}

gsa <- function(data,input, set){
  data <- as.data.frame(data)
  if(input == "dep"){
    data <- data[data$significant == T,]
  }
  data <- data[,c(1,7)]
  colnames(data) <- c("genes","log2fc")
  if(set == "case_up"){
    genes <- data$genes[data$log2fc > 0]
  } else{
    genes <- data$genes[data$log2fc < 0]
  }
  
  res_gsa <- enrichR(genes)
  return(res_gsa)
}

enrichR <- function(genes){
  dbs <- c("GO_Molecular_Function_2021", "GO_Cellular_Component_2021", "GO_Biological_Process_2021", "KEGG_2021_Human") #, "WikiPathways_2016", "Reactome_2016"
  res_er<-enrichr(genes,dbs)
  return(res_er)
}

make_myGO <- function(cate){
  dir <- "./base/"
  gmts <- c("c2.cp.kegg.v7.4.symbols.gmt","c5.go.bp.v7.4.symbols.gmt", "c5.go.cc.v7.4.symbols.gmt", "c5.go.mf.v7.4.symbols.gmt")
 
  switch(cate,
    "gobp" = {myGO <- gmtPathways(paste0(dir,gmts[2]))},
    "gocc" = {myGO <- gmtPathways(paste0(dir,gmts[3]))},
    "gomf" = {myGO <- gmtPathways(paste0(dir,gmts[4]))},
    "kegg" = {myGO <- gmtPathways(paste0(dir,gmts[1]))}
  )
  return(myGO)
}

make_geneRank <- function(rowData,stats){
  rowData <- as.data.frame(rowData)
  
  fc <- rowData[,grep("diff",colnames(rowData))]
  padj <- rowData[,grep("p.adj",colnames(rowData))]
  pval <- rowData[,grep("p.val",colnames(rowData))]
  
  if(stats == "P.adj"){
    gene_list <- as.numeric(padj)
  }else if(stats == "P.value"){
    gene_list <- as.numeric(pval)
  }else if(stats == "log2fc"){
    gene_list <- as.numeric(fc)
  }else{
    gene_list <- as.numeric(fc*-log(padj,10))  
  }
  
  gene <- as.character(rowData$name)
  names(gene_list) <- gene
  gene_list <- sort(gene_list, decreasing = T)
  
  return(gene_list)
}


gsea <- function(rowData, stats){
  #rowData <- read.csv("C:/Users/hylee4161/Downloads/iBAQ_DEP_info_data_2021-06-30.csv")
  #stats <- "log2fc"
  set.seed(1234)
  
  
  gene_list <- make_geneRank(rowData,stats)
  
  dir <- "./base/"
  gmts <- c("c2.cp.kegg.v7.4.symbols.gmt","c5.go.bp.v7.4.symbols.gmt", "c5.go.cc.v7.4.symbols.gmt", "c5.go.mf.v7.4.symbols.gmt")
  cate <- c("Kegg","GO_BP", "GO_CC","GO_MF")
  i=2
  for(i in 1:4){
    myGO <- gmtPathways(paste0(dir,gmts[i]))
    fgRes <- fgsea(pathways = myGO, stats = gene_list, minSize=2, eps=0)
    
    ## Filter FGSEA by using gage results. Must be significant and in same direction to keep 
    #gaRes <- gage(gene_list, gsets=myGO, same.dir=TRUE, set.size =c(2,Inf))
    
    #ups <- as.data.frame(gaRes$greater) %>%
      #tibble::rownames_to_column("Pathway") %>%
      #dplyr::filter(!is.na(p.geomean))#%>% dplyr::filter(q.val < 0.05)
      #dplyr::select("Pathway")
     
    #downs <- as.data.frame(gaRes$less) %>%
      #tibble::rownames_to_column("Pathway") %>%
      #dplyr::filter(!is.na(p.geomean)) %>% dplyr::filter(q.val < 0.05)
      #dplyr::select("Pathway")
    
    #print(dim(rbind(ups,downs)))
    #keepups = fgRes[fgRes$NES > 0 & !is.na(match(fgRes$pathway, ups$Pathway)), ]
    # keepdowns = fgRes[fgRes$NES < 0 & !is.na(match(fgRes$pathway, downs$Pathway)), ]
    
    keepups <- fgRes[(fgRes$NES > 0 & fgRes$padj < 0.05),]
    keepdowns <-  fgRes[(fgRes$NES < 0 & fgRes$padj < 0.05),]
    
    ### Collapse redundant pathways
    Up = fgsea::collapsePathways(keepups, pathways = myGO, stats = gene_list, pval.threshold = 0.05)
    Down = fgsea::collapsePathways(keepdowns, myGO, gene_list, pval.threshold = 0.05) 
    
    fgRes = fgRes[!is.na(match(fgRes$pathway, c( Up$mainPathways, Down$mainPathways))), ] %>% arrange(desc(NES))
    fgRes$Enrichment = ifelse(fgRes$NES > 0, "Up-regulated", "Down-regulated")
    
    # filtRes = rbind(head(fgRes, n = 10),tail(fgRes, n = 10 ))
    # 
    # g = ggplot(filtRes, aes(reorder(pathway, NES), NES)) +
    #   geom_segment(aes(reorder(pathway, NES), xend=pathway, y=0, yend=NES)) +
    #   geom_point(size=5, aes( fill = Enrichment), shape=21, stroke=2) +
    #   scale_fill_manual(values = c("Down-regulated" = "dodgerblue", "Up-regulated" = "firebrick") ) +
    #   coord_flip() +
    #   labs(x="Pathway", y="Normalized Enrichment Score", title=paste0("GSEA - ",cate[i])) +
    #   theme(text = element_text(size="15"))+
    #   theme_minimal()
    
    assign(paste0("result_",cate[i]), fgRes)
    
  }
  
  output = list("result_GO_BP" = result_GO_BP, "result_GO_CC" = result_GO_CC, "result_GO_MF" = result_GO_MF, "result_Kegg" = result_Kegg)

  return(output)
}

gsa_changePathwayID <- function(kegg) {
  #kegg <- read.csv("C:/Users/hylee4161/Downloads/GSA_EnrichR_Kegg_result_2021-06-22.csv")
  lines <- readLines("http://www.kegg.jp/kegg-bin/download_htext?htext=br08901.keg&format=htext" )
  pathways <- do.call(
    rbind,
    str_split( grep( "^[ABCD]\\s+\\d{5}\\s+.*?$", lines, value=TRUE ), "\\s{2,}" )
  )
  pathways <- as.data.frame( pathways )[-1]
  colnames( pathways )  <- c( "kegg_id", "Term" )
  
  kegg_info <- merge(pathways, kegg, by="Term")
  
  return(kegg_info)
}

gsea_changePathwayID <- function(kegg) {
  myPath <- read_excel("./base/c2.cp.kegg.v7.4.symbols.xlsx")
  
  # myGO <- cbind(kegg_id,myPath)
  # write.xlsx(myGO,"./base/c2.cp.kegg.v7.1.symbols.xlsx", row.names = F)
  
  kegg <- merge(myPath,kegg,by="pathway")
  
  lines <- readLines("http://www.kegg.jp/kegg-bin/download_htext?htext=br08901.keg&format=htext")
  pathways <- do.call(
    rbind,
    str_split( grep( "^[ABCD]\\s+\\d{5}\\s+.*?$", lines, value=TRUE ), "\\s{2,}" )
  )
  pathways <- as.data.frame( pathways )[-1]
  colnames( pathways )  <- c( "kegg_id", "pathway" )
  
  kegg_info <- merge(pathways, kegg, by="kegg_id")
  colnames(kegg_info)[2] <- "pathway"
  kegg_info <- kegg_info[,-c(3,4)]
  kegg_info <- kegg_info[,c(1:7,ncol(kegg_info))]
  kegg_info <- kegg_info[order(-abs(kegg_info$NES)),]
  return(kegg_info)
}

string_url_builder <- function(organism, gene){
  if(organism == "Homo Sapiens"){
    input_organism  <-  '9606'
  } else{
    input_organism <-  '10090'
  }
  if(length(gene)==1){
    address  <-  "network?identifier="
  } else{
    address  <-  "network?identifiers="
    gene <- paste(gene,collapse = "%0d")
  }
  URL <- paste0("http://string-db.org/api/image/",
                address,
                gene,
                "&species=",input_organism,
                "&network_flavor=confidence&hide_disconnected_nodes=1")
  return(URL)
}

string_image_download_url_builder <- function(organism, gene){
  if(organism == "Homo Sapiens"){
    input_organism  <-  '9606'
  } else{
    input_organism <-  '10090'
  }
  if(length(gene)==1){
    address  <-  "network?identifier="
  } else{
    address  <-  "network?identifiers="
    gene <- paste(gene,collapse = "%0d")
  }
  URL <- paste0("http://string-db.org/api/highres_image/",
                address,
                gene,
                "&species=",input_organism,
                "&network_flavor=confidence&hide_disconnected_nodes=1")
  return(URL)
}

string_tsv_download_url_builder <- function(organism, gene){
  if(organism == "Homo Sapiens"){
    input_organism  <-  '9606'
  } else{
    input_organism <-  '10090'
  }
  if(length(gene)==1){
    address  <-  "network?identifier="
  } else{
    address  <-  "network?identifiers="
    gene <- paste(gene,collapse = "%0d")
  }
  URL <- paste0("http://string-db.org/api/tsv/",
                address,
                gene,
                "&species=",input_organism)
  return(URL)
}
