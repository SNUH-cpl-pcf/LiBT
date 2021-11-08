ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)){
    for(i in 1:length(new.pkg)){
      install.packages(new.pkg[i], dependencies = TRUE,logical.return=TRUE)
    }
  }
  require_pkg <- sapply(pkg, require, character.only = TRUE)
  isFalse <- as.numeric(which(require_pkg == F))
  if(length(isFalse) != 0){
    if (!requireNamespace("BiocManager", quietly = TRUE)){
      install.packages("BiocManager")
    } else{
      for(j in length(isFalse)){
        BiocManager::install(pkg[isFalse[j]])
      }
    }
  }
}

#DEP, edgeR, fgsea, pathview install by BiocManager
pkg <- c("shiny","shinyjs","shinydashboard","shinydashboardPlus",
         "shinyWidgets","shinyalert","shinycssloaders","shinyjqui", "shinybusy",
         "sortable","dplyr","readr","stringr","readxl","future","promises","sparklyr",
         "DEP","sortable","ggplot2","edgeR","factoextra","enrichR","tibble",
         "fgsea","pathview","zip")
ipak(pkg)

shinyApp(ui=ui, server=server)

# BiocManager::install("DEP")
# BiocManager::install("edgeR")
# BiocManager::install("fgsea")
# BiocManager::install("pathview")
# BiocManager::install("DEP")
#install.packages("future")
# install.packages("promises")
# install.packages("sparklyr")
