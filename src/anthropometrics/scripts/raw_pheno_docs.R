
##
#
# This script installs the p
#
# Requires R 4.x.x, run in conda R_4.1.
# Execute from the root of the repo.
#
##

# Libararies

library <- "~/R/R_4.1"

# This randomly fails the first time but not the second, seems to be an error in the latest R versions

loadLibraries <- function() {
  
    library(memoise, lib.loc = library)
    library(conflicted, lib.loc = library)
    library(foreign, lib.loc = library)
    library(stringr, lib.loc = library)
    library(glue, lib.loc = library)
    library(crayon, lib.loc = library)
    library(tidyr, lib.loc = library)
    library(dplyr, lib.loc = library)
    library(janitor, lib.loc = library)
    
}

tryCatch(
    {
        
        loadLibraries()
        
    }, error = function(error_condition) {
        
        loadLibraries()
        
    }
)

# Paths from command line
args <- commandArgs(TRUE)

linkageFolder <- args[1]
mobaQuesFolder <- args[2]
tablesFolder <- args[3]
docsFolder <- args[4]
projectNumber <- args[5]

# Docs
docsFile <- file.path(docsFolder, "data.md")


# Columns that should be skipped

excludedPhenos <- c(
  paste0("PREG_ID_", projectNumber), 
  paste0("M_ID_",  projectNumber), 
  paste0("F_ID_",  projectNumber), 
  paste0("P_ID_",  projectNumber), 
  paste0("PID_MSIS_",  projectNumber), 
  paste0("PID_K_", projectNumber)
)


# Function generating the docs for a sav file

generateDocsSav <- function(
    savFolder,
    savFileName,
    savDocsFile,
    savDocsFileName,
    exportFolder
) {
    
    write(x = glue("# {savFileName}\n"), file = savDocsFile, append = F)
    
    savDF <- read.spss(
        file = file.path(savFolder, savFileName), 
        use.value.labels = T, 
        to.data.frame = T, 
        stringsAsFactors = F
    )
    
    labels <- as.data.frame(
        x = attr(savDF, "variable.labels"), 
        stringsAsFactors = F
    )
    
    for (colName in names(savDF)) {
      
      if (is.character(savDF[[colName]])) {
        
        savDF[[colName]] <- str_trim(savDF[[colName]])
        
      }
    }
    
    newName <- substr(savFileName, 1, nchar(savFileName) - 4)
    
    write.table(
        x = savDF,
        file = gzfile(file.path(exportFolder, glue("{newName}.gz"))),
        col.names = T,
        row.names = F,
        sep = "\t",
        quote = F
    )
    
    if (length(names(labels)) > 0) {
        
        names(labels) <- "description"
        
    }
    
    labels$pheno <- row.names(labels)
    row.names(labels) <- NULL
    
    write.table(
        x = labels,
        file = gzfile(file.path(exportFolder, glue("{newName}.labels.gz"))),
        col.names = T,
        row.names = F,
        sep = "\t",
        quote = F
    )
    
    labels <- labels %>% 
        filter(
            ! pheno %in% excludedPhenos
        ) %>% 
        mutate(
            sav_files = savFileName
        )
    
    if (nrow(labels) > 0) {
        
        for (phenoI in 1:nrow(labels)) {
            
            phenoName <- labels$pheno[phenoI]
            
            phenoDescription <- labels$description[phenoI]
            
            if (!is.null(phenoDescription) && length(phenoDescription) > 0 && phenoDescription != "") {
                
                write(
                    x = glue("- [{phenoName}]({savDocsFileName}#{phenoName}): {phenoDescription}"), 
                    file = savDocsFile, 
                    append = T
                )
                
            } else {
                
                write(
                    x = glue("- [{phenoName}]({savDocsFileName}#{phenoName})"), 
                    file = savDocsFile, 
                    append = T
                )
                
            }
        }
        
        write(
            x = "\n", 
            file = savDocsFile, 
            append = T
        )
        
        for (phenoI in 1:nrow(labels)) {
            
            phenoName <- labels$pheno[phenoI]
            
            phenoDescription <- labels$description[phenoI]
            
            write(
                x = glue("### {phenoName}"), 
                file = savDocsFile, 
                append = T
            )
            
            if (!is.null(phenoDescription) && length(phenoDescription) > 0 && phenoDescription != "") {
                
                write(
                    x = glue("{phenoDescription}\n"), 
                    file = savDocsFile, 
                    append = T
                )
                
            }
            
            
            write(
                x = "\n", 
                file = savDocsFile, 
                append = T
            )
            if (length(unique(savDF[[phenoName]][!is.na(savDF[[phenoName]])])) <= 20) {
                
                write(
                    x = "| Category | n |", 
                    file = savDocsFile, 
                    append = T
                )
                write(
                    x = "| -------- | - |", 
                    file = savDocsFile, 
                    append = T
                )
                
                for (category in unique(savDF[[phenoName]][!is.na(savDF[[phenoName]])])) {
                    
                    categoryN <- sum(!is.na(savDF[[phenoName]]) & savDF[[phenoName]] == category)
                    
                    write(
                        x = glue("| {category} | {categoryN} |"), 
                        file = savDocsFile, 
                        append = T
                    )
                    
                }
                
                write(
                    x = glue("| NA | {sum(is.na(savDF[[phenoName]]))} |"), 
                    file = savDocsFile, 
                    append = T
                )
                
            } else {
                
                phenoSummary <- summary(savDF[[phenoName]])
                phenoSummaryNames <- names(phenoSummary)
                
                write(
                    x = "| Key | Value |", 
                    file = savDocsFile, 
                    append = T
                )
                write(
                    x = "| --- | ----- |", 
                    file = savDocsFile, 
                    append = T
                )
                
                for (summaryI in 1:length(phenoSummary)) {
                    
                    write(
                        x = glue("| {phenoSummaryNames[summaryI]} | {phenoSummary[summaryI]} |"), 
                        file = savDocsFile, 
                        append = T
                    )
                    
                }
            }
            
            write(
                x = "\n", 
                file = savDocsFile, 
                append = T
            )
        }
    } else {
        
        write(
            x = glue("> Only identifier columns found"), 
            file = savDocsFile, 
            append = T
        )
        
    }
    
    return(labels)
    
}


# Function generating the docs for a csv file

generateDocsCsv <- function(
    csvFolder,
    csvFileName,
    csvDocsFile,
    csvDocsFileName,
    exportFolder,
    separator
) {
    
    write(x = glue("# {csvFileName}\n"), file = csvDocsFile, append = F)
    
    csvDF <- read.table(
        file = file.path(csvFolder, csvFileName), 
        header = T,
        sep = separator,
        stringsAsFactors = F,
        quote = "",
        comment.char = "",
        fill = T,
        row.names = NULL
    )
    
    if (names(csvDF)[1] == "row.names") {
        
        tempHeader <- names(csvDF)[2:length(names(csvDF))]
        
        csvDF <- csvDF[, 1:(ncol(csvDF) - 1)]
        names(csvDF) <- tempHeader
        
    }
    
    labels <- data.frame(
        pheno = names(csvDF), 
        description = names(csvDF), 
        stringsAsFactors = F
    )
    
    for (colName in names(csvDF)) {
      
      if (is.character(csvDF[[colName]])) {
        
        csvDF[[colName]] <- str_trim(csvDF[[colName]])
        
      }       
    }
    
    newName <- substr(csvFileName, 1, nchar(csvFileName) - 4)
    
    write.table(
        x = csvDF,
        file = gzfile(file.path(exportFolder, glue("{newName}.gz"))),
        col.names = T,
        row.names = F,
        sep = "\t",
        quote = F
    )
    
    write.table(
        x = labels,
        file = gzfile(file.path(exportFolder, glue("{newName}.labels.gz"))),
        col.names = T,
        row.names = F,
        sep = "\t",
        quote = F
    )
    
    labels <- labels %>% 
        filter(
            ! pheno %in% excludedPhenos
        ) %>% 
        mutate(
            csv_files = csvFileName
        )
    
    if (nrow(labels) > 0) {
        
        for (phenoI in 1:nrow(labels)) {
            
            phenoName <- labels$pheno[phenoI]
            
            phenoDescription <- labels$description[phenoI]
            
            if (!is.null(phenoDescription) && length(phenoDescription) > 0 && phenoDescription != "") {
                
                write(
                    x = glue("- [{phenoName}]({csvDocsFileName}#{phenoName}): {phenoDescription}"), 
                    file = csvDocsFile, 
                    append = T
                )
                
            } else {
                
                write(
                    x = glue("- [{phenoName}]({csvDocsFileName}#{phenoName})"), 
                    file = csvDocsFile, 
                    append = T
                )
                
            }
        }
        
        write(
            x = "\n", 
            file = csvDocsFile, 
            append = T
        )
        
        for (phenoI in 1:nrow(labels)) {
            
            phenoName <- labels$pheno[phenoI]
            
            phenoDescription <- labels$description[phenoI]
            
            write(
                x = glue("### {phenoName}"), 
                file = csvDocsFile, 
                append = T
            )
            
            if (!is.null(phenoDescription) && length(phenoDescription) > 0 && phenoDescription != "") {
                
                write(
                    x = glue("{phenoDescription}\n"), 
                    file = csvDocsFile, 
                    append = T
                )
                
            }
            
            
            write(
                x = "\n", 
                file = csvDocsFile, 
                append = T
            )
            if (length(unique(csvDF[[phenoName]][!is.na(csvDF[[phenoName]])])) <= 20) {
                
                write(
                    x = "| Category | n |", 
                    file = csvDocsFile, 
                    append = T
                )
                write(
                    x = "| -------- | - |", 
                    file = csvDocsFile, 
                    append = T
                )
                
                for (category in unique(csvDF[[phenoName]][!is.na(csvDF[[phenoName]])])) {
                    
                    categoryN <- sum(!is.na(csvDF[[phenoName]]) & csvDF[[phenoName]] == category)
                    
                    write(
                        x = glue("| {category} | {categoryN} |"), 
                        file = csvDocsFile, 
                        append = T
                    )
                    
                }
                
                write(
                    x = glue("| NA | {sum(is.na(csvDF[[phenoName]]))} |"), 
                    file = csvDocsFile, 
                    append = T
                )
                
            } else {
                
                phenoSummary <- summary(csvDF[[phenoName]])
                phenoSummaryNames <- names(phenoSummary)
                
                write(
                    x = "| Key | Value |", 
                    file = csvDocsFile, 
                    append = T
                )
                write(
                    x = "| --- | ----- |", 
                    file = csvDocsFile, 
                    append = T
                )
                
                for (summaryI in 1:length(phenoSummary)) {
                    
                    write(
                        x = glue("| {phenoSummaryNames[summaryI]} | {phenoSummary[summaryI]} |"), 
                        file = csvDocsFile, 
                        append = T
                    )
                    
                }
            }
            
            write(
                x = "\n", 
                file = csvDocsFile, 
                append = T
            )
        }
    } else {
        
        write(
            x = glue("> Only identifier columns found"), 
            file = csvDocsFile, 
            append = T
        )
        
    }
    
    return(labels)
    
}


# Set up docs

write(
    x = "# Phenotypic variables\n", 
    file = docsFile, 
    append = F
)


# MoBa ques folder

mobaQuesDocsFolder <- file.path(docsFolder, "moba_questionnaire")

if (!file.exists(mobaQuesDocsFolder)) {
    
    dir.create(
        path = mobaQuesDocsFolder,
        showWarnings = T,
        recursive = T
    )
    
}

write(
    x = "### MoBa Questionnaires", 
    file = docsFile, 
    append = T
)

for (fileName in list.files(mobaQuesFolder)) {
    
    if (endsWith(fileName, ".sav")) {
        
        print(glue("{Sys.time()} - Processing {fileName}"))
        
        shortName <- substr(fileName, 1, nchar(fileName) - 4)
        docsFileName <- paste0(shortName, ".md")
        
        write(
            x = glue("- [{fileName}](moba_questionnaire/{docsFileName})"), 
            file = docsFile, 
            append = T
        )
        
        savDocsFile = file.path(mobaQuesDocsFolder, docsFileName)
        
        tablesSubFolder <- file.path(tablesFolder, "moba_ques")
        
        if (!file.exists(tablesSubFolder)) {
          
          dir.create(
            path = tablesSubFolder,
            showWarnings = T,
            recursive = T
          )
          
        }
        
        generateDocsSav(
            savFolder = mobaQuesFolder,
            savFileName = fileName,
            savDocsFileName = docsFileName,
            savDocsFile = savDocsFile,
            exportFolder = tablesSubFolder
        )
        
    }
}


# Linkage folder

for (fileName in list.files(linkageFolder)) {
  
  if (endsWith(fileName, "sav")) {
    
    savDF <- read.spss(
      file = file.path(linkageFolder, fileName), 
      use.value.labels = T, 
      to.data.frame = T, 
      stringsAsFactors = F
    )
    
    labels <- as.data.frame(
      x = attr(savDF, "variable.labels"), 
      stringsAsFactors = F
    )
    
    for (name in names(savDF)) {
      
      if (is.character(savDF[[name]])) {
        
        savDF[[name]] <- str_trim(savDF[[name]])
        
      }
    }
    
    newName <- substr(fileName, 1, nchar(fileName) - 4)
    
    tablesSubFolder <- file.path(tablesFolder, "linkage")
    
    if (!file.exists(tablesSubFolder)) {
      
      dir.create(
        path = tablesSubFolder,
        showWarnings = T,
        recursive = T
      )
      
    }
    
    write.table(
      x = savDF,
      file = gzfile(file.path(tablesSubFolder, glue("{newName}.gz"))),
      col.names = T,
      row.names = F,
      sep = "\t",
      quote = F
    )
    
    
    if (length(names(labels)) > 0) {
      
      names(labels) <- "description"
      
    }
    
    labels$pheno <- row.names(labels)
    row.names(labels) <- NULL
    
    write.table(
      x = labels,
      file = gzfile(file.path(tablesSubFolder, glue("{newName}.labels.gz"))),
      col.names = T,
      row.names = F,
      sep = "\t",
      quote = F
    )
    
  }
}
