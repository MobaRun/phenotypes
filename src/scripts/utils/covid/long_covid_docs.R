
writeLongCovidDocs <- function(
    longCovidDocsFolder,
    longCovidDocsFile,
    phenoDF
) {
    
    
    
    roles <- c("Child", "Mother", "Father")
    
    write(x = "# Long Covid Phenotypes\n", file = longCovidDocsFile, append = F)
    
    nChildren <- sum(phenoDF$role == "Child")
    nMothers <- sum(phenoDF$role == "Mother")
    nFathers <- sum(phenoDF$role == "Father")
    write(x = glue("Genotyped samples only, PCA outliers removed ({nrow(phenoDF)} individuals, {nChildren} children, {nMothers} mothers, {nFathers} fathers)\n\n"), file = longCovidDocsFile, append = T)
    
    write(x = "## Self report vs msis\n", file = longCovidDocsFile, append = T)
    
    write(x = "### suspected or confirmed by doctor vs. msis\n", file = longCovidDocsFile, append = T)
    
    for (role in roles) {
        
        write(
            x = glue("- {role}\n\n"), 
            file = longCovidDocsFile, 
            append = T
        )
        
        n00 <- sum((is.na(phenoDF$suspected_or_confirmed_covid_doctor_past_14_days) | phenoDF$suspected_or_confirmed_covid_doctor_past_14_days == 0) & is.na(phenoDF$msis_first_registered) & phenoDF$role == role)
        n10 <- sum(!is.na(phenoDF$suspected_or_confirmed_covid_doctor_past_14_days) & phenoDF$suspected_or_confirmed_covid_doctor_past_14_days == 1 & is.na(phenoDF$msis_first_registered) & phenoDF$role == role)
        nX0 <- sum(is.na(phenoDF$suspected_or_confirmed_covid_doctor_past_14_days) & is.na(phenoDF$msis_first_registered) & phenoDF$role == role)
        n01 <- sum((is.na(phenoDF$suspected_or_confirmed_covid_doctor_past_14_days) | phenoDF$suspected_or_confirmed_covid_doctor_past_14_days == 0) & !is.na(phenoDF$msis_first_registered) & phenoDF$role == role)
        n11 <- sum(!is.na(phenoDF$suspected_or_confirmed_covid_doctor_past_14_days) & phenoDF$suspected_or_confirmed_covid_doctor_past_14_days == 1 & !is.na(phenoDF$msis_first_registered) & phenoDF$role == role)
        nX1 <- sum(is.na(phenoDF$suspected_or_confirmed_covid_doctor_past_14_days) & !is.na(phenoDF$msis_first_registered) & phenoDF$role == role)
        
        
        write(
            x = "| in msis | Suspected or confirmed covid by doctor in the past 14 days | n (%) |", 
            file = longCovidDocsFile, 
            append = T
        )
        write(
            x = "| --- | --- | --- |", 
            file = longCovidDocsFile, 
            append = T
        )
        write(
            x = glue("| No | Missing | {nX0} ({100 * nX0 / (n00 + n10 + n01 + n11)} %) |"), 
            file = longCovidDocsFile, 
            append = T
        )
        write(
            x = glue("| No | No | {n00} ({100 * n00 / (n00 + n10 + n01 + n11)} %) |"), 
            file = longCovidDocsFile, 
            append = T
        )
        write(
            x = glue("| No | Yes | {n10} ({100 * n10 / (n00 + n10 + n01 + n11)} %) |"), 
            file = longCovidDocsFile, 
            append = T
        )
        write(
            x = glue("| Yes | Missing | {nX1} ({100 * nX1 / (n00 + n10 + n01 + n11)} %) |"), 
            file = longCovidDocsFile, 
            append = T
        )
        write(
            x = glue("| Yes | No | {n01} ({100 * n01 / (n00 + n10 + n01 + n11)} %) |"), 
            file = longCovidDocsFile, 
            append = T
        )
        write(
            x = glue("| Yes | Yes | {n11} ({100 * n11 / (n00 + n10 + n01 + n11)} %) |\n\n"), 
            file = longCovidDocsFile, 
            append = T
        )
        
    }
    
    
    write(x = "### Tested positive vs. msis\n", file = longCovidDocsFile, append = T)
    
    for (role in roles) {
        
        write(
            x = glue("- {role}\n\n"), 
            file = longCovidDocsFile, 
            append = T
        )
        
        n00 <- sum((is.na(phenoDF$tested_positive) | phenoDF$tested_positive == 0) & is.na(phenoDF$msis_first_registered) & phenoDF$role == role)
        n10 <- sum(!is.na(phenoDF$tested_positive) & phenoDF$tested_positive == 1 & is.na(phenoDF$msis_first_registered) & phenoDF$role == role)
        nX0 <- sum(is.na(phenoDF$tested_positive) & is.na(phenoDF$msis_first_registered) & phenoDF$role == role)
        n01 <- sum((is.na(phenoDF$tested_positive) | phenoDF$tested_positive == 0) & !is.na(phenoDF$msis_first_registered) & phenoDF$role == role)
        n11 <- sum(!is.na(phenoDF$tested_positive) & phenoDF$tested_positive == 1 & !is.na(phenoDF$msis_first_registered) & phenoDF$role == role)
        nX1 <- sum(is.na(phenoDF$tested_positive) & !is.na(phenoDF$msis_first_registered) & phenoDF$role == role)
        
        
        write(
            x = "| in msis | Tested positive | n (%) |", 
            file = longCovidDocsFile, 
            append = T
        )
        write(
            x = "| --- | --- | --- |", 
            file = longCovidDocsFile, 
            append = T
        )
        write(
            x = glue("| No | Missing | {nX0} ({100 * nX0 / (n00 + n10 + n01 + n11)} %) |"), 
            file = longCovidDocsFile, 
            append = T
        )
        write(
            x = glue("| No | No | {n00} ({100 * n00 / (n00 + n10 + n01 + n11)} %) |"), 
            file = longCovidDocsFile, 
            append = T
        )
        write(
            x = glue("| No | Yes | {n10} ({100 * n10 / (n00 + n10 + n01 + n11)} %) |"), 
            file = longCovidDocsFile, 
            append = T
        )
        write(
            x = glue("| Yes | Missing | {nX1} ({100 * nX1 / (n00 + n10 + n01 + n11)} %) |"), 
            file = longCovidDocsFile, 
            append = T
        )
        write(
            x = glue("| Yes | No | {n01} ({100 * n01 / (n00 + n10 + n01 + n11)} %) |"), 
            file = longCovidDocsFile, 
            append = T
        )
        write(
            x = glue("| Yes | Yes | {n11} ({100 * n11 / (n00 + n10 + n01 + n11)} %) |\n\n"), 
            file = longCovidDocsFile, 
            append = T
        )
        
    }
    
    
    write(x = "### Tested positive by PCR vs. msis\n", file = longCovidDocsFile, append = T)
    
    for (role in roles) {
        
        write(
            x = glue("- {role}\n\n"), 
            file = longCovidDocsFile, 
            append = T
        )
        
        n00 <- sum((is.na(phenoDF$tested_positive_pcr) | phenoDF$tested_positive_pcr == 0) & is.na(phenoDF$msis_first_registered) & phenoDF$role == role)
        n10 <- sum(!is.na(phenoDF$tested_positive_pcr) & phenoDF$tested_positive_pcr == 1 & is.na(phenoDF$msis_first_registered) & phenoDF$role == role)
        nX0 <- sum(is.na(phenoDF$tested_positive_pcr) & is.na(phenoDF$msis_first_registered) & phenoDF$role == role)
        n01 <- sum((is.na(phenoDF$tested_positive_pcr) | phenoDF$tested_positive_pcr == 0) & !is.na(phenoDF$msis_first_registered) & phenoDF$role == role)
        n11 <- sum(!is.na(phenoDF$tested_positive_pcr) & phenoDF$tested_positive_pcr == 1 & !is.na(phenoDF$msis_first_registered) & phenoDF$role == role)
        nX1 <- sum(is.na(phenoDF$tested_positive_pcr) & !is.na(phenoDF$msis_first_registered) & phenoDF$role == role)
        
        
        write(
            x = "| in msis | Tested positive by PCR | n (%) |", 
            file = longCovidDocsFile, 
            append = T
        )
        write(
            x = "| --- | --- | --- |", 
            file = longCovidDocsFile, 
            append = T
        )
        write(
            x = glue("| No | Missing | {nX0} ({100 * nX0 / (n00 + n10 + n01 + n11 + nX0 + nX1)} %) |"), 
            file = longCovidDocsFile, 
            append = T
        )
        write(
            x = glue("| No | No | {n00} ({100 * n00 / (n00 + n10 + n01 + n11 + nX0 + nX1)} %) |"), 
            file = longCovidDocsFile, 
            append = T
        )
        write(
            x = glue("| No | Yes | {n10} ({100 * n10 / (n00 + n10 + n01 + n11 + nX0 + nX1)} %) |"), 
            file = longCovidDocsFile, 
            append = T
        )
        write(
            x = glue("| Yes | Missing | {nX1} ({100 * nX1 / (n00 + n10 + n01 + n11 + nX0 + nX1)} %) |"), 
            file = longCovidDocsFile, 
            append = T
        )
        write(
            x = glue("| Yes | No | {n01} ({100 * n01 / (n00 + n10 + n01 + n11 + nX0 + nX1)} %) |"), 
            file = longCovidDocsFile, 
            append = T
        )
        write(
            x = glue("| Yes | Yes | {n11} ({100 * n11 / (n00 + n10 + n01 + n11 + nX0 + nX1)} %) |\n\n"), 
            file = longCovidDocsFile, 
            append = T
        )
    }
    
    
    write(x = "### Tested positive by AB vs. msis\n", file = longCovidDocsFile, append = T)
    
    for (role in roles) {
        
        write(
            x = glue("- {role}\n\n"), 
            file = longCovidDocsFile, 
            append = T
        )
        
        n00 <- sum((is.na(phenoDF$tested_positive_ab) | phenoDF$tested_positive_ab == 0) & is.na(phenoDF$msis_first_registered) & phenoDF$role == role)
        n10 <- sum(!is.na(phenoDF$tested_positive_ab) & phenoDF$tested_positive_ab == 1 & is.na(phenoDF$msis_first_registered) & phenoDF$role == role)
        nX0 <- sum(is.na(phenoDF$tested_positive_ab) & is.na(phenoDF$msis_first_registered) & phenoDF$role == role)
        n01 <- sum((is.na(phenoDF$tested_positive_ab) | phenoDF$tested_positive_ab == 0) & !is.na(phenoDF$msis_first_registered) & phenoDF$role == role)
        n11 <- sum(!is.na(phenoDF$tested_positive_ab) & phenoDF$tested_positive_ab == 1 & !is.na(phenoDF$msis_first_registered) & phenoDF$role == role)
        nX1 <- sum(is.na(phenoDF$tested_positive_ab) & !is.na(phenoDF$msis_first_registered) & phenoDF$role == role)
        
        
        write(
            x = "| in msis | Tested positive by AB | n (%) |", 
            file = longCovidDocsFile, 
            append = T
        )
        write(
            x = "| --- | --- | --- |", 
            file = longCovidDocsFile, 
            append = T
        )
        write(
            x = glue("| No | Missing | {nX0} ({100 * nX0 / (n00 + n10 + n01 + n11 + nX0 + nX1)} %) |"), 
            file = longCovidDocsFile, 
            append = T
        )
        write(
            x = glue("| No | No | {n00} ({100 * n00 / (n00 + n10 + n01 + n11 + nX0 + nX1)} %) |"), 
            file = longCovidDocsFile, 
            append = T
        )
        write(
            x = glue("| No | Yes | {n10} ({100 * n10 / (n00 + n10 + n01 + n11 + nX0 + nX1)} %) |"), 
            file = longCovidDocsFile, 
            append = T
        )
        write(
            x = glue("| Yes | Missing | {nX1} ({100 * nX1 / (n00 + n10 + n01 + n11 + nX0 + nX1)} %) |"), 
            file = longCovidDocsFile, 
            append = T
        )
        write(
            x = glue("| Yes | No | {n01} ({100 * n01 / (n00 + n10 + n01 + n11 + nX0 + nX1)} %) |"), 
            file = longCovidDocsFile, 
            append = T
        )
        write(
            x = glue("| Yes | Yes | {n11} ({100 * n11 / (n00 + n10 + n01 + n11 + nX0 + nX1)} %) |\n\n"), 
            file = longCovidDocsFile, 
            append = T
        )
    }
    
    
    write(x = "## Long-covid penotypes\n", file = longCovidDocsFile, append = T)
    
    phenotypes <- list(
        reduced_smell_taste = "Reduced smell taste (kf120)", 
        brain_fog = "Brain Fog (kf480)", 
        poor_memory = "Poor Memory (kf481)", 
        dizziness = "Dizziness (kf479)", 
        heart_palpitation = "Heart Palpitation (kf476)", 
        fatigue = "Fatigue (kf468)", 
        headache = "Headache (kf484)", 
        skin_rash = "Skin rash (kf487)", 
        anxiety = "Anxiety (kf486)", 
        altered_smell_taste = "Altered smell taste (kf489)", 
        chest_pain = "Chest pain (kf475)", 
        shortness_breath = "Shortness breath (kf470)", 
        lung_function_reduced = "Lung function reduced (kf472)", 
        cough = "Cough (kf471)"
    )
    
    phenoDF$anySymptomShort <- 0 
    phenoDF$anySymptomLong <- 0 
    
    for (phenotype in names(phenotypes)) {
        
        lastPheno <- paste0(phenotype, "_last_reported")
        
        phenoDF <- phenoDF %>% 
            mutate(
                last_pheno_date = as.Date(!!sym(lastPheno)/86400, origin = "1582-10-14"),
                symptome_duration = as.numeric(last_pheno_date - msis_last_registered)
            )
        
        phenoDF$anySymptomShort[!is.na(phenoDF$msis_last_registered) & !is.na(phenoDF[[phenotype]]) & phenoDF[[phenotype]] == 1 & phenoDF$symptome_duration < 90] <- 1
        phenoDF$anySymptomLong[!is.na(phenoDF$msis_last_registered) & !is.na(phenoDF[[phenotype]]) & phenoDF[[phenotype]] == 1 & phenoDF$symptome_duration >= 90] <- 1
        
        write(x = glue("### {phenotypes[[phenotype]]}\n"), file = longCovidDocsFile, append = T)
        
        for (role in roles) {
            
            write(
                x = glue("- {role}\n\n"), 
                file = longCovidDocsFile, 
                append = T
            )
            
            tempDF <- phenoDF[phenoDF$role == role, ]
            
            n00 <- sum(is.na(tempDF$msis_last_registered) & !is.na(tempDF[[phenotype]]) & tempDF[[phenotype]] != 1)
            n01 <- sum(is.na(tempDF$msis_last_registered) & !is.na(tempDF[[phenotype]]) & tempDF[[phenotype]] == 1)
            n0X <- sum(is.na(tempDF$msis_last_registered) & is.na(tempDF[[phenotype]]))
            n10 <- sum(!is.na(tempDF$msis_last_registered) & !is.na(tempDF[[phenotype]]) & tempDF[[phenotype]] != 1)
            n11 <- sum(!is.na(tempDF$msis_last_registered) & !is.na(tempDF[[phenotype]]) & tempDF[[phenotype]] == 1 & tempDF$symptome_duration < 90)
            n12 <- sum(!is.na(tempDF$msis_last_registered) & !is.na(tempDF[[phenotype]]) & tempDF[[phenotype]] == 1 & tempDF$symptome_duration >= 90)
            n1X <- sum(!is.na(tempDF$msis_last_registered) & is.na(tempDF[[phenotype]]))
            
            
            write(
                x = glue("| in msis | {phenotype} | n (%) |"), 
                file = longCovidDocsFile, 
                append = T
            )
            write(
                x = "| --- | --- | --- |", 
                file = longCovidDocsFile, 
                append = T
            )
            write(
                x = glue("| No | Missing | {n0X} ({100 * n0X / (n00 + n01 + n0X + n10 + n11 + n12 + n1X)} %) |"), 
                file = longCovidDocsFile, 
                append = T
            )
            write(
                x = glue("| No | No | {n00} ({100 * n00 / (n00 + n01 + n0X + n10 + n11 + n12 + n1X)} %) |"), 
                file = longCovidDocsFile, 
                append = T
            )
            write(
                x = glue("| No | Yes | {n01} ({100 * n01 / (n00 + n01 + n0X + n10 + n11 + n12 + n1X)} %) |"), 
                file = longCovidDocsFile, 
                append = T
            )
            write(
                x = glue("| Yes | Missing | {n1X} ({100 * n1X / (n00 + n01 + n0X + n10 + n11 + n12 + n1X)} %) |"), 
                file = longCovidDocsFile, 
                append = T
            )
            write(
                x = glue("| Yes | No | {n10} ({100 * n10 / (n00 + n01 + n0X + n10 + n11 + n12 + n1X)} %) |"), 
                file = longCovidDocsFile, 
                append = T
            )
            write(
                x = glue("| Yes | Yes (<90 days) | {n11} ({100 * n11 / (n00 + n01 + n0X + n10 + n11 + n12 + n1X)} %) |"), 
                file = longCovidDocsFile, 
                append = T
            )
            write(
                x = glue("| Yes | Yes (≥90 days) | {n12} ({100 * n12 / (n00 + n01 + n0X + n10 + n11 + n12 + n1X)} %) |\n\n"), 
                file = longCovidDocsFile, 
                append = T
            )
            
        }
        
        phenoDF <- phenoDF %>% 
            select(
                -last_pheno_date, -symptome_duration
            )
    }
    
    write(x = glue("### Any symptom\n"), file = longCovidDocsFile, append = T)
    
    for (role in roles) {
        
        write(
            x = glue("- {role}\n\n"), 
            file = longCovidDocsFile, 
            append = T
        )
        
        tempDF <- phenoDF[phenoDF$role == role, ]
        
        nxx <- sum(is.na(tempDF$msis_last_registered))
        n00 <- sum(!is.na(tempDF$msis_last_registered) & tempDF$anySymptomShort == 0 & tempDF$anySymptomLong == 0)
        n10 <- sum(!is.na(tempDF$msis_last_registered) & tempDF$anySymptomShort == 1 & tempDF$anySymptomLong == 0)
        n11 <- sum(!is.na(tempDF$msis_last_registered) & tempDF$anySymptomShort == 1 & tempDF$anySymptomLong == 1)
        n01 <- sum(!is.na(tempDF$msis_last_registered) & tempDF$anySymptomShort == 0 & tempDF$anySymptomLong == 1)
        
        
        write(
            x = glue("| in msis | symptoms | n (%) |"), 
            file = longCovidDocsFile, 
            append = T
        )
        write(
            x = "| --- | --- | --- |", 
            file = longCovidDocsFile, 
            append = T
        )
        write(
            x = glue("| No | - | {nxx} ({100 * nxx / (nxx + n00 + n10 + n11 + n01)} %) |"), 
            file = longCovidDocsFile, 
            append = T
        )
        write(
            x = glue("| Yes | No | {n00} ({100 * n00 / (nxx + n00 + n10 + n11 + n01)} %) |"), 
            file = longCovidDocsFile, 
            append = T
        )
        write(
            x = glue("| Yes | Yes, all < 90 days | {n10} ({100 * n10 / (nxx + n00 + n10 + n11 + n01)} %) |"), 
            file = longCovidDocsFile, 
            append = T
        )
        write(
            x = glue("| Yes | Yes, some < 90 days, some ≥90 days | {n11} ({100 * n11 / (nxx + n00 + n10 + n11 + n01)} %) |"), 
            file = longCovidDocsFile, 
            append = T
        )
        write(
            x = glue("| Yes | Yes, all ≥90 days | {n01} ({100 * n01 / (nxx + n00 + n10 + n11 + n01)} %) |"), 
            file = longCovidDocsFile, 
            append = T
        )
        
    }
}