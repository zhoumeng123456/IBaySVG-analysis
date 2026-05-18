#produce the benchmark
library(spatialLIBD)#BiocManager::install("spatialLIBD")
library(dplyr)
library(readxl)


############################################################################################################################################
                                                    #1.obtain the benchmark list

############################################################################################################################################
#load the benchmark from Zeng et al.~(2012)
df <- read_excel(here::here("RealData/result_data/benchmark","model_result.xlsx"), sheet = "Table S5")
zeng_moly_benchmark=df$gene_name

#load the benchmark from Maynard et al.~(2021)
load(here::here("RealData/result_data/benchmark","Human_DLPFC_Visium_modeling_results.Rdata"))
load(here::here("RealData/result_data/benchmark","Human_DLPFC_Visium_processedData_sce_scran_sce_layer_spatialLIBD.Rdata"))

colData(sce_layer)$spatialLIBD <- TRUE
system.time(
  sig_genes <-
    sig_genes_extract_all(
      n = 150,
      modeling_results = modeling_results,
      sce_layer = sce_layer
    )
)
maynard_moly_benchmark=sig_genes@listData[["gene"]][which(sig_genes@listData[["model_type"]]=="enrichment")]
select_benchmark=Reduce(union,list(maynard_moly_benchmark,zeng_moly_benchmark))


############################################################################################################################################
                                                  #2.compare the result

############################################################################################################################################
#function for calculating metrics
benchmark_result=function(svgenelist,select_benchmark){
  method_names <- names(svgenelist)
  if (is.null(method_names)) method_names <- paste0("Method_", seq_along(svgenelist))
  
  n_scenarios <- length(svgenelist[[1]])
  scenario_names <- paste0("Scenario_", seq_len(n_scenarios))
  
  result_list <- list()
  
  for (m in c(1:5)) {
    method <- method_names[m]
    for (s in seq_len(n_scenarios)) {
      predicted <- svgenelist[[m]][[s]]
      metrics <- calc_metrics(predicted, select_benchmark)
      
      result_list[[length(result_list) + 1]] <- data.frame(
        Method = method,
        Scenario = scenario_names[s],
        Recall = metrics["Recall"],
        Precision = metrics["Precision"],
        F1 = metrics["F1"]
      )
    }
  }
  for (m in c(8:9)) {
    method <- method_names[m]
    for (s in seq_len(7)) {
      predicted <- svgenelist[[m]][[s]]
      metrics <- calc_metrics(predicted, select_benchmark)
      
      result_list[[length(result_list) + 1]] <- data.frame(
        Method = method,
        Scenario = scenario_names[s],
        Recall = metrics["Recall"],
        Precision = metrics["Precision"],
        F1 = metrics["F1"]
      )
    }
  }
  result_df <- do.call(rbind, result_list)
  rownames(result_df)=NULL
  
  result_df_summary <- result_df %>%
    mutate(Scenario = case_when(
      Scenario %in% paste0("Scenario_", 1:4) ~ "AVE",
      Scenario == "Scenario_5" ~ "Union",
      Scenario == "Scenario_6" ~ "Inter",
      Scenario == "Scenario_7" ~ "Cauchy",
      Scenario == "Scenario_8" ~ "HMP",
      Scenario == "Scenario_9" ~ "PASTE",
      TRUE ~ Scenario
    )) %>%
    group_by(Method, Scenario) %>%
    summarise(
      Recall = round(mean(Recall, na.rm = TRUE),3),
      Precision =round( mean(Precision, na.rm = TRUE),3),
      F1 = round(mean(F1, na.rm = TRUE),3),
      .groups = "drop"
    ) %>%
    arrange(Method, factor(Scenario, levels = c("AVE", "Union", "Inter", "Cauchy", "HMP", "PASTE")))
  
  result_df_summary=as.data.frame(result_df_summary)
  
  result_df_summary <- rbind(result_df_summary, c("Proposed","multi",calc_metrics(svgenelist[[6]], select_benchmark)))
  result_df_summary <- rbind(result_df_summary, c("DEspace","multi",calc_metrics(svgenelist[[7]], select_benchmark)))
  
  print(result_df_summary)
}

#Inner function of benchmark_result
calc_metrics <- function(predicted, benchmark) {
  TP <- length(intersect(predicted, benchmark))
  recall <- TP / length(benchmark)
  precision <- TP / length(predicted)
  F1 <- ifelse((precision + recall) > 0, 
               2 * precision * recall / (precision + recall), 
               0)
  recall=round(recall,3)
  precision=round(precision,3)
  F1=round(F1,3)
  return(c(Recall = recall, Precision = precision, F1 = F1))
}

#produce the result
dataset="acrossdonor"#alternative: "samedonor"
load(here::here("RealData/result_data/realdata svgene",paste0("data_dlpfc_",dataset,"_svgene_list.RData")))
benchmark_result(data_dlpfc_acrossdonor_svgene_list,select_benchmark)

















