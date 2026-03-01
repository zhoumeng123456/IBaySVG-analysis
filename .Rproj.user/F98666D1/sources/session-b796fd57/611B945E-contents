################################################################################################################################################
                                                        #comparison between ZINB and NB model for Tables35-s36


################################################################################################################################################
##1.example for using the IBaySVG model without inflation structure
#load the code. Notice:The non-spatial model is included exclusively for ablation analysis and does not constitute the main component of the IBaySVG methodology.
source(here::here("Simulations/result_data/noinflation","IBaySVG_noinflation.R"))

#load the example data
load(here::here("Simulations/result_data/noinflation","example_linear_inf5.RData"))

#produce the result
result <- NBIMSVG_noinflation(spelist = example_linear_inf5[[1]],c_alpha = example_linear_inf5[[2]],num_cores =10,max_iter=200)
tp=mean(result[[3]][1:500],na.rm=TRUE)
fp=mean(result[[3]][501:5000],na.rm=TRUE)
tp
fp


##2.load the results of comparison between ZINB and NB model 
#function for comparison between ZINB and NB model
make_table <- function(df, inf_name) {
  df %>%
    filter(inf == inf_name,
           set %in% c("F1", "FPR", "TPR")  # 只保留三种) %>%
    ) %>%
    group_by(pattern, inte_situ, method, set) %>%
    summarise(
      mean_val = round(mean(value, na.rm = TRUE),3),
      sd_val = round(sd(value, na.rm = TRUE),3),
      .groups = "drop"
    ) %>%
    mutate(mean_sd = sprintf("%.3f (%.3f)", mean_val, sd_val)) %>%
    dplyr::select(-mean_val, -sd_val) %>%
    pivot_wider(
      names_from = c(method, set),
      values_from = mean_sd
    ) %>%
    arrange(
      factor(pattern, levels = c("linear", "focal", "period")),
      factor(inte_situ, levels = c("scenario 1", "scenario 2", "scenario 3", "scenario 4"))
    )%>%
    dplyr::select(
      pattern, inte_situ,
      IBaySVG_TPR,IBaySVG_FPR, IBaySVG_F1,
      IBaySVG_noinflation_TPR, IBaySVG_noinflation_FPR, IBaySVG_noinflation_F1
    )
}

#load the result of IBaySVG
total_result=c()
for(i in c("IBaySVG")){
  for(pattern_choose in c("linear","focal","period")){
    for(inf_choose in c("inf1","inf3","inf5")){
      for(inte_situ_choose in c("scenario 1","scenario 2","scenario 3","scenario 4")){
        if(inf_choose=="inf1"){
          infname="low_dropout"
        }else if(inf_choose=="inf3"){
          infname="middle_dropout"
        }else{
          infname="high_dropout"
        }
        if(inte_situ_choose =="scenario 1"){
          setname="setting1"
        }else if(inte_situ_choose =="scenario 2"){
          setname="setting2"
        }else if(inte_situ_choose =="scenario 3"){
          setname="setting3"
        }else{
          setname="setting4"
        }
        value_result=read.csv(here::here(paste0("Simulations/result_data/simulation result/basic simulation results/",pattern_choose,"/",setname,"_",infname),paste0(i,"_basic_simu.csv")))
        TPR_value=value_result$TPR
        FPR_value=value_result$FPR
        F1_value=value_result$F1
        result_type_domain=c("TPR","FPR","F1")
        for(result_type in seq(result_type_domain)){
          temp_result=data.frame(
            method = "IBaySVG",
            set = result_type_domain[result_type],
            value = value_result[,result_type],
            inf = inf_choose,
            inte_situ = inte_situ_choose,
            pattern = pattern_choose
          )
          total_result=rbind(total_result,temp_result)
        }
      }
    }
  }
}

#load the result of IBaySVG without inflation structure
total_result_noinflation=c()
for(i in c("IBaySVG")){
  for(pattern_choose in c("linear","focal","period")){
    for(inf_choose in c("inf3","inf5")){
      for(inte_situ_choose in c("scenario 1","scenario 2","scenario 3","scenario 4")){
        if(inf_choose=="inf1"){
          infname="low_dropout"
        }else if(inf_choose=="inf3"){
          infname="middle_dropout"
        }else{
          infname="high_dropout"
        }
        if(inte_situ_choose =="scenario 1"){
          setname="setting1"
        }else if(inte_situ_choose =="scenario 2"){
          setname="setting2"
        }else if(inte_situ_choose =="scenario 3"){
          setname="setting3"
        }else{
          setname="setting4"
        }
        value_result=read.csv(here::here(paste0("Simulations/result_data/noinflation/IBaySVG_noinflation_results/",pattern_choose),paste0(i,"_noinflation_simu_",setname,"_",infname,".csv")))
        TPR_value=value_result$TPR
        FPR_value=value_result$FPR
        F1_value=value_result$F1
        result_type_domain=c("TPR","FPR","F1")
        for(result_type in seq(result_type_domain)){
          temp_result=data.frame(
            method = "IBaySVG_noinflation",
            set = result_type_domain[result_type],
            value = value_result[,result_type],
            inf = inf_choose,
            inte_situ = inte_situ_choose,
            pattern = pattern_choose
          )
          total_result_noinflation=rbind(total_result_noinflation,temp_result)
        }
      }
    }
  }
}


result=rbind(total_result_noinflation,total_result)

#plot the result
table_inf1 <- make_table(result, "inf5")
table_inf2 <- make_table(result, "inf3")
table_inf1[,-c(1:2)]
table_inf2[,-c(1:2)]
























