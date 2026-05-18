##################################################################################################################################################
                                                     #1.example for comparison between two statistics

##################################################################################################################################################
#bfdr control function in IBaySVG
BayFDR <- function(PPI, alpha){
  genenum=length(PPI)
  PPI_sorted = sort(PPI,decreasing = TRUE)
  k = 1
  fdr = 0
  while(fdr < alpha){
    fdr = mean(1 - PPI_sorted[1:k])
    k = k+1
    if(k > length(PPI_sorted)){
      k = length(PPI_sorted);
      break;
    }
  }
  if(genenum<200){
    return.value = max(PPI_sorted[k],0.95)
  }else{
    return.value = PPI_sorted[k]
  }
  return.value = ifelse(is.na(return.value), 0, return.value)
  return(return.value)
}

#load a example for changing the statistics 
#a csv file contained the u_1 and u_2 results of 5000 simulated genes produced by IBaySVG where the first 10% are set to "SV gene" and the remained "non-SV gene"
example_result=read.csv(here::here("Simulations/result_data/compare statistics","example.csv"))

#IBaysvg method procedure in simulation dataset
IBaysvg_result=apply(example_result, 1, max)
thrs=BayFDR(IBaysvg_result, (0.05 / (2*length(IBaysvg_result))))
order=sum(IBaysvg_result>thrs,na.rm = TRUE)#save the order
IBaysvg_result[IBaysvg_result > thrs] <- 1
IBaysvg_result[IBaysvg_result <= thrs] <- 0
tp=  sum(IBaysvg_result[1:500], na.rm = TRUE)/500
fp=  sum(IBaysvg_result[501:5000], na.rm = TRUE)/4500
print(c(tp,fp))

##IBaysvg method with modified statistics in simulation dataset
IBaysvg_modi_result<-sqrt(example_result[,1]^2+example_result[,2]^2)/sqrt(2)
indices <- order(IBaysvg_modi_result, decreasing = TRUE)[1:order]#choose the same order
tp=  sum(indices <= 500, na.rm = TRUE)/500
fp=  sum(indices > 500, na.rm = TRUE)/4500
print(c(tp,fp))



##################################################################################################################################################
                                                    #2.plot complete comparison between two statistics: Tables1-s3

##################################################################################################################################################
library(tidyr)
##function for comparison
make_table <- function(df, inf_name) {
  df %>%
    filter(inf == inf_name,
           set %in% c("F1", "FPR", "TPR")  ) %>%
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
      IBaySVG_TPR, IBaySVG_FPR, IBaySVG_F1,
      IBaySVG_modified_TPR,IBaySVG_modified_FPR, IBaySVG_modified_F1
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

#load the result of IBaySVG_modified
total_result_modi=c()
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
        value_result=read.csv(here::here(paste0("Simulations/result_data/compare statistics/IBaySVG_modistatis_results/",pattern_choose),paste0(i,"_modistatis_simu_",setname,"_",infname,".csv")))
        TPR_value=value_result$TPR
        FPR_value=value_result$FPR
        F1_value=value_result$F1
        result_type_domain=c("TPR","FPR","F1")
        for(result_type in seq(result_type_domain)){
          temp_result=data.frame(
            method = "IBaySVG_modified",
            set = result_type_domain[result_type],
            value = value_result[,result_type],
            inf = inf_choose,
            inte_situ = inte_situ_choose,
            pattern = pattern_choose
          )
          total_result_modi=rbind(total_result_modi,temp_result)
        }
      }
    }
  }
}

result=rbind(total_result,total_result_modi)

#plot
table_inf1 <- make_table(result, "inf1")
table_inf3 <- make_table(result, "inf3")
table_inf5 <- make_table(result, "inf5")

table_inf1
table_inf3
table_inf5


