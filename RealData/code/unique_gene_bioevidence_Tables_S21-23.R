######################################################Table S21-S23 ####################################
######################################################Table S21-S23 ####################################
######################################################Table S21-S23 ####################################
#1.obtain the unique gene sets
#DLPFC-same donor
dataset="dlpfc_samedonor"
load(here::here("RealData/result_data/realdata svgene",paste0("data_",dataset,"_svgene_list.RData")))
all_union_union_list_nogaston=Reduce(union,list(data_dlpfc_samedonor_svgene_list[[1]][[5]],data_dlpfc_samedonor_svgene_list[[3]][[5]],
                                                data_dlpfc_samedonor_svgene_list[[2]][[5]],data_dlpfc_samedonor_svgene_list[[4]][[5]],
                                                data_dlpfc_samedonor_svgene_list[[5]][[5]],data_dlpfc_samedonor_svgene_list[[8]][[5]]))
weakproblem1=setdiff(data_dlpfc_samedonor_svgene_list[[6]],all_union_union_list_nogaston)#76


#DLPFC-across donors
dataset="dlpfc_acrossdonor"
load(here::here("RealData/result_data/realdata svgene",paste0("data_",dataset,"_svgene_list.RData")))
all_union_union_list_nogaston=Reduce(union,list(data_dlpfc_acrossdonor_svgene_list[[1]][[5]],
                                                data_dlpfc_acrossdonor_svgene_list[[2]][[5]],data_dlpfc_acrossdonor_svgene_list[[4]][[5]],
                                                data_dlpfc_acrossdonor_svgene_list[[5]][[5]],data_dlpfc_acrossdonor_svgene_list[[8]][[5]]))
weakproblem2=setdiff(data_dlpfc_acrossdonor_svgene_list[[6]],all_union_union_list_nogaston)#261

#SCC
dataset="scc"
load(here::here("RealData/result_data/realdata svgene",paste0("data_",dataset,"_svgene_list.RData")))
all_union_union_list_nogaston=Reduce(union,list(data_scc_svgene_list[[1]][[4]],data_scc_svgene_list[[2]][[4]],
                                                data_scc_svgene_list[[4]][[4]],data_scc_svgene_list[[5]][[4]],
                                                data_scc_svgene_list[[8]][[4]],data_scc_svgene_list[[9]][[4]]))
weakproblem3=setdiff(data_scc_svgene_list[[6]],all_union_union_list_nogaston)
