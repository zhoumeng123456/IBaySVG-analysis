library(ComplexUpset)
library(UpSetR)

##1.load the data for dlpfc 
dir="RealData/result_data/realdata svgene"

#choose the dataset
load(here::here(dir,"data_dlpfc_acrossdonor_svgene_list.RData"))
#alternative1: data_dlpfc_samedonor_svgene_list.RData for dlpfc dataset with same donor.
#alternative2: data_scc_svgene_list.RData for scc dataset.
svgene_list=data_dlpfc_acrossdonor_svgene_list


##2.Format data
#for dlpfc dataset with same donor or across donors
#for dlpfc dataset with same donor or across donors
uni_set=list("spVC-Union"=svgene_list[[1]][[5]],"HEARTSVG-Union"=svgene_list[[2]][[5]],
             "SPARKX-Union"=svgene_list[[3]][[5]],"SPARK-Union"=svgene_list[[4]][[5]],
             "nnSVG-Union"=svgene_list[[5]][[5]],"IBaySVG"=svgene_list[[6]],
             "GASTON-Union"=svgene_list[[9]][[5]],"StarTrail-Union"=svgene_list[[8]][[5]]
)
upset_uni_set=fromList(uni_set)

inter_set=list("spVC-Inter"=svgene_list[[1]][[6]],"HEARTSVG-Inter"=svgene_list[[2]][[6]],
               "SPARKX-Inter"=svgene_list[[3]][[6]],"SPARK-Inter"=svgene_list[[4]][[6]],
               "nnSVG-Inter"=svgene_list[[5]][[6]],"IBaySVG"=svgene_list[[6]],
               "GASTON-Inter"=svgene_list[[9]][[6]],"StarTrail-Inter"=svgene_list[[8]][[6]]
)
upset_inter_set=fromList(inter_set)

integ_set=list("spVC-PASTE"=svgene_list[[1]][[9]],"HEARTSVG-PASTE"=svgene_list[[2]][[9]],
               "SPARKX-PASTE"=svgene_list[[3]][[9]],"SPARK-PASTE"=svgene_list[[4]][[9]],
               "nnSVG-PASTE"=svgene_list[[5]][[9]],"IBaySVG"=svgene_list[[6]],
               DESpace=svgene_list[[7]],
               "GASTON-PASTE"=svgene_list[[9]][[7]],"StarTrail-PASTE"=svgene_list[[8]][[7]]
)
upset_integ_set=fromList(integ_set)

cauthy_set=list("spVC-Cauchy"=svgene_list[[1]][[7]],"HEARTSVG-Cauchy"=svgene_list[[2]][[7]],
                "SPARKX-Cauchy"=svgene_list[[3]][[7]],"SPARK-Cauchy"=svgene_list[[4]][[7]],
                "nnSVG-Cauchy"=svgene_list[[5]][[7]],"IBaySVG"=svgene_list[[6]]
)
upset_cauthy_set=fromList(cauthy_set)

HMP_set=list("spVC-HMP"=svgene_list[[1]][[8]],"HEARTSVG-HMP"=svgene_list[[2]][[8]],
             "SPARKX-HMP"=svgene_list[[3]][[8]],"SPARK-HMP"=svgene_list[[4]][[8]],
             "nnSVG-HMP"=svgene_list[[5]][[8]],"IBaySVG"=svgene_list[[6]]
)
upset_HMP_set=fromList(HMP_set)

#for scc dataset
#for scc dataset
uni_set=list("nnSVG-Union"=svgene_list[[5]][[4]],"SPARK-Union"=svgene_list[[4]][[4]],
             "SPARKX-Union"=svgene_list[[3]][[4]],"spVC-Union"=svgene_list[[1]][[4]],
             "HEARTSVG-Union"=svgene_list[[2]][[4]],"StarTrail-Union"=svgene_list[[8]][[4]],
             "GASTON-Union"=svgene_list[[9]][[4]],"IBaySVG"=svgene_list[[6]]
)
upset_uni_set=fromList(uni_set)

inter_set=list("nnSVG-Inter"=svgene_list[[5]][[5]],"SPARK-Inter"=svgene_list[[4]][[5]],
               "SPARKX-Inter"=svgene_list[[3]][[5]],"spVC-Inter"=svgene_list[[1]][[5]],
               "HEARTSVG-Inter"=svgene_list[[2]][[5]],"StarTrail-Inter"=svgene_list[[8]][[5]],
               "GASTON-Inter"=svgene_list[[9]][[5]],"IBaySVG"=svgene_list[[6]]
)
upset_inter_set=fromList(inter_set)

cauthy_set=list("nnSVG-Cauchy"=svgene_list[[5]][[6]],"SPARK-Cauchy"=svgene_list[[4]][[6]],
                "SPARKX-Cauchy"=svgene_list[[3]][[6]],"spVC-Cauchy"=svgene_list[[1]][[6]],
                "HEARTSVG-Cauchy"=svgene_list[[2]][[6]],"IBaySVG"=svgene_list[[6]]
)
upset_cauthy_set=fromList(cauthy_set)

HMP_set=list("nnSVG-HMP"=svgene_list[[5]][[7]],"SPARK-HMP"=svgene_list[[4]][[7]],
             "SPARKX-HMP"=svgene_list[[3]][[7]],"spVC-HMP"=svgene_list[[1]][[7]],
             "HEARTSVG-HMP"=svgene_list[[2]][[7]],"IBaySVG"=svgene_list[[6]]
)
upset_HMP_set=fromList(HMP_set)

integ_set=list("nnSVG-PASTE"=svgene_list[[5]][[8]],"SPARK-PASTE"=svgene_list[[4]][[8]],
               "SPARKX-PASTE"=svgene_list[[3]][[8]],"spVC-PASTE"=svgene_list[[1]][[8]],
               "HEARTSVG-PASTE"=svgene_list[[2]][[8]],"StarTrail-PASTE"=svgene_list[[8]][[6]],
               "GASTON-PASTE"=svgene_list[[9]][[6]],"IBaySVG"=svgene_list[[6]],
               "DESpace"=svgene_list[[7]]
)
upset_integ_set=fromList(integ_set)


##3.plot the upset for each dataset,you need to change the "nsets", "mainbar.y.max" and "set_size.scale_max" to adapt each dataset.
#Here is the example for dlpfc datasets across donors.
dir2="RealData/result_data/upset plot"
png(here::here(dir2,"upset_cauchy_acrossdonor.png"), width = 10, height =2.8, units = "in", res = 300)
upset(upset_cauthy_set,
      nsets = 20, 
      mainbar.y.max=600,
      nintersects= 32, 
      main.bar.color = "#91CAE8", 
      matrix.color="black", 
      sets.bar.color= "#F48892", 
      set_size.show = T, 
      mb.ratio = c(0.65, 0.35),
      order.by="freq",
      set_size.scale_max = 2500,
      text.scale = c(
        1,  # intersection size title
        1,  # intersection size numbers
        1,  # set size title
        1,  # set size numbers
        1,  # set names
        1   # intersection matrix labels
      )
)
dev.off()

png(here::here(dir2,"upset_HMP_acrossdonor.png"), width = 10, height = 2.8, units = "in", res = 300)
upset(upset_HMP_set,
      nsets = 20, 
      mainbar.y.max=1500,
      nintersects= 32, 
      main.bar.color = "#91CAE8", 
      matrix.color="black", 
      sets.bar.color= "#F48892", 
      set_size.show = T, 
      mb.ratio = c(0.65, 0.35), 
      order.by="freq",
      set_size.scale_max = 4500,
      text.scale = c(
        1,  # intersection size title
        1,  # intersection size numbers
        1,  # set size title
        1,  # set size numbers
        1,  # set names
        1   # intersection matrix labels
      )
)
dev.off()

png(here::here(dir2,"upset_union_acrossdonor.png"), width = 10, height = 3.2, units = "in", res = 300)
upset(upset_uni_set,
      nsets = 20, 
      mainbar.y.max=1200,
      nintersects= 32,
      main.bar.color = "#91CAE8", 
      matrix.color="black",
      sets.bar.color= "#F48892", 
      set_size.show = T, 
      mb.ratio = c(0.6, 0.4), 
      order.by="freq",
      set_size.scale_max = 4800,
      text.scale = c(
        1,  # intersection size title
        1,  # intersection size numbers
        1,  # set size title
        1,  # set size numbers
        1,  # set names
        1   # intersection matrix labels
      )
)
dev.off()

png(here::here(dir2,"upset_inter_acrossdonor.png"), width = 10, height = 3.2, units = "in", res = 300)
upset(upset_inter_set,
      nsets = 20,
      mainbar.y.max=1000,
      nintersects= 32, 
      main.bar.color = "#91CAE8", 
      matrix.color="black",
      sets.bar.color= "#F48892", 
      set_size.show = T, 
      mb.ratio = c(0.6, 0.4), 
      order.by="freq",
      set_size.scale_max = 2000,
      text.scale = c(
        1,  # intersection size title
        1,  # intersection size numbers
        1,  # set size title
        1,  # set size numbers
        1,  # set names
        1   # intersection matrix labels
      )
)
dev.off()

png(here::here(dir2,"upset_paste_acrossdonor.png"), width = 10, height = 3.2, units = "in", res = 300)
upset(upset_integ_set,
      nsets = 20, 
      mainbar.y.max=1700,
      nintersects= 20,
      main.bar.color = "#91CAE8", 
      matrix.color="black", 
      sets.bar.color= "#F48892", 
      set_size.show = T, 
      mb.ratio = c(0.6, 0.4), 
      order.by="freq",
      set_size.scale_max = 5000,
      text.scale = c(
        1,  # intersection size title
        1,  # intersection size numbers
        1,  # set size title
        1,  # set size numbers
        1,  # set names
        1   # intersection matrix labels
      )
)
dev.off()



##4.layout and save (the result of all datasets have benn stored in the dictionary:"RealData/result_data/upset plot")
dataset="acrossdonor" #alternative:samedonor,scc
img_paths <- c(
  paste0("upset_union_",dataset,".png"),
  paste0("upset_inter_",dataset,".png"),
  paste0("upset_cauchy_",dataset,".png"),
  paste0("upset_HMP_",dataset,".png"),
  paste0("upset_paste_",dataset,".png")
)

images <- lapply(1:5, function(i) {
    img <- image_read(here::here(dir2,img_paths[i]))

    image_annotate(img, 
                   text = LETTERS[i], 
                   size = 50, 
                   color = "black",
                   gravity = "northwest",  
                   location = "+20+20")    
  })
  
combined <- image_append(do.call(c, images), stack = TRUE)
print(combined)
image_write(combined, here::here(dir2,paste0("dlp_upset_combined_",dataset,".png")))
















