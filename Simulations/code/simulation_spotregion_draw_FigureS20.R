library(ggplot2)

spot.region <- c(rep(c(rep(1,20),rep(2,12)),16),rep(c(rep(3,12),rep(1,8),rep(4,12)),4),rep(rep(3:4,each=16),12))

pal2 <- c("#AECBFA", 
          "#B2DFB2",
          "#FFF8E4", 
          "#F4C7C3")
pal2 <- setNames(pal2, c("1", "2", "3", "4"))

spot.coor <- as.data.frame(expand.grid(y=0:31,x=0:31)[,2:1])

gg <- ggplot(spot.coor, aes(x = x, y = y))
pl <- gg + geom_point(size = 3, 
                      aes(color = as.factor(spot.region))) +
  scale_color_manual(values=pal2) +
  theme(legend.text=element_text(size=15),        
        legend.position = "right", 
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(), 
        axis.ticks = element_blank(),
        panel.background = element_blank()) +
  guides(color = guide_legend(title = "Region", 
                              title.theme = element_text(size = 10),
                              override.aes = list(size = 7)))

print(pl)

ggsave("region.png", pl, width = 5, height = 5)


