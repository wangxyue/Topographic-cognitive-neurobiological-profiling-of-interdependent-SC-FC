library(ggplot2)
library(dplyr)
library(ggsignif)
library(readxl)

"%||%" <- function(a, b) {
  if (!is.null(a)) a else b
}

geom_flat_violin <- function(mapping = NULL, data = NULL, stat = "ydensity",
                             position = "dodge", trim = TRUE, scale = "area",
                             show.legend = NA, inherit.aes = TRUE, ...) {
  layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomFlatViolin,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      trim = trim,
      scale = scale,
      ...
    )
  )
}

GeomFlatViolin <-
  ggproto("GeomFlatViolin", Geom,
          setup_data = function(data, params) {
            data$width <- data$width %||%
              params$width %||% (resolution(data$x, FALSE) * 0.9)
            data %>%
              group_by(group) %>%
              mutate(ymin = min(y),
                     ymax = max(y),
                     xmin = x,
                     xmax = x + width / 2)
            
          },
          
          draw_group = function(data, panel_scales, coord) {
            
            data <- transform(data, xminv = x,
                              xmaxv = x + violinwidth * (xmax - x))
            
            newdata <- rbind(plyr::arrange(transform(data, x = xminv), y),
                             plyr::arrange(transform(data, x = xmaxv), -y))
            
            
            newdata <- rbind(newdata, newdata[1,])
            
            ggplot2:::ggname("geom_flat_violin", GeomPolygon$draw_panel(newdata, panel_scales, coord))
          },
          
          draw_key = draw_key_polygon,
          
          default_aes = aes(weight = 1, colour = "grey20", fill = "white", size = 0.5,
                            alpha = NA, linetype = "solid"),
          
          required_aes = c("x", "y")
  ) 

setwd("G:/interdependent_SCFCnet/multi_modanalysis/Rplotres/plot_reliability")
data<-read.csv("reliability.csv")
head(data)
df<-reshape2::melt(data,id='sample')
head(df)

df$variable <- factor(df$variable,levels = c('Intra-individual','Inter-individual'))
p <- ggplot(df, aes(x=variable, y=value,fill=variable,colour=variable))

p + geom_flat_violin(aes(fill=variable),
                     position=position_nudge(x=.25)) +
  geom_jitter(aes(color=variable), width=0.1,size=7) +
  geom_boxplot(width=.1,position=position_nudge(x=0.2),outlier.shape = NA,
               aes(fill=variable),size=1.5,color='black')+ 
  scale_fill_manual(values=c("Inter-individual"="#9DCD82","Intra-individual"="#F8B62D"))+
  scale_color_manual(values=c("Inter-individual"="#9DCD82","Intra-individual"="#F8B62D"))+
  scale_y_continuous(expand=c(0,0),limits=c(0,1),breaks=seq(0,1,0.2))+
  labs(x=" ",y="R")+
  ggtitle(" ")+
  theme_bw()+
  theme(panel.grid=element_blank())+
  theme(legend.title=element_blank(),legend.position ="none",
        legend.key.size=unit(0.2,'cm'),
        legend.text = element_text(family = "serif",color="black",size=5,face = "bold"),
        axis.title.x = element_text(family = "serif",color="black",size=25,face = "bold"),
        axis.title.y = element_text(family = "serif",color="black",size=25,face = "bold"),
        axis.text = element_text(color="black",size=19),
        axis.text.x= element_text(family = "serif",size = 19,face = "bold"),
        axis.text.y= element_text(family = "serif",size = 19,face = "bold"),
        axis.line=element_line(color="black",size=3.5),
        axis.ticks.x = element_line(colour = "black",size = 2.5),
        axis.ticks.y = element_line(colour = "black",size = 2.5),
        axis.ticks.length.x = unit(1.5,'cm'),
        axis.ticks.length.y = unit(1.5,'cm'))+
  theme(panel.border = element_rect(colour="white",size=.3),axis.title = element_blank())+
  geom_signif(comparisons = list(c("Inter-individual", "Intra-individual")),y_position = 0.9, 
              map_signif_level = T, textsize = 15,color = "black",tip_length = c(0, 0))+
  annotate('text',x=0.8,y=1.05,label='*** p < 0.001',size=6,color="black",fontface="bold")

ggsave("G:/interdependent_SCFCnet/multi_modanalysis/analysis_resfigs/MVreliability.tiff",width=20,height=16,dpi=300)



