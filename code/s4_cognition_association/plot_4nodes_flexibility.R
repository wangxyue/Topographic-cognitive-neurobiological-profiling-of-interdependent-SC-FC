library(ggplot2)

setwd("G:/interdependent_SCFCnet/multi_modanalysis/Rplotres/plot_flexibility")

df<-read.csv("data.csv")
head(df)
df$variable <- factor(df$variable,levels = c('Low','Moderate','Good','High'))

colors <-c(rgb(78/255, 98/255, 171/255),rgb(70/255,157/255,180/255),rgb(135/255,209/255,163/255),rgb(254/255,233/255,154/255))

ggplot(df, aes(x=variable, y=value,fill=variable)) +
  geom_violin(color = "black", width = 0.8, alpha = 0.7,trim = FALSE) +
  geom_boxplot(width = 0.3, aes(fill = variable), color = "black", lwd = 1.7)+
  geom_jitter(width=0.1,size=4)+
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors)+
  labs(x = "Neurocognitive flexibility", y = "Multilayer modular variability")+
  theme_classic()+
  theme(panel.grid=element_blank())+
  theme(legend.title=element_blank(),legend.position ="none")+
  theme(axis.ticks.length=unit(1.5,'cm'),axis.ticks=element_line(size=3),
        axis.line = element_line(size = 3))+
  theme(axis.text.x = element_text(size = 30),
        axis.text.y = element_text(size = 30),
        axis.title.x = element_text(size = 30),
        axis.title.y = element_text(size = 30))+
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2))


ggsave("G:/interdependent_SCFCnet/multi_modanalysis/analysis_resfigs/flexibility_4nodes.tiff",width=20,height=16,dpi=300)