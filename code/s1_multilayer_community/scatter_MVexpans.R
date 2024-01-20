library(ggplot2)

library(ggExtra)

library(cowplot)

setwd("G:/interdependent_SCFCnet/multi_modanalysis/Rplotres/plot_scatter/MV_corr_corexpan")

data<-read.csv("data.csv",row.names = 1,header = T,sep = ",",stringsAsFactors = F)

MV_column<-'MV'
expan_column<-'Cortical.Expansion'

X_column<-data[,MV_column]
Y_column<-data[,expan_column]

cor_R<-cor(x=X_column,y=Y_column,method='pearson')
cor_p<-cor.test(x=X_column,y=Y_column)$p.value
result<-data.frame(source_variable="MV",target_variable="Cortical.Expansion",
                   cor_R=cor_R,cor_p=cor_p)

plot_example<-ggplot(data,aes_string(x="MV",y="Cortical.Expansion"))+
  geom_point(size=10,color='#2570A4',alpha=3)+
  labs(x = "Multilayer modular variability", y = "Cortical expansion")+
  geom_smooth(method = 'lm',se=T,color='black',size=1.0)+
  theme_cowplot()+theme(axis.text=element_text(size=80),margin = margin(t = 10))+
  theme(axis.ticks.length=unit(1.5,'cm'),axis.ticks=element_line(size=3),
        axis.line = element_line(size = 3),panel.grid=element_blank(),
        text = element_text(family = "Arial"),axis.title = element_text(size = 50))+
  scale_x_continuous(limits = c(-3, 3), breaks = seq(-3, 3, by = 2))+
  scale_y_continuous(limits = c(-3, 3), breaks = seq(-3, 3, by = 2))

  
print(plot_example)

ggsave("G:/interdependent_SCFCnet/multi_modanalysis/analysis_resfigs/MVcorexpans.tiff",width=20,height=16,dpi=300)