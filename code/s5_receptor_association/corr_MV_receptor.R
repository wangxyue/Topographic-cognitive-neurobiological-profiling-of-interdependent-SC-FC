library(ggplot2)

setwd("G:/interdependent_SCFCnet/multi_modanalysis/Rplotres/plot_receptor_MV")
data<-read.csv("corrR.csv",row.names = 1,header = T,sep = ",",stringsAsFactors = F)

target_variable<-"MV"
mydata<-data[data$variable==target_variable,]
head(mydata)

mydata$receptor<-reorder(mydata$receptor,mydata$Cor)

plot<-ggplot(mydata,aes(x=Cor,y=receptor))+
  geom_segment(aes(x=0,xend=Cor,y=receptor,yend=receptor),color="black", size = 2)+
  geom_point(aes(size=abs(Cor),colour=p.value),alpha=0.8)+
  scale_colour_gradient(low="#339D5A",high="#EDB306")+
  scale_size_continuous(range=c(10,20))+
  theme_minimal()+
  scale_x_continuous(limits=c(-0.2,0.6),breaks = seq(-0.2,0.6,0.2))+
  theme(axis.line=element_line(size=1.5),panel.grid=element_blank(),
        axis.text=element_text(size=12,face="bold"),
        axis.text.y=element_text(hjust=1),
        axis.ticks.length=unit(1,'cm'),
        axis.ticks=element_line(size=1.5),
        plot.title=element_text(size=14,face="bold",hjust=0.5),
        axis.title=element_text(size=14,face="bold"))+
  theme(legend.key.size = unit(1.5, "cm"),legend.text = element_text(size = 25))+
  labs(xlab="Cor")

plot

ggsave("G:/interdependent_SCFCnet/multi_modanalysis/analysis_resfigs/corrR_MV_receptor.tiff",width=10,height=17,dpi=300)

