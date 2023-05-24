library(dplyr)
library(ggplot2)
library(viridis)

data_dleu=read.csv(file = 'data.csv')
dat.XR=data_dleu%>%filter(exp=='XR')

dat.XR.C=dat.XR%>%filter(group=='Control with anesthesia')
dat.XR.C.CD4=dat.XR.C%>%summarise(mean=mean(p.CD4), sd=sd(p.CD4), arm='Control', type='T CD4+ cells')
dat.XR.C.CD8=dat.XR.C%>%summarise(mean=mean(p.CD8), sd=sd(p.CD8), arm='Control', type='T CD8+ cells')
dat.XR.C.B=dat.XR.C%>%summarise(mean=mean(p.B.cells), sd=sd(p.B.cells), arm='Control', type='B cells')
dat.XR.C.NK=dat.XR.C%>%summarise(mean=mean(p.NK.cells), sd=sd(p.NK.cells), arm='Control', type='NK cells')
dat.XR.C.neutro=dat.XR.C%>%summarise(mean=mean(p.neutro), sd=sd(p.neutro), arm='Control', type='Neutrophils')
dat.XR.C.mono=dat.XR.C%>%summarise(mean=mean(p.mono), sd=sd(p.mono), arm='Control', type='Monocytes')

dat.XR.T=dat.XR%>%filter(RR>0)
dat.XR.T.CD4=dat.XR.T%>%summarise(mean=mean(p.CD4), sd=sd(p.CD4), arm='Irradiated', type='T CD4+ cells')
dat.XR.T.CD8=dat.XR.T%>%summarise(mean=mean(p.CD8), sd=sd(p.CD8), arm='Irradiated', type='T CD8+ cells')
dat.XR.T.B=dat.XR.T%>%summarise(mean=mean(p.B.cells), sd=sd(p.B.cells), arm='Irradiated', type='B cells')
dat.XR.T.NK=dat.XR.T%>%summarise(mean=mean(p.NK.cells), sd=sd(p.NK.cells), arm='Irradiated', type='NK cells')
dat.XR.T.neutro=dat.XR.T%>%summarise(mean=mean(p.neutro), sd=sd(p.neutro), arm='Irradiated', type='Neutrophils')
dat.XR.T.mono=dat.XR.T%>%summarise(mean=mean(p.mono), sd=sd(p.mono), arm='Irradiated', type='Monocytes')

dat.XR.sum=rbind(dat.XR.C.CD4, dat.XR.C.CD8, dat.XR.C.B, dat.XR.C.NK, dat.XR.C.neutro, dat.XR.C.mono,
                 dat.XR.T.CD4, dat.XR.T.CD8, dat.XR.T.B, dat.XR.T.NK, dat.XR.T.neutro, dat.XR.T.mono)
dat.XR.sum$type=factor(dat.XR.sum$type, levels = c('T CD4+ cells', 'T CD8+ cells', 'B cells', 'NK cells', 'Neutrophils', 'Monocytes'))
ggplot(data=dat.XR.sum, aes(x=type, y=mean, fill=arm))+ 
  geom_bar(position="dodge", stat="identity") +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                position=position_dodge(.9))+
  xlab("")+
  labs(y='Relative concentration to control group', fill=NULL)+
  scale_fill_viridis(discrete = T, option = "E")+
  theme_bw()+
  ggtitle('X-ray irradiation')

wilcox.test(dat.XR.T$p.CD4,dat.XR.C$p.CD4)
wilcox.test(dat.XR.T$p.CD8,dat.XR.C$p.CD8)
wilcox.test(dat.XR.T$p.B.cells,dat.XR.C$p.B.cells)
wilcox.test(dat.XR.T$p.NK.cells,dat.XR.C$p.NK.cells)
wilcox.test(dat.XR.T$p.neutro,dat.XR.C$p.neutro)
wilcox.test(dat.XR.T$p.mono,dat.XR.C$p.mono)


dat.PT=data_dleu%>%filter(exp=='PT')

dat.PT.C=dat.PT%>%filter(group=='Control with anesthesia')
dat.PT.C.CD4=dat.PT.C%>%summarise(mean=mean(p.CD4), sd=sd(p.CD4), arm='Control', type='T CD4+ cells')
dat.PT.C.CD8=dat.PT.C%>%summarise(mean=mean(p.CD8), sd=sd(p.CD8), arm='Control', type='T CD8+ cells')
dat.PT.C.B=dat.PT.C%>%summarise(mean=mean(p.B.cells), sd=sd(p.B.cells), arm='Control', type='B cells')
dat.PT.C.NK=dat.PT.C%>%summarise(mean=mean(p.NK.cells), sd=sd(p.NK.cells), arm='Control', type='NK cells')
dat.PT.C.neutro=dat.PT.C%>%summarise(mean=mean(p.neutro), sd=sd(p.neutro), arm='Control', type='Neutrophils')
dat.PT.C.mono=dat.PT.C%>%summarise(mean=mean(p.mono), sd=sd(p.mono), arm='Control', type='Monocytes')

dat.PT.T=dat.PT%>%filter(RR>0)
dat.PT.T.CD4=dat.PT.T%>%summarise(mean=mean(p.CD4), sd=sd(p.CD4), arm='Irradiated', type='T CD4+ cells')
dat.PT.T.CD8=dat.PT.T%>%summarise(mean=mean(p.CD8), sd=sd(p.CD8), arm='Irradiated', type='T CD8+ cells')
dat.PT.T.B=dat.PT.T%>%summarise(mean=mean(p.B.cells), sd=sd(p.B.cells), arm='Irradiated', type='B cells')
dat.PT.T.NK=dat.PT.T%>%summarise(mean=mean(p.NK.cells), sd=sd(p.NK.cells), arm='Irradiated', type='NK cells')
dat.PT.T.neutro=dat.PT.T%>%summarise(mean=mean(p.neutro), sd=sd(p.neutro), arm='Irradiated', type='Neutrophils')
dat.PT.T.mono=dat.PT.T%>%summarise(mean=mean(p.mono), sd=sd(p.mono), arm='Irradiated', type='Monocytes')

dat.PT.sum=rbind(dat.PT.C.CD4, dat.PT.C.CD8, dat.PT.C.B, dat.PT.C.NK, dat.PT.C.neutro, dat.PT.C.mono,
                 dat.PT.T.CD4, dat.PT.T.CD8, dat.PT.T.B, dat.PT.T.NK, dat.PT.T.neutro, dat.PT.T.mono)

dat.PT.sum$type=factor(dat.PT.sum$type, levels = c('T CD4+ cells', 'T CD8+ cells', 'B cells', 'NK cells', 'Neutrophils', 'Monocytes'))

ggplot(data=dat.PT.sum, aes(x=type, y=mean, fill=arm))+ 
  geom_bar(position="dodge", stat="identity") +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                position=position_dodge(.9))+
  xlab("")+
  labs(y='Relative concentration to control group', fill=NULL)+
  scale_fill_viridis(discrete = T, option = "E")+
  theme_bw()+
  ggtitle('Proton irradiation')+ylim(0,3)

wilcox.test(dat.PT.T$p.CD4,dat.PT.C$p.CD4)
wilcox.test(dat.PT.T$p.CD8,dat.PT.C$p.CD8)
wilcox.test(dat.PT.T$p.B.cells,dat.PT.C$p.B.cells)
wilcox.test(dat.PT.T$p.NK.cells,dat.PT.C$p.NK.cells)
wilcox.test(dat.PT.T$p.neutro,dat.PT.C$p.neutro)
wilcox.test(dat.PT.T$p.mono,dat.PT.C$p.mono)