library(boot)
library(bnlearn)
library(dplyr)
library(ggplot2)
rel1=as.data.frame(t(XRn(c('p.CD4','p.CD8', 'p.B.cells', 'p.NK.cells', 'p.neutro', 'p.mono'), 2)))
colnames(rel1)=c('from', 'to')
rel2=data.frame(from=rel1$to, to=rel1$from)
rel=rbind(rel1, rel2)
sum=read.csv(file = 'TPS.data.csv')

dat=read.csv(file = 'data.csv')%>%
  select(RT,TT,RR,p.CD4,p.CD8, p.B.cells, p.NK.cells, p.neutro, p.mono)%>%
  merge(sum, by=c('RT', 'TT', 'RR'))%>%
  mutate(LN=ifelse(RT=='XR', ifelse(TT=='Whole brain', 4, 2), 0))
d3=dat%>%select(p.CD4,p.CD8, p.B.cells, p.NK.cells, p.neutro, p.mono)%>%na.omit()
dag1=hc(d3)
plot(dag1)
data.dag=as.data.frame(dag1[["arcs"]])

dat.rel=NULL

for (j in 1:length(rel[,1])){
  rel1.1=rel[j,]
  dat.1=0
  for (i in 1:length(data.dag[,1])){
    a=as.numeric(ifelse(rel1.1[1]==data.dag[i,1]&rel1.1[2]==data.dag[i,2], 1,0))
    dat.1=dat.1+ a
  }
  b=data.frame(from=rel[j,1],
               to=rel[j,2],
               prob=dat.1)
  dat.rel=rbind(dat.rel,b)
}

boot.func=function(data, index){
  dag1=hc(data[index,], blacklist = bl1, whitelist = wl1)
  data.dag=as.data.frame(dag1[["arcs"]])
  dat.rel=NULL
  for (j in 1:length(rel[,1])){
    rel1.1=rel[j,]
    dat.1=0
    for (i in 1:length(data.dag[,1])){
      a=as.numeric(ifelse(rel1.1[1]==data.dag[i,1]&rel1.1[2]==data.dag[i,2], 1,0))
      dat.1=dat.1+ a
    }
    b=data.frame(from=rel[j,1],
                 to=rel[j,2],
                 prob=dat.1)
    dat.rel=rbind(dat.rel,b)
  }
  return(as.numeric(dat.rel[16,3]))
}

#X-Ray----
d3=dat%>%filter(RT=='XR')
dat.boot.XR=NULL
for (x in 1001:2000){
  for (k in 1:length(rel[,1])){
    boot.func=function(data, index){
      dag1=hc(data[index,])
      data.dag=as.data.frame(dag1[["arcs"]])
      dat.rel=NULL
      for (j in 1:length(rel[,1])){
        rel1.1=rel[j,]
        dat.1=0
        for (i in 1:length(data.dag[,1])){
          a=as.numeric(ifelse(rel1.1[1]==data.dag[i,1]&rel1.1[2]==data.dag[i,2], 1,0))
          dat.1=dat.1+ a
        }
        b=data.frame(from=rel[j,1],
                     to=rel[j,2],
                     prob=dat.1)
        dat.rel=rbind(dat.rel,b)
      }
      return(as.numeric(dat.rel[k,3]))
    }
    set.seed(x)
    bootstrap <- boot(data=d3, 
                      statistic=boot.func, R = 1)
    c=data.frame(from=dat.rel[k,1],
                 to=dat.rel[k,2],
                 prob=(bootstrap[["t"]]),
                 seed=x)
    dat.boot.XR=rbind(dat.boot.XR,c)
  }
}
bootstrap.XR=dat.boot.XR
bootstrap.XR[bootstrap.XR == 'p.CD4'] <- 'T CD4+'
bootstrap.XR[bootstrap.XR == 'p.CD8'] <- 'T CD8+'
bootstrap.XR[bootstrap.XR == 'p.B.cells'] <- 'B cells'
bootstrap.XR[bootstrap.XR == 'p.NK.cells'] <- 'NK cells'
bootstrap.XR[bootstrap.XR == 'p.neutro'] <- 'Neutrophils'
bootstrap.XR[bootstrap.XR == 'p.mono'] <- 'Monocytes'

#Cross-validation XR ----
e.XR.base = empty.graph(c("p.CD4", "p.CD8",
                          "p.B.cells", "p.NK.cells",
                          "p.neutro", "p.mono"))
arc.set.XR.base = matrix(c("p.CD4", "p.CD8", 
                           "p.neutro", "p.mono", 
                           "p.NK.cells", "p.neutro",
                           "p.CD8", "p.B.cells",
                           "p.mono","p.CD4"),
                         ncol = 2, byrow = TRUE,
                         dimnames = list(NULL, c("from", "to")))

arcs(e.XR.base) = arc.set.XR.base
plot(e.XR.base)
XR.base.fit=bn.fit(e.XR.base, data=d3)
XR.base.fit
class(XR.base.fit)
XR.base.CD4=bn.cv(data=d3,bn=e.XR.base, method = "k-fold", k = 10, runs = 1000, loss='mse',
                  loss.args=list(target='p.CD4'))
XR.base.CD8=bn.cv(data=d3,bn=e.XR.base, method = "k-fold", k = 10, runs = 1000, loss='mse',
                  loss.args=list(target='p.CD8'))
XR.base.B=bn.cv(data=d3,bn=e.XR.base, method = "k-fold", k = 10, runs = 1000, loss='mse',
                loss.args=list(target='p.B.cells'))
XR.base.NK=bn.cv(data=d3,bn=e.XR.base, method = "k-fold", k = 10, runs = 1000, loss='mse',
                 loss.args=list(target='p.NK.cells'))
XR.base.neutro=bn.cv(data=d3,bn=e.XR.base, method = "k-fold", k = 10, runs = 1000, loss='mse',
                     loss.args=list(target='p.neutro'))
XR.base.mono=bn.cv(data=d3,bn=e.XR.base, method = "k-fold", k = 10, runs = 1000, loss='mse',
                   loss.args=list(target='p.mono'))
plot(XR.base.CD4, XR.base.CD8, XR.base.B, XR.base.NK, XR.base.neutro, XR.base.mono,
     xlab=c('T CD4+', 'T CD8+', 'B cells', 'NK cells',
            'Neutrophils', 'Monocytes'))

#Proton----
d3=dat%>%filter(RT=='PT')
dat.boot.PT=NULL
for (x in 1001:2000){
  for (k in 1:length(rel[,1])){
    boot.func=function(data, index){
      dag1=hc(data[index,])
      data.dag=as.data.frame(dag1[["arcs"]])
      dat.rel=NULL
      for (j in 1:length(rel[,1])){
        rel1.1=rel[j,]
        dat.1=0
        for (i in 1:length(data.dag[,1])){
          a=as.numeric(ifelse(rel1.1[1]==data.dag[i,1]&rel1.1[2]==data.dag[i,2], 1,0))
          dat.1=dat.1+ a
        }
        b=data.frame(from=rel[j,1],
                     to=rel[j,2],
                     prob=dat.1)
        dat.rel=rbind(dat.rel,b)
      }
      return(as.numeric(dat.rel[k,3]))
    }
    set.seed(x)
    bootstrap <- boot(data=d3, 
                      statistic=boot.func, R = 1)
    c=data.frame(from=dat.rel[k,1],
                 to=dat.rel[k,2],
                 prob=(bootstrap[["t"]]),
                 seed=x)
    dat.boot.PT=rbind(dat.boot.PT,c)
  }
}
bootstrap.PT = dat.boot.PT
bootstrap.PT[bootstrap.PT == 'p.CD4'] <- 'T CD4+'
bootstrap.PT[bootstrap.PT == 'p.CD8'] <- 'T CD8+'
bootstrap.PT[bootstrap.PT == 'p.B.cells'] <- 'B cells'
bootstrap.PT[bootstrap.PT == 'p.NK.cells'] <- 'NK cells'
bootstrap.PT[bootstrap.PT == 'p.neutro'] <- 'Neutrophils'
bootstrap.PT[bootstrap.PT == 'p.mono'] <- 'Monocytes'

#Cross-validation PT ----
e.PT.base = empty.graph(c("p.CD4", "p.CD8",
                          "p.B.cells", "p.NK.cells",
                          "p.neutro", "p.mono"))
arc.set.PT.base = matrix(c("p.CD4", "p.CD8", 
                           "p.neutro", "p.mono", 
                           "p.B.cells", "p.NK.cells",
                           "p.CD8", "p.B.cells",
                           "p.CD4","p.neutro"),
                         ncol = 2, byrow = TRUE,
                         dimnames = list(NULL, c("from", "to")))

arcs(e.PT.base) = arc.set.PT.base
plot(e.PT.base)
PT.base.fit=bn.fit(e.PT.base, data=d3)
PT.base.fit
class(PT.base.fit)
PT.base.CD4=bn.cv(data=d3,bn=e.PT.base, method = "k-fold", k = 10, runs = 1000, loss='mse',
                  loss.args=list(target='p.CD4'))
PT.base.CD8=bn.cv(data=d3,bn=e.PT.base, method = "k-fold", k = 10, runs = 1000, loss='mse',
                  loss.args=list(target='p.CD8'))
PT.base.B=bn.cv(data=d3,bn=e.PT.base, method = "k-fold", k = 10, runs = 1000, loss='mse',
                loss.args=list(target='p.B.cells'))
PT.base.NK=bn.cv(data=d3,bn=e.PT.base, method = "k-fold", k = 10, runs = 1000, loss='mse',
                 loss.args=list(target='p.NK.cells'))
PT.base.neutro=bn.cv(data=d3,bn=e.PT.base, method = "k-fold", k = 10, runs = 1000, loss='mse',
                     loss.args=list(target='p.neutro'))
PT.base.mono=bn.cv(data=d3,bn=e.PT.base, method = "k-fold", k = 10, runs = 1000, loss='mse',
                   loss.args=list(target='p.mono'))
plot(PT.base.CD4, PT.base.CD8, PT.base.B, PT.base.NK, PT.base.neutro, PT.base.mono,
     xlab=c('T CD4+', 'T CD8+', 'B cells', 'NK cells',
            'Neutrophils', 'Monocytes'))

#Control----
d3=dat%>%filter(RT=='CT')
dat.boot.CT=NULL
for (x in 1001:2000){
  for (k in 1:length(rel[,1])){
    boot.func=function(data, index){
      dag1=hc(data[index,])
      data.dag=as.data.frame(dag1[["arcs"]])
      dat.rel=NULL
      for (j in 1:length(rel[,1])){
        rel1.1=rel[j,]
        dat.1=0
        for (i in 1:length(data.dag[,1])){
          a=as.numeric(ifelse(rel1.1[1]==data.dag[i,1]&rel1.1[2]==data.dag[i,2], 1,0))
          dat.1=dat.1+ a
        }
        b=data.frame(from=rel[j,1],
                     to=rel[j,2],
                     prob=dat.1)
        dat.rel=rbind(dat.rel,b)
      }
      return(as.numeric(dat.rel[k,3]))
    }
    set.seed(x)
    bootstrap <- boot(data=d3, 
                      statistic=boot.func, R = 1)
    c=data.frame(from=dat.rel[k,1],
                 to=dat.rel[k,2],
                 prob=(bootstrap[["t"]]),
                 seed=x)
    dat.boot.CT=rbind(dat.boot.CT,c)
  }
}

#Cross-validation CT----
e.CT = empty.graph(c("p.CD4", "p.CD8",
                     "p.B.cells", "p.NK.cells",
                     "p.neutro", "p.mono"))
arc.set.CT = matrix(c("p.CD4", "p.CD8", 
                      "p.neutro", "p.mono", 
                      "p.B.cells", "p.NK.cells",
                      "p.CD8", "p.B.cells"),
                    ncol = 2, byrow = TRUE,
                    dimnames = list(NULL, c("from", "to")))

arcs(e.CT) = arc.set.CT
plot(e.CT)
CT.fit=bn.fit(e.CT, data=d3)
CT.fit
class(CT.fit)
CT.CD4=bn.cv(data=d3,bn=e.CT, method = "k-fold", k = 10, runs = 1000, loss='mse',
             loss.args=list(target='p.CD4'))
CT.CD8=bn.cv(data=d3,bn=e.CT, method = "k-fold", k = 10, runs = 1000, loss='mse',
             loss.args=list(target='p.CD8'))
CT.B=bn.cv(data=d3,bn=e.CT, method = "k-fold", k = 10, runs = 1000, loss='mse',
           loss.args=list(target='p.B.cells'))
CT.NK=bn.cv(data=d3,bn=e.CT, method = "k-fold", k = 10, runs = 1000, loss='mse',
            loss.args=list(target='p.NK.cells'))
CT.neutro=bn.cv(data=d3,bn=e.CT, method = "k-fold", k = 10, runs = 1000, loss='mse',
                loss.args=list(target='p.neutro'))
CT.mono=bn.cv(data=d3,bn=e.CT, method = "k-fold", k = 10, runs = 1000, loss='mse',
              loss.args=list(target='p.mono'))
plot(CT.CD4, CT.CD8, CT.B, CT.NK, CT.neutro, CT.mono,
     xlab=c('T CD4+', 'T CD8+', 'B cells', 'NK cells',
            'Neutrophils', 'Monocytes'))
