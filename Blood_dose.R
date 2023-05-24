library(dplyr)
library(ggplot2)
#Proton----
a=500
TBV=2
func.W.2.PT <- function(N) {
  V=0.02
  pV=V/(TBV-V)
  s1 <- matrix(0, nrow = N, ncol = a)
  for (i in 1:N){
    s1[1,] <- c(1,rep(0,a-1))
    s1[2,] <- c(1,rep(0,a-1))
    if(i > 2) s1[i,] <- s1[i-1,]*(1-pV)+pV*c(0,s1[i-2,][1:a-1]);
  }
  return(s1)
}

a=500
TBV=2
func.W.1.PT <- function(N) {
  V=0.02
  pV=V/(TBV-V)
  s1 <- matrix(0, nrow = N, ncol = a)
  for (i in 1:N){
    s1[1,] <- c(1,rep(0,a-1))
    s1[2,] <- c(1,rep(0,a-1))
    if(i > 2) s1[i,] <- s1[i-1,]*(1-pV)+pV*c(0,s1[i-2,][1:a-1]);
  }
  return(s1)
}

func.H.PT <- function(N) {
  V=0.02
  pV=V/(TBV-V)
  s1 <- matrix(0, nrow = N, ncol = a)
  for (i in 1:N){
    s1[1,] <- c(1,rep(0,a-1))
    s1[2,] <- c(1,rep(0,a-1))
    if(i > 2) s1[i,] <- s1[i-1,]*(1-pV)+pV*(c(0,s1[i-2,][1:a-1])+s1[i-1,])*0.5;
  }
  return(s1)
}

W.2.PT=data.frame(dose=0, volume=1, c.volume=1, c.volume1=1, n.fraction=0, rate='2 Gy/min')
for (i in 1:60){
  dose=2.5*i
  rate=2
  n.cycle=round(dose/rate*60/4.62)
  data=data.frame(dose=2/60*4.62*c(0:(a-1)),
                  volume=((func.W.2.PT(1000)[n.cycle+3,]*0.9+func.W.2.PT(1000)[n.cycle+2,]*0.1)+
                            (func.W.2.PT(1000)[n.cycle+2,]*0.9+func.W.2.PT(1000)[n.cycle+1,]*0.1))/2)%>%
    mutate(c.volume=cumsum(volume))%>%
    mutate(c.volume1=c(0, c.volume[2:a]), n.fraction=i, rate='2 Gy/min')
  W.2.PT=rbind(data,W.2.PT)
}

W.1.PT=data.frame(dose=0, volume=1, c.volume=1, c.volume1=1, n.fraction=0, rate='1 Gy/min')
for (i in 1:60){
  dose=2.5*i
  rate=1
  n.cycle=round(dose/rate*60/4.62)
  data=data.frame(dose=rate/60*4.62*c(0:(a-1)),
                  volume=((func.W.1.PT(2000)[n.cycle+3,]*0.9+func.W.1.PT(2000)[n.cycle+2,]*0.1)+
                            (func.W.1.PT(2000)[n.cycle+2,]*0.9+func.W.1.PT(2000)[n.cycle+1,]*0.1))/2)%>%
    mutate(c.volume=cumsum(volume))%>%
    mutate(c.volume1=c(0, c.volume[2:a]), n.fraction=i, rate='1 Gy/min')
  W.1.PT=rbind(data,W.1.PT)
}

H.2.PT=data.frame(dose=0, volume=1, c.volume=1, c.volume1=1, n.fraction=0, rate='2 Gy/min')

for (i in 1:60){
  dose=2.5*i
  rate=2
  n.cycle=round(dose/rate*60/4.62)
  data=data.frame(dose=rate/60*4.62*c(0:(a-1)),
                  volume=((func.H.PT(1000)[n.cycle+3,]*0.9+func.H.PT(1000)[n.cycle+2,]*0.1)+
                            (func.H.PT(1000)[n.cycle+2,]*0.9+func.H.PT(1000)[n.cycle+1,]*0.1))/2)%>%
    mutate(c.volume=cumsum(volume))%>%
    mutate(c.volume1=c(0, c.volume[2:a]), n.fraction=i, rate='2 Gy/min')
  H.2.PT=rbind(data,H.2.PT)
}

H.1.PT=data.frame(dose=0, volume=1, c.volume=1, c.volume1=1, n.fraction=0, rate='1 Gy/min')
for (i in 1:60){
  dose=2.5*i
  rate=1
  n.cycle=round(dose/rate*60/4.62)
  data=data.frame(dose=rate/60*4.62*c(0:(a-1)),
                  volume=((func.H.PT(2000)[n.cycle+3,]*0.9+func.H.PT(2000)[n.cycle+2,]*0.1)+
                            (func.H.PT(2000)[n.cycle+2,]*0.9+func.H.PT(2000)[n.cycle+1,]*0.1))/2)%>%
    mutate(c.volume=cumsum(volume))%>%
    mutate(c.volume1=c(0, c.volume[2:a]), n.fraction=i, rate='1 Gy/min')
  H.1.PT=rbind(data,H.1.PT)
}

ggplot()+
  geom_line(data=W.2.PT[W.2.PT$n.fraction==4,], aes(x=dose, y=(1-c.volume1)*100, color=rate), size=1.2)+
  geom_line(data=W.1.PT[W.1.PT$n.fraction==4,], aes(x=dose, y=(1-c.volume1)*100, color=rate), size=1.2)+
  theme_bw()+
  xlim(0,.8)+
  labs(x="Dose (Gy)", y="Blood volume (%)", color=NULL)+ ggtitle("Whole brain - After 4 fractions")

ggplot()+
  geom_line(data=H.2.PT[H.2.PT$n.fraction==4,], aes(x=dose, y=(1-c.volume1)*100, color=rate), size=1.2)+
  geom_line(data=H.1.PT[H.1.PT$n.fraction==4,], aes(x=dose, y=(1-c.volume1)*100, color=rate), size=1.2)+
  theme_bw()+
  xlim(0,.8)+
  labs(x="Dose (Gy)", y="Blood volume (%)", color=NULL)+ ggtitle("Left hemisphere - After 4 fractions")

#X-ray----
a=500
TBV=2
func.W.2.XR <- function(N) {
  V=0.0204
  pV=V/(TBV-V)
  s1 <- matrix(0, nrow = N, ncol = a)
  for (i in 1:N){
    s1[1,] <- c(1,rep(0,a-1))
    s1[2,] <- c(1,rep(0,a-1))
    if(i > 2) s1[i,] <- s1[i-1,]*(1-pV)+pV*c(0,s1[i-2,][1:a-1]);
  }
  return(s1)
}

a=500
TBV=2
func.W.1.XR <- function(N) {
  V=0.0204
  pV=V/(TBV-V)
  s1 <- matrix(0, nrow = N, ncol = a)
  for (i in 1:N){
    s1[1,] <- c(1,rep(0,a-1))
    s1[2,] <- c(1,rep(0,a-1))
    if(i > 2) s1[i,] <- s1[i-1,]*(1-pV)+pV*c(0,s1[i-2,][1:a-1]);
  }
  return(s1)
}


func.H.XR <- function(N) {
  V=0.0204
  pV=V/(TBV-V)
  s1 <- matrix(0, nrow = N, ncol = a)
  for (i in 1:N){
    s1[1,] <- c(1,rep(0,a-1))
    s1[2,] <- c(1,rep(0,a-1))
    if(i > 2) s1[i,] <- s1[i-1,]*(1-pV)+pV*(c(0,s1[i-2,][1:a-1])+s1[i-1,])*0.5;
  }
  return(s1)
}

W.2.XR=data.frame(dose=0, volume=1, c.volume=1, c.volume1=1, n.fraction=0, rate='2 Gy/min')
for (i in 1:60){
  dose=2.5*i
  rate=2
  n.cycle=round(dose/rate*60/4.62)
  data=data.frame(dose=2/60*4.72*c(0:(a-1)),
                  volume=((func.W.2.XR(1000)[n.cycle+3,]*0.9+func.W.2.XR(1000)[n.cycle+2,]*0.1)+
                            (func.W.2.XR(1000)[n.cycle+2,]*0.9+func.W.2.XR(1000)[n.cycle+1,]*0.1))/2)%>%
    mutate(c.volume=cumsum(volume))%>%
    mutate(c.volume1=c(0, c.volume[2:a]), n.fraction=i, rate='2 Gy/min')
  W.2.XR=rbind(data,W.2.XR)
}

W.1.XR=data.frame(dose=0, volume=1, c.volume=1, c.volume1=1, n.fraction=0, rate='1 Gy/min')
for (i in 1:60){
  dose=2.5*i
  rate=1
  n.cycle=round(dose/rate*60/4.62)
  data=data.frame(dose=rate/60*4.72*c(0:(a-1)),
                  volume=((func.W.1.XR(2000)[n.cycle+3,]*0.9+func.W.1.XR(2000)[n.cycle+2,]*0.1)+
                            (func.W.1.XR(2000)[n.cycle+2,]*0.9+func.W.1.XR(2000)[n.cycle+1,]*0.1))/2)%>%
    mutate(c.volume=cumsum(volume))%>%
    mutate(c.volume1=c(0, c.volume[2:a]), n.fraction=i, rate='1 Gy/min')
  W.1.XR=rbind(data,W.1.XR)
}

H.2.XR=data.frame(dose=0, volume=1, c.volume=1, c.volume1=1, n.fraction=0, rate='2 Gy/min')
for (i in 1:60){
  dose=2.5*i
  rate=2
  n.cycle=round(dose/rate*60/4.62)
  data=data.frame(dose=rate/60*4.72*c(0:(a-1)),
                  volume=((func.H.XR(1000)[n.cycle+3,]*0.9+func.H.XR(1000)[n.cycle+2,]*0.1)+
                            (func.H.XR(1000)[n.cycle+2,]*0.9+func.H.XR(1000)[n.cycle+1,]*0.1))/2)%>%
    mutate(c.volume=cumsum(volume))%>%
    mutate(c.volume1=c(0, c.volume[2:a]), n.fraction=i, rate='2 Gy/min')
  H.2.XR=rbind(data,H.2.XR)
}

H.1.XR=data.frame(dose=0, volume=1, c.volume=1, c.volume1=1, n.fraction=0, rate='1 Gy/min')
for (i in 1:60){
  dose=2.5*i
  rate=1
  n.cycle=round(dose/rate*60/4.62)
  data=data.frame(dose=rate/60*4.72*c(0:(a-1)),
                  volume=((func.H.XR(2000)[n.cycle+3,]*0.9+func.H.XR(2000)[n.cycle+2,]*0.1)+
                            (func.H.XR(2000)[n.cycle+2,]*0.9+func.H.XR(2000)[n.cycle+1,]*0.1))/2)%>%
    mutate(c.volume=cumsum(volume))%>%
    mutate(c.volume1=c(0, c.volume[2:a]), n.fraction=i, rate='1 Gy/min')
  H.1.XR=rbind(data,H.1.XR)
}

ggplot()+
  geom_line(data=W.2.XR[W.2.XR$n.fraction==4,], aes(x=dose, y=(1-c.volume1)*100, color=rate), size=1.2)+
  geom_line(data=W.1.XR[W.1.XR$n.fraction==4,], aes(x=dose, y=(1-c.volume1)*100, color=rate), size=1.2)+
  theme_bw()+
  xlim(0,.8)+
  labs(x="Dose (Gy)", y="Blood volume (%)", color=NULL)+ ggtitle("Whole brain - After 4 fractions")

ggplot()+
  geom_line(data=H.2.XR[H.2.XR$n.fraction==4,], aes(x=dose, y=(1-c.volume1)*100, color=rate), size=1.2)+
  geom_line(data=H.1.XR[H.1.XR$n.fraction==4,], aes(x=dose, y=(1-c.volume1)*100, color=rate), size=1.2)+
  theme_bw()+
  xlim(0,.8)+
  labs(x="Dose (Gy)", y="Blood volume (%)", color=NULL)+ ggtitle("Left hemisphere - After 4 fractions")

#Parameter extraction----
H.1.PT.4=H.1.PT%>%filter(n.fraction==4)
sum.H.1.PT.4=data.frame(RT='PT', TT='Hemisphere', RR=1,
                        V0.0=1-H.1.PT.4[H.1.PT.4$dose==max(H.1.PT.4[H.1.PT.4$dose==0,]$dose),]$c.volume,
                        V0.1=1-H.1.PT.4[H.1.PT.4$dose==max(H.1.PT.4[H.1.PT.4$dose<0.1,]$dose),]$c.volume,
                        V0.2=1-H.1.PT.4[H.1.PT.4$dose==max(H.1.PT.4[H.1.PT.4$dose<0.2,]$dose),]$c.volume,
                        V0.3=1-H.1.PT.4[H.1.PT.4$dose==max(H.1.PT.4[H.1.PT.4$dose<0.3,]$dose),]$c.volume,
                        V0.4=1-H.1.PT.4[H.1.PT.4$dose==max(H.1.PT.4[H.1.PT.4$dose<0.4,]$dose),]$c.volume,
                        V0.5=1-H.1.PT.4[H.1.PT.4$dose==max(H.1.PT.4[H.1.PT.4$dose<0.5,]$dose),]$c.volume)
H.2.PT.4=H.2.PT%>%filter(n.fraction==4)
sum.H.2.PT.4=data.frame(RT='PT', TT='Hemisphere', RR=2,
                        V0.0=1-H.2.PT.4[H.2.PT.4$dose==max(H.2.PT.4[H.2.PT.4$dose==0,]$dose),]$c.volume,
                        V0.1=1-H.2.PT.4[H.2.PT.4$dose==max(H.2.PT.4[H.2.PT.4$dose<0.1,]$dose),]$c.volume,
                        V0.2=1-H.2.PT.4[H.2.PT.4$dose==max(H.2.PT.4[H.2.PT.4$dose<0.2,]$dose),]$c.volume,
                        V0.3=1-H.2.PT.4[H.2.PT.4$dose==max(H.2.PT.4[H.2.PT.4$dose<0.3,]$dose),]$c.volume,
                        V0.4=1-H.2.PT.4[H.2.PT.4$dose==max(H.2.PT.4[H.2.PT.4$dose<0.4,]$dose),]$c.volume,
                        V0.5=1-H.2.PT.4[H.2.PT.4$dose==max(H.2.PT.4[H.2.PT.4$dose<0.5,]$dose),]$c.volume)
W.2.PT.4=W.2.PT%>%filter(n.fraction==4)
sum.W.2.PT.4=data.frame(RT='PT', TT='Whole brain', RR=2,
                        V0.0=1-W.2.PT.4[W.2.PT.4$dose==max(W.2.PT.4[W.2.PT.4$dose==0,]$dose),]$c.volume,
                        V0.1=1-W.2.PT.4[W.2.PT.4$dose==max(W.2.PT.4[W.2.PT.4$dose<0.1,]$dose),]$c.volume,
                        V0.2=1-W.2.PT.4[W.2.PT.4$dose==max(W.2.PT.4[W.2.PT.4$dose<0.2,]$dose),]$c.volume,
                        V0.3=1-W.2.PT.4[W.2.PT.4$dose==max(W.2.PT.4[W.2.PT.4$dose<0.3,]$dose),]$c.volume,
                        V0.4=1-W.2.PT.4[W.2.PT.4$dose==max(W.2.PT.4[W.2.PT.4$dose<0.4,]$dose),]$c.volume,
                        V0.5=1-W.2.PT.4[W.2.PT.4$dose==max(W.2.PT.4[W.2.PT.4$dose<0.5,]$dose),]$c.volume)
W.1.PT.4=W.1.PT%>%filter(n.fraction==4)
sum.W.1.PT.4=data.frame(RT='PT', TT='Whole brain', RR=1,
                        V0.0=1-W.1.PT.4[W.1.PT.4$dose==max(W.1.PT.4[W.1.PT.4$dose==0,]$dose),]$c.volume,
                        V0.1=1-W.1.PT.4[W.1.PT.4$dose==max(W.1.PT.4[W.1.PT.4$dose<0.1,]$dose),]$c.volume,
                        V0.2=1-W.1.PT.4[W.1.PT.4$dose==max(W.1.PT.4[W.1.PT.4$dose<0.2,]$dose),]$c.volume,
                        V0.3=1-W.1.PT.4[W.1.PT.4$dose==max(W.1.PT.4[W.1.PT.4$dose<0.3,]$dose),]$c.volume,
                        V0.4=1-W.1.PT.4[W.1.PT.4$dose==max(W.1.PT.4[W.1.PT.4$dose<0.4,]$dose),]$c.volume,
                        V0.5=1-W.1.PT.4[W.1.PT.4$dose==max(W.1.PT.4[W.1.PT.4$dose<0.5,]$dose),]$c.volume)

H.1.XR.4=H.1.XR%>%filter(n.fraction==4)
sum.H.1.XR.4=data.frame(RT='XR', TT='Hemisphere', RR=1,
                        V0.0=1-H.1.XR.4[H.1.XR.4$dose==max(H.1.XR.4[H.1.XR.4$dose==0,]$dose),]$c.volume,
                        V0.1=1-H.1.XR.4[H.1.XR.4$dose==max(H.1.XR.4[H.1.XR.4$dose<0.1,]$dose),]$c.volume,
                        V0.2=1-H.1.XR.4[H.1.XR.4$dose==max(H.1.XR.4[H.1.XR.4$dose<0.2,]$dose),]$c.volume,
                        V0.3=1-H.1.XR.4[H.1.XR.4$dose==max(H.1.XR.4[H.1.XR.4$dose<0.3,]$dose),]$c.volume,
                        V0.4=1-H.1.XR.4[H.1.XR.4$dose==max(H.1.XR.4[H.1.XR.4$dose<0.4,]$dose),]$c.volume,
                        V0.5=1-H.1.XR.4[H.1.XR.4$dose==max(H.1.XR.4[H.1.XR.4$dose<0.5,]$dose),]$c.volume)
H.2.XR.4=H.2.XR%>%filter(n.fraction==4)
sum.H.2.XR.4=data.frame(RT='XR', TT='Hemisphere', RR=2,
                        V0.0=1-H.2.XR.4[H.2.XR.4$dose==max(H.2.XR.4[H.2.XR.4$dose==0,]$dose),]$c.volume,
                        V0.1=1-H.2.XR.4[H.2.XR.4$dose==max(H.2.XR.4[H.2.XR.4$dose<0.1,]$dose),]$c.volume,
                        V0.2=1-H.2.XR.4[H.2.XR.4$dose==max(H.2.XR.4[H.2.XR.4$dose<0.2,]$dose),]$c.volume,
                        V0.3=1-H.2.XR.4[H.2.XR.4$dose==max(H.2.XR.4[H.2.XR.4$dose<0.3,]$dose),]$c.volume,
                        V0.4=1-H.2.XR.4[H.2.XR.4$dose==max(H.2.XR.4[H.2.XR.4$dose<0.4,]$dose),]$c.volume,
                        V0.5=1-H.2.XR.4[H.2.XR.4$dose==max(H.2.XR.4[H.2.XR.4$dose<0.5,]$dose),]$c.volume)
W.2.XR.4=W.2.XR%>%filter(n.fraction==4)
sum.W.2.XR.4=data.frame(RT='XR', TT='Whole brain', RR=2,
                        V0.0=1-W.2.XR.4[W.2.XR.4$dose==max(W.2.XR.4[W.2.XR.4$dose==0,]$dose),]$c.volume,
                        V0.1=1-W.2.XR.4[W.2.XR.4$dose==max(W.2.XR.4[W.2.XR.4$dose<0.1,]$dose),]$c.volume,
                        V0.2=1-W.2.XR.4[W.2.XR.4$dose==max(W.2.XR.4[W.2.XR.4$dose<0.2,]$dose),]$c.volume,
                        V0.3=1-W.2.XR.4[W.2.XR.4$dose==max(W.2.XR.4[W.2.XR.4$dose<0.3,]$dose),]$c.volume,
                        V0.4=1-W.2.XR.4[W.2.XR.4$dose==max(W.2.XR.4[W.2.XR.4$dose<0.4,]$dose),]$c.volume,
                        V0.5=1-W.2.XR.4[W.2.XR.4$dose==max(W.2.XR.4[W.2.XR.4$dose<0.5,]$dose),]$c.volume)
W.1.XR.4=W.1.XR%>%filter(n.fraction==4)
sum.W.1.XR.4=data.frame(RT='XR', TT='Whole brain', RR=1,
                        V0.0=1-W.1.XR.4[W.1.XR.4$dose==max(W.1.XR.4[W.1.XR.4$dose==0,]$dose),]$c.volume,
                        V0.1=1-W.1.XR.4[W.1.XR.4$dose==max(W.1.XR.4[W.1.XR.4$dose<0.1,]$dose),]$c.volume,
                        V0.2=1-W.1.XR.4[W.1.XR.4$dose==max(W.1.XR.4[W.1.XR.4$dose<0.2,]$dose),]$c.volume,
                        V0.3=1-W.1.XR.4[W.1.XR.4$dose==max(W.1.XR.4[W.1.XR.4$dose<0.3,]$dose),]$c.volume,
                        V0.4=1-W.1.XR.4[W.1.XR.4$dose==max(W.1.XR.4[W.1.XR.4$dose<0.4,]$dose),]$c.volume,
                        V0.5=1-W.1.XR.4[W.1.XR.4$dose==max(W.1.XR.4[W.1.XR.4$dose<0.5,]$dose),]$c.volume)

sum.CT=data.frame(RT='CT', TT='Anesthesia', RR=0,
                  V0.0=0,
                  V0.1=0,
                  V0.2=0,
                  V0.3=0,
                  V0.4=0,
                  V0.5=0)
sum=rbind(sum.H.1.PT.4, sum.H.2.PT.4, sum.W.1.PT.4, sum.W.2.PT.4,
          sum.H.1.XR.4, sum.H.2.XR.4, sum.W.1.XR.4, sum.W.2.XR.4, sum.CT)
write.csv(sum, file = 'TPS.data.csv')

