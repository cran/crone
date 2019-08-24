## ----ch01----------------------------------------------------------------
library(crone)
load_structure()

## ----ch02----------------------------------------------------------------
# Make sure to type the underscore!
sdata <- load_structure("thiocyanate")

# The object returned by load_structure is a named list
class(sdata)
names(sdata)

## ----ch03----------------------------------------------------------------
sdata$a
sdata$SG

## ----ch04----------------------------------------------------------------
sdata$x0
sdata$Z
sdata$B
sdata$occ

## ----ch05,out.width='3.2in'----------------------------------------------
rtmp <- structure_gauss(sdata)
names(rtmp)
plot(rtmp$x,rtmp$rr,type="l",xlab="x",ylab=expression(rho),
     ylim=c(-1,max(rtmp$rr)))
segments(sdata$x[1],0,sdata$x[3],0,lwd=2)
points(sdata$x[1],0,pch=16,cex=3,col="yellow")
points(sdata$x[2],0,pch=16,cex=3,col="grey")
points(sdata$x[3],0,pch=16,cex=3,col="blue")

## ----ch06,out.width='3.2in'----------------------------------------------
x <- seq(0,10*sdata$a,length=1000)
rtmp <- structure_gauss(sdata,x)
plot(rtmp$x,rtmp$rr,type="l",xlab="x",ylab=expression(rho))

## ----ch07----------------------------------------------------------------
hidx <- c(0,1)
ftmp <- strufac(hidx,sdata)
names(ftmp)
ftmp$Fmod  # Structure factors' amplitudes for h=0,1
ftmp$Fpha  # Structure factors' phases for h=0,1
hidx = -2:2
ftmp <- strufac(hidx,sdata)
ftmp$Fmod  # Friedel's law |F(-h)|=|F(h)|
ftmp$Fpha  # Friedel's law phi(-h)=-phi(h)

## ----ch08,out.width='4.5in'----------------------------------------------
hidx <- 0:10
ftmp <- strufac(hidx=hidx,sdata=sdata)

# Grid
N <- 200

# Approximation with hmax=2
Fmod <- ftmp$Fmod[1:3]
Fpha <- ftmp$Fpha[1:3]
hidx <- 0:2
rtmp1 <- fousynth(sdata$a,Fmod,Fpha,hidx,N)

# Approximation with hmax=3
Fmod <- ftmp$Fmod[1:4]
Fpha <- ftmp$Fpha[1:4]
hidx <- 0:3
rtmp2 <- fousynth(sdata$a,Fmod,Fpha,hidx,N)

# Approximation with hmax=5
Fmod <- ftmp$Fmod[1:6]
Fpha <- ftmp$Fpha[1:6]
hidx <- 0:5
rtmp3 <- fousynth(sdata$a,Fmod,Fpha,hidx,N)

# Approximation with hmax=10
Fmod <- ftmp$Fmod[1:11]
Fpha <- ftmp$Fpha[1:11]
hidx <- 0:10
rtmp4 <- fousynth(sdata$a,Fmod,Fpha,hidx,N)

# Limits for plot
lims <- range(rtmp1$rr,rtmp2$rr,rtmp3$rr,rtmp4$rr,rtmp$rr)

# Comparisons
plot(rtmp1$x,rtmp1$rr,type="l",xlab="x",ylab=expression(rho),
     ylim=lims)
points(rtmp2$x,rtmp2$rr,type="l",col=2)
points(rtmp3$x,rtmp3$rr,type="l",col=3)
points(rtmp4$x,rtmp4$rr,type="l",col=4)
points(rtmp$x,rtmp$rr,type="l",col=1,lwd=3,lty=2)

## ----ch09----------------------------------------------------------------
fdata <- load_data(sname="thiocyanate")

# Names of fdata elements
ntmp <- names(fdata)
for (a in ntmp) {
  ltmp <- sprintf("%s  ",a)
  cat(ltmp)
}

# Comparison between observed and calculated amplitudes
fdata$hidx[1:4]
fdata$Fobs[1:4]
ftmp$Fmod[2:5]

# Comparison between stored phases and phases calculated
# with the Fourier synthesis
fdata$Phicalc[1:4]
Fpha[2:5]

