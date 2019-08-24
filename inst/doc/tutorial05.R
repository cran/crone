## ----ch01,out.width='3.2in',echo=FALSE-----------------------------------
library(crone)
sdata <- load_structure("thiocyanate")

# Full structure factor (h=3)
ftmp <- strufac(hidx=3,sdata=sdata)

# Contribution from Sulphur
sdataS <- sdata
sdataS$x0 <- sdata$x0[1]
sdataS$Z <- sdata$Z[1]
sdataS$B <- sdata$B[1]
sdataS$occ <- sdata$occ[1]
ftmpS <- strufac(hidx=3,sdata=sdataS)

# Contribution from carbon
sdataC <- sdata
sdataC$x0 <- sdata$x0[2]
sdataC$Z <- sdata$Z[2]
sdataC$B <- sdata$B[2]
sdataC$occ <- sdata$occ[2]
ftmpC <- strufac(hidx=3,sdata=sdataC)

# Contribution from nitrogen
sdataN <- sdata
sdataN$x0 <- sdata$x0[3]
sdataN$Z <- sdata$Z[3]
sdataN$B <- sdata$B[3]
sdataN$occ <- sdata$occ[3]
ftmpN <- strufac(hidx=3,sdata=sdataN)

# Plot
plot(0,0,type="n",xlim=c(-10,10),ylim=c(-10,10),
     xlab="Re(F)",ylab="Im(F)")
arrows(0,-10,0,10,length=0.15,lwd=2)
arrows(-10,0,10,0,length=0.15,lwd=2)
x <- ftmp$Fmod*cos(ftmp$Fpha*pi/180)
y <- ftmp$Fmod*sin(ftmp$Fpha*pi/180)
arrows(0,0,x,y,length=0.1,angle=20,lwd=2)
xS <- ftmpS$Fmod*cos(ftmpS$Fpha*pi/180)
yS <- ftmpS$Fmod*sin(ftmpS$Fpha*pi/180)
arrows(0,0,xS,yS,length=0.1,angle=20,lwd=1)
xC <- ftmpC$Fmod*cos(ftmpC$Fpha*pi/180)
yC <- ftmpC$Fmod*sin(ftmpC$Fpha*pi/180)
arrows(xS,yS,xS+xC,yS+yC,length=0.1,angle=20,lwd=1)
xN <- ftmpN$Fmod*cos(ftmpN$Fpha*pi/180)
yN <- ftmpN$Fmod*sin(ftmpN$Fpha*pi/180)
arrows(xS+xC,yS+yC,xS+xC+xN,yS+yC+yN,length=0.1,angle=20,lwd=1)
text(x=-5,y=-2,labels=expression(F[3]))
text(x=-4,y=-4,labels="S")
text(x=-5.7,y=-5.5,labels="C")
text(x=-7.7,y=-6,labels="N")

## ----ch02,out.width='3.2in',echo=FALSE-----------------------------------
plot_absorption_curves("Fe",c(0.5,4))

## ----ch03,out.width='3.2in'----------------------------------------------
# Replace the carbon of thiocyanate with Iron, as
# this scatters anomalously in a significant way
sdata$Z
sdata$Z[2] <- 26

# Reduce thermal vibration
sdata$B
sdata$B[2] <- 3

# Plot electron density
rtmp <- structure_gauss(sdata=sdata,N=1000)
plot(rtmp$x,rtmp$rr,type="l",xlab=expression(x),
     ylab=expression(rho))

## ----ch04----------------------------------------------------------------
# Wavelength used when anomalous scattering is significant
out <- fluorescent_scan("Fe")
my_lambda <- out$lambda
my_lambda

# Miller indices used
hidx <- -2:2

# No breaking of Friedel's law (lambda = 0.5). Differences
# in Bijvoet pairs are minimal
ftmp <- strufac(hidx=hidx,sdata=sdata,lbda=0.5,
         anoflag=TRUE,f1f2out=FALSE) # No message printing
ftmp$Fmod[c(2,4)]  # h=-1,1
ftmp$Fmod[c(1,5)]  # h=-2,2

# Breaking of Friedel's law 
# (lambda = 'my_lambda', roughly 1.743)
ftmp <- strufac(hidx=hidx,sdata=sdata,lbda=my_lambda,
         anoflag=TRUE,f1f2out=FALSE) # No message printing
ftmp$Fmod[c(2,4)]  # h=-1,1
ftmp$Fmod[c(1,5)]  # h=-2,2

## ----ch05----------------------------------------------------------------
ftmp <- strufac(hidx=hidx,sdata=sdata,lbda=my_lambda)
ftmp$Fmod[c(2,4)]  # h=-1,1
ftmp$Fmod[c(1,5)]  # h=-2,2

## ----ch06,out.width='3.2in'----------------------------------------------
# New structure (use 'standardise_sdata')
sdata <- standardise_sdata(a=30,
                           SG="P1",
                           x0=c(2,3.5,10,13),
                           Z=c(6,8,26,26),
                           B=c(0.5,0.5,0.5,0.5),
                           occ=c(1,1,1,1))

# Electron density
rtmp <- structure_gauss(sdata=sdata,N=1000)
plot(rtmp$x,rtmp$rr,type="l",xlab=expression(x),
     ylab=expression(rho))

## ----ch07----------------------------------------------------------------
# h from 1 to 80
hidx <- -80:80

# Structure factors with Friedel's law violated
ftmp <- strufac(hidx=hidx,sdata=sdata,anoflag=TRUE,
                lbda=my_lambda)

# F-
Fminus <- ftmp$Fmod[1:80]
Fminus <- rev(Fminus)

# F+
Fplus <- ftmp$Fmod[82:161]

# Averages
Fobs <- 0.5*(Fplus+Fminus)

# Delta_ano
Dano <- abs(Fplus-Fminus)

# Display
line1 <- sprintf("%10.3f %10.3f %10.3f %10.3f",
                Fplus[1],Fminus[1],Fobs[1],Dano[1])
line2 <- sprintf("%10.3f %10.3f %10.3f %10.3f",
                Fplus[2],Fminus[2],Fobs[2],Dano[2])
line3 <- sprintf("%10.3f %10.3f %10.3f %10.3f",
                Fplus[3],Fminus[3],Fobs[3],Dano[3])
line4 <- sprintf("%10.3f %10.3f %10.3f %10.3f",
                Fplus[4],Fminus[4],Fobs[4],Dano[4])
line5 <- sprintf("%10.3f %10.3f %10.3f %10.3f",
                Fplus[5],Fminus[5],Fobs[5],Dano[5])
cat(c(line1,line2,line3,line4,line5),sep="\n")

## ----ch08,out.width='3.2in'----------------------------------------------
# Patterson map
Fmod <- Fobs^2
Fpha <- rep(0,length=80)
Pmap <- fousynth(a=sdata$a,Fmod=Fmod,Fpha=Fpha,
                 hidx=1:80,N=1000)

# Anomalous difference map
Fmod <- (Fplus-Fminus)^2
Anomap <- fousynth(a=sdata$a,Fmod=Fmod,Fpha=Fpha,
                 hidx=1:80,N=1000)

# Scale for anomalous map
KK <- sd(Pmap$rr)/sd(Anomap$rr)

# Compare maps
ym <- min(Pmap$rr,KK*Anomap$rr)
yM <- max(Pmap$rr,KK*Anomap$rr)
plot(Pmap$x,Pmap$rr,type="l",ylim=c(ym,yM))
points(Anomap$x,KK*Anomap$rr,type="l",col=2)

## ----ch09----------------------------------------------------------------
# Patterson map peaks
idxP <- local_maxima(Pmap$rr)
jdxP <- order(Pmap$rr[idxP],decreasing=TRUE)
for (i in 1:16) {
  line <- sprintf("%10.3f\n",Pmap$rr[idxP[jdxP[i]]]) 
  cat(line)
}

## ----ch10----------------------------------------------------------------
# Anomalous map peaks
idxA <- local_maxima(Anomap$rr)
jdxA <- order(Anomap$rr[idxA],decreasing=TRUE)
for (i in 1:16) {
  line <- sprintf("%10.3f\n",Anomap$rr[idxA[jdxA[i]]]) 
  cat(line)
}

## ----ch11----------------------------------------------------------------
# First three anomalous-difference peaks
Anomap$x[idxA[jdxA[c(2,4,6)]]]

# First three Patterson peaks
Pmap$x[idxP[jdxP[c(2,4,6)]]]

## ----ch12,out.width='3.2in'----------------------------------------------
# Updated model
sdataU <- standardise_sdata(a=sdata$a,
                            SG="P1",
                            x0=c(1,4),
                            Z=c(26,26),
                            B=c(0,0),
                            occ=c(1,1))

# Calculated phases
ftmpU <- strufac(hidx=1:80,sdata=sdataU)
Fpha <- ftmpU$Fpha

# New electron density map
rtmpU <- fousynth(a=sdata$a,Fmod=Fobs,Fpha=Fpha,hidx=1:80,N=1000)
plot(rtmpU$x,rtmpU$rr,type="l",xlab=expression(x),
     ylab=expression(rho))

