## ----ch01,out.width='3.2in'----------------------------------------------
library(crone)
sdata <- load_structure("pinkerton2015")
rtmp <- structure_gauss(sdata,N=1000)
plot(rtmp$x,rtmp$rr,type="l",xlab="x",ylab=expression(rho))

## ----ch02,out.width='3.2in'----------------------------------------------
# Max resolution 1 angstrom (D=1), crystal formed of only 
# one unit cell (Ncell=1)
ltmp <- diffraction(sdata,D=1,Ncell=1)
plot(ltmp$xstar,ltmp$Imod,type="l",
     xlab=expression(paste("x"^"*")),ylab="Intensity")

# Max resolution 1 angstrom, crystal formed by 
# two unit cells (Ncell=2)
ltmp <- diffraction(sdata,D=1,Ncell=2)
plot(ltmp$xstar,ltmp$Imod,type="l",
     xlab=expression(paste("x"^"*")),ylab="Intensity")

# Max resolution 1 angstrom, crystal formed by 
# 10 unit cells (Ncell=10, default value). The number of reciprocal 
# space grid points is inreased to 1001 (n=500 -> 2*n+1=1001; 
# default value is n=100 -> 2*n+1=201)
ltmp <- diffraction(sdata,D=1,n=500)
plot(ltmp$xstar,ltmp$Imod,type="l",
     xlab=expression(paste("x"^"*")),ylab="Intensity")


## ----ch03,out.width='3.2in'----------------------------------------------
# Find all peaks
idx <- local_maxima(ltmp$Imod)

# Mean and standard deviation of electron density
M <- mean(ltmp$Imod)
S <- sd(ltmp$Imod)

# Threshold (1st attempt)
Thr <- M + 0*S
Thr

# Peaks (spots) selection
idx <- local_maxima(ltmp$Imod)
jdx <- which(ltmp$Imod[idx] > Thr)
idx <- idx[jdx]  # New index of selected peaks
length(idx)

# Too many peaks (some should be considered as noise)
plot(ltmp$xstar,ltmp$Imod,type="l",
     xlab=expression(paste("x"^"*")),ylab="Intensity")
points(ltmp$xstar[idx],ltmp$Imod[idx],pch=16,cex=0.65,col=2)

# Threshold (2nd attempt)
Thr <- M + 1*S
Thr

# Peaks (spots) selection
idx <- local_maxima(ltmp$Imod)
jdx <- which(ltmp$Imod[idx] > Thr)
idx <- idx[jdx]  # New index of selected peaks
length(idx)

# Some peaks have been missed
plot(ltmp$xstar,ltmp$Imod,type="l",
     xlab=expression(paste("x"^"*")),ylab="Intensity")
points(ltmp$xstar[idx],ltmp$Imod[idx],pch=16,cex=0.65,col=2)

## ----ch04,out.width='3.2in'----------------------------------------------
# Points for the plot
x <- c(-5,-3,-2,-1,0,1,2,3,5)
y <- ltmp$xstar[idx]
plot(x,y,pch=16,xlab=expression(h),
     ylab=expression(paste("x"^"*")))

# Least squares
model <- lm(y ~ 0+x)  # Origin included
smdl <- summary(model)
smdl

# Fit
abline(model,col=2)

# Unit cell length (approximately 10)
a = 1/smdl$coefficients[1]
a

## ----ch05,out.width='3.2in'----------------------------------------------
# Lattice
hidx <- -10:10
astar <- smdl$coefficients[1]
L <- astar * hidx

# Lattice overlapped to diffraction pattern
plot(ltmp$xstar,ltmp$Imod,type="l",
     xlab=expression(paste("x"^"*")),ylab="Intensity")
abline(v=L,col=3)

## ----ch06,out.width='3.2in'----------------------------------------------
# Miller indices
hidx <- 0:20

# Structure factors
ftmp <- strufac(hidx,sdata)

# What's in the structure factor list
names(ftmp)

# Amplitudes and phases for the Patterson
Pmod <- ftmp$Fmod^2
Ppha <- rep(0,times=length(hidx))

# Patterson as Fourier synthesis
rtmp <- fousynth(sdata$a,Pmod,Ppha,hidx,N=1000)
plot(rtmp$x,rtmp$rr,type="l",xlab="x",ylab="P")

## ----ch07----------------------------------------------------------------
# Generate the symmetry-equivalent
sdata2 <- expand_to_cell(sdata)
sdata2$x0

# Smaller inter-atomic distance
sdata2$x0[1]-(-sdata2$x0[1])

# Larger inter-atomic distance
sdata2$x0[2]-sdata2$x0[1]

# Peaks position
idx <- local_maxima(rtmp$rr)
rtmp$x[idx]

## ----ch08----------------------------------------------------------------
# Structure factors for structure with origin at x=0
hidx = 1:10
ftmp <- strufac(hidx=hidx,sdata=sdata)

# Change of origin at x=a/2=5
sdata2 <- sdata
sdata2$x0 <- sdata$x0-5+10   # Shifted back inside cell
sdata2$x0

# Structure factors for structure with origin at x=a/2
ftmp2 <- strufac(hidx=hidx,sdata=sdata2)

## ----ch09----------------------------------------------------------------
# Phases for structure with origin at x=0
ftmp$Fpha

# Phases for structure with origin at x=a/2
ftmp2$Fpha

## ----ch10----------------------------------------------------------------
# Standard structure factors
hidx <- 1:30
ftmp <- strufac(hidx=hidx,sdata=sdata)
FF <- ftmp$Fmod

# Vectors of sums of scattering factors (this structure is
# made of two carbon atoms)
ff <- sqrt(2)*scafac(h=hidx,sdata$a,sdata$Z,sdata$occ,sdata$B)

# Normalised structure factors
EE <- FF/ff

# Display
for (i in hidx) {
  line <- sprintf("%8.3f  %8.3f\n",FF[i],EE[i])
  cat(line)
}


## ----ch11,out.width='3.2in'----------------------------------------------
# Density for standard structure factors
rtmp <- fousynth(a=sdata$a,Fmod=FF,Fpha=ftmp$Fpha,
                 hidx=hidx,N=1000)

# Density for normalised structure factors
ntmp <- fousynth(a=sdata$a,Fmod=EE,Fpha=ftmp$Fpha,
                 hidx=hidx,N=1000)

# Graphical comparison
plot(rtmp$x,rtmp$rr,type="l",xlab="x",ylab=expression(rho))
points(ntmp$x,ntmp$rr,type="l",col=2)

## ----ch12----------------------------------------------------------------
# Normalised structure factors used
for (h in 1:12) {
  line <- sprintf("%5d   %10.3f\n",h,EE[h])
  cat(line)
}
# Sigma_1 relationships
hMat <- matrix(c(1:6,1:6,2,4,6,8,10,12),nrow=6,ncol=3)
colnames(hMat) <- c("h","h","2h")
PS1plus <- rep(0,length=6)
for (i in 1:6) {
  PS1plus[i] <- 0.5+0.5*tanh(EE[hMat[i,1]]*
                EE[hMat[i,2]]*EE[hMat[i,3]]/sqrt(2))
}
S1_table <- cbind(hMat,PS1plus)
idx <- order(S1_table[,4],decreasing=TRUE)
print(S1_table[idx,])

## ----ch13----------------------------------------------------------------
# Sigma_2 relationships (66 found!)
hMat <- matrix(ncol=3)
for (h in 1:11) {
  for (k in (h+1):12) {
    ll <- h-k
    if (ll >= -12 & ll <= 12) {
      hMat <- rbind(hMat,matrix(c(h,k,ll),nrow=1))
    }
  }
}
hMat <- hMat[-1,]
colnames(hMat) <- c("h","k","h-k")
PS2plus <- rep(0,length=length(hMat[,1]))
for (i in 1:length(hMat[,1])) {
  PS2plus[i] <- 0.5+0.5*tanh(EE[abs(hMat[i,1])]*
                EE[abs(hMat[i,2])]*EE[abs(hMat[i,3])]/sqrt(2))
}
S2_table <- cbind(hMat,PS2plus)
idx <- order(S2_table[,4],decreasing=TRUE)
print(S2_table[idx,])

## ----ch14,out.width='3.2in'----------------------------------------------
# Availabl set of known reflections
hidx <- c(2,3,5,6,8,10,11)
Fmod <- EE[hidx]
Fpha <- rep(0,length=length(Fmod))  # All signs are +
rtmp_test <- fousynth(sdata$a,Fmod=Fmod,
                      Fpha=Fpha,hidx=hidx,N=1000)

# Comparison: some peaks should match the two peaks for the
# structure with the origin at x=0 (rtmp) or the one with the
# origin at x=a/2 (rtmp2)
rtmp <- structure_gauss(sdata=sdata,N=1000)
rtmp2 <- structure_gauss(sdata=sdata2,N=1000)
ym <- min(rtmp$rr,rtmp2$rr,rtmp_test$rr)
yM <- max(rtmp$rr,rtmp2$rr,rtmp_test$rr)
plot(rtmp$x,rtmp$rr,type="l",xlab="x",ylab=expression(rho),
     ylim=c(ym,yM))
points(rtmp2$x,rtmp2$rr,type="l",lty=2)
points(rtmp_test$x,rtmp_test$rr,type="l",col=2)

## ----ch15,out.width='3.2in'----------------------------------------------
# Available set of known reflections
hidx <- c(2,3,5,6,8,10,11)
Fmod <- EE[hidx]
Fpha <- c(180,180,0,0,180,0,0)  # + or - signs
rtmp_test <- fousynth(sdata$a,Fmod=Fmod,
                      Fpha=Fpha,hidx=hidx,N=1000)

# Comparison
rtmp <- structure_gauss(sdata=sdata,N=1000)
rtmp2 <- structure_gauss(sdata=sdata2,N=1000)
ym <- min(rtmp$rr,rtmp2$rr,rtmp_test$rr)
yM <- max(rtmp$rr,rtmp2$rr,rtmp_test$rr)
plot(rtmp$x,rtmp$rr,type="l",xlab="x",ylab=expression(rho),
     ylim=c(ym,yM))
points(rtmp2$x,rtmp2$rr,type="l",lty=2)
points(rtmp_test$x,rtmp_test$rr,type="l",col=2)

## ----ch16----------------------------------------------------------------
idx <- local_maxima(rtmp_test$rr)
for (i in idx) {
  print(rtmp_test$rr[i])
}

# The approximate carbon position corresponds to idx[3]
rtmp_test$rr[idx[3]]
rtmp_test$x[idx[3]]

## ----ch17,out.width='3.2in'----------------------------------------------
# Correct density
rtmp <- fousynth(a=sdata$a,Fmod=ftmp$Fmod[1:12],
                 Fpha=ftmp$Fpha[1:12],hidx=1:12,N=1000)


# Initial (approximate) density
sdata0 <- sdata
sdata0$x0 <- rtmp_test$x[idx[3]]
sdata0$B <- 0
ftmp0 <- strufac(hidx=1:12,sdata=sdata0)
rtmp0 <- fousynth(a=sdata0$a,Fmod=ftmp0$Fmod,
                  Fpha=ftmp0$Fpha,hidx=1:12,N=1000)

# Compare plots
ym <- min(rtmp0$rr)
yM <- max(rtmp0$rr)
plot(rtmp$x,rtmp$rr,type="l",xlab="x",ylab=expression(rho),
     ylim=c(ym,yM))
points(rtmp0$x,rtmp0$rr,type="l",col=2)

## ----ch18----------------------------------------------------------------
# Preliminaries...
# 1) Miller indices
hidx <- 1:12

# 2) Observed amplitudes
Fo <- ftmp$Fmod[1:12]

# 3) Initial xc and Bc
xc <- sdata0$x0
Bc <- 0

# 4) Structure to be updated
sdata0 <- sdata

# 5) Calculated structure factors 
sdata0$x0 <- xc
sdata0$B <- Bc
ftmp0 <- strufac(hidx=1:12,sdata=sdata0)
Fc <- ftmp0$Fmod
Sns <- cos(ftmp0$Fpha*pi/180)

# Initial residual
RR <- 100*sum(abs(Fo-Fc))/sum(Fo)
RR

# One cycle of refinement ...
# Scattering factors
ff <- scafac(h=hidx,a=sdata$a,Zj=6,occj=1,Bj=Bc)

# y part of the regression
y <- Fo-Fc

# x1 part of the regression
x1 <- (-4*pi/a)*hidx*ff*exp(-hidx^2*Bc/(4*a^2))*
  sin(2*pi*hidx*xc/a)*Sns

# x2 part of the regression
x2 <- (-hidx^2/(2*a^2))*ff*exp(-hidx^2*Bc/(4*a^2))*
  cos(2*pi*hidx*xc/a)*Sns

# Least-Squares regression
model <- lm(y~0+x1+x2)
smm <- summary(model)

# Update model
xc <- xc+smm$coefficients[1,1]
Bc <- Bc+smm$coefficients[2,1]
xc
Bc
sdata0$x0 <- xc
sdata0$B <- Bc
ftmp0 <- strufac(hidx=1:12,sdata=sdata0)
Fc <- ftmp0$Fmod
Sns <- cos(ftmp0$Fpha*pi/180)

# Updated residual
RR <- 100*sum(abs(Fo-Fc))/sum(Fo)
RR

## ----ch19----------------------------------------------------------------
# One cycle of refinement ...
# Scattering factors
ff <- scafac(h=hidx,a=sdata$a,Zj=6,occj=1,Bj=Bc)

# y part of the regression
y <- Fo-Fc

# x1 part of the regression
x1 <- (-4*pi/a)*hidx*ff*exp(-hidx^2*Bc/(4*a^2))*
  sin(2*pi*hidx*xc/a)*Sns

# x2 part of the regression
x2 <- (-hidx^2/(2*a^2))*ff*exp(-hidx^2*Bc/(4*a^2))*
  cos(2*pi*hidx*xc/a)*Sns

# Least-Squares regression
model <- lm(y~0+x1+x2)
smm <- summary(model)

# Update model
xc <- xc+smm$coefficients[1,1]
Bc <- Bc+smm$coefficients[2,1]
xc
Bc
sdata0$x0 <- xc
sdata0$B <- Bc
ftmp0 <- strufac(hidx=1:12,sdata=sdata0)
Fc <- ftmp0$Fmod
Sns <- cos(ftmp0$Fpha*pi/180)

# Updated residual
RR <- 100*sum(abs(Fo-Fc))/sum(Fo)
RR

## ----ch20,out.width='3.2in'----------------------------------------------
# True density
ftmpT <- strufac(hidx=1:12,sdata=sdata)
rtmpT <- fousynth(a=sdata$a,Fmod=ftmpT$Fmod,Fpha=ftmpT$Fpha,
                  hidx=1:12,N=1000)

# Calculated (final) density
rtmpC <- fousynth(a=sdata$a,Fmod=ftmp0$Fmod,Fpha=ftmp0$Fpha,
                  hidx=1:12,N=1000)

# Compare densities
ym <- min(rtmpT$rr,rtmpC$rr)
yM <- max(rtmpT$rr,rtmpC$rr)
plot(rtmpT$x,rtmpT$rr,type="l",xlab="x",ylab=expression(rho),
     ylim=c(ym,yM))
points(rtmpC$x,rtmpC$rr,type="l",col=2)

