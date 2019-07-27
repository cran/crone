## ----setup, eval= TRUE, include= FALSE, cache= FALSE, echo= FALSE--------
system ("biber tutorial03")

## ----ch01----------------------------------------------------------------
library(crone)
sdata <- load_structure("cyanate")
sdata

## ----ch02,out.width='3.2in'----------------------------------------------
rtmp <- structure_gauss(sdata,N=1000) # Grid with 1000 points
plot(rtmp$x,rtmp$rr,type="l",xlab="x",ylab=expression(rho))

## ----ch03----------------------------------------------------------------
idx <- local_maxima(rtmp$rr)
idx

## ----ch04,out.width='3.2in'----------------------------------------------
plot(rtmp$x,rtmp$rr,type="l",xlab="x",ylab=expression(rho))
for (i in 1:length(idx)) {
  points(rtmp$x[idx[i]],rtmp$rr[idx[i]],pch=16,cex=1.5,col=2)
}

# Peak coordinates found and comparison with published values
for (i in 1:length(idx)) {
  line <- sprintf("%5.3f  %5.3f\n",rtmp$x[idx[i]],sdata$x0[i])
  cat(line)
}

## ----ch05,out.width='3.2in'----------------------------------------------
# Change all 3 B factors to 0
sdata$B <- c(0,0,0)

# Change all atoms to hydrogens
sdata$Z <- c(1,1,1)
sdata

# Electron density
rtmp <- structure_gauss(sdata,N=1000) # Grid with 10000 point
plot(rtmp$x,rtmp$rr,type="l",xlab="x",ylab=expression(rho))

# Local maxima
idx <- local_maxima(rtmp$rr)
idx
for (i in 1:length(idx)) {
  line <- sprintf("%5.3f  %5.3f\n",rtmp$x[idx[i]],sdata$x0[i])
  cat(line)
}

## ----ch06,out.width='3.2in'----------------------------------------------
sdata <- load_structure("cyanate")
hidx <- 0:20
ftmp <- strufac(hidx,sdata)

# This should coincide with 8(O) + 6(C) + 7(N) = 21
ftmp$Fmod[1]

# Exact analytic electron density
rtmp0 <- structure_gauss(sdata,N=1000)

# Electron density as Fourier synthesis
rtmp1 <- fousynth(sdata$a,ftmp$Fmod,ftmp$Fpha,hidx,N=1000)

# The two densities should overlap nicely
plot(rtmp0$x,rtmp0$rr,type="l",col=2,xlab="x",
     ylab=expression(rho))
points(rtmp1$x,rtmp1$rr,type="l",lty=2)

## ----ch07----------------------------------------------------------------
idx0 <- local_maxima(rtmp0$rr)  # Analytic-density case
idx1 <- local_maxima(rtmp1$rr)  # Fourier-synthesis case
for (i in 1:length(idx)) {
  line <- sprintf("%5.3f  %5.3f  %5.3f\n",
           rtmp0$x[idx0[i]],rtmp1$x[idx1[i]],sdata$x0[i])
  cat(line)
}

## ----ch08----------------------------------------------------------------
fdata <- load_data(sname="cyanate")

# Names of fdata elements
ntmp <- names(fdata)
for (a in ntmp) {
  ltmp <- sprintf("%s  ",a)
  cat(ltmp)
}

# The experimental amplitudes are available only 
# for 10 Miller indices (no h=0 term as it's not
# experimentally available)
hidx <- fdata$hidx
hidx

# Let's re-compute calculated s.f.
ftmp <- strufac(hidx,sdata)

# The observed amplitudes should be different
# from the correct amplitudes.
# The available experimental precision is to 3 decimals
for (i in 1:length(hidx)) {
  line <- sprintf("%10.3f  %10.3f\n",
                   fdata$Fobs[i],ftmp$Fmod[i])
  cat(line)
}

## ----ch09----------------------------------------------------------------
# The phases are the same (as they are calculated)
# The available precision of the phases in the data file
# is to 1 decimal. 
for (i in 1:length(hidx)) {
  line <- sprintf("%10.1f  %10.1f\n",
                   fdata$Phicalc[i],ftmp$Fpha[i])
  cat(line)
}

## ----ch10,out.width='3.2in'----------------------------------------------
hidx <- 1:10

# Correct electron density
rtmp0 <- fousynth(sdata$a,ftmp$Fmod,ftmp$Fpha,hidx,N=1000)

# Biased electron density
rtmp1 <- fousynth(sdata$a,fdata$Fobs,fdata$Phicalc,hidx,N=1000)

# Min and max to include all density in the plot
m <- min(rtmp0$rr,rtmp1$rr)
M <- max(rtmp0$rr,rtmp1$rr)

# Comparison
plot(rtmp0$x,rtmp0$rr,type="l",xlab="x",ylab=expression(rho),
     col=2,ylim=c(m,M))
points(rtmp1$x,rtmp1$rr,type="l",lty=2)

## ----ch11----------------------------------------------------------------
idx0 <- local_maxima(rtmp0$rr)  # Correct amplitudes
idx1 <- local_maxima(rtmp1$rr)  # Observed amplitudes
for (i in 1:length(idx)) {
  line <- sprintf("%5.3f  %5.3f  %5.3f\n",
           rtmp0$x[idx0[i]],rtmp1$x[idx1[i]],sdata$x0[i])
  cat(line)
}

## ----ch12,out.width='3.2in'----------------------------------------------
hidx <- 1:10

# Random errors added to phases
set.seed(8761)
pha_new <- ftmp$Fpha + rnorm(length(ftmp$Fpha),mean=30,sd=10)
idx <- which(pha_new < -180)
if (length(idx) > 0) {
  pha_new[idx] <- pha_new[idx]+360
}
idx <- which(pha_new > 180)
if (length(idx) > 0) {
  pha_new[idx] <- pha_new[idx]-360
}

# Correct electron density
rtmp0 <- fousynth(sdata$a,ftmp$Fmod,ftmp$Fpha,hidx,N=1000)

# Biased electron density
rtmp1 <- fousynth(sdata$a,ftmp$Fmod,pha_new,hidx,N=1000)

# Min and max to include all density in the plot
m <- min(rtmp0$rr,rtmp1$rr)
M <- max(rtmp0$rr,rtmp1$rr)

# Comparison
plot(rtmp0$x,rtmp0$rr,type="l",xlab="x",ylab=expression(rho),
     col=2,ylim=c(m,M))
points(rtmp1$x,rtmp1$rr,type="l",lty=2)

## ----ch13----------------------------------------------------------------
idx0 <- local_maxima(rtmp0$rr)  # Analytic-density case
idx1 <- local_maxima(rtmp1$rr)  # Fourier-synthesis case
for (i in 1:length(idx)) {
  line <- sprintf("%5.3f  %5.3f  %5.3f\n",
           rtmp0$x[idx0[i]],rtmp1$x[idx1[i]],sdata$x0[i])
  cat(line)
}

## ----ch14,out.width='3.2in'----------------------------------------------
# Correct electron density
rtmp0 <- fousynth(sdata$a,ftmp$Fmod,ftmp$Fpha,hidx,N=1000)

# Biased electron density
rtmp1 <- fousynth(sdata$a,fdata$Fobs,pha_new,hidx,N=1000)

# Min and max to include all density in the plot
m <- min(rtmp0$rr,rtmp1$rr)
M <- max(rtmp0$rr,rtmp1$rr)

# Comparison
plot(rtmp0$x,rtmp0$rr,type="l",xlab="x",ylab=expression(rho),
     col=2,ylim=c(m,M))
points(rtmp1$x,rtmp1$rr,type="l",lty=2)

