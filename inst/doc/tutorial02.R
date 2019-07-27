## ----ch01,out.width='3.2in'----------------------------------------------
library(crone)
sdata <- load_structure("carbon_dioxide")
sdata$a  # Length of unit cell
rtmp <- structure_gauss(sdata)
sdata$Z  # Only two atoms in this structure?
plot(rtmp$x,rtmp$rr,type="l",xlab="x",ylab=expression(rho),
     ylim=c(-1,max(rtmp$rr)))
segments(sdata$x[1],0,sdata$x[2],0,lwd=2)
points(sdata$x[1],0,pch=16,cex=3,col="red")
points(sdata$x[2],0,pch=16,cex=3,col="grey")
abline(v=0.5*sdata$a,lty=2,col="blue")

## ----ch02----------------------------------------------------------------
sdata$x0
sdata$Z
sdata$occ

## ----ch03----------------------------------------------------------------
xx <- sdata$a-sdata$x0[1]
xx

# The two oxygens have same distance from the carbon
sdata$x0[2]-sdata$x0[1]
xx-sdata$x0[2]

## ----ch04,out.width='3.2in'----------------------------------------------
sdata_exp <- expand_to_cell(sdata)

# We can see that there are two oxygens and one carbon now
sdata_exp$Z

plot(rtmp$x,rtmp$rr,type="l",xlab="x",ylab=expression(rho),
     ylim=c(-1,max(rtmp$rr)))
segments(sdata_exp$x[1],0,sdata_exp$x[3],0,lwd=2)
points(sdata_exp$x[1],0,pch=16,cex=3,col="red")
points(sdata_exp$x[2],0,pch=16,cex=3,col="grey")
points(sdata_exp$x[3],0,pch=16,cex=3,col="red")
abline(v=0.5*sdata_exp$a,lty=2,col="blue")

## ----ch05,out.width='3.2in'----------------------------------------------
sdata$occ
sdata$occ[2] <- 1  # Change carbon occupancy to 1
rtmp <- structure_gauss(sdata)
plot(rtmp$x,rtmp$rr,type="l",xlab="x",ylab=expression(rho),
     ylim=c(-1,max(rtmp$rr)))

## ----ch06,out.width='3.2in'----------------------------------------------
sdata$occ
sdata$occ[2] <- 0.3  # Change carbon occupancy to 0.3
rtmp <- structure_gauss(sdata)
plot(rtmp$x,rtmp$rr,type="l",xlab="x",ylab=expression(rho),
     ylim=c(-1,max(rtmp$rr)))

## ----ch07----------------------------------------------------------------
hidx = -3:3
ftmp <- strufac(hidx,sdata)

# Structure factors corresponding to -1 and 1
c(ftmp$Fmod[3],ftmp$Fmod[5])
c(ftmp$Fpha[3],ftmp$Fpha[5])

# Structure factors corresponding to -2 and 2
c(ftmp$Fmod[2],ftmp$Fmod[6])
c(ftmp$Fpha[2],ftmp$Fpha[6])

# Structure factors corresponding to -2 and 2
c(ftmp$Fmod[1],ftmp$Fmod[7])
c(ftmp$Fpha[1],ftmp$Fpha[7])

## ----ch08,out.width='3.8in'----------------------------------------------
sdata <- list()

# Unit cell length is fixed at 50 angstroms
sdata$a <- 50 

# Symmetry
sdata$SG <- "P-1"

# Populate half the cell
set.seed(9196)
sdata$x0 <- runif(5,min=0,max=0.5*sdata$a)

# Light elements
set.seed(9197)
sdata$Z <- sample(c(3:9,11:17),size=5,replace=FALSE)
sdata$Z

# Cold atoms
sdata$B <- rep(0,length=5)

# Standard occupancies
sdata$occ <- rep(1,length=5)

# Generate and display electron density
rtmp <- structure_gauss(sdata)
plot(rtmp$x,rtmp$rr,type="l",xlab="x",ylab=expression(rho))

## ----ch09----------------------------------------------------------------
names(atoms)

# Find atom names corresponding to given Z
idx <- c()
for (i in 1:5) {
  jdx <- which(atoms$Z == sdata$Z[i])
  idx <- c(idx,jdx)
}
atoms[idx,]

# These are the elements
atoms$anames[idx]

## ----ch10,out.width='3.2in'----------------------------------------------
xgrid <- seq(-0.5*sdata$a,0.5*sdata$a,length=1000)
rtmp <- structure_gauss(sdata,x=xgrid)
plot(rtmp$x,rtmp$rr,type="l",xlab="x",ylab=expression(rho))

## ----ch11,out.width='3.2in'----------------------------------------------
xgrid <- seq(-0.04*sdata$a,0.04*sdata$a,length=1000)
rtmp <- structure_gauss(sdata,x=xgrid)
plot(rtmp$x,rtmp$rr,type="l",xlab="x",ylab=expression(rho))

