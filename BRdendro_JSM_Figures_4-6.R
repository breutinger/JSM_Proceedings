# Filename: BRdendro_JSM_Figures_4-6

# Last updated by: Bailey Reutinger
# Last update date: 2024-10-02

# Purpose of this script:
# Code for reproducing JSM Proceedings Figures 4-6


### WARNING: Requires installation of dplR package ###



#### Setup ####

rm(list=ls())
options(width = 80)

filename <-"BRdendro-JSM_Figures_4-6" 

sink(paste(filename, ".log", sep = "")) # create log file to store output

# R version
R.version.string


#### Project Code ####

## 00 load libraries
library(dplR)

#set working directory
setwd ("~/Elon/BRDendro/Data")


## 01 Create Data Frame

## 01.1 Read in csv
EUF <- read.csv("EUF.csv")

## 01.2 Convert csv data to data frame
EUF.df <- as.data.frame(EUF)

## 01.3 Add year vector to data frame
years <- EUF.df$Year.Measure


## 02 Create functions

## 02.1 Calculate the p value
pCalc <- function(NumYears){
  Years67<- NumYears*0.67
  freqVal <- 1/Years67
  numerator <- 6*(cos(2*pi*freqVal)-1)^2
  denominator <- cos(2*pi*freqVal)+2
  pVal <- numerator/denominator
  pVal
}

## 02.2 Calculate spline fits

spl.fit <- function(n, p, y){
  # Fits a smoothing spline to the data for based on the CP methods. 
  
  # Note, the spline fits are for any choice of smoothing parameter. 
  
  # n = sample size                       
  # p = smoothing parameter
  # y = data sequence
  
  # Transform p into lambda
  lambda <- 1 / (2*p)
  
  # Fourier transformaiton of y
  Yomega <- rep(0,n)
  for(j in 1:n){
    Yomega[j] <- sum(y*exp(-2*pi*1i*(0:(n-1))*(j-1)/n))
  }
  
  # Fourier terms
  Wc.omega <- rep(0,n)
  for(j in 1:n){
    Wc.omega[j] <- 4/3 + 2/3*cos(2*pi*(j-1)/n) + 
      2*lambda*(6-8*cos(2*pi*(j-1)/n) + 2*cos(2*2*pi*(j-1)/n))
  }
  
  Z1.omega <- rep(0,n)
  for(j in 1:n){
    Z1.omega[j] <- -1/3 + 2*lambda*(4-exp(1i*2*pi*(j-1)/n))
  }
  
  Z2.omega <- rep(-2*lambda,n)
  
  Zn3.omega <- rep(0,n)
  for(j in 1:n){
    Zn3.omega[j] <- -2*lambda*exp(1i*2*pi*(j-1)/n)
  }
  
  Zn2.omega <- rep(0,n)
  for(j in 1:n){
    Zn2.omega[j] <- -1/3*exp(1i*2*pi*(j-1)/n)+ 2*lambda*
      (4*exp(1i*2*pi*(j-1)/n)-1)
  }
  
  Wy.omega <- rep(0,n)
  for(j in 1:n){
    Wy.omega[j] <- 2*cos(2*pi*(j-1)/n) -2
  }
  
  U0.omega <- rep(0,n)
  for(j in 1:n){
    U0.omega[j] <- 2 - exp(1i*2*pi*(j-1)/n)
  }
  
  U1.omega <- rep(-1,n)
  
  Un2.omega <- rep(0,n)
  for(j in 1:n){
    Un2.omega[j] <- -exp(1i*2*pi*(j-1)/n)
  }
  
  Un1.omega <- rep(0,n)
  for(j in 1:n){
    Un1.omega[j] <- -1 + 2*exp(1i*2*pi*(j-1)/n)
  }
  
  Vc <- matrix(rep(0,4*n), nrow = 4, ncol = n)
  for(i in 1:2){
    for(j in 1:n){
      Vc[i,j] <- exp(1i*2*pi*i*(j-1)/n) * 1/n
    }
  }
  for(j in 1:n){
    Vc[3,j] <- exp(1i*2*pi*(n-3)*(j-1)/n) * 1/n
  }
  for(j in 1:n){
    Vc[4,j] <- exp(1i*2*pi*(n-2)*(j-1)/n) * 1/n
  }
  
  Vy <- matrix(rep(0,4*n), nrow = 4, ncol = n)
  for(i in 1:2){
    for(j in 1:n){
      Vy[i,j] <- exp(1i*2*pi*(i-1)*(j-1)/n) * 1/n
    }
  }
  for(j in 1:n){
    Vy[3,j] <- exp(1i*2*pi*(n-2)*(j-1)/n) * 1/n
  }
  for(j in 1:n){
    Vy[4,j] <- exp(1i*2*pi*(n-1)*(j-1)/n) * 1/n
  }
  
  Wc <- diag(Wc.omega)
  Z <- cbind(Z1.omega, Z2.omega, Zn3.omega, Zn2.omega)
  Wy <- diag(Wy.omega)
  U <- cbind(U0.omega, U1.omega, Un2.omega, Un1.omega)
  
  M <- Wc + Z %*% Vc
  
  T.omega <- rep(0,n)
  for(j in 1:n){
    T.omega[j] <- 2*cos(2*pi*(j-1)/n) - 2
  }
  
  # Solve for Fourier transformation of c (part of cubic spline solution)
  Comega <- solve(M) %*% (Wy + U %*% Vy) %*% Yomega
  
  # Fourier transformation of the spline fit
  Aomega <- Yomega - 2*lambda * Comega * T.omega
  
  # Spline fit in time domain
  a.seq <- rep(0,n)
  for(k in 1:n){
    a.seq[k] <- 1/n*sum(Aomega*exp(2*pi*1i*(0:(n-1))*(k-1)/n))
  }
  
  return(Re(a.seq))
}

## METHOD 1: Create splines for each core and take the average

## 03 Find smoothing parameter (p) for a core

## 03.1  W core
X8.2.W.n <- length(na.omit(EUF.df$`X8.2.W`))
X8.2.W.p <- pCalc(X8.2.W.n)

## 03.2 E core
X8.2.E.n <- length(na.omit(EUF.df$`X8.2.E`))
X8.2.E.p <- pCalc(X8.2.E.n)


## 04 Calculate spline fit 

## 04.1 W core
X8.2.W.spl.Va <- spl.fit(n = X8.2.W.n, 
                         p = X8.2.W.p, 
                         y = na.omit(EUF.df$`X8.2.W`))

## 04.2 E core
X8.2.E.spl.Va <- spl.fit(n = X8.2.E.n, 
                         p = X8.2.E.p, 
                         y = na.omit(EUF.df$`X8.2.E`))


## 05 Standardize the fit

## 05.1 W core
X8.2.W.std <- na.omit(EUF.df$`X8.2.W`) / X8.2.W.spl.Va

## 05.2 E core
X8.2.E.std <- na.omit(EUF.df$`X8.2.E`) / X8.2.E.spl.Va


## METHOD 2: Create a spline of average measurements

## 06 Find average of EUF 8-2 cores
Plot.8.2.avgV2 <- rowMeans(EUF.df[, c("X8.2.W", "X8.2.E")], na.rm = T) 

## 07 Find smoothing parameter (p) for average
Plot.8.2.avgV2.n <- length(na.omit(Plot.8.2.avgV2)) 
Plot.8.2.avgV2.p <- pCalc(Plot.8.2.avgV2.n) 

## 08 Calculate spline fit 
Plot.8.2.avgV2.spl.Va <- spl.fit(n = Plot.8.2.avgV2.n, 
                                 p = Plot.8.2.avgV2.p,
                                 y = na.omit(Plot.8.2.avgV2))

## 09 Standardize the fit
Plot.8.2.avgV2.std <- na.omit(Plot.8.2.avgV2) / Plot.8.2.avgV2.spl.Va 


## 10 Create data frame with cores of interest

## 10.1
EUF.spl.df <- EUF.df[, c("Year.Measure", "X8.2.E", "X8.2.W")] 

## 10.2 for each core, add a vector to the df that contains NA for first part, then 
# adds the spline fit for the rest. The spline fit should align with the 
# raw values at the corresponding years
EUF.spl.df$Plot.8.2.E.spl <- c(rep(NA, length(EUF.spl.df$X8.2.E) - X8.2.E.n), 
                               X8.2.E.spl.Va)
EUF.spl.df$Plot.8.2.W.spl <- c(rep(NA, length(EUF.spl.df$X8.2.W) - X8.2.W.n), 
                               X8.2.W.spl.Va)


## 11 Create the average of the spline fits (stored in the df) using rowMeans
EUF.spl.df$Plot8.2.avg.spl <- 
  rowMeans(EUF.spl.df[, c("Plot.8.2.E.spl", "Plot.8.2.W.spl")], 
                                       na.rm = T)


## 12 Store average of spline fits in subset without missing data present
EUF.8.2.avg.spl <- EUF.spl.df[!is.nan(EUF.spl.df$Plot8.2.avg.spl),
                              c("Year.Measure", "Plot8.2.avg.spl")]


## 13 add the spline of averages to this new subset
EUF.8.2.avg.spl$Plot.8.2.spl.avg <- Plot.8.2.avgV2.spl.Va


## 14 Graph Figures 4-6

## 14.1 Figure 4:Raw tree ring measurements of the west (W) (black line) and 
# east (E) (gray line) cores taken from tree 8-2 in EUF
plot(EUF$X8.2.E ~ EUF$Year, type = "l", col = "gray", 
     xlim = c(1940, 2022), 
     xlab = "Year", ylab = "Ring Width (mm)")
lines(EUF$X8.2.W ~ EUF$Year)

## 14.2 Figure 5: Cubic smoothing spline (black line) created from the average 
# (gray line) of the east (E) and west (W) cores of tree 8-2 in EUF
plot(Plot.8.2.avgV2 ~ years, type = "l", col = "gray", 
     xlim = c(1940, 2022), 
     xlab = "Year", ylab = "Ring Width (mm)")
lines(Plot.8.2.avgV2.spl.Va ~ years)

## 14.3 Figure 6: Spline of the average measurements (black line) plotted
# against the average of splines (gray line) using the W and E cores of tree
# 8-2 in EUF
plot(EUF.8.2.avg.spl$Plot8.2.avg.spl ~ EUF.8.2.avg.spl$Year.Measure, 
     type = "l", col = "gray", 
     xlab = "Year", ylab = "Ring Width")
lines(EUF.8.2.avg.spl$Plot.8.2.spl.avg ~ EUF.8.2.avg.spl$Year.Measure)


## 15 Correlation of 8-2 E and W cores
EUF2 <- read.csv("EUF2.csv")
cor.test(EUF2$X3.3.E, EUF2$X3.3.W)

# close log
sink()

